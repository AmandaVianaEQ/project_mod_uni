SUBROUTINE FX(N, T, X, XDOT, F, RPAR, LRP, IPAR, LIP)
    IMPLICIT NONE
    INTEGER :: N, LRP, LIP
    DOUBLE PRECISION :: T, X(N), XDOT(N), F(N), RPAR(LRP)
    INTEGER :: IPAR(LIP)

    INTEGER :: i, NC, npairs, stats
    DOUBLE PRECISION :: Pamb, Prup, Qchauffe, Volvap, Uliq0, Kdisc
    DOUBLE PRECISION :: Uliq, Uvap, Tsys, Psys, hliq, E_deb
    DOUBLE PRECISION :: Vevap, Hvap_mel, sumx, sumy
    DOUBLE PRECISION :: Calc_Psat, Calc_Hliq, Calc_Hvap
    DOUBLE PRECISION, PARAMETER :: Tref  = 150.0D0
    DOUBLE PRECISION, PARAMETER :: R     = 8314.0D0
    DOUBLE PRECISION, PARAMETER :: Rnrtl = 1.98751D0

    DOUBLE PRECISION, ALLOCATABLE :: A_ant(:), B_ant(:), C_ant(:), D_ant(:), E_ant(:)
    DOUBLE PRECISION, ALLOCATABLE :: CPvap(:), CPliq(:), DHvap(:), Teb(:)
    DOUBLE PRECISION, ALLOCATABLE :: CNRTL(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: xliq(:), yvap(:), GAMMA(:), Psat_i(:), Ki(:)
    DOUBLE PRECISION, ALLOCATABLE :: Hvap_i(:), Hliq_i(:)

    ! ---- Récupération paramètres ----
    NC     = IPAR(1)
    stats  = IPAR(2)
    npairs = NC*(NC-1)/2

    Pamb     = RPAR(1)
    Prup     = RPAR(2)
    Qchauffe = RPAR(3)
    Volvap   = RPAR(4)
    Uliq0    = RPAR(5)
    Kdisc    = RPAR(6)

    ALLOCATE(A_ant(NC), B_ant(NC), C_ant(NC), D_ant(NC), E_ant(NC))
    ALLOCATE(CPvap(NC), CPliq(NC), DHvap(NC), Teb(NC))
    ALLOCATE(CNRTL(3, npairs))
    ALLOCATE(xliq(NC), yvap(NC), GAMMA(NC), Psat_i(NC), Ki(NC))
    ALLOCATE(Hvap_i(NC), Hliq_i(NC))

    DO i = 1, NC
        A_ant(i) = RPAR(7+i)
        B_ant(i) = RPAR(7+NC+i)
        C_ant(i) = RPAR(7+2*NC+i)
        D_ant(i) = RPAR(7+3*NC+i)
        E_ant(i) = RPAR(7+4*NC+i)
        CPvap(i) = RPAR(7+5*NC+i)
        CPliq(i) = RPAR(7+6*NC+i)
        DHvap(i) = RPAR(7+7*NC+i)
        Teb(i)   = RPAR(7+8*NC+i)
    END DO

    DO i = 1, npairs
        CNRTL(1,i) = RPAR(7+9*NC + i)
        CNRTL(2,i) = RPAR(7+9*NC + npairs + i)
        CNRTL(3,i) = RPAR(7+9*NC + 2*npairs + i)
    END DO

    ! ---- Lecture du vecteur X ----
    Uliq = X(1)
    DO i = 1, NC
        xliq(i) = X(1+i)
    END DO
    Uvap = X(NC+2)
    Tsys = X(NC+3)
    Psys = X(NC+4)
    DO i = 1, NC
        yvap(i) = X(NC+4+i)
    END DO
    hliq  = X(2*NC+5)
    E_deb = X(2*NC+6)

    ! ---- Débit d'évaporation : Vevap = -dUliq/dt ----
    Vevap = -XDOT(1)

    ! ---- Calculs thermo ----
    DO i = 1, NC
        Psat_i(i) = Calc_Psat(Tsys, A_ant(i), B_ant(i), C_ant(i), D_ant(i), E_ant(i))
    END DO
    CALL NRTL(NC, Tsys, xliq, CNRTL, Rnrtl, GAMMA)
    DO i = 1, NC
        Ki(i)     = GAMMA(i) * Psat_i(i) / Psys
        Hvap_i(i) = Calc_Hvap(Tsys, CPvap(i), CPliq(i), DHvap(i), Teb(i))
        Hliq_i(i) = Calc_Hliq(Tsys, CPliq(i))
    END DO

    Hvap_mel = 0.0D0

    DO i = 1, NC
        Hvap_mel = Hvap_mel + yvap(i) * Hvap_i(i)
    END DO

    sumx = 0.0D0
    sumy = 0.0D0

    DO i = 1, NC
        sumx = sumx + xliq(i)
        sumy = sumy + yvap(i)
    END DO

    ! ===================== RÉSIDUS =====================

    ! F(1) : bilan énergie liquide
    F(1) = hliq*XDOT(1) + Uliq*XDOT(2*NC+5) + Vevap*Hvap_mel - Qchauffe

    ! F(2..NC+1) : bilans partiels liquide
    DO i = 1, NC
        F(1+i) = Uliq*XDOT(1+i) + xliq(i)*XDOT(1) + Vevap*yvap(i)
    END DO

    ! F(NC+2) : bilan matière vapeur
    F(NC+2) = XDOT(NC+2) - Vevap + E_deb

    ! F(NC+3) : sommation
    F(NC+3) = sumy - sumx

    ! F(NC+4) : equation d'etat
    F(NC+4) = Psys*Volvap - Uvap*R*Tsys

    ! F(NC+5..2NC+4) : équilibres yi = Ki*xi
    DO i = 1, NC
        F(NC+4+i) = yvap(i) - Ki(i)*xliq(i)
    END DO

    ! F(2NC+5) : constitutive hliq
    F(2*NC+5) = hliq
    DO i = 1, NC
        F(2*NC+5) = F(2*NC+5) - xliq(i)*CPliq(i)*(Tsys - Tref)
    END DO

    ! F(2NC+6) : disque de rupture
    IF (stats == 1) THEN
        F(2*NC+6) = E_deb
    ELSE
        F(2*NC+6) = E_deb - Kdisc*SQRT(Psys - Pamb)
    END IF

    DEALLOCATE(A_ant, B_ant, C_ant, D_ant, E_ant)
    DEALLOCATE(CPvap, CPliq, DHvap, Teb, CNRTL)
    DEALLOCATE(xliq, yvap, GAMMA, Psat_i, Ki, Hvap_i, Hliq_i)

END SUBROUTINE
