PROGRAM projet_momdelisation
    USE donnees
    IMPLICIT NONE
    INTEGER :: NC
    INTEGER :: stats
    DOUBLE PRECISION :: Pamb, Prup, Qchauffe, Vvap, Uliq0, Kdisc
    DOUBLE PRECISION Calc_Psat
    DOUBLE PRECISION, ALLOCATABLE :: A_ant(:), B_ant(:), C_ant(:), D_ant(:), E_ant(:)
    DOUBLE PRECISION, ALLOCATABLE :: x0(:)
    DOUBLE PRECISION, ALLOCATABLE :: CPvap(:), CPliq(:), Delta_Hvap(:), Teb(:)
    DOUBLE PRECISION, ALLOCATABLE :: COEFNRTL(:,:)
    DOUBLE PRECISION, PARAMETER :: Tref  = 150.0D0
    DOUBLE PRECISION, PARAMETER :: R     = 8314.0D0  !8,314 x 10^3
    DOUBLE PRECISION, PARAMETER :: Rnrtl = 1.98751D0 ! R[cal/mol]

    ! ------- DECLARATION PARAMETERS ET VECTEURS DISCO -----------
    EXTERNAL FX, DFDX, DFDXDOT, GEX

    INTEGER :: N, LRW, LIW, LRP, LIP, LG, npairs

    DOUBLE PRECISION, ALLOCATABLE :: X(:), XDOT(:)
    DOUBLE PRECISION, ALLOCATABLE :: RWORK(:), RPAR(:), BOOL(:), ATOL(:)
    INTEGER, ALLOCATABLE :: IWORK(:), IPAR(:), JROOT(:)

    INTEGER :: ITOL, ITASK, ISTATE, IOPT, MF
    DOUBLE PRECISION :: T_disco, RTOL
! ------------------------------------------------------------

! ------- DECLARATION PARAMETERS DU PROJET -------------------
    INTEGER :: i, j
    DOUBLE PRECISION :: Tini, Tsys, Pini, sumKx, Keq, Psat_val
    DOUBLE PRECISION :: Uliq_ini, Uvap_ini
    DOUBLE PRECISION :: Tfin, deltat

    DOUBLE PRECISION, ALLOCATABLE :: xi_ini(:), yi_ini(:), GAMMA_ini(:)
    DOUBLE PRECISION, ALLOCATABLE :: CNRTL_l(:,:)
    DOUBLE PRECISION :: delta_T

! ------------------------------------------------------------

! ------------------- LECTURE DONNEES -----------------------
    CALL LECTURE_DONNEES(NC, npairs, Pamb, Prup, Qchauffe, Vvap,Uliq0, Kdisc, Tini, Tfin, deltat, &
                         x0, A_ant, B_ant, C_ant, D_ant, E_ant,CPvap, CPliq, Delta_Hvap, Teb, COEFNRTL)

    stats  = 1
    N = 2*NC + 6

! ------------------------------------------------------------
! ------------------- INITIALISATION DISCO -------------------
    LRP = 7+9*NC+3*npairs
    LIP = 2

    ALLOCATE(RPAR(LRP), IPAR(LIP), ATOL(N))
    IPAR(1) = NC
    IPAR(2) = stats

    RPAR(1) = Pamb
    RPAR(2) = Prup
    RPAR(3) = Qchauffe
    RPAR(4) = Vvap
    RPAR(5) = Uliq0
    RPAR(6) = Kdisc

    DO i = 1, NC
        RPAR(7+i)      = A_ant(i)
        RPAR(7+NC+i)   = B_ant(i)
        RPAR(7+2*NC+i) = C_ant(i)
        RPAR(7+3*NC+i) = D_ant(i)
        RPAR(7+4*NC+i) = E_ant(i)
        RPAR(7+5*NC+i) = CPvap(i)
        RPAR(7+6*NC+i) = CPliq(i)
        RPAR(7+7*NC+i) = Delta_Hvap(i)
        RPAR(7+8*NC+i) = Teb(i)
    END DO

    DO i = 1, npairs
        RPAR(7+9*NC + i) = COEFNRTL(1,i)
        RPAR(7+9*NC + npairs + i) = COEFNRTL(2,i)
        RPAR(7+9*NC + 2*npairs + i) = COEFNRTL(3,i)
    END DO

    LG  = 1
    LRW = 278 + 11*N + 5*N*N + 3*LG
    LIW = 251 + 36*N + 5*N*N

    ALLOCATE(X(N), XDOT(N), BOOL(N))
    ALLOCATE(RWORK(LRW), IWORK(LIW), JROOT(LG))
! ------------------------------------------------------------

! --- CONFIGURATION DISCO - TOLERANCES ----------------------
    ITOL = 2
    RTOL = 1.0D-4

    ATOL(1) = 1.0D-3                  ! Uliq
    DO i = 1, NC
        ATOL(1+i) = 1.0D-3            ! xi
    END DO
    ATOL(NC+2) = 1.0D-6               ! Uvap
    ATOL(NC+3) = 1.0D-1               ! T
    ATOL(NC+4) = 1.0D-1               ! P
    DO i = 1, NC
        ATOL(NC+4+i) = 1.0D-3         ! yi
    END DO
    ATOL(2*NC+5) = 1.0D4              ! hliq
    ATOL(2*NC+6) = 1.0D-6             ! E
! ------------------------------------------------------------

! --------- CONFIGURATION DISCO -----------------------------
    ITASK  = 1
    ISTATE = 0          ! DISCO calcule les dérivées initiales automatiquement
    IOPT   = 1
    RWORK  = 0.0D0
    IWORK  = 0
    IWORK(70) = 6       ! Messages d'erreur sur écran
    MF     = 22         ! Gear + Jacobien dense numérique
    BOOL   = 1.0D0
    T_disco = 0.0D0

! ------------------------------------------------------------

! --------- CONDITIONS INITIALES (t = 0) --------------------
    ALLOCATE(xi_ini(NC), yi_ini(NC), GAMMA_ini(NC))
    ALLOCATE(CNRTL_l(3, npairs))
    CNRTL_l = COEFNRTL

    Pini = Pamb
    Uliq_ini = Uliq0

    DO i = 1, NC
        xi_ini(i) = x0(i)
    END DO

! ------------------------------------------------------------

! ----- CALCUL TINI - TEMPERATURE DE BULLE -------------------
        Tsys = Tini   ! valeur de démarrage

        DO i = 1, 100
            CALL NRTL(NC, Tsys, xi_ini, CNRTL_l, Rnrtl, GAMMA_ini)
            sumKx = 0.0D0
            DO j = 1, NC
                Psat_val = Calc_Psat(Tsys, A_ant(j), B_ant(j), C_ant(j), D_ant(j), E_ant(j))
                sumKx = sumKx + GAMMA_ini(j) * xi_ini(j) * Psat_val / Pini
            END DO

            IF (ABS(sumKx - 1.0D0) < 1.0D-4) EXIT

            delta_T = -2.0D0 * LOG(sumKx)
            IF (delta_T >  10.0D0) delta_T =  10.0D0
            IF (delta_T < -10.0D0) delta_T = -10.0D0
            Tsys = Tsys + delta_T
        END DO

        Tini = Tsys

        WRITE(*,'(A,F8.3,A,I3,A)') "Temperature de bulle calculée : Tini = ", Tini, " K (en ", i, " iterations)"

! ------------------------------------------------------------
! ------------- PROPRIETES VAPEUR INITIALES ------------------
    sumKx = 0.0D0
    DO i = 1, NC
        Psat_val  = Calc_Psat(Tsys, A_ant(i), B_ant(i), C_ant(i), D_ant(i), E_ant(i))
        Keq       = GAMMA_ini(i) * Psat_val / Pini
        yi_ini(i) = Keq * xi_ini(i)
        sumKx     = sumKx + yi_ini(i)
    END DO

    DO i = 1, NC
        yi_ini(i) = yi_ini(i) / sumKx
    END DO

    Uvap_ini = Pini * Vvap / (R * Tsys)

! ------------------------------------------------------------

! ------DECLARATIN VECTEUR X --------------------------------

    X(1) = Uliq_ini
    DO i = 1, NC
        X(1+i) = xi_ini(i)
    END DO
    X(NC+2) = Uvap_ini
    X(NC+3) = Tsys
    X(NC+4) = Pini
    DO i = 1, NC
        X(NC+4+i) = yi_ini(i)
    END DO
    X(2*NC+5) = 0.0D0
    DO i = 1, NC
        X(2*NC+5) = X(2*NC+5) + xi_ini(i) * CPliq(i) * (Tsys - Tref)
    END DO
    X(2*NC+6) = 0.0D0
! ------------------------------------------------------------

! -------------- FICHIER DE RÉSULTATS ------------------------
    OPEN(30, FILE='resultats.csv', STATUS='REPLACE')
        WRITE(30,'(A)') 't[s];T[K];P[Pa];Uliq[kmol];Uvap[kmol];x1;x2;x3;y1;y2;y3'
        WRITE(30,'(F10.2,10(A,E14.6))') T_disco, &
            ';', X(NC+3), &
            ';', X(NC+4), &
            ';', X(1),    &
            ';', X(NC+2), &
            ';', X(2), ';', X(3), ';', X(4), &
            ';', X(NC+5), ';', X(NC+6), ';', X(NC+7)

! ------------------------------------------------------------

! ====================== BOUCLE D'INTEGRATION ================
    CALL INTEGRATION(N, NC, X, XDOT, T_disco, Tfin, deltat, &
                     ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, MF, &
                     RWORK, LRW, IWORK, LIW, &
                     LG, JROOT, BOOL, &
                     RPAR, LRP, IPAR, LIP, &
                     Pamb, Kdisc, &
                     FX, DFDX, DFDXDOT, GEX, &
                     30)
! ============================================================

! ----------- LIBERATION MEMOIRE ----------------------------
    DEALLOCATE(X, XDOT, BOOL, ATOL, RWORK, IWORK, JROOT, RPAR, IPAR)
    DEALLOCATE(xi_ini, yi_ini, GAMMA_ini, CNRTL_l)
    DEALLOCATE(A_ant, B_ant, C_ant, D_ant, E_ant)
    DEALLOCATE(x0, CPvap, CPliq, Delta_Hvap, Teb)
    DEALLOCATE(COEFNRTL)
! ------------------------------------------------------------

END PROGRAM

!   FUNCTION CALCUL DE PSAT (ANTOINE)
    DOUBLE PRECISION FUNCTION Calc_Psat(T, A, B, C, D, E)
        IMPLICIT NONE
        DOUBLE PRECISION :: T, A, B, C, D, E
        Calc_Psat = EXP(A + B/T + C*LOG(T) + D*(T**E))
        RETURN
    END FUNCTION Calc_Psat

! ------------------------------------------------------------

!   FUNCTION CALCUL Hliq
    DOUBLE PRECISION FUNCTION Calc_Hliq(T, CPliq)
        DOUBLE PRECISION :: T, CPliq
        DOUBLE PRECISION, PARAMETER :: Tref  = 150.0D0
        Calc_Hliq = CPliq*(T-Tref)
    END FUNCTION Calc_Hliq

! ------------------------------------------------------------

!   FUNCTION CALCUL Hvap
    DOUBLE PRECISION FUNCTION Calc_Hvap(T, CPvap, CPliq, DHvap, Teb)
        DOUBLE PRECISION :: T, CPvap, CPliq, DHvap, Teb
        DOUBLE PRECISION, PARAMETER :: Tref  = 150.0D0
        Calc_Hvap = CPliq*(Teb - Tref) + DHvap + CPvap*(T - Teb)
    END FUNCTION Calc_Hvap

! ------------------------------------------------------------

!   MODELE NRTL
    SUBROUTINE NRTL(NC, T, xi, CNRTL, Rnrtl, GAMMA)
        IMPLICIT NONE
        INTEGER :: NC
        DOUBLE PRECISION :: T, Rnrtl
        DOUBLE PRECISION :: xi(NC), CNRTL(3, NC*(NC-1)/2)
        DOUBLE PRECISION :: GAMMA(NC)

        INTEGER :: i, j, k, l
        DOUBLE PRECISION :: TAU(NC,NC), G(NC,NC)
        DOUBLE PRECISION :: prodGX(NC), prodTAUGX(NC), LNGAMMA(NC)

        TAU = 0.0D0
        G   = 1.0D0

        k = 0
        DO i = 1, NC-1
            DO j = i+1, NC
                k = k+1
                TAU(i,j) = CNRTL(1,k) / (Rnrtl*T)
                TAU(j,i) = CNRTL(2,k) / (Rnrtl*T)
                G(i,j)   = EXP(-CNRTL(3,k)*TAU(i,j))
                G(j,i)   = EXP(-CNRTL(3,k)*TAU(j,i))
            END DO
        END DO

        DO i = 1, NC
            TAU(i,i) = 0.0D0
            G(i,i)   = 1.0D0
        END DO

        DO i = 1, NC
            prodGX(i)   = 0.0D0
            prodTAUGX(i) = 0.0D0
            DO l = 1, NC
                prodGX(i)   = prodGX(i)   + G(l,i)*xi(l)
                prodTAUGX(i) = prodTAUGX(i) + TAU(l,i)*G(l,i)*xi(l)
            END DO
        END DO

        DO i = 1, NC
            LNGAMMA(i) = 0.0D0
            IF (prodGX(i) > 0.0D0) LNGAMMA(i) = prodTAUGX(i)/prodGX(i)
            DO j = 1, NC
                LNGAMMA(i) = LNGAMMA(i) + (xi(j)*G(i,j)/prodGX(j))*(TAU(i,j) - prodTAUGX(j)/prodGX(j))
            END DO
            GAMMA(i) = EXP(LNGAMMA(i))
        END DO

    END SUBROUTINE NRTL

! ------------------------------------------------------------

