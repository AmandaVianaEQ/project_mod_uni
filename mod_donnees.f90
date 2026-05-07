MODULE donnees
    IMPLICIT NONE
CONTAINS

    SUBROUTINE LECTURE_DONNEES(NC, npairs, Pamb, Prup, Qchauffe, Vvap, Uliq0, Kdisc, Tini, Tfin, deltat, &
                              x0, A_ant, B_ant, C_ant, D_ant, E_ant,CPvap, CPliq, Delta_Hvap, Teb, COEFNRTL)
        IMPLICIT NONE

        ! Arguments de sortie
        INTEGER :: NC, npairs
        DOUBLE PRECISION :: Pamb, Prup, Qchauffe, Vvap, Uliq0, Kdisc, Tini, Tfin, deltat
        DOUBLE PRECISION, ALLOCATABLE :: x0(:)
        DOUBLE PRECISION, ALLOCATABLE :: A_ant(:), B_ant(:), C_ant(:), D_ant(:), E_ant(:)
        DOUBLE PRECISION, ALLOCATABLE :: CPvap(:), CPliq(:), Delta_Hvap(:), Teb(:)
        DOUBLE PRECISION, ALLOCATABLE :: COEFNRTL(:,:)

        ! Variables locales
        INTEGER :: i

        ! ----- Lecture de donnees.txt -----
        OPEN(10, FILE="donnees.txt", STATUS="OLD")
            READ(10,*) Pamb
            READ(10,*) Prup
            READ(10,*) Qchauffe
            READ(10,*) Vvap
            READ(10,*) Uliq0
            READ(10,*) Kdisc
            READ(10,*) NC

            npairs = NC*(NC-1)/2

            ALLOCATE(A_ant(NC), B_ant(NC), C_ant(NC), D_ant(NC), E_ant(NC))
            ALLOCATE(x0(NC), CPvap(NC), CPliq(NC), Delta_Hvap(NC), Teb(NC))
            ALLOCATE(COEFNRTL(3, npairs))

            DO i = 1, NC
                READ(10,*) x0(i), A_ant(i), B_ant(i), C_ant(i), D_ant(i), E_ant(i)
            END DO

            READ(10,*) Tini

            DO i = 1, NC
                READ(10,*) CPvap(i), CPliq(i), Delta_Hvap(i), Teb(i)
            END DO

            READ(10,*) Tfin

            READ(10,*) deltat

        CLOSE(10)

        ! ----- Lecture de thermodat.txt -----
        OPEN(20, FILE="thermodat.txt", STATUS="OLD")
            DO i = 1, npairs
                READ(20,*) COEFNRTL(1,i), COEFNRTL(2,i), COEFNRTL(3,i)
            END DO
        CLOSE(20)

    END SUBROUTINE

    SUBROUTINE INTEGRATION(N, NC, X, XDOT, T_disco, Tfin, deltat,ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, MF, &
                       RWORK, LRW, IWORK, LIW,LG, JROOT, BOOL,RPAR, LRP, IPAR, LIP,Pamb, Kdisc, &
                       FX, DFDX, DFDXDOT, GEX, unit_csv)
        IMPLICIT NONE

        ! Arguments
        INTEGER :: N, NC, ITOL, ITASK, ISTATE, IOPT, MF
        INTEGER :: LRW, LIW, LG, LRP, LIP, unit_csv
        DOUBLE PRECISION :: T_disco, Tfin, deltat, RTOL, Pamb, Kdisc
        DOUBLE PRECISION :: X(N), XDOT(N), ATOL(N), RWORK(LRW), BOOL(N)
        DOUBLE PRECISION :: RPAR(LRP)
        INTEGER :: IWORK(LIW), JROOT(LG), IPAR(LIP)
        EXTERNAL :: FX, DFDX, DFDXDOT, GEX

        DOUBLE PRECISION :: TOUT
        INTEGER :: stats

        stats = 1   ! disque intact au démarrage

        DO WHILE (T_disco .LT. Tfin)
            TOUT = MIN(T_disco + deltat, Tfin)

            CALL DISCO(FX, DFDX, DFDXDOT, N, X, XDOT, T_disco, TOUT,ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, &
                       RWORK, LRW, IWORK, LIW, MF, GEX, LG, JROOT, BOOL, RPAR, LRP, IPAR, LIP)

            IF (ISTATE .EQ. 3) THEN
                WRITE(*,'(A,F8.2,A,E11.4,A)') &
                    "RUPTURE a t=", T_disco, "s  P=", X(NC+4), " Pa"
                IPAR(2) = 2

                IF (X(NC+4) > Pamb) THEN
                    X(2*NC+6) = Kdisc * SQRT(X(NC+4) - Pamb)
                ELSE
                    X(2*NC+6) = 0.0D0
                END IF

                ISTATE = 0    ! force DISCO à recalculer XDOT

            ELSE IF (ISTATE .LT. 0) THEN
                WRITE(*,'(A,I4,A,F10.3)') &
                    'ERREUR DISCO : ISTATE =', ISTATE, ' a t =', T_disco
                    WRITE(*,'(A,E14.6)') '  Uliq   = ', X(1)
                    WRITE(*,'(A,E14.6)') '  Uvap   = ', X(NC+2)
                    WRITE(*,'(A,E14.6)') '  T      = ', X(NC+3)
                    WRITE(*,'(A,E14.6)') '  P      = ', X(NC+4)
                    WRITE(*,'(A,E14.6)') '  hliq   = ', X(2*NC+5)
                    WRITE(*,'(A,E14.6)') '  E      = ', X(2*NC+6)
                    CLOSE(unit_csv)
                STOP
            END IF

            ! ----- Écriture des résultats -----
            WRITE(unit_csv,'(F10.2,10(A,E14.6))') T_disco, &
                ';', X(NC+3), &
                ';', X(NC+4), &
                ';', X(1),    &
                ';', X(NC+2), &
                ';', X(2), ';', X(3), ';', X(4), &
                ';', X(NC+5), ';', X(NC+6), ';', X(NC+7)
        END DO

        CLOSE(unit_csv)
    END SUBROUTINE INTEGRATION

END MODULE donnees
