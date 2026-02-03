SUBROUTINE lecture_donnees(NC,Pext,T)
    USE DONNEES
    IMPLICIT NONE

    INTEGER NC, i
    DOUBLE PRECISION Pext, T

!  LECTURE DES DONNÉES FOURNIES(donnees.txt)
    OPEN(10, FILE="donnees.txt", STATUS="OLD")
        READ(10, *) Pext
        READ(10, *) NC

        IF (ALLOCATED(A)) DEALLOCATE(A, B, C, D, E, x)
        ALLOCATE(A(NC), B(NC), C(NC), D(NC), E(NC), x(NC))

        DO i = 1, NC
            READ(10, *) x(i), A(i), B(i), C(i), D(i), E(i)
        END DO

        READ(10, *) T
    CLOSE(10)

!  LECTURE DES DONNÉES FOURNIES(thermodat.txt)
    OPEN(20, FILE="thermodat.txt", STATUS="OLD")
        IF (ALLOCATED(COEFNRTL)) DEALLOCATE(COEFNRTL)
        ALLOCATE(COEFNRTL(NC, NC))

        DO i=1,NC*(NC-1)/2
            READ(20,*) COEFNRTL(1,i),COEFNRTL(2,i),COEFNRTL(3,i)
        ENDDO
    CLOSE(20)

END SUBROUTINE
