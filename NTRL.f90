SUBROUTINE NTRL(NC,T)
    USE DONNEES
    IMPLICIT NONE
    INTEGER NC,i,j,k,l
    DOUBLE PRECISION R,T
    DOUBLE PRECISION TAU(NC,NC),G(NC,NC)
    DOUBLE PRECISION prodGX(NC),prodTAUGX(NC),LNGAMMA(NC)

    R = 1.98751

    IF (ALLOCATED(GAMMA)) DEALLOCATE(GAMMA)
    ALLOCATE(GAMMA(NC))
    GAMMA = 0.0d0

    DO i=1,NC
        DO j=1,NC
            TAU(i,j)=0
        ENDDO
    ENDDO

    k=0
    DO i=1,NC-1
        DO j=i+1,NC
            k=k+1
            TAU(i,j) = COEFNRTL(1,k)/(R*T)
            TAU(j,i) = COEFNRTL(2,k)/(R*T)
        ENDDO
    ENDDO

    DO i=1,NC
        DO j=1,NC
            G(i,j)=1.0
        END DO
    END DO

    K=0
    DO i=1,NC-1
        DO j=i+1,NC
            k=k+1
            G(i,j)=EXP(-COEFNRTL(3,k)*TAU(i,j))
            G(j,i)=EXP(-COEFNRTL(3,k)*TAU(j,i))
        ENDDO
    ENDDO

    DO i = 1, NC
        G(i,i) = 1.0d0
        TAU(i,i) = 0.0d0
    END DO

    DO i=1,NC
        prodGX(i)=0.0d0
        prodTAUGX(i)=0.0d0
        DO l=1,NC
            prodGX(i)=prodGX(i)+G(l,i)*x(l)
            prodTAUGX(i)=prodTAUGX(i)+TAU(l,i)*G(l,i)*x(l)
        END DO
    END DO

    DO i=1,NC
        LNGAMMA(i)=0.0d0
        IF(prodGX(i)>0.0d0) THEN
            LNGAMMA(i) = prodTAUGX(i)/prodGX(i)
        ENDIF

        DO j=1,NC
            LNGAMMA(i)=LNGAMMA(i)+((x(j)*G(i,j))/prodGX(j))*(TAU(i,j)-(prodTAUGX(j)/prodGX(j)))
        END DO
            GAMMA(i)=EXP(LNGAMMA(i))
             WRITE(*, '(F10.4)') GAMMA(i)
    END DO

    RETURN

END SUBROUTINE

