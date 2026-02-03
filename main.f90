! MODELISATION DES OPÉRATIONS UNITAIRES – EP8DY3
! Obj:Modéliser le comportement d’un mélange multi-constituant non idéal, placé dans une enceinte chauffée, initialement fermée.
! Amanda Gabriela Viana Fabricio et Mariana Laureano Bomfim Uchoa

MODULE DONNEES
    IMPLICIT NONE
    DOUBLE PRECISION, ALLOCATABLE :: A(:), B(:), C(:), D(:), E(:), x(:)
    DOUBLE PRECISION, ALLOCATABLE :: CPvap(:),CPliq(:),Hvap(:),Teb(:)
    DOUBLE PRECISION, ALLOCATABLE :: COEFNRTL(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: GAMMA(:)
END MODULE

PROGRAM enceinte_chauff
    USE DONNEES
    IMPLICIT NONE
    INTEGER i, j, NC

    DOUBLE PRECISION Pext, T
    DOUBLE PRECISION Calc_Psat
    DOUBLE PRECISION epsilon

    epsilon = 1.0D-4

    CALL lecture_donnees(NC,Pext,T)
    CALL NTRL(NC,T)

END PROGRAM

! DECLARATION FUNCTION Psat(T)
    DOUBLE PRECISION FUNCTION Calc_Psat(T,k)
        USE DONNEES
        INTEGER k
        DOUBLE PRECISION T
        Calc_Psat = exp(A(k) + (B(k)/T) + C(k)*log(T)+D(k)*(T**E(k)))
        RETURN
    END
