SUBROUTINE GEX(N, T, X, XDOT, G, LG, RPAR, LRP, IPAR, LIP)
    IMPLICIT NONE
    INTEGER :: N, LG, LRP, LIP
    DOUBLE PRECISION :: T
    DOUBLE PRECISION :: X(N), XDOT(N), RPAR(LRP)
    INTEGER :: IPAR(LIP)
    DOUBLE PRECISION :: G(LG)
    INTEGER :: NC, stats
    DOUBLE PRECISION :: Prup

    NC    = IPAR(1)
    stats = IPAR(2)
    Prup  = RPAR(2)

    IF (stats == 1) THEN
        G(1) = X(NC+4) - Prup    ! disque intact : événement = (P - Prup)
    ELSE
        G(1) = 1.0D0              ! disque rompu : pas de nouvel événement
    END IF
END SUBROUTINE

