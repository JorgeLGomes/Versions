    MODULE MAPOT
!>--------------------------------------------------------------------------------------------------
!> MODULE MAPOT
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>              RDPARM
!>
!> DRIVER     : INIT
!>              INITS
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : LM , LSM
    USE RDPARM  , ONLY : LP1   
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & LSL     , IXM     , IYM
!
    INTEGER(KIND=I4KIND), DIMENSION(999999)                                                     ::&
    & ISHDE
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & TSPH    , WBD     , SBD     , TLM0D   , TPH0D   , DLMD    , DPHD    , CMLD    , DP30    ,   &
    & X1P     , Y1P     , DISLP   , Z0SLP
!
    REAL   (KIND=R4KIND), DIMENSION(LSM)                                                        ::&
    & SPL     , ALSL
!
    REAL   (KIND=R4KIND), DIMENSION(999999)                                                     ::&
    & TSHDE
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ERLAM0  , CPHI0   , SPHI0

    END MODULE MAPOT
