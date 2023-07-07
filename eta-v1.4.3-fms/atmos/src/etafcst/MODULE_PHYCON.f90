    MODULE PHYCON
!>--------------------------------------------------------------------------------------------------
!> MODULE PHYCON
!>
!> USE MODULES: F77KINDS
!> 
!> DRIVER     : CLO89
!>              E1E290
!>              E290
!>              E2SPEC
!>              E3V88
!>              FST88
!>              GFDLRD
!>              HCONST
!>              LWR88
!>              RADFS
!>              SPA88
!>              SWR93
!>              TABLE
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    SAVE
!
!    REAL   (KIND=R4KIND)                                                                        ::&
!    & AMOLWT  , CSUBP   , DIFFCTR , G       , GRAVDR  , O3DIFCTR, P0      , P0XZP2  , P0XZP8  ,   &
!    & P0X2    , RADCON  , RGAS    , RGASSP  , SECPDA  ,                                           &
!    & RATCO2MW, RATH2OMW,                                                                         &
!    & RADCON1 ,                                                                                   &
!    & GINV    , P0INV   , GP0INV
!
    REAL   (KIND=R4KIND), PARAMETER :: AMOLWT   =      28.9644
    REAL   (KIND=R4KIND), PARAMETER :: CSUBP    =       1.00484E7
    REAL   (KIND=R4KIND), PARAMETER :: DIFFCTR  =       1.66
    REAL   (KIND=R4KIND), PARAMETER :: G        =     980.665
    REAL   (KIND=R4KIND), PARAMETER :: GRAVDR   =     980.0
    REAL   (KIND=R4KIND), PARAMETER :: O3DIFCTR =       1.90
    REAL   (KIND=R4KIND), PARAMETER :: P0       = 1013250.    
    REAL   (KIND=R4KIND), PARAMETER :: P0XZP2   =  202649.902
    REAL   (KIND=R4KIND), PARAMETER :: P0XZP8   =  810600.098   
    REAL   (KIND=R4KIND), PARAMETER :: P0X2     =       2.        *        1013250.
    REAL   (KIND=R4KIND), PARAMETER :: RADCON   =       8.427  
    REAL   (KIND=R4KIND), PARAMETER :: RGAS     =       8.3142E7
    REAL   (KIND=R4KIND), PARAMETER :: RGASSP   =       8.31432E7
    REAL   (KIND=R4KIND), PARAMETER :: SECPDA   =       8.64E4  
    REAL   (KIND=R4KIND), PARAMETER :: RATCO2MW =       1.519449738
    REAL   (KIND=R4KIND), PARAMETER :: RATH2OMW =        .622   
    REAL   (KIND=R4KIND), PARAMETER :: RADCON1  =       1.        /        8.427    
    REAL   (KIND=R4KIND), PARAMETER :: GINV     =       1.        /        G
    REAL   (KIND=R4KIND), PARAMETER :: P0INV    =       1.        /        P0    
    REAL   (KIND=R4KIND), PARAMETER :: GP0INV   = GINV            *        P0INV
!         
!------------
! FROM GFDLRD
!------------	
!    DATA G /980.665/
!
    END MODULE PHYCON
