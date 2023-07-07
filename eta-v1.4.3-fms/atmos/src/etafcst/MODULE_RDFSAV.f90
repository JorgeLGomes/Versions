    MODULE RDFSAV
!>--------------------------------------------------------------------------------------------------
!> MODULE RDFSAV
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : GRADFS
!>              RADFS
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & EMISP   , EMIST   , XLATT   , XLATP   , Q19001  , HP98    , H3M6    , HP75    , H6M2    ,   &
    & HP537   , H74E1   , H15E1   , Q14330  , HP2     , TWENTY  , HNINE   , DEGRAD  , HSIGMA  ,   &
    & DAYSEC  , RCO2 
!
    REAL   (KIND=R4KIND), DIMENSION(242)                                                        ::&
    & RACO2
!
    REAL   (KIND=R4KIND), DIMENSION(5)                                                          ::&
    & CAO3SW  , CAH2SW  , CBSW
!
    END MODULE RDFSAV
