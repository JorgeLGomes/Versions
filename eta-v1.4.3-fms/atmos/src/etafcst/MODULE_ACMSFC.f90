    MODULE ACMSFC
!>--------------------------------------------------------------------------------------------------
!> MODULE ACMSFC
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : CHKOUT
!>              GOSSIP
!>              INIT
!>              INITS
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>              SURFCE  
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NSRFC
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ASRFC   , TSRFC   , APHTIM 
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & SFCSHX  , SFCLHX  ,                                                                         &
    & SUBSHX  , SNOPCX  ,                                                                         &
    & SFCUVX  , SFCEVP  ,                                                                         &
    & POTEVP  , POTFLX
!
    END MODULE ACMSFC
