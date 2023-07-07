    MODULE CMICRO_CONS
!>--------------------------------------------------------------------------------------------------
!> MODULE CMICRO_CONS
!>
!> USE MODULES: F77KINDS
!> 
!> DRIVER     : GSMCOLUMN
!>              GSMCONST
!>--------------------------------------------------------------------------------------------------
!
!---------------------------------------------------- 
! CONSTANTS INITIALIZED IN GSMCONST 
!----------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ABFR              , CBFR              ,                                                     &
    & CIACW             , CIACR             ,                                                     &
    & C_N0R0            , CN0R0             ,                                                     &
    & CN0R_DMRMIN       , CN0R_DMRMAX       ,                                                     &
    & CRACW             , CRAUT             ,                                                     &
    & ESW0              , QAUT0             , RFMAX     , RHGRD            ,                      &
    & RQR_DR1           , RQR_DR2           , RQR_DR3   , RQR_DRMIN        , RQR_DRMAX        ,   &
    & RR_DR1            , RR_DR2            , RR_DR3    , RR_DRMIN         , RR_DRMAX         ,   &
    & RHGRDL            , RHGRDS            , VSNOWADST
!
    END MODULE CMICRO_CONS
