    MODULE PARA
!>--------------------------------------------------------------------------------------------------
!> MODULE PARA
!>
!> USE MODULES: F77KINDS
!> 
!> DRIVER     : COLLECT
!>              DIST
!>              MPI_FIRST
!>              QUILT
!>              SFLX
!>              SLP
!>              SLPSIG
!>              SLPSIGSPLINE
!>              UPDATE
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NUM_PROCS         , ME                 ,                                                    &
    & MY_ISD            , MY_IED             ,                                                    &
    & MY_JSD            , MY_JED             ,                                                    &
    & JSTA_I            , JEND_I             ,                                                    &
    & JSTA_IM           , JEND_IM            ,                                                    &
    & JSTA_IM2          , JEND_IM2           ,                                                    &
    & IUP               , IDN
!
    INTEGER(KIND=I4KIND), DIMENSION(0:4095)                                                     ::&
    &  JSTA             , JEND              ,                                                     &
    &  MY_IS_GLB_A      , MY_IE_GLB_A       ,                                                     &
    &  MY_JS_GLB_A      , MY_JE_GLB_A                                                                   
!
    END MODULE PARA
