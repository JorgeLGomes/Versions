    MODULE EXCH_BUF_INTEGER
!>--------------------------------------------------------------------------------------------------
!> MODULE EXCH_BUF_INTEGER
!>
!> USE MODULES: F77KINDS
!>              PARMEXCH
!> 
!> DRIVER     : MODULE_EXCHM
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMEXCH
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND), DIMENSION(IBUFEXCH)                                                   ::&
    & BUF0    , BUF1    , BUF2    , BUF3
!
    END MODULE EXCH_BUF_INTEGER


