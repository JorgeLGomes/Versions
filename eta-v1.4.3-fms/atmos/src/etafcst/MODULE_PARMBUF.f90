    MODULE PARMBUF
!>--------------------------------------------------------------------------------------------------
!> MODULE PARMBUF
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : MODULE_BUFFER
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
!    SAVE
!
    INTEGER(KIND=I4KIND), PARAMETER :: IBUFMAX = 13*2500000
    INTEGER(KIND=I4KIND), PARAMETER :: IBUFMAX = 2*2500000
    INTEGER(KIND=I4KIND), PARAMETER :: IBUFMAX = 4*2500000
!
    END MODULE PARMBUF
