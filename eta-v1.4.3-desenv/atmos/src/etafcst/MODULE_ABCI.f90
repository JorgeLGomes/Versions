    MODULE ABCI
!>--------------------------------------------------------------------------------------------------
!> MODULE ABCI
!>
!> USE MODULES: F77KINDS
!> 
!> DRIVER     : SFLX
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND), PARAMETER :: NSOLD = 20
!
    REAL   (KIND=R4KIND), DIMENSION(NSOLD)                                                      ::&
    & AI      , BI      , CI
!
    END MODULE ABCI
