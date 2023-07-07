    MODULE COMPVS0
!>--------------------------------------------------------------------------------------------------
!> MODULE COMPVS0
!>
!> USE MODULES: F77KINDS
!> 
!> DRIVER     : GSMCONST
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND), PARAMETER :: NX = 7501
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & C1XPVS0 , C2XPVS0
!
    REAL   (KIND=R4KIND), DIMENSION(NX)                                                         ::&
    & TBPVS0
!
    END MODULE COMPVS0
