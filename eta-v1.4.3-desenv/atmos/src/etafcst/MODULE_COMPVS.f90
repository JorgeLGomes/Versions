    MODULE COMPVS
!>--------------------------------------------------------------------------------------------------
!> MODULE COMPVS
!>
!> USE MODULES: COMPVS0
!>              F77KINDS
!> 
!> DRIVER     : GSMCONST
!>--------------------------------------------------------------------------------------------------
    USE COMPVS0 , ONLY : NX
    USE F77KINDS
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & C1XPVS  , C2XPVS
!
    REAL   (KIND=R4KIND), DIMENSION(NX)                                                         ::&
    & TBPVS
!
    END MODULE COMPVS
