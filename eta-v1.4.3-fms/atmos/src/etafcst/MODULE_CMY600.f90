    MODULE CMY600
!>--------------------------------------------------------------------------------------------------
!> MODULE CMY600
!>
!> USE MODULES: F77KINDS
!> 
!> DRIVER     : GSMCOLUMN
!>              GSMCONST
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND), PARAMETER :: MY_T1 =  1
    INTEGER(KIND=I4KIND), PARAMETER :: MY_T2 = 35
!
    REAL   (KIND=R4KIND), DIMENSION(MY_T1:MY_T2)                                                ::&
    & MY_GROWTH
!
    END MODULE CMY600
