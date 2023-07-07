    MODULE BOCO
!>--------------------------------------------------------------------------------------------------
!> MODULE BOCO
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : BOCOH
!>              BOCOHF
!>              BOCOV
!>              CHKOUT
!>              EXIT
!>              INIT
!>              INITS
!>              NEWsFLT
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>              WRTCOM
!>              WRTRST
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : IM, JM, LM
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND), PARAMETER :: LB = 2 * IM + JM - 3
!
    REAL   (KIND=R4KIND), DIMENSION(LB, 2)                                                      ::&
    & PDB
!
    REAL   (KIND=R4KIND), DIMENSION(LB, LM, 2)                                                  ::&
    & TB      , QB      , UB      , VB      , Q2B     , CWMB 
!
    END MODULE BOCO
