    MODULE RVELR_TABLES
!>--------------------------------------------------------------------------------------------------
!> MODULE RVELR_TABLES
!>
!> USE MODULES: F77KINDS
!>              RACCR_TABLES
!>
!> DRIVER     : GSMCOLUMN
!>              GSMCONST
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE RACCR_TABLES, ONLY : MDRMIN, MDRMAX
!
    IMPLICIT NONE
!
    SAVE
!----------------------------------------------------- 
! RVELR_TABLES - LOOKUP TABLES FOR FALL SPEEDS OF RAIN
!----------------------------------------------------- 
    REAL   (KIND=R4KIND), DIMENSION(MDRMIN:MDRMAX)                                              ::&
    & VRAIN  
!
    END MODULE RVELR_TABLES
