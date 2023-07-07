    MODULE RVENT_TABLES
!>--------------------------------------------------------------------------------------------------
!> MODULE RVENT_TABLES
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
    & VENTR1  , VENTR2
!
    END MODULE RVENT_TABLES
