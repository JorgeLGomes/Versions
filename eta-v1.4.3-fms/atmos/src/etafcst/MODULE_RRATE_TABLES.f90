    MODULE RRATE_TABLES
!>--------------------------------------------------------------------------------------------------
!> MODULE RRATE_TABLES
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
!-------------------------------------------------------------
! RRATE_TABLES - LOOKUP TABLES FOR PRECIPITATION RATES OF RAIN
!-------------------------------------------------------------
    REAL   (KIND=R4KIND), DIMENSION(MDRMIN:MDRMAX)                                              ::&
    & RRATE  
!
    END MODULE RRATE_TABLES
