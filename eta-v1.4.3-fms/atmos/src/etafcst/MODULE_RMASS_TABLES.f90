    MODULE RMASS_TABLES
!>--------------------------------------------------------------------------------------------------
!> MODULE RMASS_TABLES
!>
!> USE MODULES: F77KINDS
!>              RACCR_TABLES
!>
!> DRIVER     : GSMCOLUMN
!>              GSMCONST
!>--------------------------------------------------------------------------------------------------
!
!------------------------------------------------------
! RMASS_TABLES - LOOKUP TABLES FOR MASS CONTENT OF RAIN
!------------------------------------------------------
    USE F77KINDS
    USE RACCR_TABLES, ONLY : MDRMIN, MDRMAX
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(MDRMIN:MDRMAX)                                              ::&
    & MASSR  
!
    END MODULE RMASS_TABLES
