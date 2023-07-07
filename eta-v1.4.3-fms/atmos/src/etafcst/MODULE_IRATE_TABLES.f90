    MODULE IRATE_TABLES
!>--------------------------------------------------------------------------------------------------
!> MODULE IRATE_TABLES
!>
!> USE MODULES: F77KINDS
!>              IACCR_TABLES
!>
!> DRIVER     : GSMCOLUMN
!>              GSMCONST
!>--------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------
! IRATE_TABLES - LOOKUP TABLES FOR PRECIPITATION RATES OF ICE
!------------------------------------------------------------
    USE F77KINDS
    USE IACCR_TABLES, ONLY : MDIMIN, MDIMAX
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(MDIMIN:MDIMAX)                                              ::&
    & VSNOWI
!
    END MODULE IRATE_TABLES
