    MODULE SDENS_TABLES
!>--------------------------------------------------------------------------------------------------
!> MODULE SDENS_TABLES
!>
!> USE MODULES: F77KINDS
!>              IACCR_TABLES
!>
!> DRIVER     : GSMCONST  
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE IACCR_TABLES, ONLY : MDIMIN, MDIMAX
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(MDIMIN:MDIMAX)                                              ::&
    & SDENS
!
    END MODULE SDENS_TABLES
