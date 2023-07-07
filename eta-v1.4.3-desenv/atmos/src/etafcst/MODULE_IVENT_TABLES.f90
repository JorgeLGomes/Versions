    MODULE IVENT_TABLES
!>--------------------------------------------------------------------------------------------------
!> MODULE IVENT_TABLES
!>
!> USE MODULES: F77KINDS
!>              IACCR_TABLES
!>
!> DRIVER     : GSMCOLUMN
!>              GSMCONST
!>--------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------
! IVENT_TABLES - LOOKUP TABLES FOR VENTILATION EFFECTS OF ICE
!------------------------------------------------------------
!
!------------------------------------------------------------------------------------ 
! MASS-WEIGHTED FALL SPEED OF SNOW (LARGE ICE), USED TO CALCULATE PRECIPITATION RATES
!------------------------------------------------------------------------------------ 
    USE F77KINDS
    USE IACCR_TABLES, ONLY : MDIMIN, MDIMAX
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(MDIMIN:MDIMAX)                                              ::&
    & VENTI1  , VENTI2
!
    END MODULE IVENT_TABLES
