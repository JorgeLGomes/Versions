    MODULE IACCR_TABLES
!>--------------------------------------------------------------------------------------------------
!> MODULE IACCR_TABLES
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : GSMCOLUMN
!>              GSMCONST
!>              MODULE_IMASS_TABLES
!>              MODULE_IRATE_TABLES
!>              MODULE_IVENT_TABLES
!>              MODULE_SDENS_TABLES
!>--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------
! IACCR_TABLES - LOOKUP TABLES FOR ACCRETION RATES OF ICE
!--------------------------------------------------------
!
!------------------------------------------------------------------------
! ACCRI - INTEGRATED QUANTITY ASSOCIATED W/ CLOUD WATER COLLECTION BY ICE
!------------------------------------------------------------------------
!
!-----------------------------------------------------------------
! MEAN ICE PARTICLE DIAMETERS VARY FROM 50 MICRONS TO 1000 MICRONS
!-----------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), PARAMETER :: DMIMIN =  .05E-3
    REAL   (KIND=R4KIND), PARAMETER :: DMIMAX = 1.00E-3
    REAL   (KIND=R4KIND), PARAMETER :: DELDMI = 1.00E-6
    REAL   (KIND=R4KIND), PARAMETER :: XMIMIN = 1.00E6  * DMIMIN
    REAL   (KIND=R4KIND), PARAMETER :: XMIMAX = 1.00E6  * DMIMAX
!
    INTEGER(KIND=I4KIND), PARAMETER :: MDIMIN = XMIMIN
    INTEGER(KIND=I4KIND), PARAMETER :: MDIMAX = XMIMAX
!
!
    REAL   (KIND=R4KIND), DIMENSION(MDIMIN:MDIMAX)                                              ::&
    & ACCRI  
!
    END MODULE IACCR_TABLES
