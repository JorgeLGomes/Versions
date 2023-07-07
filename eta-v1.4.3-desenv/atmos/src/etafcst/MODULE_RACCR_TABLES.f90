    MODULE RACCR_TABLES
!>--------------------------------------------------------------------------------------------------
!> MODULE RACCR_TABLES
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
!--------------------------------------------------------
! IACCR_TABLES - LOOKUP TABLES FOR ACCRETION RATES OF ICE
!--------------------------------------------------------
!
!------------------------------------------------------------------------
! ACCRI - INTEGRATED QUANTITY ASSOCIATED W/ CLOUD WATER COLLECTION BY ICE
!------------------------------------------------------------------------
!
!-------------------------------------------------------------
! MEAN RAIN DROP DIAMETERS VARY FROM 50 MICRONS TO 450 MICRONS
!-------------------------------------------------------------
    REAL   (KIND=R4KIND), PARAMETER :: DMRMIN =  .05E-3
    REAL   (KIND=R4KIND), PARAMETER :: DMRMAX =  .45E-3
    REAL   (KIND=R4KIND), PARAMETER :: DELDMR = 1.00E-6
    REAL   (KIND=R4KIND), PARAMETER :: XMRMIN = 1.00E6 * DMRMIN 
    REAL   (KIND=R4KIND), PARAMETER :: XMRMAX = 1.00E6 * DMRMAX
!
    INTEGER(KIND=I4KIND), PARAMETER :: MDRMIN = XMRMIN
    INTEGER(KIND=I4KIND), PARAMETER :: MDRMAX = XMRMAX
!
    REAL   (KIND=R4KIND), DIMENSION(MDRMIN:MDRMAX)                                              ::&
    & ACCRR  
!
    END MODULE RACCR_TABLES
