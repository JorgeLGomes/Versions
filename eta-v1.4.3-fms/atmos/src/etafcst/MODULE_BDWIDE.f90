    MODULE BDWIDE
!>--------------------------------------------------------------------------------------------------
!> MODULE BDWIDE
!>
!> ABSTRACT:
!> MODULE BDWIDE CONTAINS RANDOM BAND PARAMETERS FOR SPECIFIC WIDE BANDS.  
!> AT PRESENT, THE INFORMATION CONSISTS OF: 
!> 1) RANDOM MODEL PARAMETERS FOR THE 15 UM BAND,560-800 CM-1; 
!> 2) THE CONTINUUM COEFFICIENT FOR THE 800-990,1070-1200 CM-1 BAND SPECIFICALLY:
!>
!> AWIDE       - RANDOM "A" PARAMETER FOR  BAND
!> BWIDE       - RANDOM "B" PARAMETER FOR  BAND
!> BETAWD      - CONTINUUM COEFFICIENTS FOR BAND
!> APWD ,BPWD  - CAPPHI COEFFICIENTS FOR  BAND
!> ATPWD,BTPWD - CAPPSI COEFFICIENTS FOR BAND
!> BDLOWD      - LOWEST FREQUENCY IN EACH  FREQ  BAND
!> BDHIWD      - HIGHEST FREQUENCY IN EACH FREQ  BAND
!> AB15WD      - THE PRODUCT ARNDM*BRNDM FOR THE ONE BAND REPRESENTING THE 15 UM BAND COMPLEX OF CO2
!> BETINW      - CONT.COEFFICIENT FOR A SPECIFIED WIDE FREQ. BAND (800-990 AND 1070-1200 CM-1).
!> SKO2D       - 1. / BETINW, USED IN SPA88 FOR CONT. COEFFS
!> SKC1R       - BETAWD/BETINW, USED FOR CONT. COEFF. FOR 15 UM BAND IN FST88
!> SKO3R       - RATIO OF CONT. COEFF. FOR 9.9 UM BAND TO BETINW USED FOR 9.6 UM CONT COEFF IN FST88
!>
!> DATA FOR AWIDE, BWIDE, APWD, BPWD, ATPWD, BTPWD, AO3WD, BO3WD ARE OBTAINED BY USING THE AFGL 1982 
!> CATALOG.
!> CONTINUUM COEFFICIENTS ARE FROM ROBERTS 1976.
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : RNDDTA
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & AWIDE   , BWIDE   ,                                                                         &
    & BETAWD  , APWD    , BPWD    , ATPWD   , BTPWD   , BDLOWD  , BDHIWD  ,                       &
    & BETINW  ,                                                                                   &
    & AB15WD  , SKO2D   ,                                                                         &
    & SKC1R   , SKO3R
!------------
! FROM GFDLRD
!------------
    DATA AWIDE  /                                                                                 &
    &  0.309801E+01 /
!
    DATA BWIDE  /                                                                                 &
    &  0.495357E-01 /
!
    DATA BETAWD /                                                                                 &
    &  0.347839E+02 /
!
    DATA APWD   /                                                                                 &
    &  0.177115E-01 /
!
    DATA BPWD   /                                                                                 &
    & -0.545226E-04 /
!
    DATA ATPWD  /                                                                                 &
    &  0.187967E-01 /
!
    DATA BTPWD  /                                                                                 &
    & -0.567449E-04 /
!
    DATA BDLOWD /                                                                                 &
    &  0.560000E+03 /
!
    DATA BDHIWD /                                                                                 &
    &  0.800000E+03 /
!
    DATA BETINW /                                                                                 &
    &  0.766811E+01 /
!
    END MODULE BDWIDE
