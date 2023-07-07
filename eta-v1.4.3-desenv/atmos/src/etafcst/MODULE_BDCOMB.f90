    MODULE BDCOMB
!>--------------------------------------------------------------------------------------------------
!> MODULE BDCOMB
!>
!> ABSTRACT:
!> MODULE BDCOMB CONTAINS RANDOM BAND PARAMETERS FOR THE LW CALCULATIONS USING COMBINED WIDE 
!> FREQUENCY BANDS BETWEEN 160 AND 1200 CM-1,AS WELL AS THE 2270-2380 BAND FOR SOURCE CALC.
!> BANDS 1-8:  COMBINED WIDE FREQUENCY BANDS FOR 160-560 CM-1
!> BANDS 9-14: FREQUENCY BANDS,AS IN BANDTA (NARROW BANDS) FOR 560-1200 CM-1
!> BAND  15:   FREQUENCY BAND 2270-2380 CM-1,USED FOR SOURCE CALCULATION ONLY THUS NBLY PRESENTLY
!>             EQUALS 15
!>
!> BANDS ARE ARRANGED IN ORDER OF INCREASING WAVENUMBER
!>
!> ACOMB        - RANDOM "A" PARAMETER FOR (NBLY) BANDS
!> BCOMB        - RANDOM "B" PARAMETER FOR (NBLY) BANDS
!> BETACM       - CONTINUUM COEFFICIENTS FOR (NBLY) BANDS
!> APCM , BPCM  - CAPPHI COEFFICIENTS FOR (NBLY) BANDS
!> ATPCM, BTPCM - CAPPSI COEFFICIENTS FOR (NBLY) BANDS
!> BDLOCM       - LOWEST FREQUENCY IN EACH OF (NBLY) FREQ. BANDS
!> BDHICM       - HIGHEST FREQUENCY IN EACH OF (NBLY) FREQ. BANDS
!> AO3CM        - RANDOM "A" PARAMETER FOR OZONE IN (3) OZONE BANDS
!> BO3CM        - RANDOM "B" PARAMETER FOR OZONE IN (3) OZONE BANDS
!> AB15CM       - THE PRODUCT ARNDM*BRNDM FOR THE TWO BANDS REPRESENTING THE 15 UM BAND COMPLEX 
!>                OF CO2
!> BETINC       - CONT.COEFFICIENT FOR A SPECIFIED WIDE FREQ.BAND (800-990 AND 1070-1200 CM-1).
!> IBAND        - INDEX NO OF THE 40 WIDE BANDS USED IN COMBINED WIDE BAND CALCULATIONS. IN OTHER
!>                WORDS, INDEX TELLING WHICH OF THE 40 WIDE BANDS BETWEEN 160-560 CM-1 ARE INCLUDED
!>                IN EACH OF THE FIRST 8 COMBINED WIDE BANDS
!>
!> DATA FOR ACOMB, BCOMB, APCM, BPCM, ATPCM, BTPCM, AO3CM, BO3CM ARE OBTAINED BY USING THE AFGL 
!> 1982 CATALOG. 
!> CONTINUUM COEFFICIENTS ARE FROM ROBERTS 1976. 
!> IBAND INDEX VALUES ARE OBTAINED BY EXPERIMENTATION.
!>
!> USE MODULES: F77KINDS
!>              RDPARM
!>
!> DRIVER     : RNDDTA
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE RDPARM  , ONLY : NBLY
!
    IMPLICIT NONE
!
    SAVE
! 
    INTEGER(KIND=I4KIND), DIMENSION(40)                                                         ::&
    & IBAND
!
    REAL   (KIND=R4KIND), DIMENSION(NBLY)                                                       ::&
    & ACOMB   , BCOMB   , BETACM  , APCM    , BPCM    , ATPCM   , BTPCM   , BDLOCM  , BDHICM
!
    REAL   (KIND=R4KIND), DIMENSION(2)                                                          ::&
    & AB15CM
!
    REAL   (KIND=R4KIND), DIMENSION(3)                                                          ::&
    & AO3CM   , BO3CM
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & BETINC
!------------
! FROM GFDLRD
!------------
    DATA IBAND  /                                                                                 &
    & 2, 1, 2, 2, 1, 2, 1, 3, 2, 2, 3, 2, 2, 4, 2, 4, 2, 3, 3, 2,                                 &
    & 4, 3, 4, 3, 7, 5, 6, 7, 6, 5, 7, 6, 7, 8, 6, 6, 8, 8, 8, 8 /
! 
    DATA ACOMB  /                                                                                 &
    &  0.152070E+05,  0.332194E+04,  0.527177E+03,  0.163124E+03,  0.268808E+03,  0.534591E+02,   &
    &  0.268071E+02,  0.123133E+02,  0.600199E+01,  0.640803E+00,  0.501549E-01,  0.167961E-01,   &
    &  0.178110E-01,  0.170166E+00,  0.537083E-02 /
!
    DATA BCOMB  /                                                                                 &
    &  0.152538E+00,  0.118677E+00,  0.103660E+00,  0.100119E+00,  0.127518E+00,  0.118409E+00,   &
    &  0.904061E-01,  0.642011E-01,  0.629660E-01,  0.643346E-01,  0.717082E-01,  0.629730E-01,   &
    &  0.875182E-01,  0.857907E-01,  0.214005E+00 /
!
    DATA BETACM /                                                                                 &
    &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.188625E+03,  0.144293E+03,   &
    &  0.174098E+03,  0.909366E+02,  0.497489E+02,  0.221212E+02,  0.113124E+02,  0.754174E+01,   &
    &  0.589554E+01,  0.495227E+01,  0.000000E+00 /
!
    DATA APCM   /                                                                                 &
    & -0.671879E-03,  0.654345E-02,  0.143657E-01,  0.923593E-02,  0.117022E-01,  0.159596E-01,   &
    &  0.181600E-01,  0.145013E-01,  0.170062E-01,  0.233303E-01,  0.256735E-01,  0.274745E-01,   &
    &  0.279259E-01,  0.197002E-01,  0.349782E-01 /
!
    DATA BPCM   /                                                                                 &
    & -0.113520E-04, -0.323965E-04, -0.448417E-04, -0.230779E-04, -0.361981E-04, -0.145117E-04,   &
    &  0.198349E-04, -0.486529E-04, -0.550050E-04, -0.684057E-04, -0.447093E-04, -0.778390E-04,   &
    & -0.982953E-04, -0.772497E-04, -0.748263E-04 /
!
    DATA ATPCM  /                                                                                 & 
    & -0.106346E-02,  0.641531E-02,  0.137362E-01,  0.922513E-02,  0.136162E-01,  0.169791E-01,   &
    &  0.206959E-01,  0.166223E-01,  0.171776E-01,  0.229724E-01,  0.275530E-01,  0.302731E-01,   &
    &  0.281662E-01,  0.199525E-01,  0.370962E-01 /
!
    DATA BTPCM  /                                                                                 &
    & -0.735731E-05, -0.294149E-04, -0.505592E-04, -0.280894E-04, -0.492972E-04, -0.341508E-04,   &
    & -0.362947E-04, -0.250487E-04, -0.521369E-04, -0.746260E-04, -0.744124E-04, -0.881905E-04,   &
    & -0.933645E-04, -0.664045E-04, -0.115290E-03 /
!
    END MODULE BDCOMB
