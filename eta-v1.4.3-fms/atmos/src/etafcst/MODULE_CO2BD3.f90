    MODULE CO2BD3
!>--------------------------------------------------------------------------------------------------
!> MODULE CO2BD3
!>
!> ABSTRACT:
!> THE FOLLOWING MODULE CONTAIN PRETABULATED CO2 TRANSMISSION FUNCTIONS, EVALUATED USING THE 
!> METHODS OF FELS AND SCHWARZKOPF (1981) AND SCHWARZKOPF AND FELS (1985), MODULE CO2BD3 CONTAINS 
!> CO2 TRANSMISSION FUNCTIONS AND TEMPERATURE AND PRESSURE DERIVATIVES FOR THE 560-800 CM-1 BAND. 
!> ALSO INCLUDED ARE THE STANDARD TEMPERATURES AND THE WEIGHTING FUNCTION. THESE DATA ARE IN BLOCK
!> DATA BD3:
!>         
!> CO251  - TRANSMISSION FCTNS. FOR T0 (STD. PROFILE) WITH P(SFC)=1013.25 MB
!> CO258  - TRANSMISSION FCTNS. FOR T0 (STD. PROFILE) WITH P(SFC)= ^810 MB
!> CDT51  - FIRST  TEMPERATURE DERIVATIVE OF CO251
!> CDT58  - FIRST  TEMPERATURE DERIVATIVE OF CO258
!> C2D51  - SECOND TEMPERATURE DERIVATIVE OF CO251
!> C2D58  - SECOND TEMPERATURE DERIVATIVE OF CO251
!> CO2M51 - TRANSMISSION FCTNS FOR T0 FOR ADJACENT PRESSURE LEVELS, WITH NO PRESSURE QUADRATURE.
!>          USED FOR NEARBY LAYER COMPUTATIONS. P(SFC)=1013.25 MB
!> CO2M58 - SAME AS CO2M51,WITH P(SFC)= ^810 MB
!> CDTM51 - FIRST  TEMPERATURE DERIVATIVE OF CO2M51
!> CDTM58 - FIRST  TEMPERATURE DERIVATIVE OF CO2M58
!> C2DM51 - SECOND TEMPERATURE DERIVATIVE OF CO2M51
!> C2DM58 - SECOND TEMPERATURE DERIVATIVE OF CO2M58
!> STEMP  - STANDARD TEMPERATURES FOR MODEL PRESSURE LEVEL STRUCTURE WITH P(SFC)=1013.25 MB
!> GTEMP  - WEIGHTING FUNCTION FOR MODEL PRESSURE LEVEL STRUCTURE WITH P(SFC)=1013.25 MB.
!> B0     - TEMP. COEFFICIENT USED FOR CO2 TRANS. FCTN. CORRECTION FOR T(K). (SEE REF. 4 AND BD3)
!> B1     - TEMP. COEFFICIENT, USED ALONG WITH B0
!> B2     - TEMP. COEFFICIENT, USED ALONG WITH B0
!> B3     - TEMP. COEFFICIENT, USED ALONG WITH B0
!>
!> USE MODULES: F77KINDS
!>              RDPARM
!>
!> DRIVER     : CO2DTA 
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE RDPARM  , ONLY : L, LP1
! 
    IMPLICIT NONE
!
    SAVE
! 
    REAL   (KIND=R4KIND), DIMENSION(LP1, LP1)                                                   ::&
    & CO251   , CO258   , CDT51   , CDT58   , C2D51   , C2D58
!
    REAL   (KIND=R4KIND), DIMENSION(L)                                                          ::&
    & CO2M51  , CO2M58  , CDTM51  , CDTM58  , C2DM51  , C2DM58
!
    REAL   (KIND=R4KIND), DIMENSION(LP1)                                                        ::&
    & STEMP   , GTEMP
!   
    REAL   (KIND=R4KIND)                                                                        ::&
    & B0      , B1      , B2      , B3 
!--------------------------------------------------------------------------------------------------
! FROM GFDLRD.f90
! B0,B1,B2,B3 ARE COEFFICIENTS USED TO CORRECT FOR THE USE OF 250K IN THE PLANCK FUNCTION USED IN 
! EVALUATING PLANCK-WEIGHTED CO2 TRANSMISSION FUNCTIONS. (SEE REF. 4)
!--------------------------------------------------------------------------------------------------
    DATA B0, B1, B2, B3 /-.51926410E-4, -.18113332E-3, -.10680132E-5, -.67303519E-7 /
!
    END MODULE CO2BD3
