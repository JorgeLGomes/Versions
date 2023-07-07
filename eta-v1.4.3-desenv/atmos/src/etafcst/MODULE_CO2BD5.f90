    MODULE CO2BD5
!>--------------------------------------------------------------------------------------------------
!> MODULE CO2BD5
!>
!> ABSTRACT:
!> MODULE CO2BD5 CONTAINS CO2 TRANSMISSION FUNCTIONS FOR THE 2270-2380 PART OF THE 4.3 UM CO2 BAND. 
!> THESE DATA ARE IN BLOCK DATA BD5.
!>
!> CO211    -  TRANSMISSION FCTNS. FOR T0 (STD. PROFILE) WITH P(SFC)=1013.25 MB
!> CO218    -  TRANSMISSION FCTNS. FOR T0 (STD. PROFILE) WITH P(SFC)= ^810 MB
!>
!> USE MODULES: F77KINDS
!>              RDPARM
!>
!> DRIVER     : CO2DTA
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE RDPARM  , ONLY : LP1
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(LP1)                                                        ::&
    & CO211   , CO218
!
    END MODULE CO2BD5
