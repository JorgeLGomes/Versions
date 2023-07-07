    MODULE CO2BD2
!>--------------------------------------------------------------------------------------------------
!> MODULE CO2BD2
!>
!> ABSTRACT:
!> MODULE CO2BD2 CONTAINS CO2 TRANSMISSION FUNCTIONS AND TEMPERATURE AND PRESSURE DERIVATIVES FOR 
!> THE 560-670 CM-1 PART OF THE 15 UM CO2 BAND. THESE DATA ARE IN BLOCK DATA BD2.
!>
!> CO231 - TRANSMISSION FCTNS. FOR T0 (STD. PROFILE) WITH P(SFC)=1013.25 MB
!> CO238 - TRANSMISSION FCTNS. FOR T0 (STD. PROFILE) WITH P(SFC)= ^810   MB
!> CDT31 - FIRST  TEMPERATURE DERIVATIVE OF CO231
!> CDT38 - FIRST  TEMPERATURE DERIVATIVE OF CO238
!> C2D31 - SECOND TEMPERATURE DERIVATIVE OF CO231
!> C2D38 - SECOND TEMPERATURE DERIVATIVE OF CO231
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
    & CO231   , CO238   ,                                                                         &
    & CDT31   , CDT38   ,                                                                         &
    & C2D31   , C2D38
!
    END MODULE CO2BD2
