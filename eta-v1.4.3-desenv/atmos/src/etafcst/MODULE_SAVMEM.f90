    MODULE SAVMEM
!>--------------------------------------------------------------------------------------------------
!> MODULE SAVMEM
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : GRADFS
!>              O3INT
!>              RADFS
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : LM
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND), PARAMETER :: LNGTH = 37 * LM
!-----------------------------   
! WINTER, SPRING, SUMMER, FALL 
!-----------------------------   
    REAL   (KIND=R4KIND), DIMENSION(37,LM)                                                      ::&
    & DDUO3N  , DDO3N2  , DDO3N3  , DDO3N4
!
    REAL   (KIND=R4KIND), DIMENSION(LNGTH)                                                      ::&
    & RAD1    , RAD2    , RAD3    , RAD4
!
    EQUIVALENCE (RAD1(1), DDUO3N(1,1)), (RAD2(1), DDO3N2(1,1))
    EQUIVALENCE (RAD3(1), DDO3N3(1,1)), (RAD4(1), DDO3N4(1,1))
!
    END MODULE SAVMEM
