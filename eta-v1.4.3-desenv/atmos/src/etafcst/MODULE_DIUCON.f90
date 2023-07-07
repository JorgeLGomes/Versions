    MODULE DIUCON
!>--------------------------------------------------------------------------------------------------
!> MODULE DIUCON
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : GFDLRD
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
! 
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & SEASON  , IXXXX   , JDNMC
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & FCSTDA  , FJDNMC  , TSLAG   , RLAG    , TIMIN   , TPI     , HPI     , YEAR    , DAY     ,   &
    & DHR
!
    INTEGER(KIND=I4KIND), DIMENSION(5)                                                          ::&
    & JTIME
!
    REAL   (KIND=R4KIND), DIMENSION(12)                                                         ::&
    & DAZ
!-------------------------------------------------------------------- 
! FROM GFDLRD
! UNUSED DATA CLEANED OUT - NOV 86 AND MAR 89 K.A. CAMPANA 
!
! FOR SEASONAL VARIATION
! SEASON=1,2,3,4 - FOR WINTER,SPRING,SUMMER,FALL ONLY (NOT ACTIVE)
! SEASON=5       - SEASONAL VARIATION(I.E.INTERPOLATE TO DAY OF FCST)
!-------------------------------------------------------------------- 
    DATA SEASON /5/
    DATA TSLAG  /45.25/
    DATA RLAG   /14.8125/
    DATA DAY    /86400./
    DATA YEAR   /365.25/
    DATA TPI    /6.283185308/
    DATA HPI    /1.570796327/
    DATA DHR    /2./
    DATA JTIME  /0, 1, 0, 0, 0/
    DATA DAZ    /0., 31., 59., 90., 120., 151., 181., 212., 243., 273., 304., 334./
!
    END MODULE DIUCON                                                                          
