    SUBROUTINE GRADFS(SIGL, KCCO2, NFILE)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE GRADFS
!>
!> SUBPROGRAM: GRADFS - ?????
!> PROGRAMMER: Q. ZHAO
!> ORG: ?????
!> DATE: 93-11-18
!>
!> ABSTRACT: 
!> ?????
!>
!> PROGRAM HISTORY LOG:
!> 93-11-18  Q. ZHAO    - ORIGINATOR
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NFILE - THE FILE NUMBER FOR O3 DATA
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> SIGL  - MIDLAYER PRESSURES IN PA (LP1=LM+1)
!> KCCO2 - = 0 (NOT USED)
!>
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>              RDFSAV
!>              RDPARM
!>              SAVMEM
!> 
!> DRIVER     : INIT
!>              INITS
!>
!> CALLS      : CONRAD
!>              HCONST
!>              O3INT
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA
    USE RDFSAV
    USE RDPARM
    USE SAVMEM
!
    IMPLICIT NONE
!
#include "sp.h"
!--------------------------------------------------------------------------------------------------
! SEASONAL CLIMATOLOGIES OF O3 (OBTAINED FROM A PREVIOUSLY RUN CODE WHICH INTERPOLATES O3 TO USER 
! VERTICAL COORDINATE).
! DEFINED AS 5 DEG LAT MEANS N.P.-> S.P.
!-------------------------------------------------------------------------------------------------- 
!
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    REAL   (KIND=R4KIND), DIMENSION(LP1)                                  , INTENT(INOUT)       ::&
    & SIGL
!
    REAL   (KIND=R4KIND), DIMENSION(5)                                                          ::&
    & XAO3SW  , XAH2SW  , XBSW
!
    REAL   (KIND=R4KIND)                                                                        ::& 
    & PI      , AVG     , A1      , B1      , B2    
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , LV       
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & NFILE
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(INOUT)       ::&
    & KCCO2 

    DATA                                                                                          &
    & XAO3SW / 0., .690, .480, .210, 0./ ,                                                        &
    & XAH2SW / 0., .690, .480, .210, 0./ ,                                                        &
    & XBSW   / 0., .035, .020, .005, 0./
!--------------------------------------------- 
! ONE TIME COMPUTATION OF NECESSARY QUANTITIES 
!--------------------------------------------- 
!
!---------------------------------------
! INITIALIZE ARRAYS,GET CONSTANTS,ETC...
!--------------------------------------- 
    PI     =   3.1415927
    Q19001 =  19.001
    HP98   =   0.98
    H3M6   =   3.0E-6
    HP537  =   0.537
    H74E1  =  74.0
    H15E1  =  15.0
    Q14330 =   1.43306E-6
    HP2    =   0.2
    TWENTY =  20.0
    HNINE  =   9.0
    DEGRAD = 180.0 / PI
    HSIGMA =   5.673E-5
    DAYSEC =   1.1574E-5
!--------------------------------------------------------------------------------------------------
! ATMOSPERIC CARBON DIOXIDE CONCENTRATION IS NOW READ BY CONRAD, BUT IT DEFAULTS TO 330 PPM 
! FOR BACKWARD COMPATIBILITY.
!--------------------------------------------------------------------------------------------------
    RCO2   =   3.3E-4
!
    CALL HCONST
!--------------------------------------------------------- 
! INTERPOLATE CLIMO O3 TO THE CURRENT VERTICAL COORDINATE.
! NEED LAYER SIGMA, GET FROM PSFC AND LAYER P FOR I=1...
!---------------------------------------------------------
    DO 3 I=1,5
        CAO3SW(I) = XAO3SW(I)
        CAH2SW(I) = XAH2SW(I)
          CBSW(I) =   XBSW(I)
  3 END DO
!----------------------------------------------- 
! CONVERT SIGL FROM PA TO MB TO BE USED IN O3INT
!----------------------------------------------- 
    DO 100 LV=1,LP1
        SIGL(LV) = 0.01 * SIGL(LV)
100 END DO
!
    CALL O3INT(SIGL)
    CALL CONRAD(NFILE)
!-------------------------------------------------------------------------------------------------------- 
! AVERAGE CLIMATOLOGICAL VALUS OF O3 FROM 5 DEG LAT MEANS, SO THAT TIME AND SPACE INTERPOLATION WILL WORK 
! (DONE ELSEWHERE IN RADFS)
!-------------------------------------------------------------------------------------------------------- 
    DO 5 I=1,LNGTH
        AVG = .25E0 * ( RAD1(I) + RAD2(I)  +  RAD3(I) + RAD4(I))
        A1  = .5E0  * ( RAD2(I) - RAD4(I))
        B1  = .5E0  * ( RAD1(I) - RAD3(I))
        B2  = .25E0 * ((RAD1(I) + RAD3(I)) - (RAD2(I) + RAD4(I)))
!
        RAD1(I) = AVG
        RAD2(I) = A1
        RAD3(I) = B1
        RAD4(I) = B2
  5 END DO
!
    EMIST =  .6E0
    EMISP =  .3E0
    XLATP = 60.E0
    XLATT = 30.E0
!
    RETURN
!
    END SUBROUTINE GRADFS
