    SUBROUTINE PRODQ2(LMHK, DTQ2, USTAR, GM, GH, EL, Q2)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE PRODQ2
!>
!> SUBROUTINE: PRODQ2 - Q2 PRODUCTION/DISSIPATION
!> PROGRAMMER: ?????
!> ORG: ?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> LEVEL 2.5 Q2 PRODUCTION/DISSIPATION 
!>
!> PROGRAM HISTORY LOG:
!> 94-??-??  ?????      - ORIGINATOR
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> LMHK  - MASS POINT MODEL SURFACE (LOWEST MODEL HEIGHT)
!> DTQ2  - PHYSICS TINE STEP
!> USTAR - SURFACE FRICTION VELOCITY
!> GM    - WIND SHEAR
!> GH    - VERTICAL GRADIENT OS POTENTIAL TEMPERATURE
!>
!> OUTPUT ARGUMENT LIST:
!> EL    -
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> Q2    - TURBULENT KINETIC ENERGY
!>
!> OUTPUT FILES:
!> NONE
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : TURBL
!>
!> CALLS      : -----             
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA
!
    IMPLICIT NONE
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: LM1 = LM - 1
!
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ2  =  0.12
    REAL   (KIND=R4KIND), PARAMETER :: EPSL   =  0.32
    REAL   (KIND=R4KIND), PARAMETER :: EPSTRB =  1.E-24
    REAL   (KIND=R4KIND), PARAMETER :: EPS1   =  1.E-12
    REAL   (KIND=R4KIND), PARAMETER :: EPS2   =  0.
!
    REAL   (KIND=R4KIND), PARAMETER :: G      =  9.8
    REAL   (KIND=R4KIND), PARAMETER :: BETA   =  1. / 270.
    REAL   (KIND=R4KIND), PARAMETER :: BTG    = BETA * G 
    REAL   (KIND=R4KIND), PARAMETER :: PRT    =  1.0 
    REAL   (KIND=R4KIND), PARAMETER :: GAM1   =  0.2222222222222222222
!---------------------------------------------------------------------
!DCR: Nakanishi and Niino, 2009
    REAL   (KIND=R4KIND), PARAMETER :: A1   =  1.18
    REAL   (KIND=R4KIND), PARAMETER :: A2   =  0.665
    REAL   (KIND=R4KIND), PARAMETER :: B1   =  24.0
    REAL   (KIND=R4KIND), PARAMETER :: B2   =  15.0
    REAL   (KIND=R4KIND), PARAMETER :: C1   =  0.137
!-------------------------------------------------------------------------
!J2001    REAL   (KIND=R4KIND), PARAMETER :: A1   =  0.659888514560862645
!J2001    REAL   (KIND=R4KIND), PARAMETER :: A2   =  0.6574209922667784586
!J2001    REAL   (KIND=R4KIND), PARAMETER :: B1   = 11.87799326209552761
!J2001    REAL   (KIND=R4KIND), PARAMETER :: B2   =  7.226971804046074028
!J2001    REAL   (KIND=R4KIND), PARAMETER :: C1   =  0.000830955950095854396
!--------------------------------------------------------- 
! N     (G=9.8,BETA=1./270.,BTG=BETA*G,
! N     PRT=1.0,GAM1=0.2222222222222222222,
! N     A1=0.3310949523016403346,A2=0.8273378704055731278,
! N     B1=5.959709141429526024,B2=3.626088092074591135,
! N     C1=-0.3330651924968952113)
!--------------------------------------------------------- 
! Y     (G=9.8,BETA=1./270.,BTG=BETA*G,
! Y     PRT=0.8,GAM1=0.2222222222222222222,
! Y     A1=0.9222222350809054114,A2=0.7350190142719400952,
! Y     B1=16.60000023145629741,B2=10.10000014082581951,
! Y     C1=0.0805318118080613468)
!--------------------------------------------------------- 
    REAL   (KIND=R4KIND), PARAMETER :: RB1  = 1. / B1    
!------------------------------------------- 
! COEFFICIENTS OF THE TERMS IN THE NUMERATOR
!-------------------------------------------
    REAL   (KIND=R4KIND), PARAMETER :: ANMM  = -3.*A1*A2*(3.*A2+3.*B2*C1+18.*A1*C1-B2)*BTG
    REAL   (KIND=R4KIND), PARAMETER :: ANMH  = -9.*A1*A2*A2*BTG*BTG 
    REAL   (KIND=R4KIND), PARAMETER :: BNMM  =  A1*(1.-3.*C1) 
    REAL   (KIND=R4KIND), PARAMETER :: BNMH  = -A2*BTG 
!---------------------------------------------                                                                                                                                     
! COEFFICIENTS OF THE TERMS IN THE DENOMINATOR
!--------------------------------------------- 
    REAL   (KIND=R4KIND), PARAMETER :: ADNM  = 18.*A1*A1*A2*(B2-3.*A2)*BTG              
    REAL   (KIND=R4KIND), PARAMETER :: ADNH  =  9.*A1*A2*A2*(12.*A1+3.*B2)*BTG*BTG 
    REAL   (KIND=R4KIND), PARAMETER :: BDNM  =  6.*A1*A1            
    REAL   (KIND=R4KIND), PARAMETER :: BDNH  =  3.*A2*(7.*A1+B2)*BTG
!-----------------------------------------                                                                                                                                        
! COEFFICIENTS OF THE EQUILIBRIUM EQUATION 
!-----------------------------------------                                                                                                                                 
    REAL   (KIND=R4KIND), PARAMETER :: AEQM  = 3.*A1*A2*B1*(3.*A2+3.*B2*C1+18.*A1*C1-B2)*BTG      &
    &                                        +18.*A1*A1*A2*(B2-3.*A2)*BTG
!
    REAL   (KIND=R4KIND), PARAMETER :: AEQH  = 9.*A1*A2*A2*B1*BTG*BTG+9.*A1*A2*A2*(12.*A1+3.*B2)  &
    &                                        *BTG*BTG
!
    REAL   (KIND=R4KIND), PARAMETER :: BEQM  =-A1*B1*(1.-3.*C1)+6.*A1*A1
    REAL   (KIND=R4KIND), PARAMETER :: BEQH  = A2*B1*BTG+3.*A2*(7.*A1+B2)*BTG
!-------------------------- 
! FORBIDDEN TURBULENCE AREA
!--------------------------                                                                                                                                       
    REAL   (KIND=R4KIND), PARAMETER :: REQU  = -AEQH/AEQM*1.02
    REAL   (KIND=R4KIND), PARAMETER :: EPSGH =1.E-9
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                   , INTENT(INOUT)       ::&
    & Q2
!
    REAL   (KIND=R4KIND), DIMENSION(LM1)                                  , INTENT(IN)          ::&
    & GM      , GH
!
    REAL   (KIND=R4KIND), DIMENSION(LM1)                                  , INTENT(INOUT)       ::&
    & EL
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LMHK
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & DTQ2    , USTAR
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & LMHM    , L        
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & GML     , GHL     , AEQU    , BEQU    , EQOL2   , ANUM    , BNUM    , ADEN    , BDEN    ,   &
    & CDEN    , ARHS    , BRHS    , CRHS    , DLOQ1   , ELOQ21  , ELOQ11  , ELOQ31  , ELOQ41  ,   &
    & ELOQ51  , RDEN1   , RHSP1   , ELOQ12  , ELOQ22  , ELOQ32  , ELOQ42  , ELOQ52  , RDEN2   ,   &
    & RHS2    , RHSP2   , RHST2   , ELOQ13  , ELOQN 
!
    LMHM = LMHK - 1
!
    DO 150 L=1,LMHM
        GML = GM(L)
        GHL = GH(L)
!-----------------------------------------
! COEFFICIENTS OF THE EQUILIBRIUM EQUATION
!-----------------------------------------
        AEQU =(AEQM * GML + AEQH * GHL) * GHL
        BEQU = BEQM * GML + BEQH * GHL
!----------------------------- 
! EQUILIBRIUM SOLUTION FOR L/Q
!-----------------------------  
        EQOL2 = -0.5 * BEQU + SQRT(BEQU * BEQU * 0.25 - AEQU)
!---------------------------------- 
! IS THERE PRODUCTION/DISSIPATION ?
!----------------------------------
        IF ((GML+GHL*GHL <= EPSTRB) .OR. (GHL >= EPSGH .AND. GML/GHL <= REQU)                     &
    &                               .OR. (EQOL2 <= EPS2))                      THEN
!-------------- 
! NO TURBULENCE
!--------------
            Q2(L) = EPSQ2
            EL(L) = EPSL
!--------------------------------
! END OF THE NO TURBULENCE BRANCH
!--------------------------------
        ELSE
!-------------------------------------------
! COEFFICIENTS OF THE TERMS IN THE NUMERATOR
!-------------------------------------------
            ANUM = (ANMM * GML + ANMH * GHL) * GHL
            BNUM =  BNMM * GML + BNMH * GHL
!---------------------------------------------
! COEFFICIENTS OF THE TERMS IN THE DENOMINATOR
!---------------------------------------------
            ADEN = (ADNM * GML + ADNH * GHL) * GHL
            BDEN =  BDNM * GML + BDNH * GHL
            CDEN =  1.
!----------------------------------------------------
! COEFFICIENTS OF THE NUMERATOR OF THE LINEARIZED EQ.
!----------------------------------------------------
            ARHS = -(ANUM * BDEN - BNUM * ADEN) * 2.
            BRHS = - ANUM * 4.
            CRHS = - BNUM * 2.
!--------------------- 
! INITIAL VALUE OF L/Q
!--------------------- 
            DLOQ1 = EL(L) / SQRT(Q2(L))
!------------------------------- 
! FIRST ITERATION FOR L/Q, RHS=0
!-------------------------------
            ELOQ21 = 1. / EQOL2
            ELOQ11 = SQRT(ELOQ21)
            ELOQ31 = ELOQ21 * ELOQ11
            ELOQ41 = ELOQ21 * ELOQ21
            ELOQ51 = ELOQ21 * ELOQ31
!--------------- 
! 1./DENOMINATOR
!--------------- 
            RDEN1 = 1. / (ADEN * ELOQ41 + BDEN * ELOQ21 + CDEN)
!--------------  
! D(RHS)/D(L/Q)
!-------------- 
            RHSP1 = (ARHS * ELOQ51 + BRHS * ELOQ31 + CRHS * ELOQ11) * RDEN1 * RDEN1
!--------------------- 
! FIRST-GUESS SOLUTION 
!---------------------
            ELOQ12 = ELOQ11 + (DLOQ1 - ELOQ11) * EXP(RHSP1 * DTQ2)
!
            ELOQ12 = AMAX1(ELOQ12, EPS1)
!-------------------------
! SECOND ITERATION FOR L/Q 
!-------------------------
            ELOQ22 = ELOQ12 * ELOQ12
            ELOQ32 = ELOQ22 * ELOQ12
            ELOQ42 = ELOQ22 * ELOQ22
            ELOQ52 = ELOQ22 * ELOQ32
!--------------- 
! 1./DENOMINATOR 
!---------------
            RDEN2 = 1. / (ADEN * ELOQ42 + BDEN * ELOQ22 + CDEN)
!
            RHS2  = -(ANUM * ELOQ42 + BNUM * ELOQ22) * RDEN2 + RB1
            RHSP2 =  (ARHS * ELOQ52 + BRHS * ELOQ32  + CRHS  * ELOQ12) * RDEN2 * RDEN2
            RHST2 =   RHS2 / RHSP2
!-------------------
! CORRECTED SOLUTION 
!-------------------
            ELOQ13 = ELOQ12 - RHST2 + (RHST2 + DLOQ1 - ELOQ12) * EXP(RHSP2 * DTQ2)
!
            ELOQ13 = AMAX1(ELOQ13, EPS1)
!---------------------------------------
! TWO ITERATIONS IS ENOUGH IN MOST CASES
!---------------------------------------
            ELOQN = ELOQ13
!
            IF (ELOQN > EPS1) THEN
                Q2(L) = EL(L) * EL(L) / (ELOQN * ELOQN)
                Q2(L) = AMAX1(Q2(L), EPSQ2)
            ELSE
                Q2(L) = EPSQ2
            END IF
!------------------------ 
! END OF TURBULENT BRANCH
!------------------------
        END IF
!-----------------------------------
! END OF PRODUCTION/DISSIPATION LOOP
!-----------------------------------
150 END DO
!-------------------------------- 
! LOWER BOUNDARY CONDITION FOR Q2
!--------------------------------
    Q2(LMHK) = AMAX1(B1 ** (2. / 3.) * USTAR * USTAR, EPSQ2)
!
    RETURN
!
    END SUBROUTINE PRODQ2
