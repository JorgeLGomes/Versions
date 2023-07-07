    SUBROUTINE MIXLEN(LMHK, LPBL, HPBL, U, V, T, Q, Q2, APE, Z, GM, GH, EL)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE MIXLEN
!>
!> SUBPROGRAM: MIXLEN - ?????
!> PROGRAMMER: ?????
!> ORG: ?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> LEVEL 2.5 MIXING LENGTH 
!>
!> JULY 1997: MODIFIED TO RESTORE AVERAGING OF LAYER VALUES OF L, ELL, A LA MESINGER 1993, RES. ACT.
!>            ATMOS. OCEAN. MOD., NO. 18, 4.36 - 4.38;
!>
!> A PROBLEM REMOVED WHICH MAY HAVE LED TO THE ÒABOVE PBLÓ SCHEME FOR PRELIMINARY EL TO BE USED 
!> INADVERTENTLY AT TIMES AND PLACES WHERE THE PBL SCHEME WAS SUPPOSED TO HAVE BEEN USED.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  JANJIC     - ORIGINATOR
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> LMVK  - MASS POINT MODEL SURFACE (LOWEST MODEL HEIGHT)
!> U     - (TRANSPOSED)
!> V     - (TRANSPOSED)
!> T     - (TRANSPOSED)
!> Q     - (TRANSPOSED)
!> Q2    - TURBULENT KINETIC ENERGY (TRANSPOSED)
!> APE   - EXNER FUNCTION
!>
!> OUTPUT ARGUMENT LIST:
!> HPBL  - HEIGHT OF THE PBL
!> GM    - WIND SHEAR
!> GH    - VERTICAL GRADIENT OS POTENTIAL TEMPERATURE
!> EL    -
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> LPBL  - THE EL VALUE OF L OF THE INTERFACE
!>
!> USE MODULES: F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : TURBL
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: LP1    = LM + 1
    INTEGER(KIND=I4KIND), PARAMETER :: LM1    = LM - 1
!
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ2  =    0.12
    REAL   (KIND=R4KIND), PARAMETER :: EPSL   =    0.32
    REAL   (KIND=R4KIND), PARAMETER :: EPSRU  =    1.E-12
    REAL   (KIND=R4KIND), PARAMETER :: EPSRS  =    1.E-12
    REAL   (KIND=R4KIND), PARAMETER :: FH     =    1.01
    REAL   (KIND=R4KIND), PARAMETER :: ALPH   =     .20 
    REAL   (KIND=R4KIND), PARAMETER :: VKRM   =     .40 
    REAL   (KIND=R4KIND), PARAMETER :: ELFC   =    0.23   
    REAL   (KIND=R4KIND), PARAMETER :: EL0MAX = 1000.
    REAL   (KIND=R4KIND), PARAMETER :: G      =    9.8
    REAL   (KIND=R4KIND), PARAMETER :: BETA   =    1.    / 270. 
    REAL   (KIND=R4KIND), PARAMETER :: BTG    = BETA     * G
    REAL   (KIND=R4KIND), PARAMETER :: PRT    =    1.0
    REAL   (KIND=R4KIND), PARAMETER :: GAM1   =    0.2222222222222222222
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
!----------------------------------------------------------
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LMHK  
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(INOUT)       ::& 
    & LPBL 
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)         ::& 
    & HPBL 
!---------------------------------------------------------- 
! N     &(G=9.8,BETA=1./270.,BTG=BETA*G
! N     &,PRT=1.0,GAM1=0.2222222222222222222
! N     &,A1=0.3310949523016403346,A2=0.8273378704055731278
! N     &,B1=5.959709141429526024,B2=3.626088092074591135
! N     &,C1=-0.3330651924968952113)
!---------------------------------------------------------- 
! Y     &(G=9.8,BETA=1./270.,BTG=BETA*G
! Y     &,PRT=0.8,GAM1=0.2222222222222222222
! Y     &,A1=0.9222222350809054114,A2=0.7350190142719400952
! Y     &,B1=16.60000023145629741,B2=10.10000014082581951
! Y     &,C1=0.0805318118080613468)
!---------------------------------------------------------- 
!
!---------------------------------------------  
! COEFFICIENTS OF THE TERMS IN THE DENOMINATOR
!---------------------------------------------  
    REAL   (KIND=R4KIND), PARAMETER                                                             ::&
    &          ADNM =18.*A1*A1*A2*(B2-3.*A2)*BTG        ,                                         &
    &          ADNH = 9.*A1*A2*A2*(12.*A1+3.*B2)*BTG*BTG,                                         &
    &          BDNM = 6.*A1*A1                          ,                                         &
    &          BDNH = 3.*A2*(7.*A1+B2)*BTG              ,                                         &
!---------------------------------------------------  
! FREE TERM IN THE EQUILIBRIUM EQUATION FOR (L/Q)**2 
!--------------------------------------------------- 
    &          AEQM = 3.*A1*A2*B1*(3.*A2+3.*B2*C1+18.*A1*C1-B2)*BTG+18.*A1*A1*A2*(B2-3.*A2)*BTG,  &
    &          AEQH = 9.*A1*A2*A2*B1*BTG*BTG+9.*A1*A2*A2*(12.*A1+3.*B2)*BTG*BTG                ,  &
!--------------------------  
! FORBIDDEN TURBULENCE AREA 
!--------------------------
    &          REQU = -AEQH/AEQM*1.02,                                                            &
    &         EPSGH = 1.E-9          ,                                                            &
    &         EPSGM = REQU*EPSGH     ,                                                            &
!------------------------------------------------------
! NEAR ISOTROPY FOR SHEAR TURBULENCE, WW/Q2 LOWER LIMIT 
!------------------------------------------------------
    &        UBRYL = (18.*REQU*A1*A1*A2*B2*C1*BTG+9.*A1*A2*A2*B2*BTG*BTG)/(REQU*ADNM+ADNH),       &
    &         UBRY = (1.+EPSRS)*UBRYL                                                     ,       &
    &        UBRY3 = 3.*UBRY                                                              ,       &
    &         AUBM = 54.*A1*A1*A2*B2*C1*BTG -ADNM*UBRY3                                   ,       &
    &         AUBH = 27.*A1*A2*A2*B2*BTG*BTG-ADNH*UBRY3                                   ,       &
    &         BUBM = 18.*A1*A1*C1           -BDNM*UBRY3                                   ,       &
    &         BUBH = (9.*A1*A2+3.*A2*B2)*BTG-BDNH*UBRY3                                   ,       &
    &         CUBR = 1.-UBRY3                                                             ,       &
    &        RCUBR = 1./CUBR
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                   , INTENT(IN)          ::&    
    & U       , V       , T       , Q       , Q2      , APE 
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                                         ::& 
    & Q1      , ELL     , THV
!
    REAL   (KIND=R4KIND), DIMENSION(LM1)                                  , INTENT(INOUT)       ::& 
    & GM      , GH      , EL      
!
    REAL   (KIND=R4KIND), DIMENSION(LM1)                                                        ::& 
    ELM
!
    REAL   (KIND=R4KIND), DIMENSION(LP1)                                  , INTENT(IN)          ::& 
    & Z
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & LMHM    , LMHP    , L       , IVI     , LPBLP   , NLEV_EPSQ2
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & SZQ     , SQ      , QDZL    , RDZ     , GML     , GHL     , EL0     , AUBR    , BUBR    ,   &
    & QOL2ST  , ELOQ2X  , ADEN    , BDEN    , QOL2UN  , VKRMZ   , SM      , DQ2
!
    LMHM = LMHK - 1
    LMHP = LMHK + 1
!--------------------------- 
! FIND THE HEIGHT OF THE PBL
!---------------------------
!    LPBL = LMHK
!
!    DO 100 IVI=1,LMHK
!        L = LMHK - IVI
!        IF (Q2(L) <= EPSQ2 * FH) THEN
!           LPBL = L
!            GOTO 110
!        END IF
!100 END DO
!
!teste JOrge - 28-01-21 - T4-NEW
          LPBL = LMHK
          NLEV_EPSQ2=0
            DO 100 L=LMHK, LMHK/2, -1   !Daniela
!             DQ2= Q2(L-1)- Q2(L)                      ! Vert Dif Q2  !
	  IF ((Q2(L)<EPSQ2*FH).and.(NLEV_EPSQ2<2)) THEN
              NLEV_EPSQ2=NLEV_EPSQ2+1
              CYCLE   ! raise L either in decreasing or in increasing tke with height
	     END IF   
              IF (Q2(L)>=EPSQ2*FH) THEN	
	     LPBL=L
                 NLEV_EPSQ2=0
                CYCLE
              ELSEIF (NLEV_EPSQ2>=2) THEN
                EXIT   ! raise L either in decreasing or in increasing tke with height
             END IF
        100 END DO
           LPBL=LPBL-1
!
!CHOU_GSM    LPBL = 1
!---------------------- 
! THE HEIGHT OF THE PBL 
!----------------------     
!110 HPBL = Z(LPBL) - Z(LMHP)
 HPBL = Z(LPBL) - Z(LMHP)
!
!
    DO 120 L=1,LMHK
        Q1(L) = 0.
        THV(L) = (0.608 * Q(L) + 1.) * T(L) * APE(L)
120 END DO
!
    DO 130 L=LPBL,LMHK
        Q1(L) = SQRT(Q2(L))
130 END DO
!
    SZQ = 0.
    SQ  = 0.
!
    DO 140 L=1,LMHM
        QDZL = (Q1(L) + Q1(L+1)) * (Z(L+1) - Z(L+2))
!
        SZQ  = (Z(L+1) + Z(L+2) - Z(LMHP) - Z(LMHP)) * QDZL + SZQ
!
        SQ = QDZL + SQ
!
        RDZ = 2. / (Z(L) - Z(L+2))
!
        GML = ((U(L) - U(L+1)) ** 2 + (V(L) - V(L+1)) ** 2) * RDZ * RDZ
!
        GM(L) = AMAX1(GML,EPSGM)
!
        GHL = (THV(L) - THV(L+1)) * RDZ
!
        IF (ABS(GHL) <= EPSGH) GHL = EPSGH
!
        GH(L) = GHL
140 END DO
!------------------------------------------------- 
! COMPUTATION OF ASYMPTOTIC L IN BLACKADAR FORMULA 
!------------------------------------------------- 
    EL0 = AMIN1(ALPH * SZQ * 0.5 / SQ, EL0MAX)
!
    DO 150 L=1,LMHM
        GML = GM(L)
        GHL = GH(L)
!
        IF (GHL >= EPSGH )THEN
            IF (GML/GHL <= REQU) THEN
                ELM(L) = EPSL
            ELSE
                AUBR = (AUBM * GML + AUBH * GHL) * GHL
                BUBR =  BUBM * GML + BUBH * GHL
!
                QOL2ST = (-0.5 * BUBR + SQRT(BUBR * BUBR * 0.25 - AUBR * CUBR)) * RCUBR
                QOL2ST = AMAX1(QOL2ST, 1.E-8)
!
                ELOQ2X = 1. / QOL2ST
                ELM(L)  = AMAX1(SQRT(ELOQ2X * Q2(L)), EPSL)
            END IF
        ELSE
            ADEN = (ADNM * GML + ADNH * GHL) * GHL
            BDEN =  BDNM * GML + BDNH * GHL
!
            QOL2UN = -0.5 * BDEN + SQRT(BDEN * BDEN * 0.25 - ADEN)
            ELOQ2X = 1. / ((1. + EPSRU) * QOL2UN)
            ELM(L) = AMAX1(SQRT(ELOQ2X * Q2(L)), EPSL)
        END IF
150 END DO
!-----------------------------------------------------------------------
! LPBLM=LPBL-1
! DO 160 L=1,LPBLM
!    EL(L)=AMIN1((Z(L)-Z(L+2))*ELFC,ELM(L))
! 60 CONTINUE
! DO 165 L=LPBL,LMHM
!    VKRMZ=(Z(L+1)-Z(LMHP))*VKRM
! 65 EL(L)=AMIN1(VKRMZ/(VKRMZ/EL0+1.),ELM(L))
!-----------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
!     NOTE: LPBL IS THE EL VALUE OF L OF THE INTERFACE STARTING WITH THE LOWEST INTERFACE AND GOING 
! UP FOR THE FIRST TIME HAVING Q2 .LE. A SPECIFIED VALUE
!--------------------------------------------------------------------------------------------------  
    DO 160 L=1,LPBL
        ELL(L) = (Z(L) - Z(L+1)) * ELFC
160 END DO
!
    IF (LPBL < LMHK) THEN
        LPBLP = LPBL + 1
        DO 165 L=LPBLP,LMHK
            VKRMZ  = (0.5 * (Z(L) + Z(L+1)) - Z(LMHP)) * VKRM
            ELL(L) = VKRMZ / (VKRMZ / EL0 + 1.)
    165 END DO
    END IF
!
    DO 260 L=1,LMHM
        EL(L) = AMIN1(0.5 * (ELL(L) + ELL(L+1)), ELM(L))
260 END DO
!
    RETURN
!
    END SUBROUTINE MIXLEN
