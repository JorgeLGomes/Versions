    SUBROUTINE FST88(HEATRA, GRNFLX, TOPFLX, QH2O  , PRESS  , P     , DELP  , DELP2 , TEMP  ,     &
    &                T     , CLDFAC, NCLDS , KTOP  , KBTM   , CAMT  , CO21  , CO2NBL, CO2SP1,     &
    &                CO2SP2, VAR1  , VAR2  , VAR3  , VAR4   , CNTVAL, TOTO3 , TPHIO3, TOTPHI,     &
    &                TOTVO2, EMX1  , EMX2  , EMPL)
!>--------------------------------------------------------------------------------------------------  
!> SUBROUTINE FST88
!> 
!> SUBPROGRAM: FST88 - IS THE MAIN COMPUTATION MODULE OF THE LONG-WAVE RADIATION CODE. 
!>
!> PROGRAMMER: ?????
!> ORG: ?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> IN IT ALL "EMISSIVITY" CALCULATIONS, INCLUDING CALLS TO TABLE LOOKUP SUBROUTINES. 
!> ALSO, AFTER CALLING SUBROUTINE "SPA88", FINAL COMBINED HEATING RATES AND GROUND FLUX ARE OBTAINED
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????  - ORIGINATOR
!> 18-01-15  LUCCI  - MODERNIZATION OF THE CODE, INCLUDING:
!>                    * F77 TO F90/F95
!>                    * INDENTATION & UNIFORMIZATION CODE
!>                    * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                    * DOCUMENTATION WITH DOXYGEN
!>                    * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> QH2O   -
!> PRESS  -
!> P      -
!> DELP   -
!> DELP2  -
!> TEMP   -
!> T      -
!> CLDFAC - 
!> NCLDS  - 
!> KTOP   -
!> KBTM   -
!> CAMT   -
!> CO21   -
!> CO2NBL -
!> CO2SP1 -
!> CO2SP2 -
!> VAR1   -
!> VAR2   -
!> VAR3   -
!> VAR4   -
!> CNTVAL -
!> TOTO3  -
!> TPHIO3 -
!> TOTPHI -
!> TOTVO2 - 
!> EMX1   -
!> EMX2   -
!> EMPL   - H2O AMOUNT,INPUT FOR E3 CALCULATION IN E3V88 (COMPUTED IN LWR88; STORED IN KDACOM.H)
!>
!> OUTPUT ARGUMENT LIST:
!> HEATRA -
!> GRNFLX -
!> TOPFLX -
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> 
!>
!> USE MODULES: F77KINDS
!>              GLB_TABLE
!>              HCON
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              PHYCON
!>              RNDDTA
!>              TABCOM
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : LWR88
!>
!> CALLS      : E1E290
!>              E290
!>              SPA88
!>              E3V88
!>              E2SPEC
!>--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
! PASSED VARIABLES:
! IN E3V88:
! EMD     =  E3 FUNCTION FOR H2O LINES (0-560,1200-2200 CM-1) COMPUTED IN E3V88
! TPL     =  TEMPERATURE INPUT FOR E3 CALCULATION IN E3V88
!
! IN E1E288:
! E1CTS1  =  E1 FUNCTION FOR THE (I+1)TH LEVEL USING THE TEMPERATURE OF THE ITH DATA LEVEL, 
!            COMPUTED OVER THE FREQUENCY RANGE 0-560, 1200-2200 CM-1. (E1CTS1-E1CTW1) IS USED IN 
!            OBTAINING THE FLUX AT THE TOP IN THE 0-160,1200-2200 CM-1 RANGE (FLX1E1).
! E1CTS2  =  E1 FUNCTION FOR THE ITH LEVEL, USING THE TEMP. OF THE ITH DATA LEVEL,COMPUTED OVER THE
!            FREQUENCY RANGE 0-560, 1200-2200 CM-1. (E1CTS2-E1CTW2) IS ALSO USED IN OBTAINING THE 
!            FLUX AT THE TOP IN THE 0-160, 1200-2200 CM-1 RANGE.
! E1FLX   =  E1 FCTN. FOR THE ITH LEVEL,USING THE TEMPERATURE AT THE TOP OF THE ATMOSPHERE. 
!            COMPUTED OVER THE FREQUENCY RANGE 0-560, 1200-2200 CM-1. USED FOR Q(APPROX) TERM. 
!            (IN COMMON BLOCK TFCOM)
! E1CTW1  =  LIKE E1CTS1, BUT COMPUTED OVER THE 160-560 CM-1 RANGE AND USED FOR 
!            Q(APPROX,CTS) CALCULATION
! E1CTW2  =  LIKE E1CTS2, BUT COMPUTED OVER THE 160-560 CM-1 RANGE AND USED FOR 
!            Q(APPROX,CTS) CALCULATION
! FXO     =  TEMPERATURE INDEX USED FOR E1 FUNCTION AND ALSO USED FOR SOURCE FUNCTION CALC. 
!            IN FST88
! DT      =  TEMP. DIFF.BETWEEN MODEL TEMPS. AND TEMPS. AT TABULAR VALUES OF E1 AND SOURCE FCTNS.
!            USED IN FST88 AND IN E1 FUNCTION CALC.
! FXOE2   =  TEMPERATURE INDEX USED FOR E2 FUNCTION
! DTE2    =  TEMP. DIFF. BETWEEN MODEL TEMP. AND TEMPS. AT TABULAR VALUES OF E2 FUNCTION.
!--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE GLB_TABLE
    USE HCON
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE PHYCON
    USE RNDDTA
    USE TABCOM
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
#include "sp.h"
!---------------------------------------------------------------------------------
! PARAMETER SETTINGS FOR THE LONGWAVE AND SHORTWAVE RADIATION CODE:
! IMAX   =  NO. POINTS ALONG THE LAT. CIRCLE USED IN CALCS.
! L      =  NO. VERTICAL LEVELS (ALSO LAYERS) IN MODEL
! NOTE: THE USER NORMALLY WILL MODIFY ONLY THE IMAX AND L PARAMETERS
! NBLW   =  NO. FREQ. BANDS FOR APPROX COMPUTATIONS. SEE BANDTA FOR DEFINITION
! NBLX   =  NO. FREQ BANDS FOR APPROX CTS COMPUTATIONS
! NBLY   =  NO. FREQ. BANDS FOR EXACT CTS COMPUTATIONS. SEE BDCOMB FOR DEFINITION
! INLTE  =  NO. LEVELS USED FOR NLTE CALCS.
! NNLTE  =  INDEX NO. OF FREQ. BAND IN NLTE CALCS. NB,KO2 ARE SHORTWAVE PARAMETERS;
!           OTHER QUANTITIES ARE DERIVED FROM THE ABOVE PARAMETERS.
!---------------------------------------------------------------------------------
    INTEGER(KIND=I4KIND), PARAMETER :: L      = LM
    INTEGER(KIND=I4KIND), PARAMETER :: NCOL   = IMAX
    INTEGER(KIND=I4KIND), PARAMETER :: NBLX   = 47
    INTEGER(KIND=I4KIND), PARAMETER :: NBLM   = NBLY  - 1
    INTEGER(KIND=I4KIND), PARAMETER :: LP2    = L     + 2
    INTEGER(KIND=I4KIND), PARAMETER :: LP3    = L     + 3
    INTEGER(KIND=I4KIND), PARAMETER :: LM1    = L     - 1
    INTEGER(KIND=I4KIND), PARAMETER :: LM2    = L     - 2
    INTEGER(KIND=I4KIND), PARAMETER :: LM3    = L     - 3
    INTEGER(KIND=I4KIND), PARAMETER :: LL     = 2     * L
    INTEGER(KIND=I4KIND), PARAMETER :: LLP1   = LL    + 1
    INTEGER(KIND=I4KIND), PARAMETER :: LLP2   = LL    + 2
    INTEGER(KIND=I4KIND), PARAMETER :: LLP3   = LL    + 3
    INTEGER(KIND=I4KIND), PARAMETER :: LLM1   = LL    - 1
    INTEGER(KIND=I4KIND), PARAMETER :: LLM2   = LL    - 2
    INTEGER(KIND=I4KIND), PARAMETER :: LLM3   = LL    - 3
    INTEGER(KIND=I4KIND), PARAMETER :: LP1M   = LP1   * LP1
    INTEGER(KIND=I4KIND), PARAMETER :: LP1M1  = LP1M  - 1
    INTEGER(KIND=I4KIND), PARAMETER :: LP121  = LP1   * NBLY
    INTEGER(KIND=I4KIND), PARAMETER :: LL3P   = 3     * L    + 2
    INTEGER(KIND=I4KIND), PARAMETER :: NB     = 12
    INTEGER(KIND=I4KIND), PARAMETER :: INLTE  = 3
    INTEGER(KIND=I4KIND), PARAMETER :: INLTEP = INLTE + 1
    INTEGER(KIND=I4KIND), PARAMETER :: NNLTE  = 56
    INTEGER(KIND=I4KIND), PARAMETER :: LP1I   = IMAX  * LP1
    INTEGER(KIND=I4KIND), PARAMETER :: LLP1I  = IMAX  * LLP1
    INTEGER(KIND=I4KIND), PARAMETER :: LL3PI  = IMAX  * LL3P
    INTEGER(KIND=I4KIND), PARAMETER :: NB1    = NB    - 1
    INTEGER(KIND=I4KIND), PARAMETER :: KO2    = 12
    INTEGER(KIND=I4KIND), PARAMETER :: KO21   = KO2   + 1
    INTEGER(KIND=I4KIND), PARAMETER :: KO2M   = KO2   - 1
!--------------------------------------------------------------------------------------------------
! COMMON BLOCK TABCOM CONTAINS QUANTITIES PRECOMPUTED IN SUBROUTINE TABLE FOR USE IN THE LONGWAVE
! RADIATION PROGRAM:
! EM1     =  E1 FUNCTION, EVALUATED OVER THE 0-560 AND 1200-2200 CM-1 INTERVALS
! EM1WDE  =  E1 FUNCTION, EVALUATED OVER THE 160-560 CM-1 INTERVAL
! TABLE1  =  E2 FUNCTION, EVALUATED OVER THE 0-560 AND 1200-2200 CM-1 INTERVALS
! TABLE2  =  TEMPERATURE DERIVATIVE OF TABLE1
! TABLE3  =  MASS DERIVATIVE OF TABLE1
! EM3     =  E3 FUNCTION, EVALUATED OVER THE 0-560 AND 1200-2200 CM-1 INTERVALS
! SOURCE  =  PLANCK FUNCTION, EVALUATED AT SPECIFIED TEMPS. FOR BANDS USED IN CTS CALCULATIONS
! DSRCE   =  TEMPERATURE DERIVATIVE OF SOURCE
! IND     =  INDEX, WITH VALUE IND(I)=I. USED IN FST88
! INDX2   =  INDEX VALUES USED IN OBTAINING "LOWER TRIANGLE" ELEMENTS OF AVEPHI,ETC., IN FST88
! KMAXV   =  INDEX VALUES USED IN OBTAINING "UPPER TRIANGLE" ELEMENTS OF AVEPHI,ETC., IN FST88
! KMAXVM  =  KMAXV(L),USED FOR DO LOOP INDICES
!--------------------------------------------------------------------------------------------------
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2)                          , INTENT(IN)          ::&
    & NCLDS 
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & ITOP    , IBOT    , INDTC
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                          , INTENT(IN)          ::&
    & EMX1    , EMX2 
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & GXCTS   , FLX1E1  , DELPTC  , PTOP    , PBOT    , FTOP    , FBOT    
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                          , INTENT(INOUT)       ::&
    & GRNFLX  , TOPFLX 
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, 2)                                             ::&
    & EMSPEC
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, L)                       , INTENT(INOUT)       ::&
    & HEATRA
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, L)                       , INTENT(IN)          ::&
    & CO2NBL  , VAR1    , VAR2    , VAR3    , VAR4     
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, L)                                             ::&
    & TO3SPC  , CTS     , EXCTS   , CTSO3       
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, L)                       , INTENT(IN)          ::&
    & DELP    , DELP2   
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, LP1)                     , INTENT(IN)          ::&
    & KTOP    , KBTM     
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                     , INTENT(IN)          ::&
    & QH2O    , PRESS   , P       , TEMP    , T       ,  CAMT   , CO2SP1  , CO2SP2  ,             &
    & CNTVAL  , TOTO3   , TPHIO3  , TOTPHI  , TOTVO2      
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & AVEPHI  , EMISS   , EMISSB  , AVPHO3  , TO31D   , CONT1D  , AVMO3   ,  OVER1D  ,            &
    & E1FLX   , CO2SP   , TO3SP   , OSS     , CSS     , SS1     , SS2     , TC       ,            &
    & DTC     , CSOUR   , AVVO2   , HEATEM
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LLP1)                    , INTENT(IN)          ::&
    & EMPL    
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LLP1)                                          ::&
    & C       , C2
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1, LP1)                , INTENT(IN)          ::&
    & CLDFAC   
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1, LP1)                , INTENT(INOUT)       ::&
    & CO21
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1, NBLY)                                     ::&
    & SORC
!------------------------------------------------------ 
! DIMENSION OF VARIABLES EQUIVALENCED TO THOSE IN VTEMP
!------------------------------------------------------
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & IXO
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & VTMP3   , DSORC   , FAC1    , DELPR1  , DELPR2  , EMISDG  , CONTDG  , TO3DG   ,             &
    & FLXNET  , VSUM1   , FLXTHK  , Z1             
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LLP1)                                          ::&
    & ALP     , CSUB    , CSUB2
!------------------------------------------------------------------------------------ 
! DIMENSION OF VARIABLES PASSED TO OTHER SUBROUTINES (AND NOT FOUND IN COMMON BLOCKS)
!------------------------------------------------------------------------------------ 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, L)                                             ::&
    & E1CTS2  , E1CTW2
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & E1CTS1  , E1CTW1
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LLP1)                                          ::&
    & EMD     , TPL
!-------------------------------------------------------------------------------------------------- 
! IT IS POSSIBLE TO EQUIVALENCE EMD,TPL TO THE ABOVE VARIABLES, AS THEY GET CALLED AT DIFFERENT 
! TIMES
!-------------------------------------------------------------------------------------------------- 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2,LP1)                                            ::&
    & FXO     , DT      , FXOE2   , DTE2                                
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2,2)                                              ::&
    & FXOSP   , DTSP
!----------------------------- 
! DIMENSION OF LOCAL VARIABLES
!----------------------------- 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2,L)                                              ::&
    & RLOG
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2,LP1)                                            ::&
    & FLX     , TOTEVV  , CNTTAU
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , K       , KP      , KLEN    , KK      , ICNT    , KCLDS   , KMIN    , KMAX    ,   &
    & J1      , J3   
!
    EQUIVALENCE (ALP,C,CSUB), (CSUB2,C2)
    EQUIVALENCE (FAC1, DSORC, OVER1D, DELPR2, FLXNET)
    EQUIVALENCE (DELPR1, HEATEM)
    EQUIVALENCE (IXO, AVVO2, FLXTHK, TO3DG)
    EQUIVALENCE (Z1, AVMO3, CONTDG)
    EQUIVALENCE (EMISDG, VSUM1, AVPHO3)
    EQUIVALENCE (EMD(IDIM1,1), E1CTS1(IDIM1,1)), (EMD(IDIM1,LP2), E1CTS2(IDIM1,1))
    EQUIVALENCE (TPL(IDIM1,1), E1CTW1(IDIM1,1)), (TPL(IDIM1,LP2), E1CTW2(IDIM1,1))
!------------------------------------------------------------------------------- 
! FIRST SECTION IS TABLE LOOKUP FOR SOURCE FUNCTION AND DERIVATIVE (B AND DB/DT)
! ALSO,THE NLTE CO2 SOURCE FUNCTION IS OBTAINED
!
! IN CALCS. BELOW, DECREMENTING THE INDEX BY 9 
! ACCOUNTS FOR THE TABLES BEGINNING AT T=100K.
! AT T=100K.
!------------------------------------------------------------------------------- 
    DO 101 K=1,LP1
        DO 101 I=MYIS,MYIE
!-----------------------------  
! TEMP. INDICES FOR E1, SOURCE
!-----------------------------  
            VTMP3(I,K) = AINT(TEMP(I,K) * HP1)
              FXO(I,K) =     VTMP3(I,K) - 9.
               DT(I,K) =      TEMP(I,K) - TEN * VTMP3(I,K)
!--------------------------------------------  
! INTEGER INDEX FOR SOURCE (USED IMMEDIATELY)
!-------------------------------------------- 
            IXO(I,K) = FXO(I,K)
101 END DO
!
    DO 103 k=1,L
        DO 103 I=MYIS,MYIE
!----------------------------------------------------------------  
! TEMP. INDICES FOR E2 (KP=1 LAYER NOT USED IN FLUX CALCULATIONS)
!---------------------------------------------------------------- 
            VTMP3(I,K) = AINT(T(I,K+1) * HP1)
            FXOE2(I,K) =  VTMP3(I,K)   - 9.
             DTE2(I,K) =      T(I,K+1) - TEN * VTMP3(I,K)
103 END DO
!----------------------------------------------------------  
! SPECIAL CASE TO HANDLE KP=LP1 LAYER AND SPECIAL E2 CALCS.
!---------------------------------------------------------- 
    DO 105 I=MYIS,MYIE
        FXOE2(I,LP1) =   FXO(I,L)
         DTE2(I,LP1) =    DT(I,L)
        FXOSP(I,1)   = FXOE2(I,LM1)
        FXOSP(I,2)   =   FXO(I,LM1)
         DTSP(I,1)   =  DTE2(I,LM1)
         DTSP(I,2)   =    DT(I,LM1)
105 END DO
!------------------------------------ 
! SOURCE FUNCTION FOR COMBINED BAND 1
!------------------------------------  
    DO 411 I=MYIS,MYIE
        DO 411 K=1,LP1
            VTMP3(I,K) = SOURCE(IXO(I,K),1)
            DSORC(I,K) =  DSRCE(IXO(I,K),1)
411 END DO
!
    DO 412 K=1,LP1
        DO 412 I=MYIS,MYIE
            SORC(I,K,1)  = VTMP3(I,K) + DT(I,K) * DSORC(I,K)
412 END DO
!------------------------------------  
! SOURCE FUNCTION FOR COMBINED BAND 2
!------------------------------------ 
    DO 421 I=MYIS,MYIE
        DO 421 K=1,LP1
            VTMP3(I,K) = SOURCE(IXO(I,K),2)
            DSORC(I,K) =  DSRCE(IXO(I,K),2)
421 END DO
!
    DO 212 K=1,LP1
        DO 212 I=MYIS,MYIE
            SORC(I,K,2)  = VTMP3(I,K) + DT(I,K) * DSORC(I,K)
212 END DO
!------------------------------------  
! SOURCE FUNCTION FOR COMBINED BAND 3
!------------------------------------  
    DO 314 I=MYIS,MYIE
        DO 314 K=1,LP1
            VTMP3(I,K) = SOURCE(IXO(I,K),3)
            DSORC(I,K) =  DSRCE(IXO(I,K),3)
314 END DO
!
    DO 312 K=1,LP1
        DO 312 I=MYIS,MYIE
            SORC(I,K,3)  = VTMP3(I,K) + DT(I,K) * DSORC(I,K)
312 END DO
!------------------------------------ 
! SOURCE FUNCTION FOR COMBINED BAND 4
!------------------------------------ 
    DO 444 I=MYIS,MYIE
        DO 444 K=1,LP1
            VTMP3(I,K) = SOURCE(IXO(I,K),4)
            DSORC(I,K) =  DSRCE(IXO(I,K),4)
444 END DO
!
    DO 442 K=1,LP1
        DO 442 I=MYIS,MYIE
            SORC(I,K,4)  = VTMP3(I,K) + DT(I,K) * DSORC(I,K)
442 END DO
!------------------------------------  
! SOURCE FUNCTION FOR COMBINED BAND 5
!------------------------------------ 
    DO 514 I=MYIS,MYIE
        DO 514 K=1,LP1
            VTMP3(I,K) = SOURCE(IXO(I,K),5)
            DSORC(I,K) =  DSRCE(IXO(I,K),5)
514 END DO
!
    DO 512 K=1,LP1
        DO 512 I=MYIS,MYIE
            SORC(I,K,5)  = VTMP3(I,K) + DT(I,K) * DSORC(I,K)
512 END DO
!------------------------------------  
! SOURCE FUNCTION FOR COMBINED BAND 6
!------------------------------------  
    DO 614 I=MYIS,MYIE
        DO 614 K=1,LP1
            VTMP3(I,K) = SOURCE(IXO(I,K),6)
            DSORC(I,K) =  DSRCE(IXO(I,K),6)
614 END DO
!
    DO 612 K=1,LP1
        DO 612 I=MYIS,MYIE
            SORC(I,K,6)  = VTMP3(I,K) + DT(I,K) * DSORC(I,K)
612 END DO
!------------------------------------
! SOURCE FUNCTION FOR COMBINED BAND 7
!------------------------------------
    DO 714 I=MYIS,MYIE
        DO 714 K=1,LP1
            VTMP3(I,K) = SOURCE(IXO(I,K),7)
            DSORC(I,K) =  DSRCE(IXO(I,K),7)
714 END DO
!
    DO 712 K=1,LP1
        DO 712 I=MYIS,MYIE
            SORC(I,K,7)  = VTMP3(I,K) + DT(I,K) * DSORC(I,K)
712 END DO
!------------------------------------
! SOURCE FUNCTION FOR COMBINED BAND 8
!------------------------------------
    DO 814 I=MYIS,MYIE
        DO 814 K=1,LP1
            VTMP3(I,K) = SOURCE(IXO(I,K),8)
            DSORC(I,K) =  DSRCE(IXO(I,K),8)
    814 END DO
!
    DO 812 K=1,LP1
        DO 812 I=MYIS,MYIE
            SORC(I,K,8)  = VTMP3(I,K) + DT(I,K) * DSORC(I,K)
812 END DO
!------------------------------------------
! SOURCE FUNCTION FOR BAND 9 (560-670 CM-1)
!------------------------------------------
    DO 914 I=MYIS,MYIE
        DO 914 K=1,LP1
            VTMP3(I,K) = SOURCE(IXO(I,K),9)
            DSORC(I,K) =  DSRCE(IXO(I,K),9)
914 END DO
!
    DO 912 K=1,LP1
        DO 912 I=MYIS,MYIE
            SORC(I,K,9)  = VTMP3(I,K) + DT(I,K) * DSORC(I,K)
912 END DO
!-------------------------------------------
! SOURCE FUNCTION FOR BAND 10 (670-800 CM-1)
!-------------------------------------------
    DO 518 I=MYIS,MYIE
        DO 518 K=1,LP1
            VTMP3(I,K) = SOURCE(IXO(I,K),10)
            DSORC(I,K) =  DSRCE(IXO(I,K),10)
518 END DO
!
    DO 516 K=1,LP1
        DO 516 I=MYIS,MYIE
            SORC(I,K,10) = VTMP3(I,K) + DT(I,K) * DSORC(I,K)
516 END DO
!-------------------------------------------
! SOURCE FUNCTION FOR BAND 11 (800-900 CM-1)
!-------------------------------------------
    DO 628 I=MYIS,MYIE
        DO 628 K=1,LP1
            VTMP3(I,K) = SOURCE(IXO(I,K),11)
            DSORC(I,K) =  DSRCE(IXO(I,K),11)
628 END DO
!
    DO 626 K=1,LP1
        DO 626 I=MYIS,MYIE
            SORC(I,K,11) = VTMP3(I,K) + DT(I,K) * DSORC(I,K)
626 END DO
!-------------------------------------------
! SOURCE FUNCTION FOR BAND 12 (900-990 CM-1)
!-------------------------------------------
    DO 528 I=MYIS,MYIE
        DO 528 K=1,LP1
            VTMP3(I,K) = SOURCE(IXO(I,K),12)
            DSORC(I,K) =  DSRCE(IXO(I,K),12)
    528 END DO
!
    DO 526 K=1,LP1
        DO 526 I=MYIS,MYIE
            SORC(I,K,12) = VTMP3(I,K) + DT(I,K) * DSORC(I,K)
    526 END DO
!--------------------------------------------
! SOURCE FUNCTION FOR BAND 13 (990-1070 CM-1)
!--------------------------------------------
    DO 428 I=MYIS,MYIE
        DO 428 K=1,LP1
            VTMP3(I,K) = SOURCE(IXO(I,K),13)
            DSORC(I,K) =  DSRCE(IXO(I,K),13)
    428 END DO
!
    DO 426 K=1,LP1
        DO 426 I=MYIS,MYIE
            SORC(I,K,13) = VTMP3(I,K) + DT(I,K) * DSORC(I,K)
426 END DO
!---------------------------------------------
! SOURCE FUNCTION FOR BAND 14 (1070-1200 CM-1)
!---------------------------------------------
    DO 328 I=MYIS,MYIE
        DO 328 K=1,LP1
            VTMP3(I,K) = SOURCE(IXO(I,K),14)
            DSORC(I,K) =  DSRCE(IXO(I,K),14)
328 END DO
!
    DO 326 K=1,LP1
        DO 326 I=MYIS,MYIE
            SORC(I,K,14) = VTMP3(I,K) + DT(I,K) * DSORC(I,K)
326 END DO
!---------------------------------------------------------------------------------------
! THE FOLLOWING SUBROUTINE OBTAINS NLTE SOURCE FUNCTION FOR CO2
!
! CALL NLTE
!
! OBTAIN SPECIAL SOURCE FUNCTIONS FOR THE 15 UM BAND (CSOUR) AND THE WINDOW REGION (SS1)
!---------------------------------------------------------------------------------------
    DO 131 K=1,LP1
        DO 131 I=MYIS,MYIE
            SS1(I,K) = SORC(I,K,11) + SORC(I,K,12) + SORC(I,K,14)
131 END DO
!
    DO 143 K=1,LP1
        DO 143 I=MYIS,MYIE
            CSOUR(I,K) = SORC(I,K,9) + SORC(I,K,10)
143 END DO
!--------------------------------------------------------------------------- 
! COMPUTE TEMP 4 (TC) AND VERTICAL TEMPERATURE DIFFERENCES (OSS,CSS,SS2,DTC) 
! ALL THESE WILL BE USED LATER IN FLUX COMPUTATIONS.
!---------------------------------------------------------------------------
    DO 901 K=1,LP1
        DO 901 I=MYIS,MYIE
            TC(I,K) = (TEMP(I,K)*TEMP(I,K)) ** 2
901 END DO
!
    DO 903 K=1,L
        DO 903 I=MYIS,MYIE
            OSS(I,K+1) =  SORC(I,K+1,13) -  SORC(I,K,13)
            CSS(I,K+1) = CSOUR(I,K+1)    - CSOUR(I,K)
            DTC(I,K+1) =    TC(I,K+1)    -    TC(I,K)
            SS2(I,K+1) =   SS1(I,K+1)    -   SS1(I,K)
903 END DO
!--------------------------------------------------------------------------------------------------
! THE FOLLOWIMG IS A DRASTIC REWRITE OF THE RADIATION CODE TO (LARGELY) ELIMINATE THREE-DIMENSIONAL
! ARRAYS.
! THE CODE WORKS ON THE FOLLOWING PRINCIPLES:
!
! LET K = FIXED FLUX LEVEL,
!    KP = VARYING FLUX LEVELTHEN FLUX(K)=SUM OVER KP: (DELTAB(KP)*TAU(KP,K)) OVER ALL KPS, FROM 1 
!         TO LP1.
!
! WE CAN BREAK DOWN THE CALCULATIONS FOR ALL KS AS FOLLOWS:
!
! FOR ALL KS K=1 TO LP1:
! FLUX(K)=SUM OVER KP: (DELTAB(KP)*TAU(KP,K))
! (1) OVER ALL KPS, FROM K+1 TO LP1 AND FOR KP FROM K+1 TO LP1: FLUX(KP) = DELTAB(K)*TAU(K,KP)              
! (2) NOW IF TAU(K,KP)=TAU(KP,K) (SYMMETRICAL ARRAYS)
!      WE CAN COMPUTE A 1-DIMENSIONAL ARRAY TAU1D(KP) FROM K+1 TO LP1, EACH TIME K IS INCREMENTED.
! EQUATIONS (1) AND (2) THEN BECOME:
!
! TAU1D(KP) = (VALUES FOR TAU(KP,K) AT THE PARTICULAR K)
! FLUX(K) = SUM OVER KP : (DELTAB(KP)*TAU1D(KP))   (3)
! FLUX(KP) = DELTAB(K)*TAU1D(KP)                   (4)
!
! THE TERMS FOR TAU (K,K) AND OTHER SPECIAL TERMS (FOR NEARBY LAYERS) MUST, OF COURSE, BE HANDLED
! SEPARATELY, AND WITH CARE.
!
! COMPUTE "UPPER TRIANGLE" TRANSMISSION FUNCTIONS FOR THE 9.6 UM BAND (TO3SP) AND THE 15 UM BAND 
! (OVER1D). ALSO, THE STAGE 1...COMPUTE O3 ,OVER TRANSMISSION FCTNS AND AVEPHI
!
! DO K=1 CALCULATION (FROM FLUX LAYER KK TO THE TOP) SEPARATELY AS VECTORIZATION IS IMPROVED,AND 
! OZONE CTS TRANSMISSIVITY MAY BE EXTRACTED HERE.
!--------------------------------------------------------------------------------------------------
    DO 121 K=1,L
        DO 121 I=MYIS,MYIE
            AVEPHI(I,K) = TOTPHI(I,K+1)
121 END DO
!--------------------------------------------------------------------------------------------------
! IN ORDER TO PROPERLY EVALUATE EMISS INTEGRATED OVER THE (LP1) LAYER, A SPECIAL EVALUATION OF 
! EMISS IS DONE THIS REQUIRES A SPECIAL COMPUTATION OF AVEPHI, AND IT IS STORED IN THE 
! (OTHERWISE VACANT) LP1TH POSITION
!--------------------------------------------------------------------------------------------------
    DO 803 I=MYIS,MYIE
        AVEPHI(I,LP1) = AVEPHI(I,LM1) + EMX1(I)
803 END DO
!----------------------- 
! COMPUTE FLUXES FOR K=1
!----------------------- 
    CALL E1E290(E1CTS1, E1CTS2, E1FLX , E1CTW1, E1CTW2, EMISS, FXO, DT,                           &
    &           FXOE2 , DTE2  , AVEPHI, TEMP  , T)
!
    DO 302 K=1,L
        DO 302 I=MYIS,MYIE
              FAC1(I,K) = BO3RND(2)       * TPHIO3(I,K+1) / TOTO3(I,K+1)
            TO3SPC(I,K) = HAF*(FAC1(I,K)  * (SQRT(ONE + (FOUR * AO3RND(2) * TOTO3(I,K+1))         &
    &                   /      FAC1(I,K)) - ONE))
!--------------------------------------------------------------------------------------------------
! FOR K=1, TO3SP IS USED INSTEAD OF TO31D (THEY ARE EQUAL IN THIS CASE); TO3SP IS PASSED TO SPA90, 
! WHILE TO31D IS A WORK-ARRAY.
!--------------------------------------------------------------------------------------------------
             TO3SP(I,K) = EXP(HM1EZ * (            TO3SPC(I,K)    + SKO3R * TOTVO2(I,K+1)))
            OVER1D(I,K) = EXP(HM1EZ * (SQRT(AB15WD*TOTPHI(I,K+1)) + SKC1R * TOTVO2(I,K+1)))
!--------------------------------------------------------------------------------------------------
! BECAUSE ALL CONTINUUM TRANSMISSIVITIES ARE OBTAINED FROM THE 2-D QUANTITY CNTTAU 
! (AND ITS RECIPROCAL TOTEVV) WE STORE BOTH OF THESE HERE. FOR K=1, CONT1D EQUALS CNTTAU
!--------------------------------------------------------------------------------------------------
            CNTTAU(I,K) = EXP(HM1EZ * TOTVO2(I,K+1))
            TOTEVV(I,K) = 1. / CNTTAU(I,K)
302 END DO
!
    DO 722 K=1,L
        DO 722 I=MYIS,MYIE
            CO2SP(I,K+1)  = OVER1D(I,K)     *   CO21(I,1,K+1)
722 END DO
!
    DO 723 K=1,L
        DO 723 I=MYIS,MYIE
            CO21(I,K+1,1) =   CO21(I,K+1,1) * OVER1D(I,K)
723 END DO
!------------------------------------------------------ 
! RLOG IS THE NBL AMOUNT FOR THE 15 UM BAND CALCULATION
!------------------------------------------------------ 
    DO 724 I=MYIS,MYIE
        RLOG(I,1) = OVER1D(I,1) * CO2NBL(I,1)
724 END DO
!-------------------------------------------------------------------------------------------------- 
! THE TERMS WHEN KP=1 FOR ALL K ARE THE PHOTON EXCHANGE WITH THE TOP OF THE ATMOSPHERE, AND ARE 
! OBTAINED DIFFERENTLY THAN THE  OTHER CALCULATIONS
!--------------------------------------------------------------------------------------------------
    DO 305 K=2,LP1
        DO 305 I=MYIS,MYIE
            FLX(I,K) = (TC(I,1) * E1FLX(I,K)   +   SS1(I,1) * CNTTAU(I,K-1)  +   SORC(I,1,13)     &
    &                *            TO3SP(I,K-1) + CSOUR(I,1) *  CO2SP(I,K))   * CLDFAC(I,1,K)
305 END DO
!
    DO 307 I=MYIS,MYIE
        FLX(I,1)     =  TC(I,1) * E1FLX(I,1)   +   SS1(I,1) +   SORC(I,1,13) + CSOUR(I,1)
307 END DO
!------------------------ 
! THE KP TERMS FOR K=1...
!------------------------ 
    DO 303 KP=2,LP1
        DO 303 I=MYIS,MYIE
            FLX(I,1) = FLX(I,1)  + (OSS(I,KP)   * TO3SP(I,KP-1) +   SS2(I,KP)    * CNTTAU(I,KP-1) &
    &                + CSS(I,KP) * CO21(I,KP,1) +   DTC(I,KP)   * EMISS(I,KP-1)) * CLDFAC(I,KP,1)
303 END DO
!-------------------------------------------------------------------------------------------------- 
! SUBROUTINE SPA88 IS CALLED TO OBTAIN EXACT CTS FOR WATER CO2 AND O3, AND APPROXIMATE CTS CO2 AND
! O3 CALCULATIONS.
!--------------------------------------------------------------------------------------------------
    CALL SPA88(EXCTS, CTSO3 , GXCTS, SORC  , CSOUR , CLDFAC, TEMP, PRESS, VAR1, VAR2, P, DELP,    &
    &          DELP2, TOTVO2, TO3SP, TO3SPC, CO2SP1, CO2SP2, CO2SP)
!-------------------------------------------------------------------------------------------------- 
! THIS SECTION COMPUTES THE EMISSIVITY CTS HEATING RATES FOR 2 EMISSIVITY BANDS: THE 0-160,
! 1200 -2200 CM-1 BAND AND THE 800-990, 1070-1200 CM-1 BAND. 
! THE REMAINING CTS COMTRIBUTIONS ARE CONTAINED IN CTSO3, COMPUTED IN SPA88.
!--------------------------------------------------------------------------------------------------
    DO 998 I=MYIS,MYIE
        VTMP3(I,1) = 1.
998 END DO
!
    DO 999 K=1,L
        DO 999 I=MYIS,MYIE
            VTMP3(I,K+1) = CNTTAU(I,K) * CLDFAC(I,K+1,1)
999 END DO
!
    DO 997 K=1,L
        DO 997 I=MYIS,MYIE
              CTS(I,K)   = RADCON * DELP(I,K) * (TC(I,K) * (E1CTW2(I,K)   * CLDFAC(I,K+1,1)       &
    &                    -                                  E1CTW1(I,K)   * CLDFAC(I,K,1))        &
    &                    +                      SS1(I,K) *  (VTMP3(I,K+1) -  VTMP3(I,K)))
997 END DO
!
    DO 996 K=1,L
        DO 996 I=MYIS,MYIE
            VTMP3(I,K)   = TC(I,K)   * (CLDFAC(I,K,1)   * (E1CTS1(I,K) -   E1CTW1(I,K))           &
    &                    -              CLDFAC(I,K+1,1) * (E1CTS2(I,K) -   E1CTW2(I,K)))
996 END DO
!
    DO 995 I=MYIS,MYIE
        FLX1E1(I)        = TC(I,LP1) *  CLDFAC(I,LP1,1) * (E1CTS1(I,LP1) - E1CTW1(I,LP1))
995 END DO
!
    DO 994 K=1,L
        DO 993 I=MYIS,MYIE
            FLX1E1(I)    = FLX1E1(I) + VTMP3(I,K)
    993 END DO
994 END DO
!--------------------------------------------------------------------------------------------------
! NOW REPEAT FLUX CALCULATIONS FOR THE K=2...LM1 CASES.
! CALCULATIONS FOR FLUX LEVEL L AND LP1 ARE DONE SEPARATELY, AS ALL EMISSIVITY AND CO2 CALCULATIONS
! ARE SPECIAL CASES OR NEARBY LAYERS.
!--------------------------------------------------------------------------------------------------
    DO 321 K=2,LM1
        KLEN = K
!    
        DO 3218 KK=1,LP1-K
            DO 3218 I=MYIS,MYIE
                AVEPHI(I,KK+K-1) = TOTPHI(I,KK+K) - TOTPHI(I,K)
   3218 END DO
!
        DO 1803 I=MYIS,MYIE
            AVEPHI(I,LP1)        = AVEPHI(I,LM1)  + EMX1(I)
   1803 END DO
!--------------------------------------------------------------------------------------------------
! COMPUTE EMISSIVITY FLUXES (E2) FOR THIS CASE. NOTE THAT WE HAVE OMITTED THE NEARBY LATER CASE 
! (EMISS(I,K,K)) AS WELL AS ALL CASES WITH K=L OR LP1.
! BUT THESE CASES HAVE ALWAYS BEEN HANDLED AS SPECIAL CASES, SO WE MAY AS WELL COMPUTE THEIR FLUXES 
! SEPARASTELY.
!--------------------------------------------------------------------------------------------------   
        CALL E290(EMISSB, EMISS, AVEPHI, KLEN, FXOE2, DTE2)
!
        DO 322 KK=1,LP1-K
            DO 322 I=MYIS,MYIE
                 AVMO3(I,KK+K-1) =  TOTO3(I,KK+K)   -  TOTO3(I,K)
                AVPHO3(I,KK+K-1) = TPHIO3(I,KK+K)   - TPHIO3(I,K)
                 AVVO2(I,KK+K-1) = TOTVO2(I,KK+K)   - TOTVO2(I,K)
                CONT1D(I,KK+K-1) = CNTTAU(I,KK+K-1) * TOTEVV(I,K-1)
    322 END DO
!    
        DO 3221 KK=1,LP1-K
            DO 3221 I=MYIS,MYIE
                  FAC1(I,K+KK-1) = BO3RND(2) * AVPHO3(I,K+KK-1)  / AVMO3(I,K+KK-1)
                 VTMP3(I,K+KK-1) = HAF * (FAC1(I,K+KK-1)                                          &
    &                            * (SQRT(ONE + (FOUR * AO3RND(2) * AVMO3(I,K+KK-1))               &
    &                            /  FAC1(I,K+KK-1)) - ONE))
!
                 TO31D(I,K+KK-1) = EXP(HM1EZ * (VTMP3(I,K+KK-1) + SKO3R * AVVO2(I,K+KK-1)))
                OVER1D(I,K+KK-1) = EXP(HM1EZ * (            SQRT(AB15WD * AVEPHI(I,K+KK-1))       &
    &                            +                                SKC1R * AVVO2(I,K+KK-1)))
!
                  CO21(I,K+KK,K) = OVER1D(I,K+KK-1) * CO21(I,K+KK,K)
   3221 END DO
!
        DO 3223 KP=K+1,LP1
            DO 3223 I=MYIS,MYIE
                CO21(I,K,KP) = OVER1D(I,KP-1) * CO21(I,K,KP)
   3223 END DO
!------------------------------------------------------ 
! RLOG IS THE NBL AMOUNT FOR THE 15 UM BAND CALCULATION
!------------------------------------------------------ 
        DO 1804 I=MYIS,MYIE
            RLOG(I,K) = OVER1D(I,K) * CO2NBL(I,K)
   1804 END DO
!-------------------------------- 
! THE KP TERMS FOR ARBIRRARY K...
!-------------------------------- 
        DO 3423 KP=K+1,LP1
            DO 3423 I=MYIS,MYIE
                FLX(I,K) =   FLX(I,K) + (OSS(I,KP) *  TO31D(I,KP-1)                               &
    &                    +               SS2(I,KP) * CONT1D(I,KP-1)                               &
    &                    +               CSS(I,KP) *   CO21(I,KP,K)                               &
    &                    +               DTC(I,KP) *  EMISS(I,KP-1))                              &
    &                    *                           CLDFAC(I,KP,K)
   3423 END DO
!
        DO 3425 KP=K+1,LP1
            DO 3425 I=MYIS,MYIE
                FLX(I,KP) = FLX(I,KP) + (OSS(I,K)    *  TO31D(I,KP-1)                             &
    &                     +              SS2(I,K)    * CONT1D(I,KP-1)                             &
    &                     +              CSS(I,K)    *   CO21(I,K,KP)                             &
    &                     +              DTC(I,K)    * EMISSB(I,KP-1))                            &
    &                     *                            CLDFAC(I,K,KP)
   3425 END DO
!
321 END DO
!--------------------------------------------------------------------------------------------------
! NOW DO K=L CASE. SINCE THE KP LOOP IS LENGTH 1, MANY SIMPLIFICATIONS OCCUR.
! ALSO, THE CO2 QUANTITIES (AS WELL AS THE EMISS QUANTITIES) ARE COMPUTED IN THE NBL SEDCTION; 
! THEREFORE, WE WANT ONLY OVER, TO3 AND CONT1D (OVER(I,L),TO31D(I,L) AND CONT1D(I,L) ACCORDING TO 
! THE NOTATION. THUS NO CALL IS MADE TO THE E290 SUBROUTINE.
! THE THIRD SECTION CALCULATES BOUNDARY LAYER AND NEARBY LAYER CORRECTIONS TO THE TRANSMISSION 
! FUNCTIONS OBTAINED ABOVE. METHODS ARE GIVEN IN REF. (4). 
! THE FOLLOWING RATIOS ARE USED IN VARIOUS NBL CALCULATIONS:
!
! THE REMAINING CALCULATIONS ARE FOR :
! 1) THE (K,K) TERMS, K=2,LM1;
! 2) THE (L,L) TERM
! 3) THE (L,LP1) TERM
! 4) THE (LP1,L) TERM
! 5) THE (LP1,LP1) TERM.
! EACH IS UNIQUELY HANDLED; DIFFERENT FLUX TERMS ARE COMPUTED DIFFERENTLY
!
! FOURTH SECTION OBTAINS WATER TRANSMISSION FUNCTIONS USED IN Q(APPROX) CALCULATIONS AND ALSO 
! MAKES NBL CORRECTIONS:
! 1) EMISS (I,J) IS THE TRANSMISSION FUNCTION MATRIX OBTAINED BY CALLING SUBROUTINE E1E288;
! 2) "NEARBY LAYER" CORRECTIONS (EMISS(I,I)) ARE OBTAINED USING SUBROUTINE E3V88;
! 3) SPECIAL VALUES AT THE SURFACE (EMISS(L,LP1),EMISS(LP1,L), EMISS(LP1,LP1)) ARE CALCULATED.
!
! OBTAIN ARGUMENTS FOR E1E288 AND E3V88:
!--------------------------------------------------------------------------------------------------
    DO 821 I=MYIS,MYIE
        TPL(I,1)    =                   TEMP(I,L)
        TPL(I,LP1)  = HAF * (T(I,LP1) + TEMP(I,L))
        TPL(I,LLP1) = HAF * (T(I,L)   + TEMP(I,L))
821 END DO
!
    DO 823 K=2,L
        DO 823 I=MYIS,MYIE
            TPL(I,K)   = T(I,K)
            TPL(I,K+L) = T(I,K)
823 END DO
!--------------------------------------------------------------------------------------------------
! E2 FUNCTIONS ARE REQUIRED IN THE NBL CALCULATIONS FOR 2 CASES, DENOTED (IN OLD CODE) AS (L,LP1) 
! AND (LP1,LP1)
!--------------------------------------------------------------------------------------------------
    DO 833 I=MYIS,MYIE
        AVEPHI(I,1) = VAR2(I,L)
        AVEPHI(I,2) = VAR2(I,L) + EMPL(I,L)
833 END DO
!
    CALL E2SPEC(EMISS, AVEPHI, FXOSP, DTSP)
!---------------------------------------- 
! CALL E3V88 FOR NBL H2O TRANSMISSIVITIES
!---------------------------------------- 
    CALL E3V88(EMD, TPL, EMPL)
!-------------------------------------------------------------------------------------------------- 
! COMPUTE NEARBY LAYER AND SPECIAL-CASE TRANSMISSIVITIES FOR EMISS USING METHODS FOR H2O GIVEN IN
! REF. (4)
!--------------------------------------------------------------------------------------------------
    DO 851 K=2,L
        DO 851 I=MYIS,MYIE
            EMISDG(I,K) = EMD(I,K+L) + EMD(I,K)
851 END DO
!-------------------------------------------------------------------  
! NOTE THAT EMX1/2 (PRESSURE SCALED PATHS) ARE NOW COMPUTED IN LWR88
!-------------------------------------------------------------------  
    DO 861 I=MYIS,MYIE
        EMSPEC(I,1)   = (EMD(I,1)   * EMPL(I,1)                                                   &
    &                 -  EMD(I,LP1) * EMPL(I,LP1))                                                &
    &                 / EMX1(I)                                                                   &
    &                 + QUARTR      * (EMISS(I,1) + EMISS(I,2))
!
        EMISDG(I,LP1) = TWO *  EMD(I,LP1)
        EMSPEC(I,2)   = TWO * (EMD(I,1)    * EMPL(I,1)                                            &
    &                 -        EMD(I,LLP1) * EMPL(I,LLP1))                                        &
    &                 / EMX2(I)
861 END DO
!
    DO 331 I=MYIS,MYIE
          FAC1(I,L) = BO3RND(2) * VAR4(I,L) / VAR3(I,L)
         VTMP3(I,L) = HAF * (FAC1(I,L)                                                            &
    &               * (SQRT(ONE + (FOUR * AO3RND(2) * VAR3(I,L)) / FAC1(I,L)) - ONE))
!
         TO31D(I,L) = EXP(HM1EZ   * (              VTMP3(I,L)  + SKO3R * CNTVAL(I,L)))
        OVER1D(I,L) = EXP(HM1EZ   * (SQRT(AB15WD *  VAR2(I,L)) + SKC1R * CNTVAL(I,L)))
        CONT1D(I,L) = CNTTAU(I,L) * TOTEVV(I,LM1)
          RLOG(I,L) = OVER1D(I,L) * CO2NBL(I,L)
331 END DO
!
    DO 618 K=1,L
        DO 618 I=MYIS,MYIE
            RLOG(I,K) = LOG(RLOG(I,K))
618 END DO
!
    DO 601 K=1,LM1
        DO 601 I=MYIS,MYIE
            DELPR1(I,K+1)     =         DELP(I,K+1)  * (PRESS(I,K+1) -     P(I,K+1))
               ALP(I,LP1+K-1) = -SQRT(DELPR1(I,K+1)) *   RLOG(I,K+1)
601 END DO
!
    DO 603 K=1,L
        DO 603 I=MYIS,MYIE
            DELPR2(I,K+1) = DELP(I,K) * (P(I,K+1)  - PRESS(I,K))
               ALP(I,K)   =   -SQRT(DELPR2(I,K+1)) *  RLOG(I,K)
603 END DO
!
    DO 625 I=MYIS,MYIE
        ALP(I,LL)   = -RLOG(I,L)
        ALP(I,LLP1) = -RLOG(I,L) * SQRT(DELP(I,L)    *     (P(I,LP1) - PRESS(I,LM1)))
625 END DO
!--------------------------------------------------------------------------------------------------
! THE FIRST COMPUTATION IS FOR THE 15 UM BAND,WITH THE FOR THE COMBINED H2O AND CO2 TRANSMISSION
! FUNCTION.
! 
! PERFORM NBL COMPUTATIONS FOR THE 15 UM BAND THE STATEMENT FUNCTION SF IN PREV. VERSIONS IS NOW 
! EXPLICITLY EVALUATED.
!-------------------------------------------------------------------------------------------------- 
    DO 631 K=1,LLP1
        DO 631 I=MYIS,MYIE
            C(I,K) = ALP(I,K) * (HMP66667 + ALP(I,K) * (QUARTR+ALP(I,K) * HM6666M2))
631 END DO
!
    DO 641 I=MYIS,MYIE
        CO21(I,LP1,LP1) = ONE + C(I,L)
        CO21(I,LP1,L)   = ONE + (DELP2(I,L)   *     C(I,LL)    - (PRESS(I,L)                      &
    &                   -            P(I,L))  *     C(I,LLM1))                                    &
    &                   /           (P(I,LP1) - PRESS(I,L))
!
        CO21(I,L,LP1)   = ONE + ((P(I,LP1) - PRESS(I,LM1)) * C(I,LLP1)                            &
    &                   -        (P(I,LP1) - PRESS(I,L))   * C(I,L))                              &
    &                   /    (PRESS(I,L)   - PRESS(I,LM1))
641 END DO
!
    DO 643 K=2,L
        DO 643 I=MYIS,MYIE
            CO21(I,K,K) = ONE + HAF * (C(I,LM1+K) + C(I,K-1))
643 END DO
!-------------------------------------------------------------------------------------------------- 
! COMPUTE NEARBY-LAYER TRANSMISSIVITIES FOR THE O3 BAND AND FOR THE ONE-BAND CONTINUUM BAND 
! (TO3 AND EMISS2). 
! THE SF2 FUNCTION IS USED. THE METHOD IS THE SAME AS DESCRIBED FOR CO2 IN REF (4).
!--------------------------------------------------------------------------------------------------
    DO 651 K=1,LM1
        DO 651 I=MYIS,MYIE
            CSUB(I,K+1)     = CNTVAL(I,K+1) * DELPR1(I,K+1)
            CSUB(I,LP1+K-1) = CNTVAL(I,K)   * DELPR2(I,K+1)
651 END DO
!---------------------------------------------------------------  
! THE SF2 FUNCTION IN PREV. VERSIONS IS NOW EXPLICITLY EVALUATED
!---------------------------------------------------------------  
    DO 655 K=1,LLM2
        DO 655 I=MYIS,MYIE
            CSUB2(I,K+1) = SKO3R       *  CSUB(I,K+1)
                C(I,K+1) = CSUB(I,K+1) * (HMP5        + CSUB(I,K+1)                               &
    &                    * (HP166666   -  CSUB(I,K+1) * H41666M2))
!
               C2(I,K+1) = CSUB2(I,K+1) * (HMP5       + CSUB2(I,K+1)                              &
    &                    * (HP166666   - CSUB2(I,K+1) * H41666M2))
655 END DO
!
    DO 661 I=MYIS,MYIE
        CONTDG(I,LP1) = 1. +  C(I,LLM1)
         TO3DG(I,LP1) = 1. + C2(I,LLM1)
661 END DO
!
    DO 663 K=2,L
        DO 663 I=MYIS,MYIE
            CONTDG(I,K) = ONE + HAF * ( C(I,K) +  C(I,LM1+K))
             TO3DG(I,K) = ONE + HAF * (C2(I,K) + C2(I,LM1+K))
663 END DO
!-------------------------------------------- 
! NOW OBTAIN FLUXES FOR THE DIAGONAL TERMS...
!--------------------------------------------  
    DO 871 K=2,LP1
        DO 871 I=MYIS,MYIE
            FLX(I,K) =    FLX(I,K) + (DTC(I,K) * EMISDG(I,K) + SS2(I,K)                           &
    &                * CONTDG(I,K) +  OSS(I,K) *  TO3DG(I,K) + CSS(I,K)                           &
    &                *   CO21(I,K,K))                                                             &
    &                * CLDFAC(I,K,K)
871 END DO
!----------------------------------  
! FOR THE TWO OFF-DIAGONAL TERMS...
!---------------------------------- 
    DO 873 I=MYIS,MYIE
        FLX(I,L)   =    FLX(I,L)  + (CSS(I,LP1) *  CO21(I,LP1,L) + DTC(I,LP1)                     &
    &              * EMSPEC(I,2)  +  OSS(I,LP1) * TO31D(I,L)     + SS2(I,LP1)                     &
    &              * CONT1D(I,L))                                                                 &
    &              * CLDFAC(I,LP1,L)
!
        FLX(I,LP1) =    FLX(I,LP1) + (CSS(I,L)  * CO21(I,L,LP1)  + OSS(I,L)                       &
    &              *  TO31D(I,L)   + SS2(I,L)   * CONT1D(I,L)    + DTC(I,L)                       &
    &              * EMSPEC(I,1))                                                                 &
    &              * CLDFAC(I,L,LP1)
873 END DO
!----------------------------------------------------------------------------------------------- 
! FINAL SECTION OBTAINS EMISSIVITY HEATING RATES, TOTAL HEATING RATES AND THE FLUX AT THE GROUND
!
! CALCULATE THE EMISSIVITY HEATING RATES
!----------------------------------------------------------------------------------------------- 
    DO 115 K=1,L
        DO 115 I=MYIS,MYIE
            HEATEM(I,K) = RADCON * (FLX(I,K+1) - FLX(I,K)) * DELP(I,K)
115 END DO
!----------------------------------  
! CALCULATE THE TOTAL HEATING RATES
!----------------------------------  
    DO 116 K=1,L
        DO 116 I=MYIS,MYIE
            HEATRA(I,K) = HEATEM(I,K) - CTS(I,K) - CTSO3(I,K) + EXCTS(I,K)
116 END DO
!--------------------------------------------------------------------------------------------------
! CALCULATE THE FLUX AT EACH FLUX LEVEL USING THE FLUX AT THE TOP (FLX1E1+GXCTS) AND THE INTEGRAL 
! OF THE HEATING RATES (VSUM1)
!--------------------------------------------------------------------------------------------------
    DO 111 K=1,L
        DO 111 I=MYIS,MYIE
            VSUM1(I,K) = HEATRA(I,K) * DELP2(I,K) * RADCON1
111 END DO
!
    DO 112 I=MYIS,MYIE
        TOPFLX(I)   = FLX1E1(I) + GXCTS(I)
        FLXNET(I,1) = TOPFLX(I)
112 END DO
!--------------------------------------------------------------------------------------------- 
! ONLY THE SURFACE VALUE OF FLUX (GRNFLX) IS NEEDED UNLESS THE THICK CLOUD SECTION IS INVOKED.
!--------------------------------------------------------------------------------------------- 
    DO 113 K=2,LP1
        DO 113 I=MYIS,MYIE
            FLXNET(I,K) = FLXNET(I,K-1) + VSUM1(I,K-1)
113 END DO
!
    DO 114 I=MYIS,MYIE
        GRNFLX(I) = FLXNET(I,LP1)
114 END DO
!-------------------------------------------------------------------------------------------------- 
! THIS IS THE THICK CLOUD SECTION.OPTIONALLY, IF THICK CLOUD FLUXES ARE TO BE "CONVECTIVELY 
! ADJUSTED", IE, DF/DP IS CONSTANT, FOR CLOUDY PART OF GRID POINT, THE FOLLOWING CODE IS EXECUTED.
! FIRST,COUNT THE NUMBER OF CLOUDS ALONG THE LAT. ROW. SKIP THE ENTIRE THICK CLOUD COMPUTATION OF 
! THERE ARE NO CLOUDS.
!--------------------------------------------------------------------------------------------------
    ICNT = 0
!
    DO 117 I=MYIS,MYIE
        ICNT = ICNT + NCLDS(I)
117 END DO
!
    IF (ICNT == 0) GO TO 6999
!------------------------------------------------------ 
! FIND THE MAXIMUM NUMBER OF CLOUDS IN THE LATITUDE ROW
!------------------------------------------------------  
    KCLDS=NCLDS(1)
    DO 118 I=MYIS,MYIE
        KCLDS = MAX(NCLDS(I),KCLDS)
118 END DO
!-------------------------------------------------------------------------------------------------- 
! OBTAIN THE PRESSURES AND FLUXES OF THE TOP AND BOTTOM OF THE NC'TH CLOUD (IT IS ASSUMED THAT ALL 
! KTOP AND KBTM'S HAVE BEEN DEFINED!).
!--------------------------------------------------------------------------------------------------
    DO 119 KK=1,KCLDS
        KMIN = LP1
        KMAX = 0
        DO 120 I=MYIS,MYIE
            J1 = KTOP(I,KK+1)
            J3 = KBTM(I,KK+1)
            IF (J3 > J1) THEN
                PTOP(I) =      P(I,J1)
                PBOT(I) =      P(I,J3+1)
                FTOP(I) = FLXNET(I,J1)
                FBOT(I) = FLXNET(I,J3+1)
!--------------------------------------------   
! OBTAIN THE "FLUX DERIVATIVE" DF/DP (DELPTC)
!--------------------------------------------  
                DELPTC(I) = (FTOP(I) - FBOT(I)) / (PTOP(I) - PBOT(I))
                KMIN = MIN(KMIN,J1)
                KMAX = MAX(KMAX,J3)
            END IF
    120 END DO
!
        KMIN = KMIN + 1
!------------------------------------------------------------------------   
! CALCULATE THE TOT. FLUX CHG. FROM THE TOP OF THE CLOUD, FOR ALL LEVELS.
!------------------------------------------------------------------------  
        DO 123 K=KMIN,KMAX
            DO 122 I=MYIS,MYIE
                IF (KTOP(I,KK+1) < K .AND. K <= KBTM(I,KK+1)) THEN
                    Z1(I,K)     = (P(I,K) - PTOP(I)) * DELPTC(I) + FTOP(I)
                    FLXNET(I,K) = Z1(I,K)
                END IF
        122 END DO
    123 END DO
119 END DO
!-------------------------------------------------------------------------------------------------- 
! USING THIS FLUX CHG. IN THE CLOUDY PART OF THE GRID BOX, OBTAIN THE NEW FLUXES, WEIGHTING THE 
! CLEAR AND CLOUDY FLUXES:AGAIN, ONLY THE FLUXES IN THICK-CLOUD LEVELS WILL EVENTUALLY BE USED.
!--------------------------------------------------------------------------------------------------
    6001 CONTINUE
    6999 CONTINUE
!-----------------------------------------------------------------------------
! HE FINAL STEP IS TO RECOMPUTE THE HEATING RATES BASED ON THE REVISED FLUXES:
!-----------------------------------------------------------------------------
    DO 124 K=1,L
        DO 124 I=MYIS,MYIE
            HEATRA(I,K) = RADCON * (FLXNET(I,K+1) - FLXNET(I,K)) * DELP(I,K)
124 END DO
!----------------------------------- 
! THE THICK CLOUD SECTION ENDS HERE.
!----------------------------------- 
    RETURN
!
    END SUBROUTINE FST88
