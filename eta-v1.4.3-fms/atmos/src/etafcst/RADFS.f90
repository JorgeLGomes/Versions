    SUBROUTINE RADFS(QS    , PP    , PPI   , QQH2O , TT    , O3QO3 , TSFC  , SLMSK , ALBEDO,      &
    &                XLAT  , CAMT  , ITYP  , KTOP  , KBTM  , NCLDS , EMCLD , RRCL  , TTCL  ,      &
    &                COSZRO, TAUDAR, IBEG  , KO3   , KALB  , SLMRF , SLYMRF, ITIMSW, ITIMLW,      &
    &                JD    , R1    , GMT   , SWH   , HLW   , FLWUP , FSWUP , FSWDN , FSWDNS,      &
    &                FSWUPS, FLWDNS, FLWUPS) 
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE RADFS
!>
!> SUBROUTINE: RADFS - RADIATION
!> PROGRAMMER:  Y.H AND K.A.C
!> ORG: ?????
!> DATE: 93-??-??
!> 
!> ABSTRACT:
!> THE INTERNAL DRIVE FOR GFDL RADIATION
!> THIS SUBROUTINE WAS FROM Y.H AND K.A.C (1993) AND MODIFIED BY Q. ZHAO FOR USE IN THE ETA MODEL 
!> 93-11-18
!>
!> UPDATE: THIS SUBROUTINE WAS MODIFIED TO USE CLOUD FRACTION ON EACH MODEL LAYER. 
!> QINGYUN ZHAO 95-3-22
!>
!> UPDATE: R1 HAS BEEN ADDED TO THE INPUTS FROM RADTN TO COMPUTE THE VARIATION OF SOLAR CONSTANT 
!>         AT THE TOP OF ATMOSPHERE WITH JULIAN DAY IN A YEAR. 
!>         QINGYUN ZHAO 96-7-23 
!>
!> PROGRAM HISTORY LOG:
!> 93-??-??  Y.H AND K.A.C    - ORIGINATOR
!> 18-01-15  LUCCI            - MODERNIZATION OF THE CODE, INCLUDING:
!>                              * F77 TO F90/F95
!>                              * INDENTATION & UNIFORMIZATION CODE
!>                              * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                              * DOCUMENTATION WITH DOXYGEN
!>                              * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> IX IS THE LENGTH OF A ROW IN THE DOMAIN
!>
!> QS(IX)          - THE SURFACE PRESSURE (PA)
!> PP(IX,L)        - THE MIDLAYER PRESSURES (PA)  (L IS THE VERT. DIMEN.)
!> PPI(IX,LP1)     - THE INTERFACE PRESSURES (PA)
!> QQH2O(IX,L)     - THE MIDLAYER WATER VAPOR MIXING RATIO (KG/KG)
!> TT(IX,L)        - THE MIDLAYER TEMPERATURE (K)
!> O3QO3(IX,L)     - THE MIDLAYER OZONE MIXING RATIO
!> TSFC(IX)        - THE SKIN TEMP. (K); NEGATIVE OVER WATER
!> SLMSK(IX)       - THE SEA MASK (LAND = 0, SEA = 1)
!> ALBEDO(IX)      - THE SURFACE ALBEDO (EXPRESSED AS A FRACTION)
!> XLAT(IX)        - THE GEODETIC LATITUDES OF EACH COLUMN IN DEGREES (N.H.> 0)
!>
!> THE FOLLOWING ARE CLOUD INFORMATION FOR EACH CLOUD LAYER
!>                 - LAYER = 1 - SURFACE
!>                 - LAYER = 2 - FIRST LAYER ABOVE GROUND, AND SO ON
!> CAMT(IX,LP1)    - CLOUD FRACTION OF EACH CLOUD LAYER
!> ITYP(IX,LP1)    - CLOUD TYPE = 1 - STRATIFORM, TYPE = 2 - CONVECTIVE
!> KTOP(IX,LP1)    - HEIGHT OF CLOUD TOP OF EACH CLOUD LAYER (IN ETA LEVEL)
!> KBTM(IX,LP1)    - BOTTOM OF EACH CLOUD LAYER
!> NCLDS(IX)       - NUMBER OF CLOUD LAYERS
!> EMCLD(IX,LP1)   - CLOUD EMISSIVITY
!> RRCL(IX,NB,LP1) - CLOUD REFLECTTANCES FOR SW SPECTRAL BANDS
!> TTCL(IX,NB,LP1) - CLOUD TRANSMITANCES FOR SW SPECTRAL BANDS
!>
!> THE ABOVE ARE CLOUD INFORMATION FOR EACH CLOUD LAYER
!> COSZRO(IX)      - THE COSINE OF THE SOLAR ZENITH ANGLE
!> TAUDAR          - = 1.0
!> IBEG            - = 1
!> KO3             - = 1 (READ IN THE QZONE DATA)
!> KALB            - = 0
!> SLMRF(LP1)      - THE INTERFACES ETA (LP1 = L + 1)
!> SLYMRF(L)       - THE MIDLAYER ETA
!> ITIMSW          - = 1 / 0 (SHORTWAVE CALC. ARE DESIRED/NOT DESIRED)
!> ITIMLW          - = 1 / 0 (LONGWAVE CALC. ARE DESIRED/NOT DESIRED)
!>
!> THE FOLLOWING ARE ADDITIONAL FOR ETA MODEL
!>
!> (JD, R1, GMT)
!>
!> JD  - JULIAN DAY IN A YEAR
!> R1  - THE NON-DIMENSIONAL SUN-EARTH DISTANCE
!> GMT - HOUR
!>
!> OUTPUT ARGUMENT LIST:
!> SWH    - ATMOSPHERIC SHORTWAVE HEATING RATES IN K/S. SWH IS A REAL ARRAY DIMENSIONED (NCOL X LM).
!> HLW    - ATMOSPHERIC LONGWAVE HEATING RATES IN K/S. HLW IS A REAL ARRAY DIMENSIONED (NCOL X LM).
!> FLWUP  - UPWARD LONGWAVE FLUX AT TOP OF THE ATMOSPHERE IN W/M**2. 
!>          FLWUP IS A REAL ARRAY DIMENSIONED (NCOL).
!> FSWUP  - UPWARD SHORTWAVE FLUX AT TOP OF THE ATMOSPHERE IN W/M**2. 
!>          FSWUP IS A REAL ARRAY DIMENSIONED (NCOL).
!> FSWDN  - DOWNWARD SHORTWAVE FLUX AT TOP OF THE ATMOSPHERE IN W/M**2. 
!>          FSWDN IS A REAL ARRAY DIMENSIONED (NCOL).
!> FSWDNS - DOWNWARD SHORTWAVE FLUX AT THE SURFACE IN W/M**2. 
!>          FSWDNS IS A REAL ARRAY DIMENSIONED (NCOL).
!> FSWUPS - UPWARD SHORTWAVE FLUX AT THE SURFACE IN W/M**2. 
!>          FSWUPS IS A REAL ARRAY DIMENSIONED (NCOL).
!> FLWDNS - DOWNWARD LONGWAVE FLUX AT THE SURFACE IN W/M**2. 
!>          FLWDNS IS A REAL ARRAY DIMENSIONED (NCOL).
!> FLWUPS - UPWARD LONGWAVE FLUX AT THE SURFACE IN W/M**2. 
!>          FLWUPS IS A REAL ARRAY DIMENSIONED (NCOL).
!>
!> USE MODULES: ASTSAV
!>              CO2DTA
!>              CTLBLK
!>              F77KINDS
!>              GLB_TABLE
!>              HCON
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              PHYCON
!>              RDFSAV
!>              RDPARM
!>              RNDDTA
!>              SAVMEM
!>              SSALB
!>              SWRSAV
!>              TABCOM
!>              TEMPCOM
!>              TOPO
!>              UPDATE_FLDS
!>
!> DRIVER     : RADTN
!>
!> CALLS      : CLO89
!>              GETDATE
!>              LWR88
!>              SWR93
!>--------------------------------------------------------------------------------------------------
    USE ASTSAV
    USE CO2DTA
    USE CTLBLK
    USE F77KINDS
    USE GLB_TABLE
    USE HCON
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE PHYCON
    USE RDFSAV
    USE RDPARM
    USE RNDDTA
    USE SAVMEM
    USE SSALB
    USE SWRSAV
    USE TABCOM
    USE TEMPCOM
    USE TOPO
    USE UPDT       , ONLY: LCO2
    USE UPDATE_FLDS
!
    IMPLICIT NONE
!
#include "sp.h"
!--------------------------------------------------------------------------------------------------
! GENERALIZED FOR PLUG-COMPATIBILITY - ORIGINAL CODE WAS CLEANED-UP GFDL CODE. K. CAMPANA MAR 89 
! EXAMPLE FOR MRF: 
! KO3    = 0  AND O3QO3 = DUMMY ARRAY. (GFDL CLIMO O3 USED)    
! KEMIS  = 0  AND HI CLD EMIS COMPUTED HERE (CEMIS=DUMMY INPUT) 
! KALB   = 0  AND SFC ALBEDO OVER OPEN WATER COMPUTED BELOW...  
! KCCO2  = 0, CO2 OBTAINED FROM BLOCK DATA = 1, CO2 COMPUTED IN HERE - NOT AVAILABLE YET 
! SLMRF  = INTERFACE (LEVELS) SIGMA    
! SLYMRF = LAYER SIGMA UPDATED FOR YUTAI HOU SIB SW RADIATION. KAC 06 MAR 92 
!          OCEAN ALBEDO FOR BEAM SET TO BULK SFCALB, SINCE COSINE ZENITH ANGLE EFFECTS ALREADY 
!          THERE (REF:PAYNE) 
! SLMSK = 0. 
! SNOW ICE ALBEDO FOR BEAM NOT ENHANCED VIA COSINE ZENITH ANGLE EITHER CAUSE VALU ALREADY HIGH 
! (WE SEE POLAR COOLING IF WE DO BEAM CALCULATION). KAC 17 MAR 92 
! ALBEDO GE .5 UPDATED TO OBTAIN CLEAR SKY FLUXES "ON THE FLY" FOR CLOUD FORCING DIAGNOSTICS 
! ELSEWHERE KAC 07 AUG 92 
! SEE ##CLR LINES. RADFS, LWR88, FST88, SPA88 UPDATED FOR USE NEW CLD SCHEME YH  DEC 92  
! INPUT CLD MAY BE AS ORIGINAL IN 3 DOMAIN (CLD, MTOP, MBOT) OR IN A VERTICAL ARRAY OF 18 MDL 
! LAYERS (CLDARY) IEMIS = 0 USE THE ORG. 
! CLD EMIS SCHEME = 1 USE TEMP DEP. 
! CLD EMIS SCHEME UPDATED TO COMPUTE CLD LAYER REFLECTTANCE AND TRANSMITTANCE
! INPUT CLD EMISSIVITY AND OPTICAL THICKNESS 'EMIS0, TAUC0' YH FEB 93         
!--------------------------------------------------------------------------------------------------
!--------------------------------  
! INPUT FROM FROM CALLING PROGRAM 
!-------------------------------- 
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2)                , INTENT (IN)                   ::&
    & NCLDS
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & ALVBR   , ALNBR   , ALVDR   , ALNDR
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                , INTENT (IN)                   ::&
    & QS      , ALBEDO  , XLAT    ,  TSFC   , SLMSK   , COSZRO  , TAUDAR  
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, L)             , INTENT (IN)                   ::&
    & PP      , QQH2O   , TT      , O3QO3    
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, LP1)           , INTENT (IN)                   ::&
    & ITYP    , KBTM     
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)           , INTENT (IN)                   ::&
    & PPI     , CAMT    , EMCLD
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, LP1)           , INTENT (IN)                   ::&
    & KTOP
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, NB, LP1)       , INTENT (INOUT)                ::&
    & RRCL    , TTCL  
!
    REAL   (KIND=R4KIND), DIMENSION(LP1)                        , INTENT (IN)                   ::&
    & SLMRF
!
    REAL   (KIND=R4KIND), DIMENSION(L)                          , INTENT (IN)                   ::&
    & SLYMRF
!-------------------------- 
! OUTPUT TO CALLING PROGRAM 
!-------------------------- 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, L)             , INTENT (INOUT)                ::&
    & SWH     , HLW    
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                , INTENT (INOUT)                ::&
    & FSWUP   , FSWUPS  , FSWDN   , FSWDNS  , FLWUP   , FLWUPS  , FLWDNS
!------------------------------------------------
! DOWNWARD SW FLUXES FOR THE SIB PARAMETERIZATION
!------------------------------------------------
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & GDFVBR  , GDFNBR  , GDFVDR  , GDFNDR
!-------------------------------------------------------
! ARRAYS NEEDED BY SWR91SIB..FOR CLEAR SKY DATA(EG.FSWL) 
!-------------------------------------------------------
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & FSWL    , HSWL    , UFL     , DFL
!------------------------------------------------- 
! ARRAYS NEEDED BY CLO88, LWR88, SWR89 OR SWR91SIB 
!------------------------------------------------- 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & EQCMT   , PRESS   , TEMP    , FSW     , HSW     , UF      , DF
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1, LP1)                                      ::&
    & CLDFAC
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, L)                                             ::&
    & RH2O    , QO3     , HEATRA
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & GRNFLX  , GRDFLX  , TOPFLX  , COSZEN  , TAUDA
!----------------------- 
! ADD PRESSURE INTERFACE
!----------------------- 
!-----------------------------------  
! VECTOR TEMPORARIES FOR CLOUD CALC. 
!----------------------------------- 
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & JJROW
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & DO3V    , DO3VP   , TTHAN
!-------------------------------------------------------------------------------------------------- 
! SEASONAL CLIMATOLOGIES OF O3 (OBTAINED FROM A PREVIOUSLY RUN CODE WHICH INTERPOLATES O3 TO USER 
! VERTICAL COORDINATE).
! DEFINED AS 5 DEG LAT MEANS N. P. -> S. P.
!--------------------------------------------------------------------------------------------------
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & ALVB    , ALNB    , ALVD    , ALND    , GDFVB   , GDFNB   , GDFVD   , GDFND   , SFCALB
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & MTOPSW  , MBTMSW
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & XAMT  
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & IBEG    , KO3     , KALB    , ITIMSW  , ITIMLW  , JD
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , IR      , K       , IQ      , JX      , N     
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & R1      , GMT     
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & YEAR    , RLAG    , TPI     , PI      , SC      , SSOLAR  , DATE    ,   &
    & RANG    , TH2     , ZEN     , DZEN    , ALB01   , ALB02   , ALBD0   , ALVD1   , ALND1   ,   &
    & RRVCO2  , RRCO2    
!
    INTEGER (KIND=I4KIND)                                                                       ::&
    & IYR     , IMO     , IDY     , IUTC    , YR      , MON     , DAY     , UTC 
!----------- 
! BEGIN HERE
!-----------
   IF (LCO2) THEN
       IYR   = IDAT(3)
       IMO   = IDAT(1)
       IDY   = IDAT(2)
       IUTC = NTSD*DT/3600
!
       CALL GETDATE(IYR, IMO, IDY, IUTC, YR, MON, DAY, UTC)
!-------------------------------------------------------------
!CO2 TRANSMISSIVITY COEFF. UPDATE FOR CLIMATE CHANGE SCENARIOS 
!-------------------------------------------------------------
       IF (NTSD.eq.1) THEN
!
           CALL FLDS_UPDATE_DRIVER('UPDATE          ','CO2COEFTRANSMISS')
!
       ELSE IF ( MOD((YR-IYR),3) == 0 .AND. MON  == 1 .AND. DAY ==1 .AND. UTC == 0) THEN
!      
          CALL FLDS_UPDATE_DRIVER('UPDATE          ','CO2COEFTRANSMISS')
!
       ENDIF
   ENDIF
!----------------------------- 
! NOTE: XLAT IS IN DEGREE HERE
!----------------------------- 
    YEAR = 365.25
    RLAG = 14.8125
    TPI  = 6.283185308
    PI   = 3.1415927
    SC   = 2.
    SOLC = SC / (R1 * R1)
!--------------------------------------------------------------------------------------------------
! SPECIAL NOTE: THE SOLAR CONSTANT IS REDUCED EXTRA 3 PERCENT TO ACCOUNT FOR THE LACK OF AEROSOLS 
! IN THE SHORTWAVE RADIATION
!               PARAMETERIZATION. Q. ZHAO 96-7-23
!--------------------------------------------------------------------------------------------------
    SSOLAR = SOLC * HP98
!   SSOLAR = SSOLAR * 0.93
    SSOLAR = SSOLAR * 0.97
!
    DATE   = JD + GMT / 24.0
    RANG   = TPI * (DATE - RLAG) / YEAR
    RSIN1  = SIN(RANG)
    RCOS1  = COS(RANG)
    RCOS2  = COS(2.0 * RANG)
!
    DO 40 I=MYIS,MYIE
        IR  = I + IBEG - 1
        TH2 = HP2 * XLAT(IR)
        JJROW(I) = Q19001 - TH2
        TTHAN(I) = (19 - JJROW(I)) - TH2
!--------------------------------------------------------------------------------------------------
! NOTE THAT THE NMC VARIABLES ARE IN MKS (THUS PRESSURE IS IN CENTIBARS) WHILE ALL GFDL VARIABLES 
! ARE IN CGS UNITS
!--------------------------------------------------------------------------------------------------
        SFCALB(I) = ALBEDO(IR)
!--------------------------------------------------------------------------------------- 
! NOW PUT SFC TEMP, PRESSURES, ZENITH ANGLE INTO SW COMMON BLOCK. ZHAO
! NOTE: ALL PRESSURES INPUT FROM THE ETA MODEL ARE IN PA THE UNIT FOR PRESS IS MICRO BAR
! SURFACE TEMPERATURE ARE NEGATIVE OVER OCEANS IN THE ETA MODEL. ZHAO
!--------------------------------------------------------------------------------------- 
        PRESS(I,LP1) = QS(IR) * 10.0
        TEMP (I,LP1) = ABS(TSFC(IR))
!
        COSZEN(I) = COSZRO(IR)
        TAUDA (I) = TAUDAR(IR)
 40 END DO
!--------------------------------------------------------------------------------------------------
! ALL GFDL VARIABLES HAVE K=1 AT THE TOP OF THE ATMOSPHERE. NMC ETA MODEL HAS THE SAME STRUCTURE. 
! ZHAO
!--------------------------------------------------------------------------------------------------
    DO 50 K=1,L
        DO 50 I=MYIS,MYIE
            IR = I + IBEG - 1
!---------------------------------------------
! NOW PUT TEMP,PRESSURES, INTO SW COMMON BLOCK 
!---------------------------------------------
            TEMP (I,K) =        TT(IR,K)
            PRESS(I,K) = 10.0 * PP(IR,K)
!----------------------------------------------
! STORE LYR MOISTURE AND ADD TO SW COMMON BLOCK
!----------------------------------------------
            RH2O(I,K) = QQH2O(IR,K)
!
            IF (RH2O(I,K) < H3M6) RH2O(I,K) = H3M6
 50 END DO
!
    IF (KO3 == 0) GOTO 65
!
    DO 60 K=1,L
        DO 60 I=MYIS,MYIE
            QO3(I,K) = O3QO3(I+IBEG-1,K)
 60 END DO
!
 65 CONTINUE
!
    IF (KALB > 0) GOTO 110
!------------------------------------------------------------------------------------------------
! THE FOLLOWING CODE GETS ALBEDO FROM PAYNE,1972 TABLES IF 1) OPEN SEA POINT (SLMSK=1); 2) KALB=0
!------------------------------------------------------------------------------------------------
    IQ = INT(TWENTY * HP537 + ONE)
!
    DO 105 I=MYIS,MYIE
        IF (COSZEN(I) > 0.0 .AND. SLMSK(I+IBEG-1) > 0.5) THEN
            ZEN = DEGRAD * ACOS(MAX(COSZEN(I),0.0))
!
            IF (ZEN >= H74E1)                   JX = INT(HAF    * (HNINETY - ZEN) + ONE  )
            IF (ZEN < H74E1 .AND. ZEN >= FIFTY) JX = INT(QUARTR * (H74E1   - ZEN) + HNINE)
            IF (ZEN < FIFTY)                    JX = INT(HP1    * (FIFTY   - ZEN) + H15E1)
            DZEN = -(ZEN - ZA(JX)) / DZA(JX)
!
            ALB01 = ALBD(IQ  ,JX)   + DZEN * (ALBD(IQ  ,JX+1) - ALBD(IQ  ,JX))
            ALB02 = ALBD(IQ+1,JX)   + DZEN * (ALBD(IQ+1,JX+1) - ALBD(IQ+1,JX))
!
            SFCALB(I) = ALB01 + TWENTY * (ALB02 - ALB01) * (HP537 - TRN(IQ))
        END IF
105 END DO
!
110 CONTINUE
!
    IF (KO3 > 0) GOTO 135
!----------------------------------------------------------------------------------------
! COMPUTE CLIMATOLOGICAL ZONAL MEAN OZONE, SEASONAL AND SPATIAL INTERPOLATION DONE BELOW.
!----------------------------------------------------------------------------------------
    DO 130 K=1,L
        DO 125 I=MYIS,MYIE
            DO3V (I) = DDUO3N(JJROW(I),K) + RSIN1 * DDO3N2(JJROW(I),K)     + RCOS1                &
    &                * DDO3N3(JJROW(I),K) + RCOS2 * DDO3N4(JJROW(I),K)
!
            DO3VP(I) = DDUO3N(JJROW(I)+1,K) + RSIN1 * DDO3N2(JJROW(I)+1,K) + RCOS1                &
    &                * DDO3N3(JJROW(I)+1,K) + RCOS2 * DDO3N4(JJROW(I)+1,K)
!------------------------------------------------------------------------------------------------
! NOW LATITUDINAL INTERPOLATION, AND CONVERT O3 INTO MASS MIXING RATIO(ORIGINAL DATA MPY BY 1.E4)
!------------------------------------------------------------------------------------------------
            QO3(I,K) = H1M4 * (DO3V(I) + TTHAN(I) * (DO3VP(I) - DO3V(I)))
    125 END DO
130 END DO
!
135 CONTINUE
!
    DO 195 I=MYIS,MYIE
!----------------------------------- 
! VISIBLE AND NEAR IR DIFFUSE ALBEDO
!----------------------------------- 
        ALVD(I) = SFCALB(I)
        ALND(I) = SFCALB(I)
!--------------------------------------- 
! VISIBLE AND NEAR IR DIRECT BEAM ALBEDO
!--------------------------------------- 
        ALVB(I) = SFCALB(I)
        ALNB(I) = SFCALB(I)
!--------------------------------------------------------------------------------------------------- 
! VISIBLE AND NEAR IR DIRECT BEAM ALBEDO,IF NOT OCEAN NOR SNOW FUNCTION OF COSINE SOLAR ZENITH ANGLE
!--------------------------------------------------------------------------------------------------- 
        IF (SLMSK(I+IBEG-1) < 0.5) THEN
            IF (SFCALB(I) <= 0.5) THEN
                ALBD0 = -18.0 * (0.5 - ACOS(COSZEN(I)) / PI)
                ALBD0 = EXP(ALBD0)
                ALVD1 = (ALVD(I) - 0.054313) / 0.945687
                ALND1 = (ALND(I) - 0.054313) / 0.945687
                ALVB(I) = ALVD1 + (1.0 - ALVD1) * ALBD0
                ALNB(I) = ALND1 + (1.0 - ALND1) * ALBD0
            END IF
        END IF
195 END DO
!-------------------------------- 
! SURFACE VALUES OF RRCL AND TTCL
!--------------------------------
    DO 200 N=1,2
        DO 200 I=MYIS,MYIE
            RRCL(I,N,1) = ALVD(I)
            TTCL(I,N,1) = ZERO
200 END DO
!
    DO 220 N=3,NB
        DO 220 I=MYIS,MYIE
            RRCL(I,N,1) = ALND(I)
            TTCL(I,N,1) = ZERO
220 END DO
!--------------------- 
! END OF CLOUD SECTION 
!--------------------- 
!
!--------------------------------------------------------------------------------------------------
! THE FOLLOWING CODE CONVERTS RRVCO2,THE VOLUME MIXING RATIO OF CO2 INTO RRCO2,THE MASS MIXING 
! RATIO.
!--------------------------------------------------------------------------------------------------
!
!-------------------------
! CO2 CONCENTRATION UPDATE
!-------------------------
!    CO21YEAR = ( 12 * 30 * 3600 * 24 ) / DT
!   IF (MOD(NTSD, CO21YEAR) == 0) THEN
    IF ( MOD((YR-IYR),1) == 0 .AND. MON  == 1 .AND. DAY ==1 .AND. UTC == 0) THEN
        CALL FLDS_UPDATE_DRIVER('UPDATE          ','CO2CONCENTRATION')
    END IF
!
    RRVCO2 = RCO2
    RRCO2  = RRVCO2 * RATCO2MW
!
    250 IF (ITIMLW == 0) GOTO 300
!-------------------- 
! LONG WAVE RADIATION 
!-------------------- 
!
!-------------------------------------------
! ACCOUNT FOR REDUCED EMISSIVITY OF ANY CLDS
!-------------------------------------------
    DO 240 K=1,LP1
        DO 240 I=MYIS,MYIE
            EQCMT(I,K) = CAMT(I,K) * EMCLD(I,K)
240 END DO
!----------------------------------- 
! GET CLD FACTOR FOR LW CALCULATIONS
!----------------------------------- 
    CALL CLO89(CLDFAC, EQCMT, NCLDS, KBTM, KTOP)
!--------------------  
! LONG WAVE RADIATION
!-------------------- 
    CALL LWR88(HEATRA, GRNFLX, TOPFLX, PRESS, TEMP, RH2O, QO3, CLDFAC, EQCMT, NCLDS, KTOP, KBTM)
!
    DO 280 I=MYIS,MYIE
        IR = I + IBEG - 1
        FLWUP(IR) = TOPFLX(I) * .001E0
        GRNFLX(I) = Q14330 * (HSIGMA * TEMP(I,LP1) ** 4 - GRNFLX(I))
!----------------------------------------------------------------
! GET LW FLUX DOWN AND UP AT GROUND(WATTS/M**2) - GRNFLX=LW DOWN.
!----------------------------------------------------------------
        FLWDNS(IR) = GRNFLX(I) / (1.43306E-06 * 1000.E0)
        FLWUPS(IR) = HSIGMA * .001E0  * TEMP(I,LP1) ** 4
280 END DO
!--------------------------------- 
! CONVERT HEATING RATES TO DEG/SEC
!--------------------------------- 
    DO 290 K=1,L
        DO 290 I=MYIS,MYIE
            HLW(I+IBEG-1,K) = HEATRA(I,K) * DAYSEC
290 END DO
!
300 CONTINUE
!
    IF (ITIMSW == 0) GOTO 350
!
    CALL SWR93(FSW   , HSW   , UF    , DF    , FSWL  , HSWL  , UFL   , DFL   , PRESS  , COSZEN,   &
    &          TAUDA , RH2O  , RRCO2 , SSOLAR, QO3   , NCLDS , KTOP  , KBTM  , CAMT   , RRCL  ,   &
    &          TTCL  , ALVB  , ALNB  , ALVD  , ALND  , GDFVB , GDFNB , GDFVD , GDFND)
!----------------------------  
! GET SW FLUXES IN WATTS/M**2
!---------------------------- 
    DO 320 I=MYIS,MYIE
        IR = I + IBEG - 1
        FSWUP (IR) = UF(I,1)   * 1.E-3
        FSWDN (IR) = DF(I,1)   * 1.E-3
        FSWUPS(IR) = UF(I,LP1) * 1.E-3
!------------------------------------------------- 
! COUPLE W/M2 DIFF, IF FSWDNS(IR)=DF(I,LP1)*1.#E-3
!------------------------------------------------- 
        FSWDNS(IR) = (GDFVB(I) + GDFNB(I) + GDFVD(I) + GDFND(I)) * 1.E-3
!------------------------------------------------------------------------- 
! DOWNWARD SFC FLUX FOR THE SIB PARAMETERATION VISIBLE AND NEAR IR DIFFUSE
!------------------------------------------------------------------------- 
        GDFVDR(IR) = GDFVD(I) * 1.E-3
        GDFNDR(IR) = GDFND(I) * 1.E-3
!--------------------------------  
! VISIBLE AND NEAR IR DIRECT BEAM
!--------------------------------  
        GDFVBR(IR) = GDFVB(I) * 1.E-3
        GDFNBR(IR) = GDFNB(I) * 1.E-3
320 END DO
!--------------------------------- 
! CONVERT HEATING RATES TO DEG/SEC
!--------------------------------- 
    DO 330 K=1,L
        DO 330 I=MYIS,MYIE
            SWH(I+IBEG-1,K) = HSW(I,K) * DAYSEC
330 END DO
!
350 CONTINUE
!
    RETURN
!
    1000 FORMAT(1H ,' YOU ARE CALLING GFDL RADIATION CODE FOR',I5,' PTS','AND',I4,                &
    &               ' LYRS,WITH KDAPRX,KO3,KCZ,KEMIS,KALB = ',5I2)
!
    END SUBROUTINE RADFS
