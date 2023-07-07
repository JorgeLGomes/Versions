    SUBROUTINE SFLX(ICE     , DT      , Z       , NSOIL   , SLDPTH  , LWDN    , SOLDN   ,         &
    &               SFCPRS  , PRCP    , SFCTMP  , TH2     , Q2      , SFCSPD  , Q2SAT   ,         &
    &               DQSDT2  , VEGTYP  , SOILTYP , SLOPETYP, SHDFAC  , PTU     , TBOT    ,         &
    &               ALB     , SNOALB  , CMC     , T1      , STC     , SMC     , SH2O    ,         &
    &               SNOWH   , SNEQV   , ALBEDO  , CH      , CM      , ETP     , ETA     ,         &
    &               H       , S       , RUNOFF1 , RUNOFF2 , Q1      , SNMAX   , SOILW   ,         &
    &               SOILM   , SMCWLT  , SMCDRY  , SMCREF  , SMCMAX )          
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SFLX
!>
!> PURPOSE: SUB-DRIVER FOR "NOAH/OSU LSM" FAMILY OF PHYSICS SUBROUTINES FOR A SOIL/VEG/SNOWPACK 
!> LAND-SURFACE MODEL TO UPDATE SOIL MOISTURE, SOIL ICE, SOIL TEMPERATURE, SKIN TEMPERATURE, 
!> SNOWPACK WATER CONTENT, SNOWDEPTH, AND ALL TERMS OF THE SURFACE ENERGY BALANCE AND SURFACE WATER
!> BALANCE (EXCLUDING INPUT ATMOSPHERIC FORCINGS OF DOWNWARD RADIATION AND PRECIP)
!>
!> VERSION 2.2 07 FEBRUARY 2001
!>
!> FROZEN GROUND VERSION 
!> ADDED STATES: SH2O(NSOIL) - UNFROZEN SOIL MOISTURE
!>               SNOWH       - SNOW DEPTH
!>
!> NOTE ON SNOW STATE VARIABLES:
!> SNOWH = ACTUAL PHYSICAL SNOW DEPTH IN M
!> SNEQV = LIQUID WATER-EQUIVALENT SNOW DEPTH IN M (TIME-DEPENDENT SNOW DENSITY IS OBTAINED FROM 
!>         SNEQV/SNOWH)
!>
!> NOTE ON ALBEDO FRACTIONS:
!> INPUT:
!> ALB    - BASELINE SNOW-FREE ALBEDO, FOR JULIAN DAY OF YEAR (USUALLY FROM TEMPORAL INTERPOLATION 
!>          OF MONTHLY MEAN VALUES) 
!>          (CALLING PROG MAY OR MAY NOT INCLUDE DIURNAL SUN ANGLE EFFECT)
!> SNOALB - UPPER BOUND ON MAXIMUM ALBEDO OVER DEEP SNOW (E.G. FROM ROBINSON AND KUKLA, 1985, J.
!>          CLIM. & APPL. METEOR.)
!>
!> OUTPUT:
!> ALBEDO - COMPUTED ALBEDO WITH SNOWCOVER EFFECTS(COMPUTED USING ALB, SNOALB, SNEQV, AND 
!>          SHDFAC->GREEN VEG FRAC)
!>
!> ARGUMENT LIST IN THE CALL TO SFLX
!>
!> 1. CALLING STATEMENT
!>
!> SUBROUTINE SFLX
!>
!> 2. INPUT (DENOTED BY "I" IN COLUMN SIX OF ARGUMENT LIST AT TOP OF ROUTINE)
!>
!> GENERAL PARAMETERS
!>
!> ICE    - SEA-ICE FLAG  (=1: SEA-ICE, =0: LAND)
!> DT     - TIMESTEP (SEC) (DT SHOULD NOT EXCEED 3600 SECS, RECOMMEND 1800 SECS OR LESS)
!> Z      - HEIGHT (M) ABOVE GROUND OF ATMOSPHERIC FORCING VARIABLES
!> NSOIL  - NUMBER OF SOIL LAYERS (AT LEAST 2, AND NOT GREATER THAN PARAMETER NSOLD SET BELOW)
!> SLDPTH - THE THICKNESS OF EACH SOIL LAYER (M)
!>
!> ATMOSPHERIC VARIABLES 
!>
!> LWDN   - LW DOWNWARD RADIATION (W M-2; POSITIVE, NOT NET LONGWAVE)
!> SOLDN  - SOLAR DOWNWARD RADIATION (W M-2; POSITIVE, NOT NET SHORTWAVE)
!> SFCPRS - PRESSURE AT HEIGHT Z ABOVE GROUND (PASCALS)
!> PRCP   - PRECIP RATE (KG M-2 S-1) (NOTE, THIS IS A RATE)
!> SFCTMP - AIR TEMPERATURE (K) AT HEIGHT Z ABOVE GROUND
!> TH2    - AIR POTENTIAL TEMPERATURE (K) AT HEIGHT Z ABOVE GROUND
!> Q2     - MIXING RATIO AT HEIGHT Z ABOVE GROUND (KG KG-1)
!> SFCSPD - WIND SPEED (M S-1) AT HEIGHT Z ABOVE GROUND
!> Q2SAT  - SAT MIXING RATIO AT HEIGHT Z ABOVE GROUND (KG KG-1)
!> DQSDT2 - SLOPE OF SAT SPECIFIC HUMIDITY CURVE AT T=SFCTMP (KG KG-1 K-1)
!>
!> CANOPY/SOIL CHARACTERISTICS 
!>
!> VEGTYP   - VEGETATION TYPE (INTEGER INDEX)
!> SOILTYP  - SOIL TYPE (INTEGER INDEX)
!> SLOPETYP - CLASS OF SFC SLOPE (INTEGER INDEX)
!> SHDFAC   - AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION (RANGE 0.0-1.0)
!> PTU      - PHOTO THERMAL UNIT (PLANT PHENOLOGY FOR ANNUALS/CROPS) (NOT YET USED, BUT PASSED TO 
!>            REDPRM FOR FUTURE USE IN VEG PARMS)
!> TBOT     - BOTTOM SOIL TEMPERATURE (LOCAL YEARLY-MEAN SFC AIR TEMPERATURE)
!> ALB      - BACKROUND SNOW-FREE SURFACE ALBEDO (FRACTION)
!> SNOALB   - ALBEDO UPPER BOUND OVER DEEP SNOW (FRACTION)
!>
!> 3. STATE VARIABLES: BOTH INPUT AND OUTPUT (NOTE: OUTPUT USUALLY MODIFIED FROM INPUT BY PHYSICS)
!>    (DENOTED BY "2" IN COLUMN SIX OF ARGUMENT LIST AT TOP OF ROUTINE)
!> 
!> STATE VARIABLES 
!>
!> CMC - CANOPY MOISTURE CONTENT (M)
!> T1  - GROUND/CANOPY/SNOWPACK) EFFECTIVE SKIN TEMPERATURE (K)
!>
!> STC(NSOIL)  - SOIL TEMP (K)
!> SMC(NSOIL)  - TOTAL SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION)
!> SH2O(NSOIL) - UNFROZEN SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION)
!> NOTE        -  FROZEN SOIL MOISTURE = SMC - SH2O
!>
!> SNOWH  - SNOW DEPTH (M)
!> SNEQV  - WATER-EQUIVALENT SNOW DEPTH (M)
!> NOTE   - SNOW DENSITY = SNEQV/SNOWH
!> ALBEDO - SURFACE ALBEDO INCLUDING SNOW EFFECT (UNITLESS FRACTION)
!> CH     - SFC EXCH COEF FOR HEAT AND MOISTURE (M S-1)
!> CM     - SFC EXCH COEF FOR MOMENTUM (M S-1)
!> NOTE   - CH AND CM ARE TECHNICALLY CONDUCTANCES SINCE THEY HAVE BEEN MULTIPLIED BY THE WIND SPEED.
!>
!> 4. OUTPUT (DENOTED BY "O" IN COLUMN SIX OF ARGUMENT LIST AT TOP OF ROUTINE)
!>
!> NOTE - SIGN CONVENTION OF SFC ENERGY FLUXES BELOW IS: NEGATIVE IF SINK OF ENERGY TO SURFACE
!>
!> ETP     - POTENTIAL EVAPORATION (W M-2)
!> ETA     - ACTUAL LATENT HEAT FLUX (W M-2: NEGATIVE, IF UP FROM SURFACE)
!> H       - SENSIBLE HEAT FLUX (W M-2: NEGATIVE, IF UPWARD FROM SURFACE)
!> S       - SOIL HEAT FLUX (W M-2: NEGATIVE, IF DOWNWARD FROM SURFACE)
!> RUNOFF1 - SURFACE RUNOFF (M S-1), NOT INFILTRATING THE SURFACE
!> RUNOFF2 - SUBSURFACE RUNOFF (M S-1), DRAINAGE OUT BOTTOM OF LAST SOIL LYR
!> Q1      - EFFECTIVE MIXING RATIO AT GRND SFC (KG KG-1)
!> NOTE    - Q1 IS NUMERICAL EXPENDIENCY FOR EXPRESSING ETA EQUIVALENTLY IN A BULK AERODYNAMIC FORM
!> SNMAX   - SNOW MELT (M) (WATER EQUIVALENT)
!> SOILW   - AVAILABLE SOIL MOISTURE IN ROOT ZONE (UNITLESS FRACTION BETWEEN SOIL SATURATION AND 
!>           WILTING POINT)
!> SOILM   - TOTAL SOIL COLUMN MOISTURE CONTENT (M) (FROZEN + UNFROZEN)
!>
!> FOR DIAGNOSTIC PURPOSES, RETURN SOME PRIMARY PARAMETERS NEXT (SET IN ROUTINE REDPRM)
!>
!> SMCWLT - WILTING POINT (VOLUMETRIC)
!> SMCDRY - DRY SOIL MOISTURE THRESHOLD WHERE DIRECT EVAP FRM TOP LYR ENDS
!> SMCREF - SOIL MOISTURE THRESHOLD WHERE TRANSPIRATION BEGINS TO STRESS
!> SMCMAX - POROSITY, I.E. SATURATED VALUE OF SOIL MOISTURE
!>
!> USE MODULES: F77KINDS
!>              RITE
!>
!> DRIVER     : SURFCE
!> 
!> CALLS      : CANRES
!>              NOPAC
!>              PENMAN
!>              REDPRM
!>              SNOPAC
!>              SNOW_NEW
!>              TDFCND
!>--------------------------------------------------------------------------------------------------    
    USE F77KINDS
    USE RITE
    USE Z0DATA
! 
    IMPLICIT NONE
!    
    INTEGER(KIND=I4KIND), PARAMETER :: NSOLD = 20
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & SNOWNG  , FRZGRA  , SATURATED
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & ICE     , NSOIL   , VEGTYP  , SOILTYP , SLOPETYP
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & K       , KZ      , NROOT  
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & ALB         
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ALBEDO  
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & B       , CFACTR
!--------------------------------------------- 
! CH IS SFC EXCHANGE COEF FOR HEAT/MOIST
! CM IS SFC MOMENTUM DRAG (NOT NEEDED IN SFLX)
!--------------------------------------------- 
    REAL   (KIND=R4KIND), DIMENSION(NSOLD)                                                      ::&
    & RTDIS   , ZSOIL
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                                      ::&
    & SLDPTH 
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & SH2O    , SMC     , STC     , Q1
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & H       , S       , RUNOFF1 , RUNOFF2 , SNMAX   , SOILW   , SOILM
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & SNOWH   , SNEQV   , ETP     , ETA           
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & DT      , Z       , LWDN    , SOLDN   , SFCPRS  , PRCP    , SFCTMP  , TH2     , Q2      ,   &
    & SFCSPD  , Q2SAT   , DQSDT2  , SHDFAC  , PTU     , TBOT    , SNOALB  , CMC     , T1      ,   &
    & CH      , CM      , SMCWLT  , SMCDRY  , SMCREF  , SMCMAX 
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & CMCMAX  , ZBOT    , CSNOW   , CSOIL   , DF1     ,                                           &
    & DKSAT   , DWSAT   , DSOIL   , DTOT    , EXPSNO  , EXPSOI  , EPSCA   ,                       &
    & EDIR1   , EC1     , ETT1    , F       , F1      , FXEXP   , FRZX    ,                       &
    & HS      , KDT     , PC      , PRCP1   , XLAI    ,                                           &
    & PSISAT  , QUARTZ  , Z0      , RCH     , REFKDT  , RR      ,                                 &
    & RGL     , RSMAX   , RC      , RCMIN   , RSNOW   , SNDENS  , SNCOND  ,                       &
    & SBETA   ,                                                                                   &
    & SNOFAC  , SN_NEW  , SLOPE   , SNUP    , SALP    ,                                           &
    & SOILWM  , SOILWW  , T1V     ,                                                               &
    & T24     , T2V     , TH2V    , TOPT            
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & CZIL

    REAL   (KIND=R4KIND), PARAMETER :: TFREEZ =  273.15
    REAL   (KIND=R4KIND), PARAMETER :: LVH2O  =    2.501000E+6
    REAL   (KIND=R4KIND), PARAMETER :: R      =  287.04
    REAL   (KIND=R4KIND), PARAMETER :: CP     = 1004.5
!-------------------------------------------------------------------------------------------------- 
! COMMON BLK "RITE" CARRIES DIAGNOSTIC QUANTITIES FOR PRINTOUT, BUT IS NOT INVOLVED IN MODEL 
! PHYSICS AND IS NOT PRESENT IN PARENT MODEL THAT CALLS SFLX
!--------------------------------------------------------------------------------------------------
!--------------- 
! INITIALIZATION
!--------------- 
    RUNOFF1 = 0.0
    RUNOFF2 = 0.0
    RUNOFF3 = 0.0
    SNMAX   = 0.0
!--------------------------------------------------- 
! THE VARIABLE "ICE" IS A FLAG DENOTING SEA-ICE CASE
!--------------------------------------------------- 
    IF (ICE == 1) THEN
!-------------------------------------------------------  
! SEA-ICE LAYERS ARE EQUAL THICKNESS AND SUM TO 3 METERS
!-------------------------------------------------------                 
        DO KZ = 1,NSOIL
            ZSOIL(KZ) = -3. * FLOAT(KZ) / FLOAT(NSOIL)
        END DO
    ELSE
!---------------------------------------------------------------------------------------- 
! CALCULATE DEPTH (NEGATIVE) BELOW GROUND FROM TOP SKIN SFC TO BOTTOM OF EACH SOIL LAYER.
! NOTE: SIGN OF ZSOIL IS NEGATIVE (DENOTING BELOW GROUND)
!---------------------------------------------------------------------------------------- 
        ZSOIL(1) = -SLDPTH(1)
!                
        DO KZ=2,NSOIL
            ZSOIL(KZ) = - SLDPTH(KZ) + ZSOIL(KZ-1)
        END DO
!
    END IF
!--------------------------------------------------------------------------------------------------
! NEXT IS CRUCIAL CALL TO SET THE LAND-SURFACE PARAMETERS, INCLUDING SOIL-TYPE AND VEG-TYPE 
! DEPENDENT PARAMETERS.
!--------------------------------------------------------------------------------------------------
    CALL REDPRM(VEGTYP  , SOILTYP , SLOPETYP, CFACTR  , CMCMAX  , RSMAX   , TOPT    , REFKDT  ,   &
    &           KDT     , SBETA   , SHDFAC  , RCMIN   , RGL     , HS      , ZBOT    , FRZX    ,   &
    &           PSISAT  , SLOPE   , SNUP    , SALP    , B       , DKSAT   , DWSAT   , SMCMAX  ,   &
    &           SMCWLT  , SMCREF  , SMCDRY  , F1      , QUARTZ  , FXEXP   , RTDIS   , SLDPTH  ,   &
    &           ZSOIL   , NROOT   , NSOIL   , Z0      , CZIL    , XLAI    , CSOIL   , PTU)
!--------------------------------------------------------------------------------------------------
! NEXT CALL ROUTINE SFCDIF TO CALCULATE
! THE SFC EXCHANGE COEF (CH) FOR HEAT AND MOISTURE
!
! NOTE  NOTE
!
! COMMENT OUT CALL SFCDIF, IF SFCDIF ALREADY CALLED IN CALLING PROGRAM (SUCH AS IN COUPLED 
! ATMOSPHERIC MODEL)
!
! NOTE: DO NOT CALL SFCDIF UNTIL AFTER ABOVE CALL TO REDPRM, IN CASE ALTERNATIVE VALUES OF 
!       ROUGHNESS LENGTH (Z0) AND ZILINTINKEVICH  COEF (CZIL) ARE SET THERE VIA NAMELIST I/O
!
! NOTE: ROUTINE SFCDIF RETURNS A CH THAT REPRESENTS THE WIND SPD TIMES THE "ORIGINAL" 
!       NONDIMENSIONAL "Ch" TYPICAL IN LITERATURE. HENCE THE CH RETURNED FROM SFCDIF HAS UNITS 
!       OF M/S. THE IMPORTANT COMPANION COEFFICIENT OF CH, CARRIED HERE AS "RCH", IS THE CH FROM
!       SFCDIF TIMES AIR DENSITY AND PARAMETER "CP" "RCH" IS COMPUTED IN "CALL PENMAN". RCH RATHER
!       THAN CH IS THE COEFF USUALLY INVOKED LATER IN EQNS.
!
! NOTE: SFCDIF ALSO RETURNS THE SURFACE EXCHANGE COEFFICIENT FOR MOMENTUM, CM, ALSO KNOWN AS THE 
!       SURFACE DRAGE COEFFICIENT, BUT CM IS NOT USED HERE
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------- 
! CALC VIRTUAL TEMPS AND VIRTUAL POTENTIAL TEMPS NEEDED BY SUBROUTINES SFCDIF AND PENMAN
!--------------------------------------------------------------------------------------- 
    T2V = SFCTMP * (1.0 + 0.61 * Q2)
!
    SNOWNG = .FALSE. 
    FRZGRA = .FALSE. 
!--------------------------------------------------------  
! IF SEA-ICE CASE, ASSIGN DEFAULT WATER-EQUIV SNOW ON TOP
!--------------------------------------------------------  
    IF (ICE == 1) THEN
        SNEQV = 0.01
        SNOWH = 0.05
    END IF
!-------------------------------------------------------------------------------------------------- 
! IF INPUT SNOWPACK IS NONZERO, THEN COMPUTE SNOW DENSITY "SNDENS" AND SNOW THERMAL CONDUCTIVITY 
! "SNCOND" (NOTE THAT CSNOW IS A FUNCTION SUBROUTINE)
!-------------------------------------------------------------------------------------------------- 
    IF (SNEQV == 0.0) THEN
        SNDENS = 0.0
        SNOWH  = 0.0
        SNCOND = 1.0
    ELSE
        SNDENS = SNEQV / SNOWH
        SNCOND = CSNOW(SNDENS)
    END IF
! -------------------------------------------------------------------------------------------------
! DETERMINE IF ITS PRECIPITATING AND WHAT KIND OF PRECIP IT IS.
! IF IT'S PRCPING AND THE AIR TEMP IS COLDER THAN 0 C, IT'S SNOWING.
! IF ITS PRCPING AND THE AIR TEMP IS WARMER THAN 0 C, BUT THE GRND TEMP IS COLDER THAN 0 C, 
! FREEZING RAIN IS PRESUMED TO BE FALLING.
! -------------------------------------------------------------------------------------------------
    IF (PRCP > 0.0) THEN
        IF (SFCTMP <= TFREEZ) THEN
            SNOWNG = .TRUE. 
        ELSE
            IF ( T1 <= TFREEZ ) FRZGRA = .TRUE. 
        END IF
    END IF
!--------------------------------------------------------------------------------------------------
! IF EITHER PRCP FLAG IS SET, DETERMINE NEW SNOWFALL (CONVERTING PRCP RATE FROM KG M-2 S-1 TO A 
! LIQUID EQUIV SNOW DEPTH IN METERS) AND ADD IT TO THE EXISTING SNOWPACK.
! NOTE THAT SINCE ALL PRECIP IS ADDED TO SNOWPACK, NO PRECIP INFILTRATES INTO THE SOIL SO THAT 
! PRCP1 IS SET TO ZERO.
!--------------------------------------------------------------------------------------------------
    IF ((SNOWNG) .OR. (FRZGRA)) THEN
        SN_NEW = PRCP  * DT * 0.001
        SNEQV  = SNEQV + SN_NEW
        PRCP1  = 0.0
!-------------------------------------------------------------------  
! UPDATE SNOW DENSITY BASED ON NEW SNOWFALL, USING OLD AND NEW SNOW.
!------------------------------------------------------------------- 
        CALL SNOW_NEW(SFCTMP, SN_NEW, SNOWH, SNDENS)
!---------------------------------   
! UPDATE SNOW THERMAL CONDUCTIVITY
!--------------------------------- 
        SNCOND = CSNOW(SNDENS)
!
    ELSE
!--------------------------------------------------------------------------------------------------    
! PRECIP IS LIQUID (RAIN), HENCE SAVE IN THE PRECIP VARIABLE THAT LATER CAN WHOLELY OR PARTIALLY 
! INFILTRATE THE SOIL (ALONG WITH ANY CANOPY "DRIP" ADDED TO THIS LATER)
!--------------------------------------------------------------------------------------------------     
        PRCP1 = PRCP
!
    END IF
!--------------------------------------  
! THERMAL CONDUCTIVITY FOR SEA-ICE CASE
!-------------------------------------- 
    IF (ICE == 1) THEN
        DF1 = 2.2
    ELSE
!--------------------------------------------------------------------------------------------------   
! NEXT CALCULATE THE SUBSURFACE HEAT FLUX, WHICH FIRST REQUIRES CALCULATION OF THE THERMAL 
! DIFFUSIVITY. TREATMENT OF THE LATTER FOLLOWS THAT ON PAGES 148-149 FROM HEAT TRANSFER IN COLD 
! CLIMATES, BY V. J. LUNARDINI (PUBLISHED IN 1981 BY VAN NOSTRAND REINHOLD) I.E. TREATMENT OF TWO 
! CONTIGUOUS "PLANE PARALLEL" MEDIUMS (NAMELY HERE THE FIRST SOIL LAYER AND THE SNOWPACK LAYER, IF 
! ANY). THIS DIFFUSIVITY TREATMENT BEHAVES WELL FOR BOTH ZERO AND NONZERO SNOWPACK, INCLUDING THE 
! LIMIT OF VERY THIN SNOWPACK. THIS TREATMENT ALSO ELIMINATES THE NEED TO IMPOSE AN ARBITRARY UPPER
! BOUND ON SUBSURFACE HEAT FLUX WHEN THE SNOWPACK BECOMES EXTREMELY THIN.
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
! FIRST CALCULATE THERMAL DIFFUSIVITY OF TOP SOIL LAYER, USING BOTH THE FROZEN AND LIQUID SOIL 
! MOISTURE, FOLLOWING THE SOIL THERMAL DIFFUSIVITY FUNCTION OF PETERS-LIDARD ET AL. (1998,JAS, 
! VOL 55, 1209-1224), WHICH REQUIRES THE SPECIFYING THE QUARTZ CONTENT OF THE GIVEN SOIL CLASS 
! (SEE ROUTINE REDPRM)
! -------------------------------------------------------------------------------------------------   
        CALL TDFCND(DF1, SMC(1), QUARTZ, SMCMAX, SH2O(1))
!--------------------------------------------------------------------------------------------------
! NEXT ADD SUBSURFACE HEAT FLUX REDUCTION EFFECT FROM THE OVERLYING GREEN CANOPY, ADAPTED FROM
! SECTION 2.1.2 OF PETERS-LIDARD ET AL.(1997, JGR, VOL 102(D4))
!--------------------------------------------------------------------------------------------------   
        DF1 = DF1 * EXP(SBETA * SHDFAC)
    END IF
!--------------------------------------------------------------------------------------------------
! FINALLY "PLANE PARALLEL" SNOWPACK EFFECT FOLLOWING V.J. LINARDINI REFERENCE CITED ABOVE. NOTE 
! THAT DTOT IS COMBINED DEPTH OF SNOWDEPTH AND THICKNESS OF FIRST SOIL LAYER
!--------------------------------------------------------------------------------------------------
    DSOIL = - (0.5 * ZSOIL(1))
!
    IF (SNEQV == 0.) THEN
        S = DF1 * (T1 - STC(1)) / DSOIL
    ELSE
        DTOT   = SNOWH + DSOIL
        EXPSNO = SNOWH / DTOT
        EXPSOI = DSOIL / DTOT
!----------------------------------------------------------- 
! 1. HARMONIC MEAN (SERIES FLOW)
!    DF1 = (SNCOND * DF1) / (EXPSOI * SNCOND + EXPSNO * DF1)
! 2. ARITHMETIC MEAN (PARALLEL FLOW)
!    DF1 = EXPSNO * SNCOND + EXPSOI * DF1
! 3. GEOMETRIC MEAN (INTERMEDIATE BETWEEN
!    HARMONIC AND ARITHMETIC MEAN)
!----------------------------------------------------------- 
        DF1 = (SNCOND ** EXPSNO) * (DF1 ** EXPSOI)
!--------------------------------------------------------------------------------------------------
! CALCULATE SUBSURFACE HEAT FLUX, S, FROM FINAL THERMAL DIFFUSIVITY OF SURFACE MEDIUMS, DF1 ABOVE, 
! AND SKIN TEMPERATURE AND TOP MID-LAYER SOIL TEMPERATURE
!--------------------------------------------------------------------------------------------------
        S = DF1 * (T1 - STC(1)) / DTOT
    END IF
!-----------------------------------  
! UPDATE ALBEDO, EXCEPT OVER SEA-ICE
!----------------------------------- 
    IF (ICE == 0) THEN
!-------------------------------------------------------------------------------------------------- 
! NEXT IS TIME-DEPENDENT SURFACE ALBEDO MODIFICATION DUE TO TIME-DEPENDENT SNOWDEPTH STATE AND 
! TIME-DEPENDENT CANOPY GREENNESS
!--------------------------------------------------------------------------------------------------
        IF (SNEQV == 0.0) THEN
            ALBEDO = ALB
        ELSE
!--------------------------------------------------------------------------------------------------
! SNUP IS VEG-CLASS DEPENDENT SNOWDEPTH THRESHHOLD (SET IN ROUTINE REDPRM) WHERE MAX SNOW ALBEDO 
! EFFECT IS FIRST ATTAINED
!--------------------------------------------------------------------------------------------------
            IF (SNEQV < SNUP) THEN
                RSNOW = SNEQV / SNUP
                SNOFAC = 1.0 - (EXP(- SALP * RSNOW) - RSNOW * EXP(- SALP))
            ELSE
                SNOFAC = 1.0
            END IF
!--------------------------------------------------------------------------------------------------
! SNOALB IS ARGUMENT REPRESENTING MAXIMUM ALBEDO OVER DEEP SNOW, AS PASSED INTO SFLX, AND ADAPTED 
! FROM THE SATELLITE-BASED MAXIMUM SNOW ALBEDO FIELDS PROVIDED BY D. ROBINSON AND G. KUKLA 
! (1985, JCAM, VOL 24, 402-411)
!--------------------------------------------------------------------------------------------------
            ALBEDO = ALB + (1.0 - SHDFAC) * SNOFAC * (SNOALB - ALB)
!
            IF (ALBEDO > SNOALB) ALBEDO = SNOALB
        END IF
!
    ELSE
!-------------------- 
! ALBEDO OVER SEA-ICE
!-------------------- 
        ALBEDO = 0.60
        SNOFAC = 1.0
    END IF
!--------------------------------------------------------------------------------------------------
! CALCULATE TOTAL DOWNWARD RADIATION (SOLAR PLUS LONGWAVE) NEEDED IN PENMAN EP SUBROUTINE THAT 
! FOLLOWS
!--------------------------------------------------------------------------------------------------              
    F = SOLDN * (1.0 - ALBEDO) + LWDN
!--------------------------------------------------------------------------------------------------
! CALL PENMAN SUBROUTINE TO CALCULATE POTENTIAL EVAPORATION (ETP) (AND OTHER PARTIAL PRODUCTS AND 
! SUMS SAVE IN COMMON/RITE FOR LATER CALCULATIONS)
!--------------------------------------------------------------------------------------------------
    CALL PENMAN (SFCTMP  , SFCPRS  , CH      , T2V     , TH2     , PRCP    , F       , T24     ,  &
    &            S       , Q2      , Q2SAT   , ETP     , RCH     , EPSCA   , RR      , SNOWNG  ,  &
    &            FRZGRA  , DQSDT2)
!------------------------------------- 
! FOLLOWING OLD CONSTRAINT IS DISABLED
! IF(SATURATED) ETP = 0.0
!------------------------------------- 
!
!--------------------------------------------------------------------------------------------------
! CALL CANRES TO CALCULATE THE CANOPY RESISTANCE AND CONVERT IT INTO PC IF MORE THAN TRACE AMOUNT 
! OF CANOPY GREENNESS FRACTION
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------  
! IF(SHDFAC .GT. 1.E-6) THEN MAKE THIS THRESHOLD CONSISTENT WITH THE ONE IN SMFLX FOR TRANSP 
! AND EC(ANOPY)
!-------------------------------------------------------------------------------------------------- 
    IF (SHDFAC > 0.) THEN
!--------------------------------------------------------------------------------------------------              
! FROZEN GROUND EXTENSION: TOTAL SOIL WATER "SMC" WAS REPLACED BY UNFROZEN SOIL WATER "SH2O" IN 
! CALL TO CANRES BELOW
!--------------------------------------------------------------------------------------------------    
        CALL CANRES(SOLDN   , CH      , SFCTMP  , Q2      , SFCPRS  , SH2O    , ZSOIL   ,         &
    &               NSOIL   , SMCWLT  , SMCREF  , RCMIN   , RC      , PC      , NROOT   ,         &
    &               Q2SAT   , DQSDT2  , TOPT    , RSMAX   , RGL     , HS      , XLAI)
!
    END IF
!------------------------------------------------------------------------------------ 
! NOW DECIDE MAJOR PATHWAY BRANCH TO TAKE DEPENDING ON WHETHER SNOWPACK EXISTS OR NOT
!------------------------------------------------------------------------------------ 
    IF (SNEQV == 0.0) THEN

        CALL NOPAC (ETP     , ETA     , PRCP    , SMC     , SMCMAX  , SMCWLT  , SMCREF  ,         &
    &               SMCDRY  , CMC     , CMCMAX  , NSOIL   , DT      , SHDFAC  , SBETA   ,         &
    &               Q1      , Q2      , T1      , SFCTMP  , T24     , TH2     , F       ,         &
    &               F1      , S       , STC     , EPSCA   , B       , PC      , RCH     ,         &
    &               RR      , CFACTR  , SH2O    , SLOPE   , KDT     , FRZX    , PSISAT  ,         &
    &               ZSOIL   , DKSAT   , DWSAT   , TBOT    , ZBOT    , RUNOFF1 , RUNOFF2 ,         &
    &               RUNOFF3 , EDIR1   , EC1     , ETT1    , NROOT   , ICE     , RTDIS   ,         &
    &               QUARTZ  , FXEXP   , CSOIL)
!
    ELSE
!
        CALL SNOPAC (ETP     , ETA     , PRCP    , PRCP1   , SNOWNG  , SMC     , SMCMAX  ,        &
    &                SMCWLT  , SMCREF  , SMCDRY  , CMC     , CMCMAX  , NSOIL   , DT      ,        &
    &                SBETA   , Q1      , DF1     , Q2      , T1      , SFCTMP  , T24     ,        &
    &                TH2     , F       , F1      , S       , STC     , EPSCA   , SFCPRS  ,        &
    &                B       , PC      , RCH     , RR      , CFACTR  , SNOFAC  , SNEQV   ,        &
    &                SNDENS  , SNOWH   , SH2O    , SLOPE   , KDT     , FRZX    , PSISAT  ,        &
    &                SNUP    , ZSOIL   , DWSAT   , DKSAT   , TBOT    , ZBOT    , SHDFAC  ,        &
    &                RUNOFF1 , RUNOFF2 , RUNOFF3 , EDIR1   , EC1     , ETT1    , NROOT   ,        &
    &                SNMAX   , ICE     , RTDIS   , QUARTZ  , FXEXP   , CSOIL)
!                
    END IF
!-----------------------------------------------------  
! PREPARE SENSIBLE HEAT (H) FOR RETURN TO PARENT MODEL
!----------------------------------------------------- 
    H = - (CH * CP * SFCPRS) / (R * T2V) * (TH2 - T1)             
! -------------------------------------------------------------------------------------------------
! CONVERT UNITS AND/OR SIGN OF TOTAL EVAP (ETA), POTENTIAL EVAP (ETP), SUBSURFACE HEAT FLUX (S), 
! AND RUNOFFS FOR WHAT PARENT MODEL EXPECTS
! -------------------------------------------------------------------------------------------------
!
!-------------------------------------  
! CONVERT ETA FROM KG M-2 S-1 TO W M-2
!------------------------------------- 
    ETA = ETA * LVH2O
    ETP = ETP * LVH2O
!--------------------------------------------- 
! CONVERT THE SIGN OF SOIL HEAT FLUX SO THAT:
! S > 0 : WARM THE SURFACE  (NIGHT TIME)
! S < 0 : COOL THE SURFACE  (DAY TIME)
!---------------------------------------------
    S = -1.0 * S
!--------------------------------------------------------------------------------------------------
! CONVERT RUNOFF3 (INTERNAL LAYER RUNOFF FROM SUPERSAT) FROM M TO M S-1 AND ADD TO SUBSURFACE 
! RUNOFF/DRAINAGE/BASEFLOW
!--------------------------------------------------------------------------------------------------
    RUNOFF3 = RUNOFF3 / DT
    RUNOFF2 = RUNOFF2 + RUNOFF3
!--------------------------------------------------------------------------------------------------
! TOTAL COLUMN SOIL MOISTURE IN METERS (SOILM) AND ROOT-ZONE SOIL MOISTURE AVAILABILITY (FRACTION) 
!RELATIVE TO POROSITY/SATURATION
!--------------------------------------------------------------------------------------------------
    SOILM = - 1.0 * SMC(1) * ZSOIL(1)
!          
    DO K = 2,NSOIL
        SOILM = SOILM + SMC(K) * (ZSOIL(K-1) - ZSOIL(K))
    END DO
!
    SOILWM = - 1.0 * (SMCMAX - SMCWLT) * ZSOIL(1)
    SOILWW = - 1.0 * (SMC(1) - SMCWLT) * ZSOIL(1)
!
    DO K = 2,NROOT
        SOILWM = SOILWM + (SMCMAX - SMCWLT) * (ZSOIL(K-1) - ZSOIL(K))
        SOILWW = SOILWW + (SMC(K) - SMCWLT) * (ZSOIL(K-1) - ZSOIL(K))
    END DO
!
    SOILW = SOILWW / SOILWM
!
    RETURN
!
    END SUBROUTINE SFLX
!
!
!
    SUBROUTINE CANRES(SOLAR   , CH      , SFCTMP  , Q2      , SFCPRS  , SMC     , ZSOIL   ,       &
    &                 NSOIL   , SMCWLT  , SMCREF  , RCMIN   , RC      , PC      , NROOT   ,       &
    &                 Q2SAT   , DQSDT2  , TOPT    , RSMAX   , RGL     , HS      , XLAI)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE CANRES
!>
!> SUBPROGRAM: CANRES - ?????
!> PROGRAMMER: JARVIS
!> ORG: ????? 
!> DATE: 76-??-??
!>
!> ABSTRACT:
!> THIS ROUTINE CALCULATES THE CANOPY RESISTANCE WHICH DEPENDS ON INCOMING SOLAR RADIATION, AIR 
!> TEMPERATURE, ATMOSPHERIC WATER VAPOR  PRESSURE DEFICIT AT THE LOWEST MODEL LEVEL, AND SOIL 
!> MOISTURE (PREFERABLY UNFROZEN SOIL MOISTURE RATHER THAN TOTAL).
!>
!> PROGRAM HISTORY LOG:
!> 76-??-??  JARVIS              - ORIGINATOR
!> 90-??-??  JACQUEMIN & NOILHAN - BLM
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> SOLAR  - INCOMING SOLAR RADIATION
!> CH     - SURFACE EXCHANGE COEFFICIENT FOR HEAT AND MOISTURE
!> SFCTMP - AIR TEMPERATURE AT 1ST LEVEL ABOVE GROUND
!> Q2     - AIR HUMIDITY AT 1ST LEVEL ABOVE GROUND
!> Q2SAT  - SATURATION AIR HUMIDITY AT 1ST LEVEL ABOVE GROUND
!> DQSDT2 - SLOPE OF SATURATION HUMIDITY FUNCTION WRT TEMP
!> SFCPRS - SURFACE PRESSURE
!> SMC    - VOLUMETRIC SOIL MOISTURE
!> ZSOIL  - SOIL DEPTH (NEGATIVE SIGN, AS IT IS BELOW GROUND)
!> NSOIL  - NO. OF SOIL LAYERS
!> NROOT  - NO. OF SOIL LAYERS IN ROOT ZONE (1.LE.NROOT.LE.NSOIL)
!> XLAI   - LEAF AREA INDEX
!> SMCWLT - WILTING POINT
!> SMCREF - REFERENCE SOIL MOISTURE (WHERE SOIL WATER DEFICIT STRESS SETS IN)
!>
!> OUTPUT ARGUMENT LIST:
!> PC     - PLANT COEFFICIENT
!> RC     - CANOPY RESISTANCE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE 
!>
!> USE MODULES: F77KINDS
!>  
!> DRIVER     : SFLX
!>
!> CALLS      : ----- 
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND), PARAMETER :: NSOLD = 20
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & SMC     , ZSOIL
!
    REAL   (KIND=R4KIND), DIMENSION(NSOLD)                                                      ::&
    & PART
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & NSOIL   , NROOT
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & K          
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & SOLAR   , CH      , SFCTMP  , Q2      , SFCPRS  , SMCWLT  , SMCREF  , RCMIN   , Q2SAT   ,   &
    & DQSDT2  , TOPT    , RSMAX   , RGL     , HS      , XLAI    
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & RC      , PC         
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & RCS     , RCT     , RCQ     , RCSOIL  , FF      , P       , QS      , GX      , TAIR4   ,   &
    & ST1     , SLVCP   , RR      , DELTA
!
    REAL   (KIND=R4KIND), PARAMETER :: SIGMA =    5.67E-8
    REAL   (KIND=R4KIND), PARAMETER :: RD    =  287.04
    REAL   (KIND=R4KIND), PARAMETER :: CP    = 1004.5
    REAL   (KIND=R4KIND), PARAMETER :: SLV   =    2.501000E6
!
    RCS    = 0.0
    RCT    = 0.0
    RCQ    = 0.0
    RCSOIL = 0.0
    RC     = 0.0
!---------------------------------------------  
! CONTRIBUTION DUE TO INCOMING SOLAR RADIATION
!--------------------------------------------- 
    FF  = 0.55   * 2.0   * SOLAR  / (RGL * XLAI)
    RCS = (FF    + RCMIN / RSMAX) / (1.0 + FF)
    RCS = MAX(RCS, 0.0001)
!---------------------------------------------------------------------- 
! CONTRIBUTION DUE TO AIR TEMPERATURE AT FIRST MODEL LEVEL ABOVE GROUND
!---------------------------------------------------------------------- 
    RCT = 1.0    - 0.0016 * ((TOPT - SFCTMP) ** 2.0)
    RCT = MAX(RCT, 0.0001)
!----------------------------------------------------------------- 
! CONTRIBUTION DUE TO VAPOR PRESSURE DEFICIT AT FIRST MODEL LEVEL.
!-----------------------------------------------------------------  
    QS = Q2SAT
!------------------------- 
! RCQ EXPRESSION FROM SSIB
!------------------------- 
    RCQ = 1.0    /(1.0 + HS * (QS - Q2))
    RCQ = MAX(RCQ, 0.01)
!---------------------------------------------------------------  
! CONTRIBUTION DUE TO SOIL MOISTURE AVAILABILITY.
! DETERMINE CONTRIBUTION FROM EACH SOIL LAYER, THEN ADD THEM UP.
!---------------------------------------------------------------  
    GX = (SMC(1) - SMCWLT) / (SMCREF - SMCWLT)
!
    IF (GX > 1.) GX = 1.0
    IF (GX < 0.) GX = 0.0
!------------------------------------- 
! USING SOIL DEPTH AS WEIGHTING FACTOR
!------------------------------------- 
    PART(1) = (ZSOIL(1) / ZSOIL(NROOT)) * GX
!-------------------------------------------- 
! USING ROOT DISTRIBUTION AS WEIGHTING FACTOR
!-------------------------------------------- 
    DO K=2,NROOT
        GX = (SMC(K) - SMCWLT) / (SMCREF - SMCWLT)
        IF (GX > 1.) GX = 1.0
        IF (GX < 0.) GX = 0.0
!-------------------------------------  
! USING SOIL DEPTH AS WEIGHTING FACTOR
!------------------------------------- 
        PART(K) = ((ZSOIL(K) - ZSOIL(K-1)) / ZSOIL(NROOT)) * GX
!-------------------------------------------- 
! USING ROOT DISTRIBUTION AS WEIGHTING FACTOR
!--------------------------------------------           
    END DO
!
    DO K=1,NROOT
        RCSOIL = RCSOIL + PART(K)
    END DO
!
    RCSOIL = MAX(RCSOIL, 0.0001)
!----------------------------------------------------------  
! DETERMINE CANOPY RESISTANCE DUE TO ALL FACTORS.
! CONVERT CANOPY RESISTANCE (RC) TO PLANT COEFFICIENT (PC).
!----------------------------------------------------------  
    RC    = RCMIN / (XLAI * RCS * RCT * RCQ * RCSOIL)
!               
    TAIR4 = SFCTMP ** 4.
    ST1   = (4.    * SIGMA  * RD)     / CP
    SLVCP = SLV    / CP
    RR    = ST1    * TAIR4  / (SFCPRS * CH) + 1.0
    DELTA = SLVCP  * DQSDT2
!          
    PC    = (RR    + DELTA) / (RR     * (1. + RC * CH) + DELTA)
!          
    RETURN
!
    END SUBROUTINE CANRES
!
!
!
    FUNCTION CSNOW(DSNOW)
!--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & DSNOW
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & C       , CSNOW
!
    REAL   (KIND=R4KIND), PARAMETER :: UNIT = 0.11631
!-------------------------------------------------------------------------------                                            
! SIMULATION OF TERMAL SNOW CONDUCTIVITY
! SIMULATION UNITS OF CSNOW IS CAL/(CM*HR* C) AND IT WILL BE RETURND IN W/(M* C)
! BASIC VERSION IS DYACHKOVA EQUATION
!------------------------------------------------------------------------------- 
!
!--------------------------------------------- 
! DYACHKOVA EQUATION (1960), FOR RANGE 0.1-0.4
!--------------------------------------------- 
    C     = 0.328 * 10 ** (2.25 * DSNOW)
    CSNOW = UNIT  * C
!------------------------------------------ 
! DE VAUX EQUATION (1933), IN RANGE 0.1-0.6
! CSNOW = 0.0293 * (1. + 100. * DSNOW ** 2)          
! E. ANDERSEN FROM FLERCHINGER
! CSNOW = 0.021 + 2.51 * DSNOW ** 2
!------------------------------------------        
    RETURN
    END FUNCTION CSNOW
!
!
!
    FUNCTION DEVAP (ETP1    , SMC     , ZSOIL   , SHDFAC  , SMCMAX  , B       , DKSAT   ,         &
    &               DWSAT   , SMCDRY  , SMCREF  , SMCWLT  , FXEXP)
!--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!--------------------------------------------------------
! NAME: DIRECT EVAPORATION (DEVAP) FUNCTION  VERSION: N/A
!--------------------------------------------------------
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & ETP1    , SMC     , ZSOIL   , SHDFAC  , SMCMAX  , B       , DKSAT   , DWSAT   , SMCDRY ,    &
    & SMCREF  , SMCWLT  , FXEXP   
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & DEVAP   ,   FX 
!
    FX = (ABS((SMC - SMCDRY) / (SMCMAX - SMCDRY))) ** FXEXP
!--------------------------------- 
! FX > 1 REPRESENTS DEMAND CONTROL
! FX < 1 REPRESENTS FLUX CONTROL
!---------------------------------
    FX = MAX(MIN(FX, 1.0), 0.0)
!--------------------------------------------------- 
! ALLOW FOR THE DIRECT-EVAP-REDUCING EFFECT OF SHADE
!--------------------------------------------------- 
    DEVAP = FX * (1.0 - SHDFAC) * ETP1
!
    RETURN
    END FUNCTION DEVAP
!
!
!
    FUNCTION FRH2O(TKELV, SMC, SH2O, SMCMAX, B, PSIS)
!>--------------------------------------------------------------------------------------------------
!> FUNCTION FRH2O
!>
!> FUNCTION: FRH2O - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> CALCULATE AMOUNT OF SUPERCOOLED LIQUID SOIL WATER CONTENT IF TEMPERATURE IS BELOW 273.15K (T0).
!> REQUIRES NEWTON-TYPE ITERATION TO SOLVE THE NONLINEAR IMPLICIT EQUATION GIVEN IN EQN 17 OF KOREN
!> ET AL. (1999, JGR, VOL 104(D16), 19569-19585).
!>
!> NEW VERSION (FEB 2001): MUCH FASTER AND MORE ACCURATE NEWTON ITERATION ACHIEVED BY FIRST TAKING 
!> LOG OF EQN CITED ABOVE LESS THAN 4 (TYPICALLY 1 OR 2) ITERATIONS ACHIEVES CONVERGENCE. ALSO, 
!> EXPLICIT 1-STEP SOLUTION OPTION FOR SPECIAL CASE OF PARAMETER CK=0, WHICH REDUCES THE ORIGINAL 
!> IMPLICIT EQUATION TO A SIMPLER EXPLICIT FORM, KNOWN AS THE "FLERCHINGER EQN". 
!> IMPROVED HANDLING OF SOLUTION IN THE LIMIT OF FREEZING POINT TEMPERATURE T0.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> TKELV  - TEMPERATURE (KELVIN)
!> SMC    - TOTAL SOIL MOISTURE CONTENT (VOLUMETRIC)
!> SH2O   - LIQUID SOIL MOISTURE CONTENT (VOLUMETRIC)
!> SMCMAX - SATURATION SOIL MOISTURE CONTENT (FROM REDPRM)
!> B      - SOIL TYPE "B" PARAMETER (FROM REDPRM)
!> PSIS   - SATURATED SOIL MATRIC POTENTIAL (FROM REDPRM)
!>
!> OUTPUT ARGUMENT LIST:
!> FRH2O  - SUPERCOOLED LIQUID WATER CONTENT.
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & TKELV   , SMC     , SH2O    , SMCMAX  , B       , PSIS    
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & BX      , DENOM   , DF      , DSWL    , FK      , FRH2O   , SWL     , SWLK    
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NLOG    , KCOUNT
!
    REAL   (KIND=R4KIND), PARAMETER :: CK    =    8.0
    REAL   (KIND=R4KIND), PARAMETER :: BLIM  =    5.5
    REAL   (KIND=R4KIND), PARAMETER :: ERROR =    0.005
    REAL   (KIND=R4KIND), PARAMETER :: HLICE =    3.335E5
    REAL   (KIND=R4KIND), PARAMETER :: GS    =    9.81
    REAL   (KIND=R4KIND), PARAMETER :: DICE  =  920.0
    REAL   (KIND=R4KIND), PARAMETER :: DH2O  = 1000.0
    REAL   (KIND=R4KIND), PARAMETER :: T0    =  273.15
!--------------------------------------------------------------------------------------------------
! LIMITS ON PARAMETER B: B < 5.5  (USE PARAMETER BLIM) 
! SIMULATIONS SHOWED IF B > 5.5 UNFROZEN WATER CONTENT IS NON-REALISTICALLY HIGH AT VERY LOW 
! TEMPERATURES 
!--------------------------------------------------------------------------------------------------
    BX = B
    IF (B > BLIM) BX = BLIM
!------------------------------------------------------------- 
! INITIALIZING ITERATIONS COUNTER AND ITERATIVE SOLUTION FLAG.
!------------------------------------------------------------- 
    NLOG   = 0.0
    KCOUNT = 0.0
!------------------------------------------------------------------ 
!  IF TEMPERATURE NOT SIGNIFICANTLY BELOW FREEZING (T0), SH2O = SMC
!------------------------------------------------------------------ 
    IF (TKELV > (T0 - 1.E-3)) THEN
!
        FRH2O = SMC
!
    ELSE
!
        IF (CK /= 0.0) THEN
!----------------------------------------------------------------------------- 
! OPTION 1: ITERATED SOLUTION FOR NONZERO CK IN KOREN ET AL, JGR, 1999, EQN 17
!----------------------------------------------------------------------------- 
!
!---------------------------------------     
! INITIAL GUESS FOR SWL (FROZEN CONTENT)
!--------------------------------------- 
            SWL = SMC - SH2O
!--------------------   
! KEEP WITHIN BOUNDS.
!--------------------   
            IF (SWL > (SMC   - 0.02)) THEN
                SWL = MAX(0.0, SMC - 0.02)
            END IF
!--------------------   
! START OF ITERATIONS
!--------------------  
            DO WHILE (NLOG < 10 .AND. KCOUNT == 0)
                NLOG  = NLOG       + 1
                DF    = ALOG((PSIS * GS / HLICE) * ((1.  + CK     * SWL) ** 2.)                   &
    &                 * (SMCMAX / (SMC- SWL)) ** BX) - ALOG( - (TKELV - T0)  /TKELV)
!            
                DENOM = 2.         * CK / (1.    + CK * SWL)   + BX   / (SMC  - SWL)
                SWLK  = SWL        - DF / DENOM
                DSWL  = ABS(SWLK   - SWL)
                SWL   = SWLK
!--------------------  
! KEEP WITHIN BOUNDS.
!--------------------  
                IF (SWL > (SMC-0.02)) THEN
                    SWL = MAX(0., SMC - 0.02)
                END IF
!----------------------------------------------------------------- 
! IF MORE THAN 10 ITERATIONS, USE EXPLICIT METHOD (CK = 0 APPROX.)
! WHEN DSWL LESS OR EQ. ERROR, NO MORE ITERATIONS REQUIRED.
!----------------------------------------------------------------- 
                IF (DSWL <= ERROR)  THEN
                    KCOUNT = KCOUNT + 1
                END IF
            END DO
!------------------   
! END OF ITERATIONS
!------------------ 
            FRH2O = SMC - SWL
!-------------       
! END OPTION 1 
!------------- 
        END IF
!
        IF (KCOUNT == 0) THEN
            PRINT*, 'FLERCHINGER USED. ITERATIONS=', NLOG
!-------------------------------------------------------------------------------------------------- 
! OPTION 2: EXPLICIT SOLUTION FOR FLERCHINGER EQ. i.e. CK=0 IN KOREN ET AL., JGR, 1999, EQN 17 
!--------------------------------------------------------------------------------------------------         
            FK = (((HLICE / (GS * (- PSIS))) * ((TKELV - T0) / TKELV)) ** (- 1 / BX)) * SMCMAX
            IF (FK < 0.02) FK = 0.02
            FRH2O = MIN(FK, SMC)
!-------------          
! END OPTION 2 
!-------------   
        END IF
!
    END IF
!
    RETURN
    END FUNCTION FRH2O
!
!
!
    SUBROUTINE HRT (RHSTS   , STC     , SMC     , SMCMAX  , NSOIL   , ZSOIL   , YY      ,         &
    &               ZZ1     , TBOT    , ZBOT    , PSISAT  , SH2O    , DT      , B       ,         &
    &               F1      , DF1     , QUARTZ  , CSOIL)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE HRT
!>
!> SUBROUTINE: HRT - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> TO CALCULATE THE RIGHT HAND SIDE OF THE TIME TENDENCY TERM OF THE SOIL THERMAL DIFFUSION
!> EQUATION. ALSO TO COMPUTE (PREPARE) THE MATRIX COEFFICIENTS FOR THE TRI-DIAGONAL MATRIX 
!> OF THE IMPLICIT TIME SCHEME.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> STC    -
!> SMC    -
!> SMCMAX -
!> NSOIL  -
!> ZSOIL  -
!> YY     -
!> ZZ1    -
!> TBOT   -
!> ZBOT   -
!> PSISAT -
!> SH2O   -
!> DT     -
!> B      -
!> F1     -
!> DF1    -
!> QUARTZ -
!> CSOIL  -
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> RHTS - 
!>
!> USE MODULES: ABCI
!>              F77KINDS
!>
!> DRIVER     : SFLX
!>
!> CALLS      : TDFCND
!>              TBND
!>-------------------------------------------------------------------------------------------------
    USE ABCI
    USE F77KINDS
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & NSOIL
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , K       
!--------------------------------- 
! DECLARE SPECIFIC HEAT CAPACITIES
!--------------------------------- 
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & SMCMAX  , YY      , ZZ1     , TBOT    , ZBOT    , PSISAT  , DT      , B       , F1      ,   &
    & DF1     , QUARTZ  , CSOIL   
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & DDZ     , DDZ2    , DENOM   , DF1N    ,                                                     &
    & DF1K    , DTSDZ   , DTSDZ2  , HCPCT   ,  QTOT    , S       ,                                &
    & TSURF   ,                                                                                   &
    & SICE    , TBK     , TSNSR   , TBK1    , SNKSRC
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(INOUT)       ::&
    & RHSTS 
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & STC     , SMC     , ZSOIL   , SH2O        
!
    REAL   (KIND=R4KIND), PARAMETER :: T0 = 273.15
!-------------------------------------------------------------- 
! SET SPECIFIC HEAT CAPACITIES OF AIR, WATER, ICE, SOIL MINERAL
!-------------------------------------------------------------- 
    REAL   (KIND=R4KIND), PARAMETER :: CAIR = 1004.0
    REAL   (KIND=R4KIND), PARAMETER :: CH2O =    4.2E6
    REAL   (KIND=R4KIND), PARAMETER :: CICE =    2.106E6
!---------------------------------------------------- 
! NOTE: CSOIL NOW SET IN ROUTINE REDPRM AND PASSED IN
!---------------------------------------------------- 
!
!---------------------------------  
! BEGIN SECTION FOR TOP SOIL LAYER 
!--------------------------------- 
!
!--------------------------------------------- 
! CALC THE HEAT CAPACITY OF THE TOP SOIL LAYER
!---------------------------------------------  
    HCPCT = SH2O(1)  * CH2O + (1.0 - SMCMAX) * CSOIL + (SMCMAX - SMC(1)) * CAIR + (SMC(1)         &
    &     - SH2O(1)) * CICE
!--------------------------------------------------------------  
! CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER
!-------------------------------------------------------------- 
    DDZ   = 1.0     / (- 0.5 *             ZSOIL(2))
    AI(1) = 0.0
    CI(1) = (DF1    * DDZ)   / (           ZSOIL(1) * HCPCT)
    BI(1) = - CI(1) + DF1    / (0.5      * ZSOIL(1) * ZSOIL(1) * HCPCT * ZZ1)
!--------------------------------------------------------------------------------------------------
! CALC THE VERTICAL SOIL TEMP GRADIENT BTWN THE 1ST AND 2ND SOIL LAYERS.
! THEN CALCULATE THE SUBSURFACE HEAT FLUX. USE THE TEMP GRADIENT AND SUBSFC HEAT FLUX TO CALC 
! "RIGHT-HAND SIDE TENDENCY TERMS", OR "RHSTS", FOR TOP SOIL LAYER.
!--------------------------------------------------------------------------------------------------
    DTSDZ    = (STC(1) -  STC(2)) / (-  0.5       * ZSOIL(2))
    S        =  DF1    * (STC(1)  - YY)    / (0.5 * ZSOIL(1) * ZZ1)
    RHSTS(1) = (DF1    * DTSDZ    - S)     /       (ZSOIL(1) * HCPCT)
!-------------------------------------------------------------------------------------------------- 
! NEXT, SET TEMP "TSURF" AT TOP OF SOIL COLUMN (FOR USE IN FREEZING SOIL PHYSICS LATER IN FUNCTION
! SUBROUTINE SNKSRC). 
! IF SNOWPACK CONTENT IS ZERO, THEN EXPRESSION BELOW GIVES TSURF = SKIN TEMP.
! IF SNOWPACK IS NONZERO (HENCE ARGUMENT ZZ1=1), THEN EXPRESSION BELOW YIELDS SOIL COLUMN TOP 
! TEMPERATURE UNDER SNOWPACK.
!-------------------------------------------------------------------------------------------------- 
    TSURF = (YY + (ZZ1 - 1) * STC(1)) / ZZ1
!-------------------------------------------------------------------------------------------------- 
! NEXT CAPTURE THE VERTICAL DIFFERENCE OF THE HEAT FLUX AT TOP AND BOTTOM OF FIRST SOIL LAYER FOR
! USE IN HEAT FLUX CONSTRAINT APPLIED TO POTENTIAL SOIL FREEZING/THAWING IN ROUTINE SNKSRC
!-------------------------------------------------------------------------------------------------- 
    QTOT = S - DF1 * DTSDZ
!--------------------------------------------------------------------------------------------------
! CALCULATE TEMPERATURE AT BOTTOM INTERFACE OF 1ST SOIL LAYER FOR USE LATER IN FCN SUBROUTINE 
! SNKSRC
!--------------------------------------------------------------------------------------------------
    CALL TBND(STC(1), STC(2), ZSOIL, ZBOT, 1, NSOIL, TBK)
!-------------------------------------------------- 
! CALCULATE FROZEN WATER CONTENT IN 1ST SOIL LAYER.
!-------------------------------------------------- 
    SICE = SMC(1) - SH2O(1)
!--------------------------------------------------------------------------------------------------
! IF FROZEN WATER PRESENT OR ANY OF LAYER-1 MID-POINT OR BOUNDING INTERFACE TEMPERATURES BELOW 
! FREEZING, THEN CALL SNKSRC TOCOMPUTE HEAT SOURCE/SINK (AND CHANGE IN FROZEN WATER CONTENT) DUE
! TO POSSIBLE SOIL WATER PHASE CHANGE.
!--------------------------------------------------------------------------------------------------
    IF ((SICE > 0.) .OR. (TSURF < T0) .OR. (STC(1) < T0) .OR. (TBK < T0)) THEN
!         
        TSNSR    = SNKSRC(TSURF, STC(1), TBK, SMC(1), SH2O(1), ZSOIL, NSOIL, SMCMAX, PSISAT, B,   & 
    &                     DT   , 1     , QTOT)
!
        RHSTS(1) = RHSTS(1) - TSNSR / (ZSOIL(1) * HCPCT)
!
    END IF
!-------------------------------------    
! THIS ENDS SECTION FOR TOP SOIL LAYER 
!-------------------------------------                  
!
!---------------- 
! INITIALIZE DDZ2
!---------------- 
    DDZ2 = 0.0
!--------------------------------------------------------------------------------------------------
! LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABOVE PROCESS EXCEPT SUBSFC OR "GROUND" HEAT 
! FLUX NOT REPEATED IN LOWER LAYERS)
!--------------------------------------------------------------------------------------------------
    DF1K = DF1
!
    DO K=2,NSOIL
!------------------------------------ 
! CALC THIS SOIL LAYERS HEAT CAPACITY
!------------------------------------
        HCPCT = SH2O(K)  * CH2O + (1.0 - SMCMAX) * CSOIL + (SMCMAX - SMC(K)) * CAIR + (SMC(K)     &
    &         - SH2O(K)) * CICE
!    
        IF (K /= NSOIL) THEN
!-------------------------------------------------------- 
! THIS SECTION FOR LAYER 2 OR GREATER, BUT NOT LAST LAYER 
!-------------------------------------------------------- 
!
!---------------------------------------------  
! CALCULATE THERMAL DIFFUSIVITY FOR THIS LAYER
!--------------------------------------------- 
            CALL TDFCND(DF1N, SMC(K), QUARTZ, SMCMAX, SH2O(K))
!----------------------------------------------------- 
! CALC THE VERTICAL SOIL TEMP GRADIENT THRU THIS LAYER
!----------------------------------------------------- 
            DENOM  = 0.5     * (ZSOIL(K-1)  - ZSOIL(K+1))
            DTSDZ2 = (STC(K) -    STC(K+1)) / DENOM
!----------------------------------------------------------- 
! CALC THE MATRIX COEF, CI, AFTER CALCNG ITS PARTIAL PRODUCT
!----------------------------------------------------------- 
            DDZ2  = 2.    / (ZSOIL(K-1) -   ZSOIL(K+1))
            CI(K) = -DF1N * DDZ2        / ((ZSOIL(K-1) - ZSOIL(K)) * HCPCT)
!---------------------------------- 
! CALCULATE TEMP AT BOTTOM OF LAYER
!---------------------------------- 
            CALL TBND (STC(K), STC(K+1), ZSOIL, ZBOT, K, NSOIL, TBK1)
!
        ELSE
!---------------------------------- 
! SPECIAL CASE OF BOTTOM SOIL LAYER 
!---------------------------------- 
!
!--------------------------------------------- 
! CALCULATE THERMAL DIFFUSIVITY FOR THIS LAYER
!--------------------------------------------- 
            CALL TDFCND(DF1N, SMC(K), QUARTZ, SMCMAX, SH2O(K))
!----------------------------------------------------- 
! CALC THE VERTICAL SOIL TEMP GRADIENT THRU THIS LAYER
!----------------------------------------------------- 
            DENOM  = .5      * (ZSOIL(K-1) + ZSOIL(K)) - ZBOT
            DTSDZ2 = (STC(K) - TBOT)       / DENOM
!-------------------------------------------- 
! SET MATRIX COEF, CI TO ZERO IF BOTTOM LAYER
!-------------------------------------------- 
            CI(K) = 0.0 
!--------------------------------------- 
! CALCULATE TEMP AT BOTTOM OF LAST LAYER
!--------------------------------------- 
            CALL TBND(STC(K), TBOT, ZSOIL, ZBOT, K, NSOIL, TBK1)
!
        END IF
!---------------------------------------- 
! THIS ENDS SPECIAL CODE FOR BOTTOM LAYER 
!---------------------------------------- 
!
!--------------------------------------------------------- 
! CALC RHSTS FOR THIS LAYER AFTER CALCNG A PARTIAL PRODUCT
!--------------------------------------------------------- 
        DENOM    = (ZSOIL(K) - ZSOIL(K-1)) * HCPCT
        RHSTS(K) = (DF1N     * DTSDZ2      - DF1K  * DTSDZ) / DENOM
!
        QTOT = - 1.0  * DENOM * RHSTS(K)
!
        SICE = SMC(K) - SH2O(K)
!
        IF ((SICE > 0.) .OR. (TBK < T0) .OR. (STC(K) < T0) .OR. (TBK1 < T0)) THEN
!
            TSNSR    = SNKSRC( TBK, STC(K), TBK1, SMC(K), SH2O(K), ZSOIL, NSOIL, SMCMAX, PSISAT, &
    &                          B  , DT    , K   , QTOT)
!
            RHSTS(K) = RHSTS(K) - TSNSR / DENOM
!
        END IF
!---------------------------------------------- 
! CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER.
!---------------------------------------------- 
        AI(K) = -DF1 * DDZ / ((ZSOIL(K-1) - ZSOIL(K)) * HCPCT)
        BI(K) = -(AI(K) + CI(K))
!------------------------------------------------------------------- 
! RESET VALUES OF DF1, DTSDZ, DDZ, AND TBK FOR LOOP TO NEXT SOIL LYR
!------------------------------------------------------------------- 
        TBK   = TBK1
        DF1K  = DF1N
        DTSDZ = DTSDZ2
        DDZ   = DDZ2
!    
    END DO
!
    RETURN
!
    END SUBROUTINE HRT
!
!
!
    SUBROUTINE HRTICE (RHSTS, STC, NSOIL, ZSOIL, YY, ZZ1, DF1)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE HRTICE
!>
!> SUBROUTINE: HRTICE - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> TO CALCULATE THE RIGHT HAND SIDE OF THE TIME TENDENCY TERM OF THE SOIL THERMAL DIFFUSION EQUATION
!> IN THE CASE OF SEA-ICE PACK. ALSO TO COMPUTE (PREPARE) THE MATRIX COEFFICIENTS FOR THE 
!> TRI-DIAGONAL MATRIX OF THE IMPLICIT TIME SCHEME.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> STC   -
!> NSOIL - 
!> ZSOIL - 
!> YY    -
!> ZZ1   -
!> DF1   -
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> RHSTS - 
!>
!> USE MODULES: ABCI
!>              F77KINDS
!>
!> DRIVER     : SFLX
!>
!> CALLS      : -----
!>-------------------------------------------------------------------------------------------------
    USE ABCI
    USE F77KINDS
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & NSOIL 
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & K
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(INOUT)       ::&
    & RHSTS
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & STC     , ZSOIL
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & YY      , ZZ1     , DF1
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & DDZ     , DDZ2    , DENOM   , DTSDZ   , DTSDZ2  , HCPCT   , S       , TBOT    , ZBOT  
!--------------------------------------------------------------------------------------------------
! THE INPUT ARGUMENT DF1 A UNIVERSALLY CONSTANT VALUE OF SEA-ICE THERMAL DIFFUSIVITY, SET IN 
! ROUTINE SNOPAC AS DF1 = 2.2
!
! SET LOWER BOUNDARY DEPTH AND BOUNDARY TEMPERATURE OF UNFROZEN SEA WATER AT BOTTOM OF SEA ICE 
! PACK. ASSUME ICE PACK IS OF NSOIL LAYERS SPANNING A UNIFORM CONSTANT ICE PACK THICKNESS AS 
! DEFINED IN ROUTINE SFLX
!--------------------------------------------------------------------------------------------------
    ZBOT = ZSOIL(NSOIL)
    TBOT = 271.16
!--------------------------------------------------------------------  
! SET A NOMINAL UNIVERSAL VALUE OF THE SEA-ICE SPECIFIC HEAT CAPACITY
!--------------------------------------------------------------------         
    HCPCT = 1880.0 * 917.0
!-------------------------------------------------------------- 
! CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER
!-------------------------------------------------------------- 
    DDZ   = 1.0 / (-0.5 * ZSOIL(2))
    AI(1) = 0.0
    CI(1) = (DF1*DDZ) / (ZSOIL(1) * HCPCT)
    BI(1) = -CI(1) + DF1 / (0.5 * ZSOIL(1) * ZSOIL(1) * HCPCT * ZZ1)
!--------------------------------------------------------------------------------------------------
! CALC THE VERTICAL SOIL TEMP GRADIENT BTWN THE TOP AND 2ND SOIL LAYERS.
! RECALC/ADJUST THE SOIL HEAT FLUX. USE THE GRADIENT AND FLUX TO CALC RHSTS FOR THE TOP SOIL LAYER.
!--------------------------------------------------------------------------------------------------
    DTSDZ    = (STC(1) - STC(2))  / (- 0.5 * ZSOIL(2))
    S        =  DF1     * (STC(1) - YY)    / (0.5      * ZSOIL(1) * ZZ1)
    RHSTS(1) = (DF1     * DTSDZ   - S)     / (ZSOIL(1) * HCPCT)
!----------------  
! INITIALIZE DDZ2
!---------------- 
    DDZ2 = 0.0
!----------------------------------------------------------------- 
! LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABOVE PROCESS
!----------------------------------------------------------------- 
    DO K=2,NSOIL
!
        IF (K /= NSOIL) THEN
!------------------------------------------------------  
! CALC THE VERTICAL SOIL TEMP GRADIENT THRU THIS LAYER.
!------------------------------------------------------ 
            DENOM  = 0.5     * (ZSOIL(K-1)  - ZSOIL(K+1))
            DTSDZ2 = (STC(K) -    STC(K+1)) / DENOM
!-----------------------------------------------------------  
! CALC THE MATRIX COEF, CI, AFTER CALCNG ITS PARTIAL PRODUCT
!----------------------------------------------------------- 
            DDZ2  = 2.    / (ZSOIL(K-1) -   ZSOIL(K+1))
            CI(K) = - DF1 * DDZ2        / ((ZSOIL(K-1) - ZSOIL(K)) * HCPCT)
        ELSE
!----------------------------------------------------------- 
! CALC THE VERTICAL SOIL TEMP GRADIENT THRU THE LOWEST LAYER
!-----------------------------------------------------------  
            DTSDZ2 = (STC(K) - TBOT) / (.5 * (ZSOIL(K-1) + ZSOIL(K)) - ZBOT)
!---------------------------- 
! SET MATRIX COEF, CI TO ZERO
!----------------------------  
            CI(K) = 0.0
        END IF
!--------------------------------------------------------- 
! CALC RHSTS FOR THIS LAYER AFTER CALCNG A PARTIAL PRODUCT
!--------------------------------------------------------- 
        DENOM    = (ZSOIL(K) - ZSOIL(K-1)) * HCPCT
        RHSTS(K) = (DF1      * DTSDZ2      - DF1 * DTSDZ) / DENOM
!----------------------------------------------  
! CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER.
!---------------------------------------------- 
        AI(K) = - DF1    * DDZ  / ((ZSOIL(K-1) - ZSOIL(K)) * HCPCT)
        BI(K) = - (AI(K) + CI(K))
!-------------------------------------------------------- 
! RESET VALUES OF DTSDZ AND DDZ FOR LOOP TO NEXT SOIL LYR
!-------------------------------------------------------- 
        DTSDZ = DTSDZ2
        DDZ   = DDZ2
!
    END DO
!
    RETURN
!
    END SUBROUTINE HRTICE
!
!
!
    SUBROUTINE HSTEP (STCOUT, STCIN, RHSTS, DT, NSOIL)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE HSTEP
!>
!> SUBROUTINE: HSTEP - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> TO CALCULATE/UPDATE THE SOIL TEMPERATURE FIELD.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> STCIN -
!> DT    - 
!> NSOIL -
!>
!> OUTPUT ARGUMENT LIST:
!> STCOUT - 
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> RHSTS - 
!>
!> USE MODULES: ABCI
!>              F77KINDS
!>
!> DRIVER     : SFLX
!>
!> CALLS      : ROSR12
!>-------------------------------------------------------------------------------------------------
    USE ABCI
    USE F77KINDS
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & NSOIL
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & K
!
    REAL   (KIND=R4KIND), DIMENSION(NSOLD)                                                      ::&
    & CIIN
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & DT
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & STCIN
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(INOUT)       ::&
    & STCOUT
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(INOUT)       ::&
    & RHSTS
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                                      ::&
    & RHSTSIN
!----------------------------------------------------------  
! CREATE FINITE DIFFERENCE VALUES FOR USE IN ROSR12 ROUTINE
!---------------------------------------------------------- 
    DO K=1,NSOIL
        RHSTS(K) =   RHSTS(K) * DT
           AI(K) =      AI(K) * DT
           BI(K) = 1. + BI(K) * DT
           CI(K) =      CI(K) * DT
    END DO
!------------------------------------------------------  
! COPY VALUES FOR INPUT VARIABLES BEFORE CALL TO ROSR12
!------------------------------------------------------ 
    DO K=1,NSOIL
        RHSTSIN(K) = RHSTS(K)
    END DO
!
    DO K=1,NSOLD
        CIIN(K) = CI(K)
    END DO
!---------------------------------------  
! SOLVE THE TRI-DIAGONAL MATRIX EQUATION
!--------------------------------------- 
    CALL ROSR12(                                                                                  &
    & CI      , AI      , BI      , CIIN    , RHSTSIN , RHSTS   , NSOIL)
!------------------------------------------------- 
! CALC/UPDATE THE SOIL TEMPS USING MATRIX SOLUTION
!-------------------------------------------------  
    DO K=1,NSOIL
        STCOUT(K) = STCIN(K) + CI(K)
    END DO
!
    RETURN
!    
    END SUBROUTINE HSTEP
!
!
!
    SUBROUTINE NOPAC(ETP     , ETA     , PRCP    , SMC     , SMCMAX  , SMCWLT  , SMCREF  ,        &
    &                SMCDRY  , CMC     , CMCMAX  , NSOIL   , DT      , SHDFAC  , SBETA   ,        &
    &                Q1      , Q2      , T1      , SFCTMP  , T24     , TH2     , F       ,        &
    &                F1      , S       , STC     , EPSCA   , B       , PC      , RCH     ,        &
    &                RR      , CFACTR  , SH2O    , SLOPE   , KDT     , FRZFACT , PSISAT  ,        &
    &                ZSOIL   , DKSAT   , DWSAT   , TBOT    , ZBOT    , RUNOFF1 , RUNOFF2 ,        &
    &                RUNOFF3 , EDIR1   , EC1     , ETT1    , NROOT   , ICE     , RTDIS   ,        &
    &                QUARTZ  , FXEXP   , CSOIL)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE NOPAC
!>
!> SUBROUTINE: NOPAC - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> TO CALCULATE SOIL MOISTURE AND HEAT FLUX VALUES AND UPDATE SOIL MOISTURE CONTENT AND 
!> SOIL HEAT CONTENT VALUES FOR THE CASE WHEN NO SNOW PACK IS PRESENT.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> ETP     -
!> PRCP    - 
!> NSOIL   -
!> SMC     - 
!> SMCMAX  - 
!> SMCWLT  -
!> SMCREF  -
!> SMCDRY  -
!> CMC     - 
!> CMCMAX  -
!> DT      - 
!> SHDFAC  - 
!> SBETA   - 
!> Q1      - 
!> Q2      -  
!> T1      -  
!> SMCREF  - 
!> T24     - 
!> TH2     -
!> F       - 
!> F1      -
!> S       -
!> STC     - 
!> EPSCA   -
!> B       -
!> PC      -
!> RCH     -
!> RR      - 
!> CFACTR  - 
!> SH2O    -
!> SLOPE   -
!> KDT     - 
!> PSISAT  -
!> ZSOIL   - 
!> DKSAT   -
!> DWSAT   -
!> TBOT    -
!> ZBOT    -
!> RUNOFF1 - 
!> RUNOFF2 -
!> RUNOFF3 - 
!> EDIR1   -
!> EC1     - 
!> ETT1    - 
!> NROOT   -
!> ICE     - 
!> RTDIS   -
!> QUARTZ  -
!> FXEXP   - 
!> CSOIL   -
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> FRZFACT -  
!>
!> USE MODULES: F77KINDS
!>              RITE
!>
!> DRIVER     : SFLX
!>
!> CALLS      : SHFLX
!>              SMFLX
!>              TDFCND 
!>-------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE RITE, RUNOXX3 => RUNOFF3
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ICE     , NROOT   , NSOIL
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ZZ1     , YYNUM   , YY      , PRCP1   , ETP1    , ETA1    , DF1  
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & FRZFACT
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & B       , CFACTR  , CMC     , ETT1    , CMCMAX  , CSOIL   , DKSAT   , DT      ,             &
    & DWSAT   , EPSCA   , F       , F1      , FXEXP   ,                                           &
    & KDT     , PC      , PRCP    , Q2      , RCH     , RR      , S       , SBETA   ,             &
    & SFCTMP  , SHDFAC  , SMCDRY  , QUARTZ  , SMCMAX  , SMCREF  , SMCWLT  , T1      , T24     ,   &
    & TBOT    , ZBOT    , TH2     , Q1      , SLOPE   ,                                           &
    & PSISAT  , RUNOFF1 , RUNOFF2 , RUNOFF3 , EDIR1   , EC1     , ETP   
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & ETA     
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & RTDIS   , SMC     , SH2O    , STC     , ZSOIL
!
    REAL   (KIND=R4KIND), PARAMETER :: CP    = 1004.5
    REAL   (KIND=R4KIND), PARAMETER :: SIGMA =    5.67E-8
!----------------------------- 
! EXECUTABLE CODE BEGINS HERE.
!-----------------------------
!
!--------------------------------------------------------
! CONVERT ETP FROM KG M-2 S-1 TO MS-1 AND INITIALIZE DEW.
!-------------------------------------------------------- 
    PRCP1 = PRCP * 0.001
    ETP1  = ETP  * 0.001
    DEW   = 0.0
!
    IF (ETP > 0.0) THEN
!----------------------------------------- 
! CONVERT PRCP FROM  KG M-2 S-1  TO  M S-1
!-----------------------------------------  
        CALL SMFLX (ETA1    , SMC     , NSOIL   , CMC     , ETP1    , DT      , PRCP1   ,         &
    &               ZSOIL   , SH2O    , SLOPE   , KDT     , FRZFACT , SMCMAX  , B       ,         &
    &               PC      , SMCWLT  , DKSAT   , DWSAT   , SMCREF  , SHDFAC  , CMCMAX  ,         &
    &               SMCDRY  , CFACTR  , RUNOFF1 , RUNOFF2 , RUNOFF3 , EDIR1   , EC1     ,         &
    &               ETT1    , SFCTMP  , Q2      , NROOT   , RTDIS   , FXEXP)
!-------------------------------------------------------------  
! CONVERT MODELED EVAPOTRANSPIRATION FM  M S-1  TO  KG M-2 S-1
!------------------------------------------------------------- 
        ETA = ETA1 * 1000.0
    ELSE
!------------------------------------------------------------------------------------- 
! IF ETP < 0, ASSUME DEW FORMS (TRANSFORM ETP1 INTO DEW AND REINITIALIZE ETP1 TO ZERO)
!-------------------------------------------------------------------------------------  
        DEW  = -ETP1
        ETP1 = 0.0
!----------------------------------------------------------  
! CONVERT PRCP FROM  KG M-2 S-1  TO  M S-1  AND ADD DEW AMT
!----------------------------------------------------------  
        PRCP1 = PRCP1 + DEW
!    
        CALL SMFLX (ETA1    , SMC     , NSOIL   , CMC     , ETP1    , DT      , PRCP1   ,         &
    &               ZSOIL   , SH2O    , SLOPE   , KDT     , FRZFACT , SMCMAX  , B       ,         &
    &               PC      , SMCWLT  , DKSAT   , DWSAT   , SMCREF  , SHDFAC  , CMCMAX  ,         &
    &               SMCDRY  , CFACTR  , RUNOFF1 , RUNOFF2 , RUNOFF3 , EDIR1   , EC1     ,         &
    &               ETT1    , SFCTMP  , Q2      , NROOT   , RTDIS   , FXEXP)
!-------------------------------------------------------------  
! CONVERT MODELED EVAPOTRANSPIRATION FM  M S-1  TO  KG M-2 S-1
!-------------------------------------------------------------  
        ETA  = ETA1 * 1000.0
!
    END IF
!------------------------------------------  
! BASED ON ETP AND E VALUES, DETERMINE BETA
!------------------------------------------  
    IF (ETP <= 0.0) THEN
        BETA = 0.0
        IF (ETP < 0.0) THEN
            BETA = 1.0
            ETA  = ETP
        END IF
    ELSE
        BETA = ETA / ETP
    END IF
!--------------------------------------------------------------------------------------------------
! GET SOIL THERMAL DIFFUXIVITY/CONDUCTIVITY FOR TOP SOIL LYR, CALC. ADJUSTED TOP LYR SOIL TEMP AND 
! ADJUSTED SOIL FLUX, THEN CALL SHFLX TO COMPUTE/UPDATE SOIL HEAT FLUX AND SOIL TEMPS.
!-------------------------------------------------------------------------------------------------- 
    CALL TDFCND (DF1, SMC(1), QUARTZ, SMCMAX, SH2O(1))
!--------------------------------------------------------------------------------------------------
! VEGETATION GREENNESS FRACTION REDUCTION IN SUBSURFACE HEAT FLUX VIA REDUCTION FACTOR, WHICH IS 
! CONVENIENT TO APPLY HERE TO THERMAL DIFFUSIVITY THAT IS LATER USED IN HRT TO COMPUTE SUB SFC HEAT 
! FLUX (SEE ADDITIONAL COMMENTS ON VEG EFFECT SUB-SFC HEAT FLX IN ROUTINE SFLX)
!--------------------------------------------------------------------------------------------------
    DF1 = DF1 * EXP(SBETA * SHDFAC)
!--------------------------------------------------------------------------------------------------
! COMPUTE INTERMEDIATE TERMS PASSED TO ROUTINE HRT (VIA ROUTINE SHFLX BELOW) FOR USE IN COMPUTING 
! SUBSURFACE HEAT FLUX IN HRT
!--------------------------------------------------------------------------------------------------
    YYNUM = F      - SIGMA  * T24
    YY    = SFCTMP + (YYNUM / RCH      + TH2 - SFCTMP - BETA * EPSCA) / RR
    ZZ1   = DF1    / (-0.5  * ZSOIL(1) * RCH * RR)    + 1.0
!
    CALL SHFLX (S       , STC     , SMC     , SMCMAX  , NSOIL   , T1      , DT      , YY      ,   &
    & ZZ1     , ZSOIL   , TBOT    , ZBOT    , SMCWLT  , PSISAT  , SH2O    , B       , F1      ,   &
    & DF1     , ICE     , QUARTZ  , CSOIL)
!-------------------------------------------------------------------------------------------------- 
! SET FLX1, AND FLX3 TO ZERO SINCE THEY ARE NOT USED. FLX2 WAS SIMILARLY INITIALIZED IN THE PENMAN
! ROUTINE.
!--------------------------------------------------------------------------------------------------
    FLX1 = 0.0
    FLX3 = 0.0
!
    RETURN
!
    END SUBROUTINE NOPAC
!
!
!
    SUBROUTINE PENMAN (SFCTMP  , SFCPRS  , CH      , T2V     , TH2     , PRCP    , F       ,      &
    &                  T24     , S       , Q2      , Q2SAT   , ETP     , RCH     , EPSCA   ,      &
    &                  RR      , SNOWNG  , FRZGRA  , DQSDT2)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE PENMAN
!>
!> SUBROUTINE: PENMAN - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> TO CALCULATE POTENTIAL EVAPORATION FOR THE CURRENT POINT. VARIOUS PARTIAL SUMS / PRODUCTS
!> ARE ALSO CALCULATED AND PASSED BACK TO THE CALLING ROUTINE FOR LATER USE.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> SFCTMP - 
!> SFCPRS - 
!> CH     -
!> T2V    -
!> TH2    - 
!> PRCP   - 
!> F      -
!> S      -
!> Q2     - 
!> Q2SAT  -
!> SNOWNG - 
!> FRZGRA -
!> DQSDT2 -
!>
!> OUTPUT ARGUMENT LIST:
!> T24    -
!> ETP    - 
!> RCH    - 
!> EPSCA  - 
!> RR     - 
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> ETA     - 
!>
!> USE MODULES: F77KINDS
!>              RITE
!>
!> DRIVER     : SFLX
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE RITE, RUNOXX3 => RUNOFF3
!
    IMPLICIT NONE
!
    LOGICAL(KIND=L4KIND)                                                  , INTENT(IN)          ::&
    & SNOWNG  , FRZGRA
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & T24     , ETP     , RCH     , EPSCA   , RR

    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & CH      , F       , PRCP    , Q2      , Q2SAT   , S       , SFCPRS  , SFCTMP  , T2V     ,   &
    & TH2     , DQSDT2
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & A       , RHO     , RAD     , FNET    , DELTA     
!
    REAL   (KIND=R4KIND), PARAMETER :: CP    = 1004.6
    REAL   (KIND=R4KIND), PARAMETER :: CPH2O =    4.218E+3
    REAL   (KIND=R4KIND), PARAMETER :: CPICE =    2.106E+3
    REAL   (KIND=R4KIND), PARAMETER :: R     =  287.04 
    REAL   (KIND=R4KIND), PARAMETER :: ELCP  =    2.4888E+3
    REAL   (KIND=R4KIND), PARAMETER :: LSUBF =    3.335E+5
    REAL   (KIND=R4KIND), PARAMETER :: LSUBC =    2.501000E+6
    REAL   (KIND=R4KIND), PARAMETER :: SIGMA =    5.67E-8
!---------------------------- 
! EXECUTABLE CODE BEGINS HERE 
!---------------------------- 
    FLX2 = 0.0
!------------------------------------------------ 
! PREPARE PARTIAL QUANTITIES FOR PENMAN EQUATION.
!------------------------------------------------ 
    DELTA = ELCP   * DQSDT2
    T24   = SFCTMP * SFCTMP  * SFCTMP  * SFCTMP
    RR    = T24    * 6.48E-8 / (SFCPRS * CH)   + 1.0
    RHO   = SFCPRS / (R      * T2V)
    RCH   = RHO    * CP      * CH
!--------------------------------------------------------------------------------------------------
! ADJUST THE PARTIAL SUMS / PRODUCTS WITH THE LATENT HEAT EFFECTS CAUSED BY FALLING PRECIPITATION.
!-------------------------------------------------------------------------------------------------- 
    IF (.NOT. SNOWNG) THEN
        IF (PRCP > 0.0) RR = RR + CPH2O * PRCP / RCH
    ELSE
        RR = RR + CPICE * PRCP / RCH
    END IF
!
    FNET = F - SIGMA * T24 - S
!--------------------------------------------------------------------------------------------------
! INCLUDE THE LATENT HEAT EFFECTS OF FRZNG RAIN CONVERTING TO ICE ON IMPACT IN THE CALCULATION OF 
! FLX2 AND FNET.
!-------------------------------------------------------------------------------------------------- 
    IF (FRZGRA) THEN
        FLX2 = - LSUBF * PRCP
        FNET =   FNET  - FLX2
    END IF
!------------------------------------- 
! FINISH PENMAN EQUATION CALCULATIONS.
!------------------------------------- 
    RAD   = FNET  / RCH    + TH2  - SFCTMP
    A     = ELCP  * (Q2SAT - Q2)
    EPSCA = (A    * RR     + RAD  * DELTA) / (DELTA + RR)
    ETP   = EPSCA * RCH    / LSUBC
!
    RETURN
!
    END SUBROUTINE PENMAN
!
!
!
    SUBROUTINE REDPRM (VEGTYP  , SOILTYP , SLOPETYP, CFACTR  , CMCMAX  , RSMAX   , TOPT    ,      &
    &                  REFKDT  , KDT     , SBETA   , SHDFAC  , RCMIN   , RGL     , HS      ,      &
    &                  ZBOT    , FRZX    , PSISAT  , SLOPE   , SNUP    , SALP    , B       ,      &
    &                  DKSAT   , DWSAT   , SMCMAX  , SMCWLT  , SMCREF  , SMCDRY  , F1      ,      &
    &                  QUARTZ  , FXEXP   , RTDIS   , SLDPTH  , ZSOIL   , NROOT   , NSOIL   ,      &
    &                  Z0      , CZIL    , LAI     , CSOIL   , PTU)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE REDPRM
!>
!> SUBROUTINE: REDPRM - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> THIS SUBROUTINE INTERNALLY SETS (DEFAULTS), OR OPTIONALLY READS-IN VIA NAMELIST I/O, ALL THE SOIL 
!> AND VEGETATION PARAMETERS REQUIRED FOR THE EXECUSION OF THE NOAH - LSM
!>
!> OPTIONAL NON-DEFAULT PARAMETERS CAN BE READ IN, ACCOMMODATING UP TO 30 SOIL, VEG, OR SLOPE 
!> CLASSES, IF THE DEFAULT MAX NUMBER OF SOIL, VEG, AND/OR SLOPE TYPES IS RESET.
!>
!> FUTURE UPGRADES OF ROUTINE REDPRM MUST EXPAND TO INCORPORATE SOME OF THE EMPIRICAL PARAMETERS OF 
!> THE FROZEN SOIL AND SNOWPACK PHYSICS (SUCH AS IN ROUTINES FRH2O, SNOWPACK, AND SNOW_NEW) NOT YET 
!> SET IN THIS REDPRM ROUTINE, BUT RATHER SET IN LOWER LEVEL SUBROUTINES 
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> VEGTYP   -
!> SOILTYP  -
!> SLOPETYP -
!> SLDPTH   -
!> ZSOIL    -
!> NSOIL    -
!> PTU      -
!>
!> OUTPUT ARGUMENT LIST:
!> CFACTR   - 
!> CMCMAX   -
!> RSMAX    -
!> TOPT     - 
!> REFKDT   -
!> KDT      -
!> SBETA    -
!> SHDFAC   -
!> RCMIN    - 
!> RGL      - 
!> HS       - 
!> ZBOT     -
!> FRZX     - 
!> PSISAT   -
!> SLOPE    - 
!> SNUP     -
!> SALP     - 
!> B        -
!> SALP     - 
!> DWSAT    -
!> SMCMAX   -
!> SMCWLT   - 
!> SMCREF   -
!> SMCDRY   -
!> F1       -
!> QUARTZ   -
!> FXEXP    - 
!> RTDIS    -
!> Z0       - 
!> CZIL     - 
!> LAI      -
!> CSOIL    -
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> 
!>
!> USE MODULES: F77KINDS
!>              NSOILTYPE
!>
!> DRIVER     : SFLX
!>
!> CALLS      : -----
!>-------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE NSOILTYPE
    USE Z0DATA  
!
    IMPLICIT NONE
!------------------------------------------------------------------ 
! SET MAXIMUM NUMBER OF SOIL-, VEG-, AND SLOPETYP IN DATA STATEMENT
!------------------------------------------------------------------ 
    INTEGER(KIND=I4KIND), PARAMETER :: MAX_SOILTYP  = 30
!    INTEGER(KIND=I4KIND), PARAMETER :: MAX_VEGTYP   = 30
    INTEGER(KIND=I4KIND), PARAMETER :: MAX_SLOPETYP = 30
!--------------------------------------------------- 
!  NUMBER OF DEFINED SOIL-, VEG-, AND SLOPETYPS USED
!--------------------------------------------------- 
    INTEGER(KIND=I4KIND)                                                                        ::&
    & DEFINED_VEG       , DEFINED_SOIL      , DEFINED_SLOPE
!------------------------- 
! ANIEL:DEFINED SOIL DOYLE
!------------------------- 
    INTEGER(KIND=I4KIND)                                                                        ::&
    & DEFINED_SOILD
!---------------------------  
! ANIEL DEFINED SOIL ZOGBLER
!--------------------------- 
    INTEGER(KIND=I4KIND)                                                                        ::&
    & DEFINED_SOILZ
!---- 
! ROL
!---- 
    INTEGER(KIND=I4KIND)                                                                        ::&
    & DEFINED_SOILG
!
!IP    DATA DEFINED_VEG / 14 /
    DATA DEFINED_VEG / 16 /
!------  
! ANIEL
!------  
    DATA DEFINED_SOILD / 18 /
    DATA DEFINED_SOILZ /  9 /
!---- 
! ROL
!---- 
    DATA DEFINED_SOILG / 15 /
!------ 
! ANIEL
!------ 
    DATA DEFINED_SLOPE / 9 /
!--------------------------------------------------------------------------------------------------
! SET-UP SOIL PARAMETERS FOR GIVEN SOIL TYPE
! INPUT : SOLTYP: SOIL TYPE (INTEGER INDEX)
! OUTPUT: SOIL PARAMETERS
!
!    MAXSMC - MAX SOIL MOISTURE CONTENT (POROSITY)
!    REFSMC - REFERENCE SOIL MOISTURE (ONSET OF SOIL MOISTURE STRESS IN TRANSPIRATION)
!    WLTSMC - WILTING PT SOIL MOISTURE CONTENTS
!    DRYSMC - AIR DRY SOIL MOIST CONTENT LIMITS
!    SATPSI - SATURATED SOIL POTENTIAL
!    SATDK  - SATURATED SOIL HYDRAULIC CONDUCTIVITY
!    BB     - THE 'B' PARAMETER
!    SATDW  - SATURATED SOIL DIFFUSIVITY
!    F11    - USED TO COMPUTE SOIL DIFFUSIVITY/CONDUCTIVITY
!    QUARTZ - SOIL QUARTZ CONTENT
!
! SOIL TYPES   ZOBLER (1986)     COSBY ET AL (1984) (QUARTZ CONT.(1))
!     1        COARSE            LOAMY SAND         (0.82)
!     2        MEDIUM            SILTY CLAY LOAM    (0.10)
!     3        FINE              LIGHT CLAY         (0.25)
!     4        COARSE-MEDIUM     SANDY LOAM         (0.60)
!     5        COARSE-FINE       SANDY CLAY         (0.52)
!     6        MEDIUM-FINE       CLAY LOAM          (0.35)
!     7        COARSE-MED-FINE   SANDY CLAY LOAM    (0.60)
!     8        ORGANIC           LOAM               (0.40)
!     9        GLACIAL LAND ICE  LOAMY SAND         (NA USING 0.82)
!--------------------------------------------------------------------------------------------------
    REAL   (KIND=R4KIND), DIMENSION(MAX_SOILTYP)                                                ::&
    & BB      , DRYSMC  , F11     , MAXSMC  , REFSMC  , SATPSI  , SATDK   , SATDW   , WLTSMC  ,   &
    & QTZ   
!------------------------  
! ANIEL: PARAMETROS DOYLE
!------------------------  
    REAL   (KIND=R4KIND), DIMENSION(MAX_SOILTYP)                                                ::&
    & BBD     , DRYSMCD , F11D    , MAXSMCD , REFSMCD , SATPSID , SATDKD  , SATDWD  , WLTSMCD ,   &
    & QTZD   
!--------------------------  
! ANIEL: PARAMETROS ZOGBLER
!-------------------------- 
    REAL   (KIND=R4KIND), DIMENSION(MAX_SOILTYP)                                                ::&
    & BBZ     , DRYSMCZ , F11Z    , MAXSMCZ , REFSMCZ , SATPSIZ , SATDKZ  , SATDWZ  , WLTSMCZ ,   &
    & QTZZ   
!---------------------------------------- 
! ROL - PARAMETROS DO MAPA DE SOLO GLOBAL
!---------------------------------------- 
    REAL   (KIND=R4KIND), DIMENSION(MAX_SOILTYP)                                                ::&
    & BBG     , DRYSMCG , F11G    , MAXSMCG , REFSMCG , SATPSIG , SATDKG  , SATDWG  , WLTSMCG ,   &
    & QTZG
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & B       , DKSAT   , DWSAT   , SMCMAX  , SMCWLT  , SMCREF  , SMCDRY  ,  F1      , QUARTZ  
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & REFSMC1 , WLTSMC1
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & PTU     
!--------------------------------------------------------------------------------------------------
! DANIEL: DOYLE
! THE M. EK PARAMETERS WAS CHANGED FOR THE MOIRA'S ONES (DOYLE@CIMA.FCEN.UBA.AR) RELATED TO THE 
! 0.25 X 0.25 SOIL TYPE
!--------------------------------------------------------------------------------------------------
    DATA MAXSMCD /                                                                                &
    & 0.240   , 0.354   , 0.514   , 0.402   , 0.477   , 0.564   , 0.440   , 0.566   , 0.631   ,   &
    & 0.443   , 0.517   , 0.566   , 0.478   , 0.562   , 0.675   , 0.505   , 0.534   , 0.580   ,   &  
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!
    DATA SATPSID /                                                                                &
    & 0.265   , 0.103   , 0.252   , 0.106   , 0.088   , 0.052   , 0.109   , 0.042   , 0.018   ,   &
    & 0.424   , 0.400   , 0.220   , 0.601   , 0.500   , 0.132   , 0.929   , 1.294   , 0.850   ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!
    DATA SATDKD /                                                                                 &
    & 1.47E-6 , 1.23E-6 , 6.72E-6 , 1.15E-5 , 7.57E-6 , 4.87E-6 , 4.46E-5 , 9.57E-6 , 7.40E-6 ,   &
    & 4.80E-5 , 1.61E-5 , 7.31E-6 , 6.73E-5 , 2.30E-5 , 4.85E-5 , 3.08E-5 , 6.03E-5 , 3.30E-5 ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!
    DATA BBD    /                                                                                 &
    & 2.02    , 2.78    , 3.03    , 2.25    , 2.62    , 3.02    , 1.97    , 2.53    , 2.98    ,   &
    & 2.13    , 2.44    , 2.90    , 2.01    , 2.42    , 2.97    , 2.81    , 1.83    , 2.29    ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!
    DATA QTZD   /                                                                                 &
    & 0.60    , 0.25    , 0.25    , 0.60    , 0.45    , 0.25    , 0.82    , 0.25    , 0.25    ,   &
    & 0.60    , 0.35    , 0.25    , 0.60    , 0.10    , 0.10    , 0.10    , 0.25    , 0.10    ,   &
    & 0.00    , 0.00    , 0.00    ,0.00     , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!------------------------------------ 
! DANIEL: PARAMETROS ORIGINAL ZOGBLER
!------------------------------------ 
    DATA MAXSMCZ /                                                                                &
    & 0.421   , 0.464   , 0.468   , 0.434   , 0.406   , 0.465   , 0.404   , 0.439   , 0.421   ,   &
    & 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000/
!
    DATA SATPSIZ /                                                                                &
    & 0.04    , 0.62    , 0.47    , 0.14    , 0.10    , 0.26    , 0.14    , 0.36    , 0.04    ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!
    DATA SATDKZ /                                                                                 &
    & 1.41E-5 , 0.20E-5 , 0.10E-5 , 0.52E-5 , 0.72E-5 , 0.25E-5 , 0.45E-5 , 0.34E-5 , 1.41E-5 ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!
    DATA BBZ    /                                                                                 &
    & 4.26    , 8.72    , 11.55   , 4.74    , 10.73   , 8.17    , 6.77    , 5.25    , 4.26    ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!
    DATA QTZZ   /                                                                                 &
    & 0.82    , 0.10    , 0.25    , 0.60    , 0.52    , 0.35    , 0.60    , 0.40    , 0.82    ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!----------------------- 
! ROL: PARAMETROS GLOBAL
!----------------------- 
    DATA MAXSMCG /                                                                                &
    & 0.339   , 0.421   , 0.434   , 0.476   , 0.476   , 0.439   , 0.404   , 0.464   , 0.465   ,   &
    & 0.406   , 0.468   , 0.468   , 0.439   , 0.200   , 0.421   , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!
    DATA SATPSIG /                                                                                &
    & 0.069   , 0.036   , 0.141   , 0.759   , 0.759   , 0.355   , 0.135   , 0.617   , 0.263   ,   &
    & 0.098   , 0.324   , 0.468   , 0.355   , 0.069   , 0.036   , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!
    DATA SATDKG /                                                                                 &
    & 1.07E-6 , 1.41E-5 , 5.23E-6 , 2.81E-6 , 2.81E-6 , 3.38E-6 , 4.45E-6 , 2.04E-6 , 2.45E-6 ,   &
    & 7.22E-6 , 1.34E-6 , 9.74E-7 , 3.38E-6 , 1.41E-4 , 1.41E-5 , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!
    DATA BBG    /                                                                                 &
    & 2.79    , 4.26    , 4.74    , 5.33    , 5.33    , 5.25    , 6.66    , 8.72    , 8.17    ,   &
    & 10.73   , 10.39   , 11.55   , 5.25    , 2.79    , 4.26    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!
    DATA QTZG   /                                                                                 &
    & 0.92    , 0.82    , 0.60    , 0.25    , 0.10    , 0.40    , 0.60    , 0.10    , 0.35    ,   &
    & 0.52    , 0.10    , 0.25    , 0.05    , 0.07    , 0.25    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!--------------------------------------------------------------------------------------------------
! THE FOLLOWING 5 PARAMETERS ARE DERIVED LATER IN REDPRM.F FROM THE SOIL DATA, AND ARE JUST GIVEN
! HERE FOR REFERENCE AND TO FORCE STATIC STORAGE ALLOCATION. -DAG LOHMANN, FEB. 2001
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
! DANIEL: DOYLE
! THE M. EK PARAMETERS WAS CHANGED FOR THE MOIRA'S ONES (DOYLE@CIMA.FCEN.UBA.AR) RELATED TO THE 
! 0.25X0.25 SOIL TYPE
!--------------------------------------------------------------------------------------------------

    DATA REFSMCD /                                                                                &
    & 0.117   , 0.207   , 0.389   , 0.159   , 0.245   , 0.330   , 0.114   , 0.237   , 0.318   ,   &
    & 0.214   , 0.302   , 0.359   , 0.227   , 0.325   , 0.386   , 0.354   , 0.296   , 0.357   ,   &
    & 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000/
!
    DATA WLTSMCD /                                                                                &
    & 0.084   , 0.163   , 0.322   , 0.111   , 0.188   , 0.269   , 0.070   , 0.185   , 0.267   ,   &
    & 0.121   , 0.197   , 0.263   , 0.111   , 0.195   , 0.275   , 0.196   , 0.117   , 0.194   ,   &
    & 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000/
!
    DATA DRYSMCD /                                                                                &
    & 0.084   , 0.163   , 0.322   , 0.111   , 0.188   , 0.269   , 0.070   , 0.185   , 0.267   ,   &
    & 0.121   , 0.197   , 0.263   , 0.111   , 0.195   , 0.275   , 0.196   , 0.117   , 0.194   ,   &
    & 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000/
!
    DATA SATDWD  /                                                                                &
    & 1.48E-06, 1.21E-06, 6.68E-06, 1.16E-05, 7.60E-06, 4.90E-06, 3.85E-05, 9.35E-06, 7.32E-06,   &
    & 4.85E-05, 1.61E-05, 7.22E-06, 6.77E-05, 2.30E-05, 4.95E-05, 3.06E-05, 6.05E-05, 3.28E-05,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!
    DATA F11D    /                                                                                &
    & -0.209  , -0.996  , -0.773  , -0.174  , -0.466  , -0.880  , 0.128   , -0.448  , -1.053  ,   &
    &  0.583  ,  0.398  , -0.170  ,  0.907  ,  0.648  , -0.051  , 0.537   ,  1.412  ,  0.977  ,   &
    &  0.000  ,  0.000  ,  0.000  ,  0.000  ,  0.000  ,  0.000  , 0.000   ,  0.000  ,  0.000  ,   &
    &  0.000  ,  0.000  ,  0.000/
!--------------------------------------------------------------------------------------------------
! DANIEL: PARAMETROS ORIGINAL ZOGBLER
! THE FOLLOWING 5 PARAMETERS ARE DERIVED LATER IN REDPRM.F FROM THE SOIL DATA, AND ARE JUST GIVEN 
! HERE FOR REFERENCE AND TO FORCE STATIC STORAGE ALLOCATION DAG LOHMANN, FEB. 2001
!--------------------------------------------------------------------------------------------------
    DATA REFSMCZ /                                                                                &
    & 0.283   , 0.387   , 0.412   , 0.312   , 0.338   , 0.382   , 0.315   , 0.329   , 0.283   ,   &
    & 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000/
!
    DATA WLTSMCZ /                                                                                &
    & 0.029   , 0.119   , 0.139   , 0.047   , 0.100   , 0.103   , 0.069   , 0.066   , 0.029   ,   &
    & 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000/
!
    DATA DRYSMCZ /                                                                                &
    & 0.029   , 0.119   , 0.139   , 0.047   , 0.100   , 0.103   , 0.069   , 0.066   , 0.029   ,   &
    & 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000/
!
    DATA SATDWZ  /                                                                                &
    & 5.71E-6 , 2.33E-5 , 1.16E-5 , 7.95E-6 , 1.90E-5 , 1.14E-5 , 1.06E-5 , 1.46E-5 , 5.71E-6 ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    ,   &
    & 0.00    , 0.00    , 0.00/
!
    DATA F11Z    /                                                                                &
    & -0.999  , -1.116  , -2.137  , -0.572  , -3.201  , -1.302  , -1.519  , -0.329  , -0.999  ,   &
    &  0.000  ,  0.000  ,  0.000  ,  0.000  ,  0.000  ,  0.000  ,  0.000  ,  0.000  ,  0.000  ,   &
    &  0.000  ,  0.000  ,  0.000  ,  0.000  ,  0.000  ,  0.000  ,  0.000  ,  0.000  ,  0.000  ,   &
    &  0.000  ,  0.000  ,  0.000/
!--------------------------- 
! ROL - PARAMETROS DO GLOBAL
!--------------------------- 
    DATA REFSMCG /                                                                                &
    & 0.236   , 0.383   , 0.383   , 0.360   , 0.383   , 0.329   , 0.314   , 0.387   , 0.382   ,   &
    & 0.338   , 0.404   , 0.412   , 0.329   , 0.170   , 0.283   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000/
!
    DATA WLTSMCG /                                                                                &
    & 0.010   , 0.028   , 0.047   , 0.084   , 0.084   , 0.066   , 0.067   , 0.120   , 0.103   ,   &
    & 0.100   , 0.126   , 0.138   , 0.066   , 0.006   , 0.028   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000/
!
    DATA DRYSMCG /                                                                                &
    & 0.010   , 0.028   , 0.047   , 0.084   , 0.084   , 0.066   , 0.067   , 0.120   , 0.103   ,   &
    & 0.100   , 0.126   , 0.138   , 0.066   , 0.006   , 0.028   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   , 0.000   ,   &
    & 0.000   , 0.000   , 0.000/
!
    DATA SATDWG  /                                                                                &
    & 0.608E-06         , 0.514E-05         , 0.805E-05         , 0.239E-04         ,             &
    & 0.239E-04         , 0.143E-04         , 0.990E-05         , 0.237E-04         ,             &
    & 0.113E-04         , 0.187E-04         , 0.964E-05         , 0.112E-04         ,             &
    & 0.143E-04         , 0.136E-03         , 0.514E-05         , 0.00              ,             &
    & 0.00              , 0.00              , 0.00              , 0.00              , 0.00    ,   &
    & 0.00              , 0.00              , 0.00              , 0.00              , 0.00    ,   &
    & 0.00              , 0.00              , 0.00              , 0.00/
!
    DATA F11G    /                                                                                &
    & -0.472   , -1.044   , -0.569   ,  0.162   ,  0.162   , -0.327   , -1.491   , -1.118   ,     &
    & -1.297   , -3.209   , -1.916   , -2.138   , -0.327   , -1.111   , -1.044   ,  0.000   ,     &
    &  0.000   ,  0.000   ,  0.000   ,  0.000   ,  0.000   ,  0.000   ,  0.000   ,  0.000   ,     &
    &  0.000   ,  0.000   ,  0.000   ,  0.000   ,  0.000   ,  0.000/
!--------------------------------------------------------------------------------------------------
! SET-UP VEGETATION PARAMETERS FOR A GIVEN VEGETAION TYPE
!
! INPUT: VEGTYP - VEGETATION TYPE (INTEGER INDEX)
! OUPUT: VEGETATION PARAMETERS
!         SHDFAC - VEGETATION GREENNESS FRACTION
!         RCMIN  - MIMIMUM STOMATAL RESISTANCE
!         RGL    - PARAMETER USED IN SOLAR RAD TERM OF CANOPY RESISTANCE FUNCTION
!         HS     - PARAMETER USED IN VAPOR PRESSURE DEFICIT TERM OF CANOPY RESISTANCE FUNCTION
!         SNUP   - THRESHOLD SNOW DEPTH (IN WATER EQUIVALENT M) THAT IMPLIES 100% SNOW COVER
!
! SSIB VEGETATION TYPES (DORMAN AND SELLERS, 1989; JAM)
!
!  1:   BROADLEAF-EVERGREEN TREES  (TROPICAL FOREST)
!  2:   BROADLEAF-DECIDUOUS TREES
!  3:   BROADLEAF AND NEEDLELEAF TREES (MIXED FOREST)
!  4:   NEEDLELEAF-EVERGREEN TREES
!  5:   NEEDLELEAF-DECIDUOUS TREES (LARCH)
!  6:   BROADLEAF TREES WITH GROUNDCOVER (SAVANNA)
!  7:   GROUNDCOVER ONLY (PERENNIAL)
!  8:   BROADLEAF SHRUBS WITH PERENNIAL GROUNDCOVER
!  9:   BROADLEAF SHRUBS WITH BARE SOIL
! 10:   DWARF TREES AND SHRUBS WITH GROUNDCOVER (TUNDRA)
! 11:   BARE SOIL
! 12:   CULTIVATIONS (THE SAME PARAMETERS AS FOR TYPE 7)
! 13:   GLACIAL (THE SAME PARAMETERS AS FOR TYPE 11)
! 14:   URBAN (PARAMETERS CAROL)
! 15:
! 16:
!--------------------------------------------------------------------------------------------------
    INTEGER(KIND=I4KIND), DIMENSION(MAX_VEGTYP)                                                 ::&
    & NROOT_DATA
!
    REAL   (KIND=R4KIND), DIMENSION(MAX_VEGTYP)                                                 ::&
    & RSMTBL  , RGLTBL  , HSTBL   , SNUPX   , LAI_DATA
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(INOUT)       ::&
    & NROOT
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & FRZFACT 
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & SHDFAC  , RCMIN   , RGL     , HS      , PSISAT  , SNUP    , Z0      , LAI 
!
    DATA NROOT_DATA /8  , 4       , 4       , 4       , 4       , 7       , 3       , 3       ,   &
    &                3  , 2       , 3       , 3       , 2       , 1       , 5       , 4       ,   &
    &                0  , 0       , 0       , 0       , 0       , 0       , 0       , 0       ,   &
    &                0  , 0       , 0       , 0       , 0       , 0/
!
    DATA RSMTBL /150.0  , 100.0   , 125.0   , 150.0   , 100.0   ,  70.0   ,  40.0   , 300.0   ,   &
    &            400.0  , 150.0   , 400.0   ,  40.0   , 150.0   , 200.0   , 150.0   ,  40.0   ,   &
    &              0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   &
    &              0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0/
!
    DATA RGLTBL / 30.0  ,  30.0   ,  30.0   ,  30.0   ,  30.0   ,  65.0   , 100.0   , 100.0   ,   &
    &            100.0  , 100.0   , 100.0   , 100.0   , 100.0   , 100.0   ,  30.0   , 100.0   ,   &
    &              0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   &
    &              0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0/
!
    DATA HSTBL  / 41.69 ,  54.53  ,  51.93  ,  47.35  ,  47.35  ,  54.53  ,  36.35  ,  42.00  ,   &
    &             42.00 ,  42.00  ,  42.00  ,  36.35  ,  42.00  ,  42.00  ,  41.69  ,  36.35  ,   &
    &              0.00 ,   0.00  ,   0.00  ,   0.00  ,   0.00  ,   0.00  ,   0.00  ,   0.00  ,   &
    &              0.00 ,   0.00  ,   0.00  ,   0.00  ,   0.00  ,   0.00/
!
    DATA SNUPX  /  0.080,   0.080 ,   0.080 ,   0.080 ,   0.080 ,   0.080 ,   0.040 ,   0.040 ,   &
    &              0.040,   0.040 ,   0.025 ,   0.040 ,   0.025 ,   0.020 ,   0.080 ,   0.040 ,   &
    &              0.000,   0.000 ,   0.000 ,   0.000 ,   0.000 ,   0.000 ,   0.000 ,   0.000 ,   &
    &              0.000,   0.000 ,   0.000 ,   0.000 ,   0.000 ,   0.000/
!
!    DATA Z0_DATA / 2.653,   0.826 ,   0.563 ,   1.089 ,   0.854 ,   0.856 ,   0.035 ,   0.238 ,   &
!    &              0.065,   0.076 ,   0.011 ,   0.035 ,   0.011 ,   1.000 ,   1.800 ,   0.070 ,   &
!    &              0.000,   0.000 ,   0.000 ,   0.000 ,   0.000 ,   0.000 ,   0.000 ,   0.000 ,   &
!    &              0.000,   0.000 ,   0.000 ,   0.000 ,   0.000 ,   0.000/
!
    DATA LAI_DATA /4.0  ,   4.0   ,   4.0   ,   4.0   ,   4.0   ,   4.0   ,   4.0   ,   4.0   ,   &
    &              4.0  ,   4.0   ,   4.0   ,   4.0   ,   4.0   ,   4.0   ,   4.0   ,   4.0   ,   &
    &              0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   &
    &              0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0/
!--------------------------------------------------------------------------------------------------
! CLASS PARAMETER 'SLOPETYP' WAS INCLUDED TO ESTIMATE LINEAR RESERVOIR COEFFICIENT 'SLOPE' TO THE 
! BASEFLOW RUNOFF OUT OF THE BOTTOM LAYER. LOWEST CLASS (SLOPETYP=0) MEANS HIGHEST SLOPE 
! PARAMETER = 1 DEFINITION OF SLOPETYP FROM 'ZOBLER' SLOPE TYPE
! SLOPE CLASS      PERCENT SLOPE
! 1                0-8
! 2                8-30
! 3                > 30
! 4                0-30
! 5                0-8 & > 30
! 6                8-30 & > 30
! 7                0-8, 8-30, > 30
! 9                GLACIAL ICE
! BLANK            OCEAN/SEA
! NOTE:  CLASS 9 FROM 'ZOBLER' FILE SHOULD BE REPLACED BY 8 AND 'BLANK'  9
!--------------------------------------------------------------------------------------------------
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & SLOPE
!
    REAL   (KIND=R4KIND), DIMENSION(MAX_SLOPETYP)                                               ::&
    & SLOPE_DATA
!
    DATA SLOPE_DATA /                                                                             &
    & 0.1     , 0.6     , 1.0     , 0.35    , 0.55    , 0.8     , 0.63    , 0.0     , 0.0     ,   &
    & 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     ,   &
    & 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     ,   &
    & 0.0     , 0.0     , 0.0/
!----------------------- 
! SET NAMELIST FILE NAME
!-----------------------
    CHARACTER(50) :: NAMELIST_NAME
!------------------------------------------------------------------
! SET UNIVERSAL PARAMETERS (NOT DEPENDENT ON SOIL, VEG, SLOPE TYPE)
!------------------------------------------------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & VEGTYP  , SOILTYP , SLOPETYP, NSOIL
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , BARE
!
    DATA    BARE /11/
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & LPARAM
!
    DATA    LPARAM / .TRUE. /
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & LFIRST
!
    DATA    LFIRST / .FALSE. /
!------------------------------------------------------ 
!  PARAMETER USED TO CALCULATE ROUGHNESS LENGTH OF HEAT
!------------------------------------------------------ 
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & CZIL
!
!GSM    REAL   (KIND=R4KIND)                                                                        ::&
!GSM    & CZIL_DATA
!
!GSM    DATA CZIL_DATA /0.2/
!------------------------------------------------------------------
!  PARAMETER USED TO CALUCULATE VEGETATION EFFECT ON SOIL HEAT FLUX
!------------------------------------------------------------------
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & SBETA
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & SBETA_DATA
!
    DATA SBETA_DATA /-2.0/
!--------------------------------------------- 
! BARE SOIL EVAPORATION EXPONENT USED IN DEVAP
!--------------------------------------------- 
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & FXEXP
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & FXEXP_DATA
!
    DATA FXEXP_DATA /2.0/
!----------------------------- 
! SOIL HEAT CAPACITY [J/M^3/K]
!----------------------------- 
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & CSOIL
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & CSOIL_DATA
!
    DATA CSOIL_DATA /1.26E+6/
!--------------------------------------------------------------------------------------------------
! SPECIFY SNOW DISTRIBUTION SHAPE PARAMETER SALP - SHAPE PARAMETER OF DISTRIBUTION FUNCTION OF SNOW
! COVER.FROM ANDERSONS DATA (HYDRO-17) BEST FIT IS WHEN SALP = 2.6
!--------------------------------------------------------------------------------------------------
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & SALP
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & SALP_DATA
!
    DATA SALP_DATA /2.6/
!--------------------------------------------------------------------------------------------------
! KDT IS DEFINED BY REFERENCE REFKDT AND DKSAT REFDK=2.E-6 IS THE SAT. DK. 
! VALUE FOR THE SOIL TYPE 2
!--------------------------------------------------------------------------------------------------
    REAL   (KIND=R4KIND)                                                                        ::&
    & REFDK   , REFDK_DATA
!
    DATA REFDK_DATA /2.0E-6/
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & REFKDT
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & REFKDT_DATA
!
    DATA REFKDT_DATA /3.0/
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & KDT     , FRZX
!--------------------------------------------------------------------------------------------------
! FROZEN GROUND PARAMETER, FRZK, DEFINITION FRZK IS ICE CONTENT THRESHOLD ABOVE WHICH FROZEN SOIL 
! IS IMPERMEABLE REFERENCE VALUE OF THIS PARAMETER FOR THE LIGHT CLAY SOIL (TYPE=3) FRZK = 0.15 M
!--------------------------------------------------------------------------------------------------
    REAL   (KIND=R4KIND)                                                                        ::&
    & FRZK    , FRZK_DATA
!
    DATA FRZK_DATA /0.15/
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & SLDPTH  , ZSOIL
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(INOUT)       ::&
    & RTDIS   
!-------------------------------- 
! SET TWO CANOPY WATER PARAMETERS
!--------------------------------
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & CFACTR  , CMCMAX
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & CFACTR_DATA       ,                                                               &
    & CMCMAX_DATA
!
    DATA CFACTR_DATA /0.5/
    DATA CMCMAX_DATA /0.5E-3/
!-----------------------------
! SET MAX. STOMATAL RESISTANCE
!-----------------------------
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & RSMAX
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & RSMAX_DATA
!
    DATA RSMAX_DATA /5000.0/
!------------------------------------------
! SET OPTIMUM TRANSPIRATION AIR TEMPERATURE
!------------------------------------------
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & TOPT
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & TOPT_DATA
!
    DATA TOPT_DATA /298.0/
!----------------------------------------------------
! SPECIFY DEPTH[M] OF LOWER BOUNDARY SOIL TEMPERATURE
!----------------------------------------------------
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & ZBOT    
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ZBOT_DATA
!
    DATA ZBOT_DATA /-15.0/
!-------------------- 
! NAMELIST DEFINITION
!--------------------
!GSM    NAMELIST /SOIL_VEG/ SLOPE_DATA   , RSMTBL       , RGLTBL       , HSTBL        , SNUPX        ,&
!GSM    &                   BB           , DRYSMC       , F11          , MAXSMC       , REFSMC       ,&
!GSM    &                   SATPSI       , SATDK        , SATDW        , WLTSMC       , QTZ          ,&
!GSM    &                   LPARAM       , ZBOT_DATA    , SALP_DATA    , CFACTR_DATA  , CMCMAX_DATA  ,&
!GSM    &                   SBETA_DATA   , RSMAX_DATA   , TOPT_DATA    , REFDK_DATA   , FRZK_DATA    ,&
!GSM    &                   BARE         , DEFINED_VEG  , DEFINED_SOIL , DEFINED_SLOPE, FXEXP_DATA   ,&
!GSM    &                   NROOT_DATA   , REFKDT_DATA  , Z0_DATA      , CZIL_DATA    , LAI_DATA     ,&
!GSM    &                   CSOIL_DATA
!------ 
! ANIEL
!------ 
    IF (NSOTYP == 18) THEN
        DEFINED_SOIL = DEFINED_SOILD
!
        DO I=1,DEFINED_SOIL
                BB(I) =     BBD(I)
            DRYSMC(I) = DRYSMCD(I)
               F11(I) =    F11D(I)
            MAXSMC(I) = MAXSMCD(I)
            REFSMC(I) = REFSMCD(I)
            SATPSI(I) = SATPSID(I)
             SATDK(I) =  SATDKD(I)
             SATDW(I) =  SATDWD(I)
            WLTSMC(I) = WLTSMCD(I)
               QTZ(I) =    QTZD(I)
        END DO
    END IF
!------ 
! ANIEL
!------ 
!
!----------------
! RIGINAL ZOGBLER
!----------------
    IF (NSOTYP == 9) THEN
        DEFINED_SOIL = DEFINED_SOILZ
!
        DO I=1,DEFINED_SOIL
                BB(I) =     BBZ(I)
            DRYSMC(I) = DRYSMCZ(I)
               F11(I) =    F11Z(I)
            MAXSMC(I) = MAXSMCZ(I)
            REFSMC(I) = REFSMCZ(I)
            SATPSI(I) = SATPSIZ(I)
             SATDK(I) =  SATDKZ(I)
             SATDW(I) =  SATDWZ(I)
            WLTSMC(I) = WLTSMCZ(I)
               QTZ(I) =    QTZZ(I)
        END DO
    END IF
!------------- 
! ROL - GLOBAL
!-------------
    IF (NSOTYP == 15) THEN
        DEFINED_SOIL = DEFINED_SOILG
!        		
        DO I=1,DEFINED_SOIL
                BB(I) =     BBG(I)
            DRYSMC(I) = DRYSMCG(I)
               F11(I) =    F11G(I)
            MAXSMC(I) = MAXSMCG(I)
            REFSMC(I) = REFSMCG(I)
            SATPSI(I) = SATPSIG(I)
             SATDK(I) =  SATDKG(I)
             SATDW(I) =  SATDWG(I)
            WLTSMC(I) = WLTSMCG(I)
               QTZ(I) =    QTZG(I)
        END DO
    END IF
!---- 
! ROL
!---- 
!
!-------------------------------------------------------------    
! READ NAMELIST FILE TO OVERRIDE DEFAULT PARAMETERS ONLY ONCE.
!------------------------------------------------------------- 
!    IF (LFIRST) THEN
!        OPEN (58, FILE = 'NAMELIST_FILENAME.txt')
!---------------------------------------------
! NAMELIST_NAME MUST BE 50 CHARACTERS OR LESS.
!---------------------------------------------
!        READ (58,'(A)') NAMELIST_NAME
!        CLOSE (58)
!
!        WRITE(*,*) 'NAMELIST FILENAME IS', NAMELIST_NAME
!
         50 CONTINUE
!
        100 CONTINUE
!
        LFIRST = .FALSE. 
!
        IF (DEFINED_SOIL > MAX_SOILTYP) THEN
            PRINT*,DEFINED_SOIL
            WRITE(*,*) 'WARNING: DEFINED_SOIL TOO LARGE IN NAMELIST'
            STOP 222
        END IF
!
        IF (DEFINED_VEG > MAX_VEGTYP) THEN
            WRITE(*,*) 'WARNING: DEFINED_VEG TOO LARGE IN NAMELIST'
            STOP 222
        END IF
!
        IF (DEFINED_SLOPE > MAX_SLOPETYP) THEN
            WRITE(*,*) 'WARNING: DEFINED_SLOPE TOO LARGE IN NAMELIST'
            STOP 222
        END IF
!
        DO I = 1, DEFINED_SOIL
            SATDW(I)  = BB(I) * SATDK(I) * (SATPSI(I) / MAXSMC(I))
            F11(I)    = ALOG10(SATPSI(I)) + BB(I) * ALOG10(MAXSMC(I)) + 2.0
            REFSMC1   = MAXSMC(I) * (5.79E-9 / SATDK(I)) ** (1.0 / (2.0 * BB(I) + 3.0))
            REFSMC(I) = REFSMC1 + (MAXSMC(I) - REFSMC1) / 3.0
            WLTSMC1   = MAXSMC(I) * (200.0 / SATPSI(I)) ** (-1.0 / BB(I))
            WLTSMC(I) = WLTSMC1 - 0.5 * WLTSMC1
!------------------------------------------------------------------
! CURRENT VERSION DRYSMC VALUES THAT EQUATE TO WLTSMC
! FUTURE VERSION COULD LET DRYSMC BE INDEPENDENTLY SET VIA NAMELIST
!------------------------------------------------------------------
            DRYSMC(I) = WLTSMC(I)
        END DO
!
!    END IF
!        
    IF (SOILTYP > DEFINED_SOIL) THEN
        WRITE(*,*) 'WARNING: TOO MANY SOIL TYPES'
        STOP 333
    END IF
!
    IF (VEGTYP > DEFINED_VEG) THEN
        WRITE(*,*) 'WARNING: TOO MANY VEG TYPES'
        STOP 333
    END IF
!
    IF (SLOPETYP > DEFINED_SLOPE) THEN
        WRITE(*,*) 'WARNING: TOO MANY SLOPE TYPES'
        STOP 333
    END IF
!---------------------------------------------------------------------------
! SET-UP UNIVERSAL PARAMETERS (NOT DEPENDENT ON SOILTYP, VEGTYP OR SLOPETYP)
!---------------------------------------------------------------------------
    ZBOT   =   ZBOT_DATA
    SALP   =   SALP_DATA
    CFACTR = CFACTR_DATA
    CMCMAX = CMCMAX_DATA
    SBETA  =  SBETA_DATA
    RSMAX  =  RSMAX_DATA
    TOPT   =   TOPT_DATA
    REFDK  =  REFDK_DATA
    FRZK   =   FRZK_DATA
    FXEXP  =  FXEXP_DATA
    REFKDT = REFKDT_DATA
    CZIL   =   CZIL_DATA(VEGTYP)
    CSOIL  =  CSOIL_DATA
!-----------------------
! SET-UP SOIL PARAMETERS
!-----------------------
    B       =     BB(SOILTYP)
    SMCDRY  = DRYSMC(SOILTYP)
    F1      =    F11(SOILTYP)
    SMCMAX  = MAXSMC(SOILTYP)
    SMCREF  = REFSMC(SOILTYP)
    PSISAT  = SATPSI(SOILTYP)
    DKSAT   =  SATDK(SOILTYP)
    DWSAT   =  SATDW(SOILTYP)
    SMCWLT  = WLTSMC(SOILTYP)
    QUARTZ  =    QTZ(SOILTYP)
!
    FRZFACT = (SMCMAX / SMCREF) * (0.412 / 0.468)
    KDT     = REFKDT * DKSAT / REFDK
!-------------------------------------------------------------
! TO ADJUST FRZK PARAMETER TO ACTUAL SOIL TYPE: FRZK * FRZFACT
!-------------------------------------------------------------
    FRZX = FRZK * FRZFACT
!----------------------------- 
! SET-UP VEGETATION PARAMETERS
!-----------------------------
    NROOT = NROOT_DATA(VEGTYP)
    SNUP  =      SNUPX(VEGTYP)
    RCMIN =     RSMTBL(VEGTYP)
    RGL   =     RGLTBL(VEGTYP)
    HS    =      HSTBL(VEGTYP)
    Z0    =    Z0_DATA(VEGTYP)
    LAI   =   LAI_DATA(VEGTYP)
!
    IF (VEGTYP == BARE) SHDFAC = 0.0
!----
! ROL
!----
    IF (VEGTYP == 14) SHDFAC = 0.1
!    						
    IF (NROOT > NSOIL) THEN
        WRITE(*,*) 'WARNING: TOO MANY ROOT LAYERS '
        STOP 333
    END IF
!------------------------------------------------------------------
! CALCULATE ROOT DISTRIBUTION
! PRESENT VERSION ASSUMES UNIFORM DISTRIBUTION BASED ON SOIL LAYERS
!------------------------------------------------------------------
    DO I=1,NROOT
        RTDIS(I) = -SLDPTH(I) / ZSOIL(NROOT)
    END DO
!----------------------- 
! SET-UP SLOPE PARAMETER
!----------------------- 
    SLOPE = SLOPE_DATA(SLOPETYP)
!
    RETURN
    END SUBROUTINE REDPRM
!
!
!
    SUBROUTINE ROSR12(P, A, B, C, D, DELTA, NSOIL)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE ROSR12
!>
!> SUBROUTINE: ROSR12 - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> TO INVERT (SOLVE) THE TRI-DIAGONAL MATRIX PROBLEM SHOWN
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NSOIL -
!> A     -
!> B     -
!> D     -
!>
!> OUTPUT ARGUMENT LIST:
!> DELTA -
!> P     -
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> C     -
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : SFLX
!>
!> CALLS      : -----
!>-------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & K       , KK      , NSOIL
!          
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(INOUT)       ::&
    & P       , DELTA
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & A       , B       , D
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(INOUT)       ::&
    & C
!------------------------------------------------- 
! INITIALIZE EQN COEF C FOR THE LOWEST SOIL LAYER.
!------------------------------------------------- 
    C(NSOIL) = 0.0
!--------------------------------------- 
! SOLVE THE COEFS FOR THE 1ST SOIL LAYER
!---------------------------------------  
        P(1) = -C(1) / B(1)
    DELTA(1) =  D(1) / B(1)
!--------------------------------------------- 
! SOLVE THE COEFS FOR SOIL LAYERS 2 THRU NSOIL
!--------------------------------------------- 
    DO K=2,NSOIL
            P(K) = -C(K) * (1.0 / (B(K) + A(K) * P(K-1)))
        DELTA(K) = (D(K) - A(K) * DELTA(K-1)) * (1.0 / (B(K) + A(K) * P(K-1)))
    END DO 
!-------------------------------------- 
! SET P TO DELTA FOR LOWEST SOIL LAYER.
!-------------------------------------- 
    P(NSOIL) = DELTA(NSOIL)
!-------------------------------------- 
! ADJUST P FOR SOIL LAYERS 2 THRU NSOIL
!-------------------------------------- 
    DO K=2,NSOIL
        KK = NSOIL - K + 1
        P(KK) = P(KK) * P(KK+1) + DELTA(KK)
    END DO
!
    RETURN
    END SUBROUTINE ROSR12
!
!
!
    SUBROUTINE SHFLX(S     , STC   , SMC   , SMCMAX, NSOIL , T1    , DT    , YY    , ZZ1   ,      &
    &                ZSOIL , TBOT  , ZBOT  , SMCWLT, PSISAT, SH2O  , B     , F1    , DF1   ,      &
    &                ICE   , QUARTZ, CSOIL)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SHFLX
!>
!> SUBROUTINE: SHFLX - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> UPDATE THE TEMPERATURE STATE OF THE SOIL COLUMN BASED ON THE THERMAL DIFFUSION EQUATION AND 
!> UPDATE THE FROZEN SOIL MOISTURE CONTENT BASED ON THE TEMPERATURE.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> S      -
!> SMC    -
!> SMCMAX -
!> NSOIL  -
!> DT     - 
!> YY     -
!> ZZ1    -
!> ZSOIL  -
!> TBOT   -
!> ZBOT   -
!> SMCWLT - 
!> PSISAT - 
!> SH2O   -
!> B      - 
!> F1     -
!> DF1    -
!> ICE    -
!> QUARTZ - 
!> CSOIL  -
!>
!> OUTPUT ARGUMENT LIST:
!> T1     -
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> STC    - 
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : SFLX
!>
!> CALLS      : HRT
!>              HRTICE
!>              HSTEP
!>-------------------------------------------------------------------------------------------------   
    USE F77KINDS
!   
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND), PARAMETER :: NSOLD =  20
!
    REAL   (KIND=R4KIND), PARAMETER :: T0    = 273.15
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , ICE     , IFRZ    
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & NSOIL
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & B       , DF1     , CSOIL   , DT      , F1      , PSISAT  , QUARTZ  , S       , SMCMAX  ,   &
    & SMCWLT  , TBOT    , ZBOT    , YY      , ZZ1
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & T1
!
    REAL   (KIND=R4KIND), DIMENSION(NSOLD)                                                      ::&
    & RHSTS   , STCF
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & SH2O    , ZSOIL
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & SMC  
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(INOUT)       ::&
    & STC
!--------------------------------------------------------------- 
! HRT ROUTINE CALCS THE RIGHT HAND SIDE OF THE SOIL TEMP DIF EQN
!--------------------------------------------------------------- 
    IF (ICE == 1) THEN
!------------- 
! SEA-ICE CASE
!-------------
        CALL HRTICE(RHSTS, STC, NSOIL, ZSOIL, YY, ZZ1, DF1)
!
        CALL HSTEP (STCF, STC, RHSTS, DT, NSOIL)
 !                
    ELSE
!---------------
! LAND-MASS CASE
!---------------
        CALL HRT(RHSTS, STC, SMC, SMCMAX, NSOIL, ZSOIL, YY, ZZ1, TBOT, ZBOT, PSISAT, SH2O, DT, B, &
    &            F1   , DF1, QUARTZ, CSOIL)
 !                
        CALL HSTEP(STCF, STC, RHSTS, DT, NSOIL)
!
    END IF
!
    DO I=1,NSOIL
        STC(I) = STCF(I)
    END DO
!--------------------------------------------------------------------------------------------------
! IN THE NO SNOWPACK CASE (VIA ROUTINE NOPAC BRANCH,) UPDATE THE GRND (SKIN) TEMPERATURE HERE IN 
! RESPONSE TO THE UPDATED SOIL TEMPERATURE PROFILE ABOVE. 
!(NOTE: INSPECTION OF ROUTINE SNOPAC SHOWS THAT T1 BELOW IS A DUMMY VARIABLE ONLY, AS SKIN 
! TEMPERATURE IS UPDATED DIFFERENTLY IN ROUTINE SNOPAC)
!--------------------------------------------------------------------------------------------------
    T1 = (YY + (ZZ1 - 1.0) * STC(1)) / ZZ1
!---------------------------- 
! CALC THE SFC SOIL HEAT FLUX
!---------------------------- 

    RETURN
!
    END SUBROUTINE SHFLX
!
!
!
    SUBROUTINE SMFLX (ETA1    , SMC    , NSOIL  , CMC    , ETP1   , DT     , PRCP1  , ZSOIL  ,    &
    &                 SH2O    , SLOPE  , KDT    , FRZFACT, SMCMAX , B      , PC     , SMCWLT ,    &
    &                 DKSAT   , DWSAT  , SMCREF , SHDFAC , CMCMAX , SMCDRY , CFACTR , RUNOFF1,    &
    &                 RUNOFF2 , RUNOFF3, EDIR1  , EC1    , ETT1   , SFCTMP , Q2     , NROOT  ,    &
    &                 RTDIS   , FXEXP)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SMFLX
!>
!> SUBROUTINE: SMFLX - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> TO CALCULATE SOIL MOISTURE FLUX. THE SOIL MOISTURE CONTENT (SMC - A PER UNIT VOLUME MEASUREMENT)
!> IS A DEPENDENT VARIABLE THAT IS UPDATED WITH PROGNOSTIC EQNS. 
!> THE CANOPY MOISTURE CONTENT (CMC) IS ALSO UPDATED.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> SMC     -
!> NSOIL   -
!> CMC     -
!> ETP1    -
!> DT      -
!> PRCP1   -
!> ZSOIL   -
!> SH2O    -
!> SLOPE   -
!> KDT     -
!> FRZFACT -
!> SMCMAX  -
!> B       -
!> PC      -
!> SMCWLT  -
!> DKSAT   -
!> DWSAT   -
!> SMCREF  -
!> SHDFAC  -
!> CMCMAX  -
!> SMCDRY  -
!> CFACTR  -
!> RUNOFF1 -
!> RUNOFF2 -
!> FUNOFF3 -
!> EDIR1   -
!> EC1     -
!> ETT1    -
!> SFCTMP  -
!> Q2      -
!> NROOT   -
!> RTDIS   -
!> FXEXP   -
!>
!> OUTPUT ARGUMENT LIST:
!> ETA1    -
!> EDIR1   - 
!> EC1     -
!> ETT1    -
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: F77KINDS
!>              RITE
!>
!> DRIVER     : SFLX
!>
!> CALLS      : SRT
!>              SSTEP
!>              TRANSP
!>------------------------------------------------------------------------------------------------- 
    USE F77KINDS
    USE RITE    , RUNOXX3 => RUNOFF3
!
    IMPLICIT NONE
!----------------------------------------------------------------------------------------  
! FROZEN GROUND VERSION
! NEW STATES ADDED: SH2O, AND FROZEN GROUD CORRECTION FACTOR, FRZFACT AND PARAMETER SLOPE
!---------------------------------------------------------------------------------------- 
    INTEGER(KIND=I4KIND), PARAMETER :: NSOLD  =  20
    REAL   (KIND=R4KIND), PARAMETER :: TFREEZ = 273.15
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & NSOIL   , NROOT
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & CMC     , ETP1    , DT      , PRCP1   , SLOPE   , KDT     , FRZFACT , SMCMAX  ,             &
    & B       , PC      , SMCWLT  , DKSAT   , DWSAT   , SMCREF  , SHDFAC  , CMCMAX  , SMCDRY  ,   &
    & CFACTR  , RUNOFF1 , RUNOFF2 , RUNOFF3 , SFCTMP  , Q2      , FXEXP
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & ETA1    , EDIR1   , EC1     , ETT1 
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & SMC     , ZSOIL   , SH2O    , RTDIS
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & K
!
    REAL   (KIND=R4KIND), DIMENSION(NSOLD)                                                      ::&
    & ET      , RHSTT
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & EXCESS  , PCPDRP  , RHSCT   , RUNOFF
!---------------------- 
! FROZEN GROUND VERSION 
!---------------------- 
    REAL   (KIND=R4KIND), DIMENSION(NSOLD)                                                      ::&
    & SICE    , SH2OA   , SH2OFG         
! 
    REAL   (KIND=R4KIND)                                                                        ::&
    & TRHSCT
!-----------------------------------------------------------------------------
! TEMPERATURE CRITERIA FOR SNOWFALL TFREEZ SHOULD HAVE SAME VALUE AS IN SFLX.F
!-----------------------------------------------------------------------------
    REAL   (KIND=R4KIND)                                                                        ::&
    & DUMMY   , CMC2MS  , DEVAP
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I
!---------------------------------------------------------------------------------------
! EXECUTABLE CODE BEGINS HERE. IF THE POTENTIAL EVAPOTRANSPIRATION IS GREATER THAN ZERO.
!---------------------------------------------------------------------------------------
    DUMMY = 0.
    EDIR  = 0.
    EC    = 0.
    ETT   = 0.
!
    DO K=1,NSOIL
        ET(K) = 0.
    END DO
!
    IF (ETP1 > 0.0) THEN
!---------------------------------------------- 
! RETRIEVE DIRECT EVAPORATION FROM SOIL SURFACE
!----------------------------------------------
!
!-------------------------------------------------- 
! CALL THIS FUNCTION ONLY IF VEG COVER NOT COMPLETE
! FROZEN GROUND VERSION
! SMC STATES WERE REPLACED BY SH2O STATES
!-------------------------------------------------- 
           IF (SHDFAC < 1.) THEN
            EDIR = DEVAP(ETP1  , SH2O(1), ZSOIL(1), SHDFAC, SMCMAX, B, DKSAT, DWSAT, SMCDRY,      &
    &                    SMCREF, SMCWLT , FXEXP)
           END IF
!--------------------------------------------------------------------------------------------------
! INITIALIZE PLANT TOTAL TRANSPIRATION, RETRIEVE PLANT TRANSPIRATION, AND ACCUMULATE IT FOR ALL 
! SOIL LAYERS.
!--------------------------------------------------------------------------------------------------
        IF (SHDFAC > 0.0) THEN
!---------------------------------------- 
! FROZEN GROUND VERSION 
! SMC STATES WERE REPLACED BY SH2O STATES
!----------------------------------------         
            CALL TRANSP(ET    , NSOIL , ETP1  , SH2O  , CMC   , ZSOIL , SHDFAC, SMCWLT, CMCMAX,   &
    &                   PC    , CFACTR, SMCREF, SFCTMP, Q2    , NROOT , RTDIS)
!                      
            DO K=1,NSOIL
                ETT = ETT + ET(K)
            END DO
!----------------------------------------------------------------- 
! MOVE THIS ENDIF AFTER CANOPY EVAP CALCS SINCE CMC=0 FOR SHDFAC=0
!----------------------------------------------------------------- 
!
!----------------------------- 
! CALCULATE CANOPY EVAPORATION
!----------------------------- 
!
!------------------------------------------------------------- 
! IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR CMC=ZERO
!-------------------------------------------------------------
            IF (CMC > 0.0) THEN
                EC = SHDFAC * ((CMC / CMCMAX) ** CFACTR) * ETP1
            ELSE
                EC = 0.0
            END IF
!--------------------------------------------------------------------------------------------------
! EC SHOULD BE LIMITED BY THE TOTAL AMOUNT OF AVAILABLE WATER ON THE CANOPY. 
! MODIFIED BY F.CHEN ON 10/18/94
!--------------------------------------------------------------------------------------------------
            CMC2MS = CMC / DT
            EC = MIN(CMC2MS, EC)
        END IF
    END IF
!------------------------------------------------------------ 
! TOTAL UP EVAP AND TRANSP TYPES TO OBTAIN ACTUAL EVAPOTRANSP
!------------------------------------------------------------ 
    EDIR1 = EDIR
    EC1   = EC
    ETT1  = ETT
!          
    ETA1 = EDIR + ETT + EC
!----------------------------------------------------------- 
! COMPUTE THE RIGHT HAND SIDE OF THE CANOPY EQN TERM (RHSCT)
!-----------------------------------------------------------
    RHSCT = SHDFAC * PRCP1 - EC
!--------------------------------------------------------------------------------------------------
! CONVERT RHSCT (A RATE) TO TRHSCT (AN AMT) AND ADD IT TO EXISTING CMC. 
! IF RESULTING AMT EXCEEDS MAX CAPACITY, IT BECOMES DRIP AND WILL FALL TO THE GRND.
!--------------------------------------------------------------------------------------------------
    DRIP   = 0.
    TRHSCT = DT  * RHSCT
    EXCESS = CMC + TRHSCT
!
    IF (EXCESS > CMCMAX) DRIP = EXCESS - CMCMAX
!------------------------------------------------------------------------- 
! PCPDRP IS THE COMBINED PRCP1 AND DRIP (FROM CMC) THAT GOES INTO THE SOIL
!------------------------------------------------------------------------- 
    PCPDRP = (1. - SHDFAC) * PRCP1 + DRIP / DT

!---------------------------------------------------------------- 
! FROZEN GROUND VERSION
! STORE ICE CONTENT AT EACH SOIL LAYER BEFORE CALLING SRT & SSTEP
!---------------------------------------------------------------- 
    DO I=1,NSOIL
        SICE(I) = SMC(I) - SH2O(I)
    END DO
!--------------------------------------------------------------------------------------------------   
! CALL SUBROUTINES SRT AND SSTEP TO SOLVE THE SOIL MOISTURE TENDENCY EQUATIONS.
!
! IF THE INFILTRATING PRECIP RATE IS NONTRIVIAL, (WE CONSIDER NONTRIVIAL TO BE A PRECIP TOTAL OVER 
! THE TIME STEP EXCEEDING ONE ONE-THOUSANDTH OF THE WATER HOLDING CAPACITY OF THE FIRST SOIL LAYER)
!
! THEN CALL THE SRT/SSTEP SUBROUTINE PAIR TWICE IN THE MANNER OF TIME SCHEME "F" (IMPLICIT STATE, 
! AVERAGED COEFFICIENT) OF SECTION 2 OF KALNAY AND KANAMITSU (1988, MWR, VOL 116, PAGES 1945-1958) 
! TO MINIMIZE 2-DELTA-T OSCILLATIONS IN THE SOIL MOISTURE VALUE OF THE TOP SOIL LAYER THAT CAN 
! ARISE BECAUSE OF THE EXTREME NONLINEAR DEPENDENCE OF THE SOIL HYDRAULIC DIFFUSIVITY COEFFICIENT 
! AND THE HYDRAULIC CONDUCTIVITY ON THE SOIL MOISTURE STATE
!
! OTHERWISE CALL THE SRT/SSTEP SUBROUTINE PAIR ONCE IN THE MANNER OF TIME SCHEME "D" 
! (IMPLICIT STATE, EXPLICIT COEFFICIENT) OF SECTION 2 OF KALNAY AND KANAMITSU
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------- 
! PCPDRP IS UNITS OF KG/M**2/S OR MM/S, ZSOIL IS NEGATIVE DEPTH IN M IF (PCPDRP .GT. 0.0) THEN
!---------------------------------------------------------------------------------------------
    IF ((PCPDRP * DT) > (0.001 * 1000.0 * (-ZSOIL(1)) * SMCMAX)) THEN
!--------------------------------------------------------------------------------------------------
! FROZEN GROUND VERSION 
! SMC STATES REPLACED BY SH2O STATES IN SRT SUBR.
! SH2O & SICE STATES INCLUDED IN SSTEP SUBR.
! FROZEN GROUND CORRECTION FACTOR, FRZFACT, ADDED ALL WATER BALANCE CALCULATIONS USING 
! UNFROZEN WATER
!--------------------------------------------------------------------------------------------------
        CALL SRT (RHSTT  , RUNOFF , EDIR   , ET     , SH2O   , SH2O   , NSOIL  , PCPDRP , ZSOIL  ,&
    &             DWSAT  , DKSAT  , SMCMAX , B      , RUNOFF1, RUNOFF2, DT     , SMCWLT , SLOPE  ,&
    &             KDT    , FRZFACT, SICE)
!                 
        CALL SSTEP (SH2OFG , SH2O   , DUMMY  , RHSTT  , RHSCT  , DT     , NSOIL  , SMCMAX ,       &
    &               CMCMAX , RUNOFF3, ZSOIL  , SMC    , SICE)
!                 
        DO K=1,NSOIL
            SH2OA(K) = (SH2O(K) + SH2OFG(K)) * 0.5
        END DO
!                
        CALL SRT (RHSTT  , RUNOFF , EDIR   , ET     , SH2O   , SH2OA  , NSOIL  , PCPDRP , ZSOIL  ,&
    &             DWSAT  , DKSAT  , SMCMAX , B      , RUNOFF1, RUNOFF2, DT     , SMCWLT , SLOPE  ,&
    &             KDT    , FRZFACT, SICE)
!                 
        CALL SSTEP (SH2O   , SH2O   , CMC    , RHSTT  , RHSCT  , DT     , NSOIL  , SMCMAX ,       &
    &               CMCMAX , RUNOFF3, ZSOIL  , SMC    , SICE)
!                
    ELSE
!                 
        CALL SRT (RHSTT  , RUNOFF , EDIR   , ET     , SH2O   , SH2O   , NSOIL  , PCPDRP , ZSOIL  ,&
    &             DWSAT  , DKSAT  , SMCMAX , B      , RUNOFF1, RUNOFF2, DT     , SMCWLT , SLOPE  ,&
    &             KDT    , FRZFACT, SICE)
!                
        CALL SSTEP (SH2O   , SH2O   , CMC    , RHSTT  , RHSCT  , DT     , NSOIL  , SMCMAX ,       &
    &               CMCMAX , RUNOFF3, ZSOIL  , SMC    , SICE)
!                 
    END IF
!          
    RUNOF = RUNOFF
!
    RETURN
!
    END SUBROUTINE SMFLX
!
!
!
    FUNCTION SNKSRC(TUP, TM, TDN, SMC, SH2O, ZSOIL, NSOIL, SMCMAX, PSISAT, B, DT, K, QTOT)
!>--------------------------------------------------------------------------------------------------
!> FUNCTION SNKSRC
!>
!> FUNCTION: SNKSRC - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> TO CALCULATE SINK/SOURCE TERM OF THE TERMAL DIFFUSION EQUATION. (SH2O) IS AVAILABLE 
!> LIQUED WATER.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> TUP    -
!> TM     -
!> TDN    -
!> SMC    -
!> ZSOIL  -
!> NSOIL  -
!> SMCMAX -
!> PSISAT -
!> B      -
!> DT     -
!> K      - 
!> QTOT   -
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> SH2O   -
!>
!> USE MODULES: F77KINDS
!>         
!>
!> DRIVER     : SFLX
!>------------------------------------------------------------------------------------------------- 
    USE F77KINDS
!       
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND), PARAMETER :: HLICE = 3.3350E5
    REAL   (KIND=R4KIND), PARAMETER :: DH2O  = 1.0000E3
    REAL   (KIND=R4KIND), PARAMETER :: T0    = 2.7315E2 

    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & K       , NSOIL
!       
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & TUP     , TM      , TDN     , SMC     , SMCMAX  , PSISAT  , B       , DT      , QTOT
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & SH2O    
! 
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & ZSOIL
!  
    REAL   (KIND=R4KIND)                                                                        ::&
    & DF      , DFH2O   , DFICE   , DZ      , DZH     , FREE    , FRH2O   , SNKSRC  , TAVG    ,   &
    & TZ      , X0      , XDN     , XH2O    , XUP      
!      
    IF (K == 1) THEN
        DZ = -ZSOIL(1)
    ELSE
        DZ = ZSOIL(K-1) - ZSOIL(K)
    END IF
!------------------------------------------------------  
! CALCULATE POTENTIAL REDUCTION OF LIQUED WATER CONTENT
!------------------------------------------------------    
    XH2O = QTOT * DT / (DH2O * HLICE * DZ) + SH2O
!-------------------------------------------------------------------------------------------------  
! ESTIMATE UNFROZEN WATER AT TEMPERATURE TAVG, AND CHECK IF CALCULATED WATER CONTENT IS REASONABLE
!-------------------------------------------------------------------------------------------------  
!            
!--------------------------------------------------------------------------------------------------
! NEW CALCULATION OF AVERAGE TEMPERATURE (TAVG) IN FREEZING/THAWING LAYER USING UP, DOWN, AND 
! MIDDLE LAYER TEMPERATURES (TUP,TDN,TM)
!--------------------------------------------------------------------------------------------------   
    DZH = DZ * 0.5
!
    IF (TUP < T0) THEN
!
        IF (TM < T0) THEN
!
            IF (TDN < T0) THEN
!------------------ 
! TUP, TM, TDN < T0
!------------------
                TAVG = (TUP + 2.0 * TM + TDN) / 4.0
!                            
            ELSE
!--------------------------
! TUP & TM < T0,  TDN >= T0 
!--------------------------
                X0   = (T0 - TM) * DZH / (TDN - TM)
                TAVG = 0.5 * (TUP * DZH + TM * (DZH + X0) + T0 * (2. * DZH - X0)) / DZ
!                                       
            END IF
!
        ELSE
!                    
            IF (TDN < T0) THEN
!-----------------------------
! TUP < T0, TM >= T0, TDN < T0 
!-----------------------------
                 XUP = (T0 - TUP) * DZH / (TM - TUP)
                 XDN = DZH - (T0 - TM) * DZH / (TDN - TM)
                TAVG = 0.5 * (TUP * XUP + T0 * (2. * DZ - XUP - XDN) + TDN * XDN) / DZ
!
            ELSE
!------------------------------
! TUP < T0, TM >= T0, TDN >= T0 
!------------------------------
                XUP  = (T0 - TUP) * DZH / (TM - TUP)
                TAVG = 0.5 * (TUP * XUP + T0 * (2. * DZ - XUP)) / DZ
!                                      
            END IF
!                    
        END IF
!
    ELSE
!
        IF (TM < T0) THEN
!
            IF (TDN < T0) THEN
!-----------------------------
! TUP >= T0, TM < T0, TDN < T0 
!-----------------------------
                XUP  = DZH - (T0 - TUP) * DZH / (TM - TUP)
                TAVG = 0.5 * (T0 * (DZ - XUP) + TM * (DZH + XUP) + TDN * DZH) / DZ
!                                      
            ELSE
!------------------------------
! TUP >= T0, TM < T0, TDN >= T0 
!------------------------------
                XUP  = DZH - (T0 - TUP) * DZH / (TM - TUP)
                XDN  = (T0 - TM) * DZH / (TDN - TM)
                TAVG = 0.5 * (T0 * (2. * DZ - XUP - XDN) + TM * (XUP + XDN)) / DZ
!                                                   
            END IF
!
        ELSE
!
            IF (TDN < T0) THEN
!------------------------------
! TUP >= T0, TM >= T0, TDN < T0 
!------------------------------
                XDN  = DZH -  (T0 - TM) * DZH / (TDN - TM)
                TAVG = (T0 * (DZ - XDN) + 0.5 * (T0 + TDN) * XDN) / DZ
!                                 
            ELSE
!-------------------------------
! TUP >= T0, TM >= T0, TDN >= T0 
!-------------------------------
                TAVG = (TUP + 2.0 * TM + TDN) / 4.0
!                                      
            END IF
!
        END IF
!
    END IF
!
    FREE = FRH2O(TAVG, SMC, SH2O, SMCMAX, B, PSISAT)
!
    IF (XH2O < SH2O .AND. XH2O < FREE) THEN
        IF (FREE > SH2O) THEN
            XH2O = SH2O
        ELSE
            XH2O = FREE
        END IF
    END IF
!                  
    IF (XH2O > SH2O .AND. XH2O > FREE)  THEN
        IF (FREE < SH2O) THEN
            XH2O = SH2O
        ELSE
            XH2O = FREE
        END IF
    END IF
!
    IF (XH2O < 0. ) XH2O = 0.
    IF (XH2O > SMC) XH2O = SMC
!--------------------------------------------------------------
! CALCULATE SINK/SOURCE TERM AND REPLACE PREVIOUS WATER CONTENT
!--------------------------------------------------------------       
    SNKSRC = -DH2O * HLICE * DZ * (XH2O - SH2O) / DT
!
    SH2O = XH2O
!          
    77 RETURN
!
    END FUNCTION SNKSRC
!
!
!
    SUBROUTINE SNOPAC(ETP    , ETA    , PRCP   , PRCP1  , SNOWNG , SMC    , SMCMAX , SMCWLT ,     &
    &                 SMCREF , SMCDRY , CMC    , CMCMAX , NSOIL  , DT     , SBETA  , Q1     ,     &
    &                 DF1    , Q2     , T1     , SFCTMP , T24    , TH2    , F      , F1     ,     &
    &                 S      , STC    , EPSCA  , SFCPRS , B      , PC     , RCH    , RR     ,     &
    &                 CFACTR , SNCOVER, ESD    , SNDENS , SNOWH  , SH2O   , SLOPE  , KDT    ,     &
    &                 FRZFACT, PSISAT , SNUP   , ZSOIL  , DWSAT  , DKSAT  , TBOT   , ZBOT   ,     &
    &                 SHDFAC , RUNOFF1, RUNOFF2, RUNOFF3, EDIR1  , EC1    , ETT1   , NROOT  ,     &
    &                 SNMAX  , ICE    , RTDIS  , QUARTZ , FXEXP  , CSOIL)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SNOPAC
!>
!> SUBROUTINE: SNOPAC - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> TO CALCULATE SOIL MOISTURE AND HEAT FLUX VALUES & UPDATE SOIL MOISTURE CONTENT AND SOIL
!> HEAT CONTENT VALUES FOR THE CASE WHEN A SNOW PACK IS PRESENT.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NSOIL    -
!> NROOT    -
!> ICE      -
!> SNOWNG   -
!> PRCP     -
!> SMCMAX   -
!> SMCWLT   -
!> SMCREF   -
!> SMCDRY   -
!> CMC      -
!> CMCMAX   -
!> DT       -
!> SBETA    -
!> Q1       -
!> DF1      -
!> Q2       -
!> SFCTMP   -
!> T24      -
!> TH2      -
!> F        -
!> F1       -
!> EPSCA    -
!> SFCPRS   -
!> B        -
!> PC       -
!> RCH      -
!> RR       -
!> CFACTR   -
!> SNCOVER  -
!> SLOPE    -
!> KDT      -
!> FRZFACT  -
!> PSISAT   -
!> SNUP     -
!> DWSAT    - 
!> DKSAT    -
!> TBOT     -
!> ZBOT     -
!> SHDFAC   -
!> RUNOFF1  -
!> RUNOFF2  -
!> RUNOFF3  -
!> EDIR1    -
!> EC1      -
!> ETT1     -
!> QUARTZ   -
!> FXEXP    -
!> CSOIL    -
!> SMC      -
!> STC      -
!> SH2O     -
!> ZSOIL    -
!> RTDIS    -
!>
!> OUTPUT ARGUMENT LIST:
!> ETA      -
!> S        -
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> ETP      -
!> PRCP1    - 
!> T1       -
!> ESD      -
!> SNDENS   -
!> SNOWH    -
!>
!> USE MODULES: F77KINDS
!>              RITE
!>
!> DRIVER     : SFLX
!>
!> CALLS      : SHFLX
!>              SMFLX
!>              SNOWPACK 
!>------------------------------------------------------------------------------------------------- 
    USE F77KINDS 
    USE RITE    , RUNOXX3 => RUNOFF3
! 
     IMPLICIT NONE
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & NSOIL   , NROOT   , ICE
!
    LOGICAL(KIND=L4KIND)                                                  , INTENT(IN)          ::&
    & SNOWNG
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & ETP     , PRCP1   , T1      , ESD     , SNDENS  , SNOWH
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & ETA     , S
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & PRCP    , SMCMAX  , SMCWLT  , SMCREF  , SMCDRY  , CMC     , CMCMAX  , DT      , SBETA   ,   &
    & Q1      , DF1     , Q2      , SFCTMP  , T24     , TH2     , F       , F1      , EPSCA   ,   &
    & SFCPRS  , B       , PC      , RCH     , RR      , CFACTR  , SNCOVER , SLOPE   , KDT     ,   &
    & FRZFACT , PSISAT  , SNUP    , DWSAT   , DKSAT   , TBOT    , ZBOT    , SHDFAC  , RUNOFF1 ,   &
    & RUNOFF2 , RUNOFF3 , EDIR1   , EC1     , ETT1    , QUARTZ  , FXEXP   , CSOIL
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & SMC     , STC     , SH2O    , ZSOIL   , RTDIS
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & DENOM   , DSOIL   , DTOT    , EXPSNO  , EXPSOI  , ETA1    , ETP1    , ETP2    , EX      ,   &
    & EXPFAC  , RUNOFF  , S1      , SNMAX   , T11     , T12     , T12A    , T12B    , YY      ,   &
    & ZZ1     , SALP    , SNCOND  , RSNOW   , QSAT    , ETP3    , SEH     , T14     , CSNOW
!         
    REAL   (KIND=R4KIND), PARAMETER :: CP     = 1004.5
    REAL   (KIND=R4KIND), PARAMETER :: CPH2O  =    4.218E+3
    REAL   (KIND=R4KIND), PARAMETER :: CPICE  =    2.106E+3
    REAL   (KIND=R4KIND), PARAMETER :: LSUBF  =    3.335E+5
    REAL   (KIND=R4KIND), PARAMETER :: LSUBC  =    2.501000E+6
    REAL   (KIND=R4KIND), PARAMETER :: LSUBS  =    2.83E+6
    REAL   (KIND=R4KIND), PARAMETER :: SIGMA  =    5.67E-8
    REAL   (KIND=R4KIND), PARAMETER :: TFREEZ =  273.15
!--------------------------------------------------------------------------------------------------
! EXECUTABLE CODE BEGINS HERE.
! CONVERT POTENTIAL EVAP (ETP) FROM KG M-2 S-1 TO M S-1 AND THEN TO AN AMOUNT (M) GIVEN TIMESTEP 
! (DT) AND CALL IT AN EFFECTIVE SNOWPACK REDUCTION AMOUNT, ETP2 (M). THIS IS THE AMOUNT THE 
! SNOWPACK WOULD BE REDUCED DUE TO EVAPORATION FROM THE SNOW SFC DURING THE TIMESTEP. 
! EVAPORATION WILL PROCEED AT THE POTENTIAL RATE UNLESS THE SNOW DEPTH IS LESS THAN THE EXPECTED 
! SNOWPACK REDUCTION.
! IF SEAICE (ICE=1), BETA REMAINS=1.
!--------------------------------------------------------------------------------------------------
    PRCP1 = PRCP1 * 0.001
!
    ETP2 = ETP * 0.001 * DT
    BETA = 1.0
!
    IF (ICE /= 1) THEN
        IF (ESD < ETP2) THEN
            BETA = ESD / ETP2
        END IF
    END IF
! ----------------------------------------------------------- 
! IF ETP<0 (DOWNWARD) THEN DEWFALL (=FROSTFALL IN THIS CASE).
! ----------------------------------------------------------- 
    DEW = 0.0
    IF (ETP < 0.0) THEN
        DEW = -ETP * 0.001
    END IF
!--------------------------------------------------------------------------------------------------
! IF PRECIP IS FALLING, CALCULATE HEAT FLUX FROM SNOW SFC TO NEWLY ACCUMULATING PRECIP. NOTE THAT 
!THIS REFLECTS THE FLUX APPROPRIATE FOR THE NOT-YET-UPDATED SKIN TEMPERATURE (T1). 
!ASSUMES TEMPERATURE OF THE SNOWFALL STRIKING THE GOUND IS =SFCTMP (LOWEST MODEL LEVEL AIR TEMP).
!--------------------------------------------------------------------------------------------------
    FLX1 = 0.0
    IF (SNOWNG) THEN
        FLX1 = CPICE * PRCP * (T1 - SFCTMP)
    ELSE
        IF (PRCP > 0.0) FLX1 = CPH2O * PRCP * (T1 - SFCTMP)
    END IF
!
    DSOIL = -(0.5 * ZSOIL(1))
    DTOT  =  SNOWH + DSOIL
!--------------------------------------------------------------------------------------------------
! CALCULATE AN 'EFFECTIVE SNOW-GRND SFC TEMP' (T12) BASED ON HEAT FLUXES BETWEEN THE SNOW PACK AND
! THE SOIL AND ON NET RADIATION.
! INCLUDE FLX1 (PRECIP-SNOW SFC) AND FLX2 (FREEZING RAIN LATENT HEAT) FLUXES. FLX1 FROM ABOVE, 
! FLX2 BROUGHT IN VIA COMMOM BLOCK RITE FLX2 REFLECTS FREEZING RAIN LATENT HEAT FLUX USING T1 
! CALCULATED IN PENMAN.
!--------------------------------------------------------------------------------------------------
    DENOM = 1.0 + DF1 / (DTOT * RR * RCH)
    T12A  = ((F - FLX1 - FLX2 - SIGMA * T24) / RCH + TH2 - SFCTMP - BETA * EPSCA) / RR
    T12B  = DF1 * STC(1) / (DTOT * RR * RCH)
    T12   = (SFCTMP + T12A + T12B) / DENOM
!--------------------------------------------------------------------------------------------------
! IF THE 'EFFECTIVE SNOW-GRND SFC TEMP' IS AT OR BELOW FREEZING, NO SNOW MELT WILL OCCUR. SET THE 
! SKIN TEMP TO THIS EFFECTIVE TEMP AND SET THE EFFECTIVE PRECIP TO ZERO.
!--------------------------------------------------------------------------------------------------
    IF (T12 <= TFREEZ) THEN
        ESD = MAX(0.0, ESD - ETP2)
!---------------------- 
! GG UPDATE SNOW DEPTH.
!----------------------    
        SNOWH = ESD / SNDENS
           T1 = T12
!----------------------------------------------------------
! UPDATE SOIL HEAT FLUX (S) USING NEW SKIN TEMPERATURE (T1)
!----------------------------------------------------------
        S     = DF1 * (T1 - STC(1)) / (DTOT)
        FLX3  = 0.0
        EX    = 0.0
        SNMAX = 0.0
!--------------------------------------------------------------------------------------------------
! IF THE 'EFFECTIVE SNOW-GRND SFC TEMP' IS ABOVE FREEZING, SNOW MELT WILL OCCUR.  CALL THE SNOW 
! MELT RATE,EX AND AMT, SNMAX. REVISE THE EFFECTIVE SNOW DEPTH. REVISE THE SKIN TEMP BECAUSE IT 
! WOULD HAVE CHGD DUE TO THE LATENT HEAT RELEASED BY THE MELTING. CALC THE LATENT HEAT RELEASED, 
! FLX3. SET THE EFFECTIVE PRECIP, PRCP1 TO THE SNOW MELT RATE, EX FOR USE IN SMFLX.
! ADJUSTMENT TO T1 TO ACCOUNT FOR SNOW PATCHES.
!--------------------------------------------------------------------------------------------------
    ELSE
        T1   = TFREEZ * SNCOVER + T12 * (1.0 - SNCOVER)
        QSAT = (0.622 * 6.11E2) / (SFCPRS - 0.378 * 6.11E2)
        ETP  = RCH * (QSAT - Q2) / CP
        ETP2 = ETP * 0.001 * DT
        BETA = 1.0     	
!----------------------------------------------------------------------- 
! IF POTENTIAL EVAP (SUBLIMATION) GREATER THAN DEPTH OF SNOWPACK. BETA<1
!----------------------------------------------------------------------- 
        IF (ESD <= ETP2) THEN
            BETA = ESD / ETP2
            ESD  = 0.0
!----------------------------------------------- 
! GG SNOW PACK HAS SUBLIMATED, SET DEPTH TO ZERO
!-----------------------------------------------
            SNOWH = 0.0
!
            SNMAX = 0.0
            EX    = 0.0
!----------------------------------------------------------
! UPDATE SOIL HEAT FLUX (S) USING NEW SKIN TEMPERATURE (T1)
!----------------------------------------------------------
            S = DF1 * (T1 - STC(1)) / (DTOT)             
!--------------------------------------------------------------------------------------------------
! POTENTIAL EVAP (SUBLIMATION) LESS THAN DEPTH OF SNOWPACK, BETA=1.
! SNOWPACK (ESD) REDUCED BY POT EVAP RATE ETP3 (CONVERT TO FLUX) UPDATE SOIL HEAT FLUX BECAUSE T1 
! PREVIOUSLY CHANGED.
! SNOWMELT REDUCTION DEPENDING ON SNOW COVER IF SNOW COVER LESS THAN 5% NO SNOWMELT REDUCTION
!--------------------------------------------------------------------------------------------------
        ELSE
            ESD = ESD - ETP2
!------------------------------------------------------------ 
! GG      SNOW PACK REDUCED BY SUBLIMATION, REDUCE SNOW DEPTH
!------------------------------------------------------------ 
            SNOWH = ESD / SNDENS
!
            ETP3 = ETP * LSUBC
            S    = DF1 * (T1 - STC(1)) / (DTOT)
            SEH  = RCH * (T1 - TH2)
            T14  =  T1 * T1
            T14  = T14 * T14
            FLX3 = F - FLX1 - FLX2 - SIGMA * T14 - S - SEH - ETP3
!
            IF (FLX3 <= 0.0) FLX3 = 0.0
            EX = FLX3 * 0.001 / LSUBF
!-------------------------------------------------------------
! DOES BELOW FAIL TO MATCH THE MELT WATER WITH THE MELT ENERGY
!-------------------------------------------------------------
            IF (SNCOVER > 0.05) EX = EX * SNCOVER
            SNMAX = EX * DT
        END IF
!--------------------------------------------------------------------------------------------------
! THE 1.E-6 VALUE REPRESENTS A SNOWPACK DEPTH THRESHOLD VALUE (0.1 MM) BELOW WHICH WE CHOOSE NOT TO 
! RETAIN ANY SNOWPACK, AND INSTEAD INCLUDE IT IN SNOWMELT.
!--------------------------------------------------------------------------------------------------
        IF (SNMAX < ESD-1.E-6) THEN
            ESD = ESD - SNMAX
!------------------------------------------------------- 
! GG      SNOW MELT REDUCED SNOW PACK, REDUCE SNOW DEPTH
!------------------------------------------------------- 
            SNOWH = ESD / SNDENS
        ELSE
            EX    = ESD / DT
            SNMAX = ESD
            ESD   = 0.0
!------------------------------------- 
! GG      SNOW MELT EXCEEDS SNOW DEPTH
!------------------------------------- 
            SNOWH = 0.0
            FLX3  = EX * 1000.0 * LSUBF
        END IF
!
        PRCP1 = PRCP1 + EX
!
    END IF
!--------------------------------------------------------------------------------------------------
! SET THE EFFECTIVE POTNL EVAPOTRANSP (ETP1) TO ZERO SINCE SNOW CASE SO SURFACE EVAP NOT CALCULATED
! FROM EDIR, EC, OR ETT IN SMFLX (BELOW).
! IF SEAICE (ICE=1) SKIP CALL TO SMFLX.
! SMFLX RETURNS SOIL MOISTURE VALUES AND PRELIMINARY VALUES OF EVAPOTRANSPIRATION. IN THIS, THE 
! SNOW PACK CASE, THE PRELIM VALUES (ETA1) ARE NOT USED IN SUBSEQUENT CALCULATION OF EVAP.
! NEW STATES ADDED: SH2O, AND FROZEN GROUND CORRECTION FACTOR EVAP EQUALS POTENTIAL EVAP UNLESS
! BETA<1.
!--------------------------------------------------------------------------------------------------
    ETP1 = 0.0
!
    IF (ICE /= 1) THEN
        CALL SMFLX (ETA1   , SMC    , NSOIL  , CMC    , ETP1   , DT     , PRCP1  , ZSOIL  ,       &
    &               SH2O   , SLOPE  , KDT    , FRZFACT, SMCMAX , B      , PC     , SMCWLT ,       &
    &               DKSAT  , DWSAT  , SMCREF , SHDFAC , CMCMAX , SMCDRY , CFACTR , RUNOFF1,       &
    &               RUNOFF2, RUNOFF3, EDIR1  , EC1    , ETT1   , SFCTMP , Q2     , NROOT  ,       &
    &               RTDIS  , FXEXP)
    END IF
!
    ETA = BETA * ETP
!--------------------------------------------------------------------------------------------------
! THE 'ADJUSTED TOP SOIL LYR TEMP' (YY) AND THE ADJUSTED SOIL HEAT FLUX (ZZ1) ARE SET TO THE TOP 
! SOIL LYR TEMP, AND 1, RESPECTIVELY.
! THESE ARE CLOSE-ENOUGH APPROXIMATIONS BECAUSE THE SFC HEAT FLUX TO BE COMPUTED IN SHFLX WILL 
! EFFECTIVELY BE THE FLUX AT THE SNOW TOP SURFACE. 
! T11 IS A DUMMY ARGUEMENT SINCE WE WILL NOT USE ITS VALUE AS REVISED BY SHFLX.
!--------------------------------------------------------------------------------------------------
    ZZ1 = 1.0
    YY  = STC(1) - 0.5 * S * ZSOIL(1) * ZZ1 / DF1
    T11 = T1
!--------------------------------------------------------------------------------------------------
! SHFLX WILL CALC/UPDATE THE SOIL TEMPS.  NOTE:  THE SUB-SFC HEAT FLUX (S1) AND THE SKIN TEMP (T11) 
! OUTPUT FROM THIS SHFLX CALL ARE NOT USED IN ANY SUBSEQUENT CALCULATIONS. 
! RATHER, THEY ARE DUMMY VARIABLES HERE IN THE SNOPAC CASE, SINCE THE SKIN TEMP AND SUB-SFC HEAT 
! FLUX ARE UPDATED INSTEAD NEAR THE BEGINNING OF THE CALL TO SNOPAC.
!--------------------------------------------------------------------------------------------------
    CALL SHFLX(S1    ,  STC  , SMC , SMCMAX, NSOIL, T11, DT , YY    ,ZZ1   , ZSOIL, TBOT, ZBOT,   &
    &          SMCWLT, PSISAT, SH2O, B     , F1   , DF1, ICE, QUARTZ, CSOIL)
!----------------------------------------------------------------------- 
! SNOW DEPTH AND DENSITY ADJUSTMENT BASED ON SNOW COMPACTION.
! YY IS ASSUMED TO BE THE SOIL TEMPERTURE AT THE TOP OF THE SOIL COLUMN.
! ----------------------------------------------------------------------
    IF (ESD > 0.) THEN
        CALL SNOWPACK(ESD, DT, SNOWH, SNDENS, T1, YY)
    ELSE
           ESD = 0.
         SNOWH = 0.
        SNDENS = 0.
        SNCOND = 1.
    END IF
!
    RETURN
!
    END SUBROUTINE SNOPAC
!
!
!
    SUBROUTINE SNOWPACK(W, DTS, HC, DS, TSNOW, TSOIL)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SNOWPACK
!>
!> SUBROUTINE: SNOWPACK - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> SUBROUTINE TO CALCULATE COMPACTION OF SNOWPACK UNDER CONDITIONS OF INCREASING SNOW DENSITY, AS 
!> OBTAINED FROM AN APPROXIMATE SOLUTION OF E. ANDERSONS DIFFERENTIAL EQUATION (3.29), NOAA 
!> TECHNICAL REPORT NWS 19.
!>
!> PROGRAM HISTORY LOG:
!> 95-03-25  VICTOR KOREN        - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> W     - WATER EQUIVALENT OF SNOW, IN M
!> DTS   - TIME STEP, IN SEC 
!> TSNOW - SNOW SURFACE TEMPERATURE, K 
!> TSOIL - SOIL SURFACE TEMPERATURE, K
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> HC    - SNOW DEPTH, IN M 
!> DS    - SNOW DENSITY, IN G/CM3
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : SFLX
!> 
!> CALLS      : ----- 
!>-------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IPOL    , J
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & HC      , DS     
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & W       , DTS     , TSNOW   , TSOIL   
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & H       , WX      , DT      , TSNOWX  , TSOILX  , TAVG    , B       , DSX     , DW      ,   &
    & PEXP    , WXX
!
    REAL   (KIND=R4KIND), PARAMETER :: C1  =  0.01
    REAL   (KIND=R4KIND), PARAMETER :: C2  = 21.0
!--------------------------------- 
! CONVERSION INTO SIMULATION UNITS 
!--------------------------------- 
    H  = HC * 100.
    WX = W  * 100.
!
    DT = DTS / 3600.
!
    TSNOWX = TSNOW - 273.15
    TSOILX = TSOIL - 273.15
!------------------------------------------------ 
! CALCULATING OF AVERAGE TEMPERATURE OF SNOW PACK 
!------------------------------------------------
    TAVG = 0.5 * (TSNOWX + TSOILX)
!------------------------------------------------------------------------
! CALCULATING OF SNOW DEPTH AND DENSITY AS A RESULT OF COMPACTION
! DS = DS0 * (EXP(B * W) - 1.) / (B * W)
! B  = DT  * C1 * EXP(0.08 * TAVG - C2 * DS0)
! NOTE: B*W IN DS EQN ABOVE HAS TO BE CAREFULLY TREATED NUMERICALLY BELOW
! C1 IS THE FRACTIONAL INCREASE IN DENSITY (1/(CM*HR))
! C2 IS A CONSTANT (CM3/G) KOJIMA ESTIMATED AS 21 CMS/G
!------------------------------------------------------------------------
    IF (WX > 1.E-2) THEN
        WXX = WX
    ELSE
        WXX = 1.E-2
    END IF
!
    B = DT * C1 * EXP(0.08 * TAVG - C2 * DS)
!--------------------------------------------------------------------------------------------------
! THE FUNCTION OF THE FORM (E**X-1)/X IMBEDDED IN ABOVE EXPRESSION FOR DSX WAS CAUSING NUMERICAL 
! DIFFICULTIES WHEN THE DENOMINATOR "X" (I.E. B*WX) BECAME ZERO OR APPROACHED ZERO (DESPITE THE 
! FACT THAT THE ANALYTICAL FUNCTION (E**X-1)/X HAS A WELL DEFINED LIMIT AS "X" APPROACHES ZERO), 
! HENCE BELOW WE REPLACE THE (E**X-1)/X EXPRESSION WITH AN EQUIVALENT, NUMERICALLY WELL-BEHAVED 
! POLYNOMIAL EXPANSION.
!
! NUMBER OF TERMS OF POLYNOMIAL EXPANSION, AND HENCE ITS ACCURACY, IS GOVERNED BY ITERATION LIMIT 
! "IPOL".
! IPOL GREATER THAN 9 ONLY MAKES A DIFFERENCE ON DOUBLE PRECISION (RELATIVE ERRORS GIVEN IN 
! PERCENT %).
! IPOL=9, FOR REL.ERROR <~ 1.6 E-6 % (8 SIGNIFICANT DIGITS)
! IPOL=8, FOR REL.ERROR <~ 1.8 E-5 % (7 SIGNIFICANT DIGITS)
! IPOL=7, FOR REL.ERROR <~ 1.8 E-4 % ..
!--------------------------------------------------------------------------------------------------
    IPOL = 4
    PEXP = 0.
!
    DO J=IPOL,1,-1
        PEXP = (1. + PEXP) * B * WXX / REAL(J+1)
    END DO
!
    PEXP = PEXP + 1.
!
    DSX = DS * (PEXP)
!---------------------------------------- 
! ABOVE LINE ENDS POLYNOMIAL SUBSTITUTION
!---------------------------------------- 
    IF (DSX > 0.40) DSX = 0.40
!-------------------------------------------------------------------
! MBEK - APRIL 2001
! SET LOWER LIMIT ON SNOW DENSITY, RATHER THAN JUST PREVIOUS VALUE.
! IF(DSX .LT. 0.05) DSX=DS
!-------------------------------------------------------------------
    IF (DSX < 0.05) DSX = 0.05
!
    DS = DSX
!--------------------------------------------------------------------------------------------------
! UPDATE OF SNOW DEPTH AND DENSITY DEPENDING ON LIQUID WATER DURING SNOWMELT. ASSUMED THAT 13% OF 
! LIQUID WATER CAN BE STORED IN SNOW PER DAY DURING SNOWMELT TILL SNOW DENSITY 0.40
!--------------------------------------------------------------------------------------------------
    IF (TSNOWX >= 0.) THEN
        DW = 0.13 * DT / 24.
        DS = DS * (1. - DW) + DW
!
        IF (DS > 0.40) DS = 0.40
    END IF
!-----------------------------------------------------------------------
! CALCULATE SNOW DEPTH (CM) FROM SNOW WATER EQUIVALENT AND SNOW DENSITY.
!-----------------------------------------------------------------------
    H = WX / DS
!---------------------------------- 
! CHANGE SNOW DEPTH UNITS TO METERS
!---------------------------------- 
    HC = H * 0.01
!
    RETURN
!
    END SUBROUTINE SNOWPACK
!
!
!
    SUBROUTINE SNOW_NEW(T, P, HC, DS)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SNOW_NEW
!>
!> SUBROUTINE: SNOW_NEW - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> SUBROUTINE TO CALCULATING SNOW DEPTH AND DENSITITY TO ACCOUNT FOR THE NEW SNOWFALL.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> T     - AIR TEMPERATURE, K
!> P     - NEW SNOWFALL, M 
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> HC    - SNOW DEPTH, IN M 
!> DS    - SNOW DENSITY, IN G/CM3
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : SFLX
!> 
!> CALLS      : -----
!>-------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE         
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & HC      , DS  
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & T       , P
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & H       , PX      , TX      , DS0     , HNEW    , ESD
!--------------------------------- 
! CONVERSION INTO SIMULATION UNITS
!--------------------------------- 
    H  = HC * 100.
    PX = P  * 100.
    TX = T  - 273.15         
!--------------------------------------------------------------------------------------------------
! CALCULATING NEW SNOWFALL DENSITY DEPENDING ON TEMPERATURE EQUATION FROM GOTTLIB L. 'A GENERAL 
! RUNOFF MODEL FOR SNOWCOVERED' 'AND GLACIERIZED BASIN', 6TH NORDIC HYDROLOGICAL CONFERENCE, 
! VEMADOLEN, SWEDEN, 1980, 172-177PP.
!--------------------------------------------------------------------------------------------------
    IF (TX <= -15.) THEN
        DS0 = 0.05
    ELSE
        DS0 = 0.05 + 0.0017 * (TX + 15.) ** 1.5
    END IF
!----------------------------------------------------- 
! ADJUSTMENT OF SNOW DENSITY DEPENDING ON NEW SNOWFALL
!----------------------------------------------------- 
    HNEW = PX / DS0
    DS   = (H * DS + HNEW * DS0) / (H + HNEW)
    H    = H + HNEW
    HC   = H * 0.01
!
    RETURN
!
    END SUBROUTINE SNOW_NEW
!
!
!
    SUBROUTINE SRT (RHSTT  , RUNOFF , EDIR   , ET     , SH2O   , SH2OA  , NSOIL  , PCPDRP ,       &
    &               ZSOIL  , DWSAT  , DKSAT  , SMCMAX , B      , RUNOFF1, RUNOFF2, DT     ,       &
    &               SMCWLT , SLOPE  , KDT    , FRZX   , SICE)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SRT
!>
!> SUBROUTINE: SRT - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> TO CALCULATE THE RIGHT HAND SIDE OF THE TIME TENDENCY TERM OF THE SOIL WATER DIFFUSION EQUATION.
!> ALSO TO COMPUTE (PREPARE) THE MATRIX COEFFICIENTS FOR THE TRI-DIAGONAL MATRIX OF THE IMPLICIT 
!> TIME SCHEME.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NSOIL   -
!> ET      -
!> SH2O    -
!> SH2OA   -
!> ZSOIL   -
!> SICE    -
!> RUNOFF  -
!> EDIR    -
!> PCPDRP  -
!> DWSAT   -
!> DKSAT   -
!> SMCMAX  -
!> B       -
!> DT      -
!> SMCWLT  -
!> SLOPE   -
!> KDT     -
!> FRZX    -
!>
!> OUTPUT ARGUMENT LIST:
!> RHSTT   - 
!> RUNOFF1 -
!> RUNOFF2 -
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: ABCI
!>              F77KINDS
!>
!> DRIVER     : SFLX
!>
!> CALLS      : WDFCND
!>-------------------------------------------------------------------------------------------------
    USE ABCI
    USE F77KINDS
!
    IMPLICIT NONE
!---------------------------------------------------- 
!   CURRENT LOGIC DOESNT ALLOW CVFRZ BE BIGGER THAN 3
!----------------------------------------------------
    INTEGER(KIND=I4KIND), PARAMETER :: CVFRZ = 3
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & NSOIL   
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(INOUT)       ::&
    & RHSTT
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & ET      , SH2O    , SH2OA   , ZSOIL   , SICE
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & RUNOFF  , EDIR    , PCPDRP  , DWSAT   , DKSAT   , SMCMAX  , B       , DT      , SMCWLT  ,   &
    & SLOPE   , KDT     , FRZX     
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & RUNOFF1 , RUNOFF2  
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IALP1   , IOHINF  , J       , JJ      , K       , KS
!
    REAL   (KIND=R4KIND), DIMENSION(NSOLD)                                                      ::&
    & DMAX
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & DDZ     , DDZ2    , DENOM   , DENOM2  , DSMDZ   , DSMDZ2  , INFMAX  , MXSMC   , MXSMC2  ,   &
    & NUMER   , PDDUM   , SICEMAX , WCND    , WCND2   , WDF     , WDF2    , DT1     , SMCAV   ,   &
    & DICE    , DD      , VAL     , DDT     , PX      , FCR     , ACRT    , SUM     , SSTT    ,   &
    & SLOPX
!--------------------------------------------------------------------------------------------------
! FROZEN GROUND VERSION 
! REFERENCE FROZEN GROUND PARAMETER, CVFRZ, IS A SHAPE PARAMETER OF AREAL DISTRIBUTION FUNCTION OF 
! SOIL ICE CONTENT WHICH EQUALS 1/CV.
! CV IS A COEFFICIENT OF SPATIAL VARIATION OF SOIL ICE CONTENT. BASED ON FIELD DATA CV DEPENDS ON 
! AREAL MEAN OF FROZEN DEPTH, AND IT CLOSE TO CONSTANT = 0.6 IF AREAL MEAN FROZEN DEPTH IS ABOVE 
! 20 CM. THAT IS WHY PARAMETER CVFRZ = 3 (INT{1/0.6*0.6})
!-------------------------------------------------------------------------------------------------- 
!
!------------------------------------------------ 
! DETERMINE RAINFALL INFILTRATION RATE AND RUNOFF
!------------------------------------------------ 
!
!-------------------------------------------------------------- 
! INCLUDE THE INFILTRATION FORMULE FROM SCHAAKE AND KOREN MODEL
!-------------------------------------------------------------- 
!
!-------------------  
! MODIFIED BY Q DUAN
!-------------------
    IOHINF = 1
!------------------------------------------------------------------------------
! LET SICEMAX BE THE GREATEST, IF ANY, FROZEN WATER CONTENT WITHIN SOIL LAYERS.
!------------------------------------------------------------------------------
    SICEMAX = 0.0
!
    DO KS=1,NSOIL
        IF (SICE(KS) > SICEMAX) SICEMAX = SICE(KS)
    END DO
!------------------------------------------------ 
! DETERMINE RAINFALL INFILTRATION RATE AND RUNOFF
!------------------------------------------------
    PDDUM   = PCPDRP
    RUNOFF1 = 0.0
!
    IF (PCPDRP /= 0.0) THEN
!---------------------------------- 
! C++  MODIFIED BY Q. DUAN, 5/16/94
!---------------------------------- 
        DT1    =  DT / 86400.
        SMCAV  =  SMCMAX - SMCWLT
        DMAX(1)= -ZSOIL(1) * SMCAV
!---------------------- 
! FROZEN GROUND VERSION 
!-----------------------    
        DICE = -ZSOIL(1) * SICE(1)
!
        DMAX(1) = DMAX(1) * (1.0 - (SH2OA(1) + SICE(1) - SMCWLT) / SMCAV)
        DD      = DMAX(1)
!
        DO KS=2,NSOIL
!----------------------                       
! FROZEN GROUND VERSION 
!----------------------       
            DICE = DICE + (ZSOIL(KS-1) - ZSOIL(KS)) * SICE(KS)
!
            DMAX(KS) = (ZSOIL(KS-1) - ZSOIL(KS)) * SMCAV
            DMAX(KS) =  DMAX(KS) * (1.0 - (SH2OA(KS) + SICE(KS) - SMCWLT) / SMCAV)
            DD=DD+DMAX(KS)
        END DO
!----------------------------------- 
! IN BELOW, REMOVE THE SQRT IN ABOVE
!-----------------------------------
        VAL = (1. -EXP(-KDT * DT1))
        DDT = DD * VAL
        PX  = PCPDRP * DT
!
        IF (PX < 0.0) PX = 0.0
        INFMAX = (PX * (DDT / (PX + DDT))) / DT
!------------------------------------------------------------                 
! FROZEN GROUND VERSION  
! REDUCTION OF INFILTRATION BASED ON FROZEN GROUND PARAMETERS
!------------------------------------------------------------      
        FCR = 1.
        IF (DICE > 1.E-2) THEN
            ACRT  = CVFRZ * FRZX / DICE
            SUM   = 1.
            IALP1 = CVFRZ - 1
!
            DO J=1,IALP1
                K=1
                DO JJ=J+1,IALP1
                    K = K * JJ
                END DO
                SUM = SUM + (ACRT ** (CVFRZ - J)) / FLOAT(K)
            END DO
            FCR = 1. - EXP(-ACRT) * SUM
        END IF
!
        INFMAX = INFMAX * FCR
!--------------------------------------------------------------------------------------
! CORRECTION OF INFILTRATION LIMITATION
! IF INFMAX .LE. HYDROLIC CONDUCTIVITY ASSIGN INFMAX THE VALUE OF HYDROLIC CONDUCTIVITY
!--------------------------------------------------------------------------------------
        MXSMC = SH2OA(1)
!
        CALL WDFCND(WDF, WCND, MXSMC, SMCMAX, B, DKSAT, DWSAT, SICEMAX)
!
        INFMAX = MAX(INFMAX, WCND)
        INFMAX = MIN(INFMAX, PX)
!
        IF (PCPDRP > INFMAX) THEN
            RUNOFF1 = PCPDRP - INFMAX
            PDDUM   = INFMAX
        END IF
!
    END IF
!--------------------------------------------------------------------------------------------------
! TO AVOID SPURIOUS DRAINAGE BEHAVIOR IDENTIFIED BY P. GRUNMANN, FORMER APPROACH IN LINE BELOW 
! REPLACED WITH NEW APPROACH IN 2ND  LINE.
!--------------------------------------------------------------------------------------------------
    MXSMC = SH2OA(1)
!
    CALL WDFCND (WDF, WCND, MXSMC, SMCMAX, B, DKSAT, DWSAT, SICEMAX)
!-------------------------------------------------------------- 
! CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER
!--------------------------------------------------------------
    DDZ = 1. / (-.5 * ZSOIL(2))
!
    AI(1) =  0.0
    BI(1) =  WDF * DDZ / (-ZSOIL(1))
    CI(1) = -BI(1)
!--------------------------------------------------------------------------------------------------
! CALC RHSTT FOR THE TOP LAYER AFTER CALCING THE VERTICAL SOIL MOISTURE GRADIENT BTWN THE TOP AND 
! NEXT TO TOP LAYERS.
!--------------------------------------------------------------------------------------------------
    DSMDZ    = (SH2O(1) - SH2O(2)) / (-.5 * ZSOIL(2))
    RHSTT(1) = (WDF * DSMDZ + WCND - PDDUM + EDIR + ET(1)) / ZSOIL(1)
    SSTT     =  WDF * DSMDZ + WCND + EDIR         + ET(1)
!---------------- 
! INITIALIZE DDZ2
!---------------- 
    DDZ2 = 0.0
!--------------------------------------------------------------- 
! LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABV PROCESS
!--------------------------------------------------------------- 
    DO K=2,NSOIL
        DENOM2 = (ZSOIL(K-1) - ZSOIL(K))
        IF (K /= NSOIL) THEN
            SLOPX = 1.
!--------------------------------------------------------------------------------------------------       
! AGAIN, TO AVOID SPURIOUS DRAINAGE BEHAVIOR IDENTIFIED BY P. GRUNMANN, FORMER APPROACH IN LINE 
! BELOW REPLACED WITH NEW APPROACH IN 2ND LINE
!--------------------------------------------------------------------------------------------------
            MXSMC2 = SH2OA(K)
!
            CALL WDFCND(WDF2, WCND2, MXSMC2, SMCMAX, B, DKSAT, DWSAT, SICEMAX)
!--------------------------------------------------------- 
! CALC SOME PARTIAL PRODUCTS FOR LATER USE IN CALCNG RHSTT
!--------------------------------------------------------- 
            DENOM  = (ZSOIL(K-1) - ZSOIL(K+1))
            DSMDZ2 = ( SH2O(K  ) -  SH2O(K+1)) / (DENOM * 0.5)
!----------------------------------------------------------- 
! CALC THE MATRIX COEF, CI, AFTER CALCNG ITS PARTIAL PRODUCT
!----------------------------------------------------------- 
            DDZ2  =  2.0 / DENOM
            CI(K) = -WDF2 * DDZ2 / DENOM2
        ELSE
!------------------------------------ 
! SLOPE OF BOTTOM LAYER IS INTRODUCED 
!------------------------------------         
            SLOPX = SLOPE
!------------------------------------------------------------------------------ 
! RETRIEVE THE SOIL WATER DIFFUSIVITY AND HYDRAULIC CONDUCTIVITY FOR THIS LAYER
!------------------------------------------------------------------------------ 
            CALL WDFCND (WDF2, WCND2, SH2OA(NSOIL), SMCMAX, B, DKSAT, DWSAT, SICEMAX)
!----------------------------------------------------- 
! CALC A PARTIAL PRODUCT FOR LATER USE IN CALCNG RHSTT
!----------------------------------------------------- 
            DSMDZ2 = 0.0
!--------------------------- 
! SET MATRIX COEF CI TO ZERO
!--------------------------- 
            CI(K) = 0.0
        END IF
!----------------------------------------------------- 
! CALC RHSTT FOR THIS LAYER AFTER CALCNG ITS NUMERATOR
!----------------------------------------------------- 
        NUMER    = (WDF2 * DSMDZ2) + SLOPX * WCND2 - (WDF * DSMDZ) - WCND + ET(K)
        RHSTT(K) = NUMER / (-DENOM2)
!--------------------------------------------- 
! CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER
!--------------------------------------------- 
        AI(K) = -WDF * DDZ / DENOM2
        BI(K) = -(AI(K) + CI(K))
!--------------------------------------------------------------- 
! RESET VALUES OF WDF, WCND, DSMDZ, AND DDZ FOR LOOP TO NEXT LYR
!--------------------------------------------------------------- 
        IF (K == NSOIL) THEN
!----------------------------- 
! RUNOFF2: GROUND WATER RUNOFF 
!-----------------------------
            RUNOFF2 = SLOPX * WCND2
        END IF
!
        IF (K /= NSOIL) THEN
            WDF   = WDF2
            WCND  = WCND2
            DSMDZ = DSMDZ2
            DDZ   = DDZ2
        END IF
    END DO
!
    RETURN
!
    END SUBROUTINE SRT
!
!
!
    SUBROUTINE SSTEP (SH2OOUT, SH2OIN , CMC    , RHSTT  , RHSCT  , DT     , NSOIL  , SMCMAX ,     &
    &                 CMCMAX , RUNOFF3, ZSOIL  , SMC    , SICE)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SSTEP
!>
!> SUBROUTINE: SSTEP - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> TO CALCULATE/UPDATE THE SOIL MOISTURE CONTENT VALUES AND THE CANOPY MOISTURE CONTENT VALUES.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NSOIL   -
!> SH2OIN  - 
!> ZSOIL   -
!> SICE    -
!> RHSCT   -
!> DT      -
!> SMCMAX  -
!> CMCMAX  -
!>
!> OUTPUT ARGUMENT LIST:
!> SH2OOUT -
!> SMC     -
!> RUNOFF3 -
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> RHSTT   -
!> CMC     -
!>
!> USE MODULES: ABCI
!>              F77KINDS
!>
!> DRIVER     : SFLX
!>
!> CALLS      : ROSR12
!>-------------------------------------------------------------------------------------------------
    USE ABCI
    USE F77KINDS
!
    IMPLICIT NONE     
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & NSOIL   
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(INOUT)       ::&
    & SH2OOUT , SMC
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & SH2OIN  , ZSOIL   , SICE
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(INOUT)       ::&
    & RHSTT
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & CMC
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & RHSCT   , DT      , SMCMAX  , CMCMAX
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & RUNOFF3
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , K       , KK11
!
    REAL   (KIND=R4KIND), DIMENSION(NSOLD)                                                      ::&
    & CIIN
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                                      ::&
    & RHSTTIN
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & RUNOFS  , WPLUS   , DDZ     , STOT    , WFREE   , DPLUS
!------------------------------------------------------------------------------------ 
! CREATE 'AMOUNT' VALUES OF VARIABLES TO BE INPUT TO THE TRI-DIAGONAL MATRIX ROUTINE.
!------------------------------------------------------------------------------------ 
    DO K=1,NSOIL
        RHSTT(K) =   RHSTT(K) * DT
           AI(K) =      AI(K) * DT
           BI(K) = 1. + BI(K) * DT
           CI(K) =      CI(K) * DT
    END DO
!------------------------------------------------------ 
! COPY VALUES FOR INPUT VARIABLES BEFORE CALL TO ROSR12
!------------------------------------------------------ 
    DO K=1,NSOIL
        RHSTTIN(K) = RHSTT(K)
    END DO
!
    DO K=1,NSOLD
        CIIN(K) = CI(K)
    END DO
!--------------------------------------------- 
! CALL ROSR12 TO SOLVE THE TRI-DIAGONAL MATRIX
!--------------------------------------------- 
    CALL ROSR12(CI, AI, BI, CIIN, RHSTTIN, RHSTT, NSOIL)
!--------------------------------------------------------------------------------------------------
! SUM THE PREVIOUS SMC VALUE AND THE MATRIX SOLUTION TO GET A NEW VALUE. MIN ALLOWABLE VALUE OF SMC
! WILL BE 0.02.
!--------------------------------------------------------------------------------------------------
!
!----------------------------------- 
! RUNOFF3: RUNOFF WITHIN SOIL LAYERS
!-----------------------------------
    RUNOFS  = 0.0
    WPLUS   = 0.0
    RUNOFF3 = 0.
!
    DDZ = -ZSOIL(1)
!          
    DO K=1,NSOIL
        IF (K /= 1) DDZ = ZSOIL(K-1) - ZSOIL(K)
!
        SH2OOUT(K) = SH2OIN(K) + CI(K) + WPLUS / DDZ
!
        STOT = SH2OOUT(K) + SICE(K)
!
        IF (STOT > SMCMAX) THEN
            IF (K == 1) THEN
                DDZ = -ZSOIL(1)
            ELSE
                KK11=K-1
                DDZ = -ZSOIL(K) + ZSOIL(KK11)
            END IF
!
            WPLUS = (STOT - SMCMAX) * DDZ
        ELSE
            WPLUS = 0.
        END IF
!
            SMC(K) = MAX (MIN(STOT, SMCMAX), 0.02)
        SH2OOUT(K) = MAX((SMC(K) - SICE(K)), 0.0)
    END DO
!-------------------------------------------------
! V. KOREN   9/01/98 WATER BALANCE CHECKING UPWARD
!-------------------------------------------------
    IF (WPLUS > 0.) THEN
        DO I=NSOIL-1,1,-1
            IF (I == 1) THEN
                DDZ = -ZSOIL(1)
            ELSE
                DDZ = -ZSOIL(I) + ZSOIL(I-1)
            END IF
!
            WFREE = (SMCMAX - SH2OOUT(I) - SICE(I)) * DDZ
            DPLUS = WFREE - WPLUS
!
            IF (DPLUS >= 0.) THEN
                SH2OOUT(I) = SH2OOUT(I) + WPLUS / DDZ
                    SMC(I) = SH2OOUT(I) + SICE(I)
                WPLUS      = 0.              
            ELSE
                SH2OOUT(I) = SH2OOUT(I) + WFREE / DDZ
                    SMC(I) = SH2OOUT(I) + SICE(I)
                WPLUS      = -DPLUS
            END IF
        END DO
!
     30 RUNOFF3 = WPLUS
!
    END IF
!--------------------------------------------------------------------------------------------------
! UPDATE CANOPY WATER CONTENT/INTERCEPTION (CMC). 
! CONVERT RHSCT TO AN 'AMOUNT' VALUE AND ADD TO PREVIOUS CMC VALUE TO GET NEW CMC.
!--------------------------------------------------------------------------------------------------
    CMC = CMC + DT * RHSCT
!
    IF (CMC < 1.E-20) CMC = 0.0
    CMC = MIN(CMC, CMCMAX)
!
    RETURN
!
    END SUBROUTINE SSTEP
!
!
!
    SUBROUTINE TBND(TU, TB, ZSOIL, ZBOT, K, NSOIL, TBND1)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE TBND
!>
!> SUBROUTINE: TBND - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> CALCULATE TEMPERATURE ON THE BOUNDARY OF THE LAYER BY INTERPOLATION OF THE MIDDLE LAYER 
!> TEMPERATURES
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> K        -
!> NSOIL    -
!> ZSOIL    -
!> TU       -
!> TB       -
!> ZBOT     -
!>
!> OUTPUT ARGUMENT LIST:
!> TBND1    -
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: F77KINDS        
!>
!> DRIVER     : SFLX
!>
!> CALLS      : -----
!>-------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND), PARAMETER :: T0 = 273.15
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & K       , NSOIL
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & ZSOIL
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & TU      , TB      , ZBOT
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & TBND1
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ZB      , ZUP
!------------------------------------------------------ 
! USE SURFACE TEMPERATURE ON THE TOP OF THE FIRST LAYER
!------------------------------------------------------      
    IF (K == 1) THEN
        ZUP = 0.
    ELSE
        ZUP = ZSOIL(K-1)
    END IF
!------------------------------------------------------------------------------------------------------- 
! USE DEPTH OF THE CONSTANT BOTTOM TEMPERATURE WHEN INTERPOLATE TEMPERATURE INTO THE LAST LAYER BOUNDARY
!-------------------------------------------------------------------------------------------------------       
    IF (K == NSOIL) THEN
        ZB = 2. * ZBOT - ZSOIL(K  )
    ELSE
        ZB =             ZSOIL(K+1)
    END IF
!------------------------------------------------------------ 
! LINEAR INTERPOLATION BETWEEN THE AVERAGE LAYER TEMPERATURES
!------------------------------------------------------------         
    TBND1 = TU + (TB - TU) * (ZUP - ZSOIL(K)) / (ZUP - ZB)
!          
    RETURN
!
    END SUBROUTINE TBND
!
!
!
    SUBROUTINE TDFCND(DF, SMC, Q, SMCMAX, SH2O)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE TDFCND
!>
!> SUBROUTINE: TDFCND - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> TO CALCULATE THERMAL DIFFUSIVITY AND CONDUCTIVITY OF THE SOIL FOR A GIVEN POINT AND TIME.
!> VERSION: PETERS-LIDARD APPROACH (PETERS-LIDARD ET AL., 1998)
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> SMC    -
!> Q      - QUARTZ CONTENT (SOIL TYPE DEPENDENT)
!> SMCMAX - POROSITY (= SMCMAX)
!> SH2O   -
!>
!> OUTPUT ARGUMENT LIST:
!> DF     - 
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: F77KINDS        
!>
!> DRIVER     : SFLX
!>
!> CALLS      : -----
!>-------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & DF
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & SMC     , Q       , SMCMAX  , SH2O
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & GAMMD   , THKDRY  , AKE     , THKICE  , THKO    , THKQTZ  , THKSAT  , THKS    , THKW    ,   &
    & SATRATIO, XU      , XUNFROZ
!------------------------------------------------------------------ 
! WE NOW GET QUARTZ AS AN INPUT ARGUMENT (SET IN ROUTINE REDPRM):
! ATA QUARTZ /0.82, 0.10, 0.25, 0.60, 0.52, 0.35, 0.60, 0.40, 0.82/
!------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
! IF THE SOIL HAS ANY MOISTURE CONTENT COMPUTE A PARTIAL SUM/PRODUCT OTHERWISE USE A CONSTANT VALUE 
! WHICH WORKS WELL WITH MOST SOILS
!--------------------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------- 
! THKW   - WATER THERMAL CONDUCTIVITY
! THKQTZ - THERMAL CONDUCTIVITY FOR QUARTZ
! THKO   - THERMAL CONDUCTIVITY FOR OTHER SOIL COMPONENTS
! THKS   - THERMAL CONDUCTIVITY FOR THE SOLIDS COMBINED(QUARTZ+OTHER)
! THKICE - ICE THERMAL CONDUCTIVITY
!----------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
! USE AS IN PETERS-LIDARD, 1998 (MODIF. FROM JOHANSEN, 1975).
!
! PABLO GRUNMANN, 08/17/98
! REFS.:
!      FAROUKI, O.T.,1986: THERMAL PROPERTIES OF SOILS. SERIES ON ROCK AND SOIL MECHANICS, VOL. 11,
!      TRANS TECH, 136 PP. JOHANSEN, O. 1975: THERMAL CONDUCTIVITY OF SOILS. PH.D. THESIS, 
!      UNIVERSITY OF TRONDHEIM, PETERS-LIDARD, C. D., ET AL., 1998: THE EFFECT OF SOIL THERMAL 
!      CONDUCTIVITY PARAMETERIZATION ON SURFACE ENERGY FLUXES AND TEMPERATURES. JOURNAL OF THE 
!      ATMOSPHERIC SCIENCES, VOL. 55, PP. 1209-1224.
!--------------------------------------------------------------------------------------------------
    SATRATIO = SMC / SMCMAX
    THKICE   = 2.2
    THKW     = 0.57
    THKO     = 2.0
    THKQTZ   = 7.7
!--------------------- 
!  SOLIDS CONDUCTIVITY
!--------------------- 
    THKS = (THKQTZ ** Q) * (THKO ** (1. - Q))
!--------------------------------------------------------------------- 
! UNFROZEN FRACTION (FROM 1.0, I.E., 100%LIQUID, TO 0.0 (100% FROZEN))
!---------------------------------------------------------------------
    XUNFROZ = (SH2O + 1.E-9)/(SMC + 1.E-9)
!-------------------------------------------------- 
! UNFROZEN VOLUME FOR SATURATION (POROSITY*XUNFROZ)
!--------------------------------------------------
    XU = XUNFROZ * SMCMAX
!------------------------------- 
! SATURATED THERMAL CONDUCTIVITY
!------------------------------- 
    THKSAT = THKS ** (1. - SMCMAX) * THKICE ** (SMCMAX - XU) * THKW ** (XU)
!--------------------- 
! DRY DENSITY IN KG/M3
!---------------------
    GAMMD = (1. - SMCMAX) * 2700.
!--------------------------------------
! DRY THERMAL CONDUCTIVITY IN W.M-1.K-1
!--------------------------------------
    THKDRY = (0.135 * GAMMD + 64.7) / (2700. - 0.947 * GAMMD)
!-----------------------------------------
! RANGE OF VALIDITY FOR THE KERSTEN NUMBER
!-----------------------------------------
    IF (SATRATIO > 0.1) THEN
!-----------------------------------------------------------------
! KERSTEN NUMBER (FINE FORMULA, AT LEAST 5% OF PARTICLES<(2.E-6)M)
!-----------------------------------------------------------------
        IF ((XUNFROZ + 0.0005) < SMC) THEN
!------- 
! FROZEN
!-------
            AKE = SATRATIO
        ELSE
!---------
! UNFROZEN
!---------
            AKE = LOG10(SATRATIO) + 1.0
        END IF
!
    ELSE
!-------------                
! USE K = KDRY
!-------------
        AKE = 0.0
!
    END IF
!---------------------
! THERMAL CONDUCTIVITY
!---------------------
    DF = AKE * (THKSAT - THKDRY) + THKDRY
!
    RETURN
!
    END SUBROUTINE TDFCND
!
!
!
    SUBROUTINE TRANSP (ET    , NSOIL , ETP1  , SMC   , CMC   , ZSOIL , SHDFAC, SMCWLT, CMCMAX,    &
    &                  PC    , CFACTR, SMCREF, SFCTMP, Q2    , NROOT , RTDIS)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE TRANSP
!>
!> SUBROUTINE: TRANSP - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> TO CALCULATE TRANSPIRATION FROM THE VEGTYP FOR THIS PT.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> SMC     -
!> ZSOIL   - 
!> RTDIS   -
!> ETP1    -
!> CMC     -
!> SHDFAC  -
!> SMCWLT  -
!> CMCMAX  -
!> PC      -
!> CFACTR  -
!> SMCREF  -
!> SFCTMP  -
!> Q2      -
!> NSOIL   -
!> NROOT   -
!>
!> OUTPUT ARGUMENT LIST:
!> ET      -
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: F77KINDS        
!>
!> DRIVER     : SFLX
!>
!> CALLS      : -----
!>-------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(INOUT)       ::&
    & ET
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                , INTENT(IN)          ::&
    & SMC     , ZSOIL   , RTDIS
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & ETP1    , CMC     , SHDFAC  , SMCWLT  , CMCMAX  , PC      , CFACTR  , SMCREF  , SFCTMP  ,   &
    & Q2      
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & NSOIL   , NROOT
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , K
! 
    REAL   (KIND=R4KIND)                                                                        ::&
    & ETP1A
!
    REAL   (KIND=R4KIND), DIMENSION(7)                                                          ::&
    & GX
!------------------ 
! REAL PART (NSOIL)
!------------------ 
    REAL   (KIND=R4KIND)                                                                        ::&
    & SGX     , DENOM   , RTX
!----------------------------------------------------- 
! INITIALIZE PLANT TRANSP TO ZERO FOR ALL SOIL LAYERS.
!----------------------------------------------------- 
    DO K=1,NSOIL
        ET(K) = 0.
    END DO
!----------------------------------------  
! CALC AN 'ADJUSTED' POTNTL TRANSPIRATION
!----------------------------------------  
!
!--------------------------------------------------------- 
! IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO
!---------------------------------------------------------
    IF (CMC /= 0.0) THEN
        ETP1A = SHDFAC * PC * ETP1 * (1.0 - (CMC / CMCMAX) ** CFACTR)
    ELSE
        ETP1A = SHDFAC * PC * ETP1
    END IF
!          
    SGX = 0.0
!
    DO I=1,NROOT
        GX(I) = (SMC(I) - SMCWLT) / (SMCREF - SMCWLT)
        GX(I) = MAX(MIN(GX(I), 1.), 0.)
!
        SGX = SGX + GX(I)
    END DO
!
    SGX = SGX / NROOT
!          
    DENOM = 0.
!
    DO I=1,NROOT
        RTX   = RTDIS(I)+ GX(I) - SGX
        GX(I) = GX(I)   * MAX(RTX,0.)
        DENOM = DENOM   + GX(I)
    END DO
!
    IF (DENOM <= 0.0) DENOM = 1.
!          
    DO I=1,NROOT
        ET(I) = ETP1A * GX(I) / DENOM
    END DO
!--------------------------------------------------------------------------------------------------
! ABOVE CODE ASSUMES A VERTICALLY UNIFORM ROOT DISTRIBUTION
!
! CODE BELOW TESTS A VARIABLE ROOT DISTRIBUTION
!
! ET(1)=(ZSOIL(1)/ZSOIL(NROOT))*GX*ETP1A
! ET(1)=(ZSOIL(1)/ZSOIL(NROOT))*ETP1A
!
! USING ROOT DISTRIBUTION AS WEIGHTING FACTOR
! ET(1)=RTDIS(1)*ETP1A
! ET(1)=ETP1A*PART(1)
! 
! LOOP DOWN THRU THE SOIL LAYERS REPEATING THE OPERATION ABOVE, BUT USING THE THICKNESS OF THE SOIL
! LAYER (RATHER THAN THE ABSOLUTE DEPTH OF EACH LAYER) IN THE FINAL CALCULATION.
!--------------------------------------------------------------------------------------------------   
    RETURN
!
    END SUBROUTINE TRANSP
!
!
!
    SUBROUTINE WDFCND(WDF, WCND, SMC, SMCMAX, B, DKSAT, DWSAT, SICEMAX)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE WDFCND
!>
!> SUBROUTINE: WDFCND - ?????
!> PROGRAMMER: ?????
!> ORG: ????? 
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> TO CALCULATE SOIL WATER DIFFUSIVITY AND SOIL HYDRAULIC CONDUCTIVITY.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????               - ORIGINATOR
!> 18-01-15  LUCCI               - MODERNIZATION OF THE CODE, INCLUDING:
!>                                 * F77 TO F90/F95
!>                                 * INDENTATION & UNIFORMIZATION CODE
!>                                 * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                 * DOCUMENTATION WITH DOXYGEN
!>                                 * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> B       -
!> DKSAT   -
!> DWSAT   -
!> SICEMAX -
!>
!> OUTPUT ARGUMENT LIST:
!> WDF     -
!> WCND    -
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> SMC     - 
!> SMCMAX  -
!>
!> USE MODULES: F77KINDS        
!>
!> DRIVER     : SFLX
!>
!> CALLS      : -----
!>------------------------------------------------------------------------------------------------- 
    USE F77KINDS
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & B       , DKSAT   , DWSAT   , SICEMAX
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & WDF     , WCND
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & SMC     , SMCMAX
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & EXPON   , FACTR1  , FACTR2  , VKWGT
!--------------------------------------------------------------  
! CALC THE RATIO OF THE ACTUAL TO THE MAX PSBL SOIL H2O CONTENT
!-------------------------------------------------------------- 
    SMC    = SMC
    SMCMAX = SMCMAX
    FACTR1 = 0.2 / SMCMAX
    FACTR2 = SMC / SMCMAX
!--------------------------------------------------------  
! PREP AN EXPNTL COEF AND CALC THE SOIL WATER DIFFUSIVITY
!--------------------------------------------------------  
    EXPON = B + 2.0
    WDF   = DWSAT * FACTR2 ** EXPON
!--------------------------------------------------------------------------------------------------
! FROZEN SOIL HYDRAULIC DIFFUSIVITY. VERY SENSITIVE TO THE VERTICAL GRADIENT OF UNFROZEN WATER. 
! THE LATTER GRADIENT CAN BECOME VERY EXTREME IN FREEZING/THAWING SITUATIONS, AND GIVEN THE 
! RELATIVELY FEW AND THICK SOIL LAYERS, THIS GRADIENT SUFFERES SERIOUS TRUNCTION ERRORS YIELDING 
! ERRONEOUSLY HIGH VERTICAL TRANSPORTS OF UNFROZEN WATER IN BOTH DIRECTIONS FROM HUGE HYDRAULIC 
! DIFFUSIVITY. THEREFORE, WE FOUND WE HAD TO ARBITRARILY CONSTRAIN WDF
! VERSION D_10CM:........................ FACTR1=0.2/SMCMAX
! WEIGHTED APPROACH...................... PABLO GRUNMANN, 09/28/99.
!--------------------------------------------------------------------------------------------------
    IF (SICEMAX > 0.0)  THEN
        VKWGT = 1. / (1. + (500. * SICEMAX) ** 3.)
        WDF   = VKWGT * WDF + (1. - VKWGT) * DWSAT * FACTR1 ** EXPON
    END IF
!---------------------------------------------------------- 
! RESET THE EXPNTL COEF AND CALC THE HYDRAULIC CONDUCTIVITY
!---------------------------------------------------------- 
    EXPON = (2.0 * B) + 3.0
    WCND  = DKSAT * FACTR2 ** EXPON
!
    RETURN
!
    END SUBROUTINE WDFCND
