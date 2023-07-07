    SUBROUTINE GSMCOLUMN (ARAIN    , ASNOW , DTPH  , I_INDEX, J_INDEX, LSFC     ,                 &
    &                     P_COL    , QI_COL, QR_COL, QV_COL , QW_COL , RIMEF_COL, T_COL,          &
    &                     THICK_COL, WC_COL)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE GSMCOLUMN
!>
!> SUBPROGRAM: GSMCOLUMN - GRID-SCALE MICROPHYSICAL PROCESSES - CONDENSATION AND PRECIPITATION
!> PROGRAMMER: FERRIER
!> ORG: W/NP22
!> DATE: 01-08-??
!>
!> ABSTRACT: 
!> MERGES ORIGINAL GSCOND AND PRECPD SUBROUTINES.
!> CODE HAS BEEN SUBSTANTIALLY STREAMLINED AND RESTRUCTURED.
!> EXCHANGE BETWEEN WATER VAPOR & SMALL CLOUD CONDENSATE IS CALCULATED USING THE ORIGINAL ASAI 
!> (1965, J. JAPAN) ALGORITHM.
!> SEE ALSO REFERENCES TO YAU AND AUSTIN (1979, JAS), RUTLEDGE AND HOBBS (1983, JAS), AND TAO 
!> ET AL. (1989, MWR).  
!> THIS ALGORITHM REPLACES THE SUNDQVIST ET AL. (1989, MWR) PARAMETERIZATION.
!>
!> VERSION OF MICROPHYSICS DESIGNED FOR HIGHER RESOLUTION MESO ETA MODEL
!> (1) REPRESENTS SEDIMENTATION BY PRESERVING A PORTION OF THE PRECIPITATION THROUGH TOP-DOWN 
!>     INTEGRATION FROM CLOUD-TOP. MODIFIED PROCEDURE TO ZHAO AND CARR (1997).
!>
!> (2) MICROPHYSICAL EQUATIONS ARE MODIFIED TO BE LESS SENSITIVE TO TIME STEPS BY USE OF CLAUSIUS-
!>     CLAPEYRON EQUATION TO ACCOUNT FOR CHANGES IN SATURATION MIXING RATIOS IN RESPONSE TO LATENT 
!>     HEATING/COOLING.
!>
!> (3) PREVENT SPURIOUS TEMPERATURE OSCILLATIONS ACROSS 0C DUE TO MICROPHYSICS.
!>
!> (4) USES LOOKUP TABLES FOR: CALCULATING TWO DIFFERENT VENTILATION COEFFICIENTS IN CONDENSATION 
!>     AND DEPOSITION PROCESSES;
!>     ACCRETION OF CLOUD WATER BY PRECIPITATION; PRECIPITATION MASS; 
!>     PRECIPITATION RATE (AND MASS-WEIGHTED PRECIPITATION FALL SPEEDS).
!>
!> (5) ASSUMES TEMPERATURE-DEPENDENT VARIATION IN MEAN DIAMETER OF LARGE ICE 
!>     (HOUZE ET AL., 1979; RYAN ET AL., 1996).
!>     8/22/01: THIS RELATIONSHIP HAS BEEN EXTENDED TO COLDER TEMPERATURES TO PARAMETERIZE SMALLER 
!>              LARGE-ICE PARTICLES DOWN TO MEAN SIZES OF MDIMIN, WHICH IS 50 MICRONS REACHED AT 
!>              -55.9C.
!>
!> (6) ATTEMPTS TO DIFFERENTIATE GROWTH OF LARGE AND SMALL ICE, MAINLY FOR IMPROVED TRANSITION FROM
!>     THIN CIRRUS TO THICK, PRECIPITATING ICE ANVILS.
!>     8/22/01: THIS FEATURE HAS BEEN DIMINISHED BY EFFECTIVELY ADJUSTING TO ICE SATURATION DURING 
!>              DEPOSITIONAL GROWTH AT TEMPERATURES COLDER THAN -10C.  ICE SUBLIMATION IS CALCULATED
!>              MORE EXPLICITLY. THE LOGIC IS THAT SOURCES OF ARE EITHER POORLY UNDERSTOOD 
!>              (E.G., NUCLEATION FOR NWP) OR ARE NOT REPRESENTED IN THE ETA MODEL 
!>              (E.G., DETRAINMENT OF ICE FROM CONVECTION). OTHERWISE THE MODEL IS TOO WET COMPARED 
!>              TO THE RADIOSONDE OBSERVATIONS BASED ON 1 FEB - 18 MARCH 2001 RETROSPECTIVE RUNS.
!>
!> (7) TOP-DOWN INTEGRATION ALSO ATTEMPTS TO TREAT MIXED-PHASE PROCESSES, ALLOWING A MIXTURE OF ICE 
!>     AND WATER.  BASED ON NUMEROUS OBSERVATIONAL STUDIES, ICE GROWTH IS BASED ON NUCLEATION AT 
!>     CLOUD TOP SUBSEQUENT GROWTH BY VAPOR DEPOSITION AND RIMING AS THE ICE PARTICLES FALL THROUGH 
!>     THE CLOUD. EFFECTIVE NUCLEATION RATES ARE A FUNCTION OF ICE SUPERSATURATION FOLLOWING MEYERS
!>     ET AL. (JAM, 1992).
!>     8/22/01: THE SIMULATED RELATIVE HUMIDITIES WERE FAR TOO MOIST COMPARED TO THE RAWINSONDE 
!>              OBSERVATIONS. THIS FEATURE HAS BEEN SUBSTANTIALLY DIMINISHED, LIMITED TO A MUCH 
!>              NARROWER TEMPERATURE RANGE OF 0 TO -10C.
!>
!> (8) DEPOSITIONAL GROWTH OF NEWLY NUCLEATED ICE IS CALCULATED FOR LARGE TIME STEPS USING FIG. 8 
!>     OF MILLER AND YOUNG (JAS, 1979), AT 1 DEG INTERVALS USING THEIR ICE CRYSTAL MASSES CALCULATED
!>     AFTER 600 S OF GROWTH IN WATER SATURATED CONDITIONS. THE GROWTH RATES ARE NORMALIZED BY TIME
!>     STEP ASSUMING 3D GROWTH WITH TIME**1.5 FOLLOWING EQ. (6.3) IN YOUNG (1993).
!>     8/22/01: THIS FEATURE HAS BEEN EFFECTIVELY LIMITED TO 0 TO -10C.
!>
!> (9) ICE PRECIPITATION RATES CAN INCREASE DUE TO INCREASE IN RESPONSE TO CLOUD WATER RIMING DUE TO
!>     (A) INCREASED DENSITY AND MASS OF THE RIMED ICE, AND (B) INCREASED FALL SPEEDS OF RIMED ICE.
!>     8/22/01: THIS FEATURE HAS BEEN EFFECTIVELY LIMITED TO 0 TO -10C.
!>
!> NOTE:  CODE IS CURRENTLY SET UP W/O THREADING
!> 
!> PROGRAM HISTORY LOG:
!> 01-08-??  FERRIER    - ORIGINATOR
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> DTPH       - PHYSICS TIME STEP (S)
!> I_INDEX    - I INDEX
!> J_INDEX    - J INDEX
!> LSFC       - ETA LEVEL OF LEVEL ABOVE SURFACE, GROUND
!> P_COL      - VERTICAL COLUMN OF MODEL PRESSURE (PA)
!> THICK_COL  - VERTICAL COLUMN OF MODEL MASS THICKNESS (DENSITY*HEIGHT INCREMENT)
!>
!> OUTPUT ARGUMENT LIST:
!> ARAIN      - ACCUMULATED RAINFALL AT THE SURFACE (KG)
!> ASNOW      - ACCUMULATED SNOWFALL AT THE SURFACE (KG)
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> QI_COL     - VERTICAL COLUMN OF MODEL ICE MIXING RATIO (KG/KG)
!> QR_COL     - VERTICAL COLUMN OF MODEL RAIN RATIO (KG/KG)
!> QV_COL     - VERTICAL COLUMN OF MODEL WATER VAPOR SPECIFIC HUMIDITY (KG/KG)
!> QW_COL     - VERTICAL COLUMN OF MODEL CLOUD WATER MIXING RATIO (KG/KG)
!> RIMEF_COL  - VERTICAL COLUMN OF RIME FACTOR FOR ICE IN MODEL (RATIO, DEFINED BELOW)
!> T_COL      - VERTICAL COLUMN OF MODEL TEMPERATURE (DEG K)
!> WC_COL     - VERTICAL COLUMN OF MODEL MIXING RATIO OF TOTAL CONDENSATE (KG/KG)
!>
!> USE MODULES: CMICRO_CONS
!>              CMICRO_STATS
!>              CMY600
!>              F77KINDS
!>              GLB_TABLE
!>              IACCR_TABLES 
!>              IMASS_TABLES
!>              IRATE_TABLES
!>              IRIME_TABLES
!>              IVENT_TABLES
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              RACCR_TABLES
!>              RRATE_TABLES
!>              RVELR_TABLES
!>              RVENT_TABLES
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : GSMDRIVE
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE CMICRO_CONS
    USE CMICRO_STATS
    USE CMY600
    USE F77KINDS
    USE GLB_TABLE
    USE IACCR_TABLES
    USE IMASS_TABLES
    USE IRATE_TABLES
    USE IRIME_TABLES
    USE IVENT_TABLES
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE RACCR_TABLES
    USE RRATE_TABLES
    USE RVELR_TABLES
    USE RVENT_TABLES
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!-------------------------------------- 
! ARRAYS AND CONSTANTS IN ARGUMENT LIST
!--------------------------------------
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)         ::&
    & ARAIN   , ASNOW      
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                   , INTENT(IN)          ::& 
    & P_COL   , THICK_COL 
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                   , INTENT(INOUT)       ::& 
    & QI_COL  , QR_COL  , QV_COL  , QW_COL  , RIMEF_COL, T_COL  , WC_COL
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & I_INDEX , J_INDEX , LSFC
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & DTPH     
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & TK       , TC      , PP        , QV       , WV       , ESW     , FPVS0   , QSW       ,      &
    & WS       , ESI     , FPVS      , QSI      , QSWGRD   , QSIGRD  , WSGRD   , RHO       ,      &
    & DTRHO    , BLDTRH  , THICK     , ARAINNEW , ASNOWNEW , QI      , WC      , RRHO      ,      &
    & QINEW    , QR      , QRNEW     , QW       , QWNEW    , PCOND   , PIDEP   , PIACW     ,      &
    & PIACWI   , PIACWR  , PIACR     , PICND    , PIEVP    , PIMLT   , PRAUT   , PRACW     ,      &
    & PREVP    , XLV     , XLF       , XLV1     , XLF1     , TK2     , XLV2    , DENOMW    ,      &
    & TFACTOR  , DYNVIS  , THERM_COND, DIFFUS   , GAMMAS   , GAMMAR  , QLICE   , QTICE     ,      &
    & WVQW     , FLARGE  , FSMALL    , XSIMASS  , TOT_ICE  , PILOSS  , RIMEF1  , TOT_ICENEW,      &
    & VRIMEF   , VEL_INC , VSNOW     , XLIMASS  , FLIMASS  , DUM     , DUM1    , DUM2      ,      &
    & XRF      , XLI     , DLI       , QW0      , FWS      , DEPOSIT , DWVI    , DENOMI    ,      &
    & PIDEP_MAX, SFACTOR , ABI       , VENTIL   , VENTIS   , DIDEP   , CONDENSE, DENOMWI   ,      &
    & DENOMF   , TCC     , AIEVP     , DIEVP    , QSW0     , DWV0    , TOT_RAIN, VRAIN1    ,      &
    & QTRAIN   , PRLOSS  , RQR       , RR       , DWVR     , RFACTOR ,                            &
    & ABW      , VENTR   , CREVP     , FWR      , FIR      , DELW    , DELT    , TNEW      ,      &
    & DELV     , WVNEW   , DELI      , RIMEF    , VRAIN2   , WCNEW   , QT      , QTNEW     ,      &
    & BUDGET   , QSWNEW  , QSINEW    , WSNEW    , DELR     , TOT_RAINNEW
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & L        , LBEF    , IPASS     , IXS      , IXRF     , INDEX_MY, IDR     , ITDX      ,      &
    & INDEXR1  , INDEXR  , INDEXS   
!--------------------------------------------------------------------------------------------------
! KEY PARAMETERS, LOCAL VARIABLES, AND IMPORTANT COMMENTS 
!
! KEY PARAMETERS:
! COMMENTS ON 2 AUGUST 2000
! EPSQ=1.E-20 IS THE LOWER LIMIT FOR SPECIFIC HUMIDITY AND CLOUD CONDENSATE. 
! THE VALUE OF EPSQ WILL NEED TO BE CHANGED IN THE OTHER SUBROUTINES IN ORDER TO MAKE IT CONSISTENT
! THROUGHOUT THE ETA CODE.
! - NLIMAX - MAXIMUM NUMBER CONCENTRATION OF LARGE ICE CRYSTALS (20,000 /M**3, 20 PER LITER)
! - NLIMIN - MINIMUM NUMBER CONCENTRATION OF LARGE ICE CRYSTALS (100 /M**3, 0.1 PER LITER)
!--------------------------------------------------------------------------------------------------
    REAL   (KIND=R4KIND), PARAMETER :: CP         = 1004.6
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ       =    1.E-20      
    REAL   (KIND=R4KIND), PARAMETER :: EPSW       =    1.E-12    
    REAL   (KIND=R4KIND), PARAMETER :: RD         =  287.04    
    REAL   (KIND=R4KIND), PARAMETER :: RHOL       = 1000.   
    REAL   (KIND=R4KIND), PARAMETER :: RV         =  461.5 
    REAL   (KIND=R4KIND), PARAMETER :: T0C        =  273.15        
    REAL   (KIND=R4KIND), PARAMETER :: XLS        =    2.834E6    
    REAL   (KIND=R4KIND), PARAMETER :: EPS        = RD         / RV     
    REAL   (KIND=R4KIND), PARAMETER :: NLIMAX     =   20.E3   
    REAL   (KIND=R4KIND), PARAMETER :: NLIMIN     =  100.    
    REAL   (KIND=R4KIND), PARAMETER :: T_ICE      =  -10.        
    REAL   (KIND=R4KIND), PARAMETER :: T_ICE_INIT =   -5. 
    REAL   (KIND=R4KIND), PARAMETER :: TOLER      =    5.E-7   
    REAL   (KIND=R4KIND), PARAMETER :: CLIMIT     =   10.      * EPSQ      
    REAL   (KIND=R4KIND), PARAMETER :: CLIMIT1    = -CLIMIT 
    REAL   (KIND=R4KIND), PARAMETER :: EPS1       = RV         / RD - 1.     
    REAL   (KIND=R4KIND), PARAMETER :: RCP        =    1.      / CP      
    REAL   (KIND=R4KIND), PARAMETER :: RCPRV      = RCP        / RV  
    REAL   (KIND=R4KIND), PARAMETER :: RRHOL      =    1.      / RHOL  
    REAL   (KIND=R4KIND), PARAMETER :: XLS1       = XLS        * RCP   
    REAL   (KIND=R4KIND), PARAMETER :: XLS2       = XLS        * XLS * RCPRV 
    REAL   (KIND=R4KIND), PARAMETER :: XLS3       = XLS        * XLS / RV
    REAL   (KIND=R4KIND), PARAMETER :: C1         =    1.      / 3.      
    REAL   (KIND=R4KIND), PARAMETER :: C2         =    1.      / 6.       
    REAL   (KIND=R4KIND), PARAMETER :: C3         =    3.31    / 6.     
    REAL   (KIND=R4KIND), PARAMETER :: DMR1       =     .1E-3        
    REAL   (KIND=R4KIND), PARAMETER :: DMR2       =     .2E-3     
    REAL   (KIND=R4KIND), PARAMETER :: DMR3       =     .32E-3   
    REAL   (KIND=R4KIND), PARAMETER :: N0R0       =    8.E6      
    REAL   (KIND=R4KIND), PARAMETER :: N0RMIN     =    1.E4    
    REAL   (KIND=R4KIND), PARAMETER :: N0S0       =    4.E6         
    REAL   (KIND=R4KIND), PARAMETER :: RHO0       =    1.194    
    REAL   (KIND=R4KIND), PARAMETER :: XMR1       =    1.E6    * DMR1 
    REAL   (KIND=R4KIND), PARAMETER :: XMR2       =    1.E6    * DMR2 
    REAL   (KIND=R4KIND), PARAMETER :: XMR3       =    1.E6    * DMR3 
    REAL   (KIND=R4KIND), PARAMETER :: XRATIO     =     .025
!
    INTEGER(KIND=I4KIND), PARAMETER :: MDR1       = XMR1
    INTEGER(KIND=I4KIND), PARAMETER :: MDR2       = XMR2
    INTEGER(KIND=I4KIND), PARAMETER :: MDR3       = XMR3
!-------------------------------------------------------------------------------------------------- 
! IF BLEND=1:
! PRECIPITATION (LARGE) ICE AMOUNTS ARE ESTIMATED AT EACH LEVEL AS A BLEND OF ICE FALLING FROM THE 
! GRID POINT ABOVE AND THE PRECIP ICE PRESENT AT THE START OF THE TIME STEP (SEE TOT_ICE BELOW).
! IF BLEND=0:
! PRECIPITATION (LARGE) ICE AMOUNTS ARE ESTIMATED TO BE THE PRECIP ICE PRESENT AT THE START OF THE
! TIME STEP.
!
! EXTENDED TO INCLUDE SEDIMENTATION OF RAIN ON 2/5/01
!-------------------------------------------------------------------------------------------------- 
    REAL   (KIND=R4KIND), PARAMETER :: BLEND = 1.
!----------------------------------------------------  
! THIS VARIABLE IS FOR DEBUGGING PURPOSES (IF .TRUE.)
!---------------------------------------------------- 
    LOGICAL(KIND=L4KIND), PARAMETER :: PRINT_DIAG = .TRUE. 
!----------------  
! LOCAL VARIABLES
!----------------
    REAL   (KIND=R4KIND)                                                                        ::&
    & EMAIRI  , N0R     , NLICE   , NSMICE
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & CLEAR   , ICE_LOGICAL       , DBG_LOGICAL           , DIAG_PRINT    , RAIN_LOGICAL
!---------------- 
! BEGIN EXECUTION  
!----------------
    ARAIN      = 0. ! ACCUMULATED RAINFALL INTO GRID BOX FROM ABOVE (KG/M**2)
    ASNOW      = 0. ! ACCUMULATED SNOWFALL INTO GRID BOX FROM ABOVE (KG/M**2)
!
    DIAG_PRINT = .FALSE. 
!---------------------------------------- 
! LOOP FROM TOP (L=1) TO SURFACE (L=LSFC) 
!----------------------------------------
    DO 10 L=1,LSFC
!-------------------------------------------------------------------------------------------------
! SKIP THIS LEVEL AND GO TO THE NEXT LOWER LEVEL IF NO CONDENSATE AND VERY LOW SPECIFIC HUMIDITIES
!-------------------------------------------------------------------------------------------------
        IF (QV_COL(L) <= EPSQ .AND. WC_COL(L) <= EPSQ) GOTO 10
!--------------------------------------------- 
! PROCEED WITH CLOUD MICROPHYSICS CALCULATIONS 
!--------------------------------------------- 
        TK = T_COL(L)     ! TEMPERATURE (DEG K)
        TC = TK - T0C     ! TEMPERATURE (DEG C)
        PP = P_COL(L)     ! PRESSURE (PA)
        QV = QV_COL(L)    ! SPECIFIC HUMIDITY OF WATER VAPOR (KG/KG)
        WV = QV / (1.-QV) ! WATER VAPOR MIXING RATIO (KG/KG)
        WC = WC_COL(L)    ! GRID-SCALE MIXING RATIO OF TOTAL CONDENSATE (WATER OR ICE; KG/KG)   
!-------------------------------------------------------------------- 
! MOISTURE VARIABLES BELOW ARE MIXING RATIOS & NOT SPECIFC HUMIDITIES
!--------------------------------------------------------------------   
        CLEAR = .TRUE. 
!-------------------------------------------------------------------------------    
! THIS CHECK IS TO DETERMINE GRID-SCALE SATURATION WHEN NO CONDENSATE IS PRESENT
!-------------------------------------------------------------------------------   
        ESW = 1000. * FPVS0(TK)          ! SATURATION VAPOR PRESSURE W/R/T WATER
        QSW = EPS   * ESW / (PP-ESW)     ! SATURATION MIXING RATIO W/R/T WATER
        WS  = QSW                        ! GENERAL SATURATION MIXING RATIO (WATER/ICE)
!
        IF (TC < 0.) THEN
            ESI = 1000. * FPVS(TK)       ! SATURATION VAPOR PRESSURE W/R/T ICE
            QSI = EPS   * ESI / (PP-ESI) ! SATURATION MIXING RATIO W/R/T WATER
            WS  = QSI                    ! GENERAL SATURATION MIXING RATIO (WATER/ICE)
        END IF
!----------------------------------------------     
! EFFECTIVE GRID-SCALE SATURATION MIXING RATIOS
!-----------------------------------------------------------------------
        QSWGRD = RHGRD * QSW
        QSIGRD = RHGRD * QSI
        WSGRD  = RHGRD * WS
!------------------------------------------------     
! CHECK IF AIR IS SUBSATURATED AND W/O CONDENSATE
!------------------------------------------------   
        IF (WV > WSGRD .OR. WC > EPSW) CLEAR = .FALSE. 
!---------------------------------------------------     
! CHECK IF ANY RAIN IS FALLING INTO LAYER FROM ABOVE
!---------------------------------------------------    
        IF (ARAIN > CLIMIT) THEN
            CLEAR = .FALSE. 
        ELSE
            ARAIN = 0.
        END IF
!-----------------------------------------------------------------------------------------    
! CHECK IF ANY ICE IS FALLING INTO LAYER FROM ABOVE
!    
! NOTE THAT "SNOW" IN VARIABLE NAMES IS SYNONOMOUS WITH LARGE, PRECIPITATION ICE PARTICLES
!----------------------------------------------------------------------------------------- 
        IF (ASNOW > CLIMIT) THEN
            CLEAR = .FALSE. 
        ELSE
            ASNOW = 0.
        END IF
!----------------------------------------------------------------- 
! LOOP TO THE END IF IN CLEAR, SUBSATURATED AIR FREE OF CONDENSATE 
!----------------------------------------------------------------- 
        IF (CLEAR) GOTO 10
!-------------------------------------------------------------------------------------------------- 
! INITIALIZE RHO, THICK AND MICROPHYSICAL PROCESSES 
!
! VIRTUAL TEMPERATURE, TV=T*(1./EPS-1)*Q, Q IS SPECIFIC HUMIDITY; 
! (SEE PP. 63-65 IN FLEAGLE & BUSINGER, 1963)
!--------------------------------------------------------------------------------------------------
        RHO    = PP / (RD * TK * (1. + EPS1 * QV)) ! AIR DENSITY (KG/M**3)
        RRHO   = 1. / RHO                          ! RECIPROCAL OF AIR DENSITY
        DTRHO  = DTPH  * RHO                       ! TIME STEP * AIR DENSITY
        BLDTRH = BLEND * DTRHO                 ! BLEND PARAMETER * TIME STEP * AIR DENSITY
        THICK  = THICK_COL(L)                      ! LAYER THICKNESS = RHO*DZ = -DP/G = 
                                                   ! (PSFC-PTOP)*D_ETA/(G*ETA_SFC)
!    
        ARAINNEW = 0.                      ! UPDATED ACCUMULATED RAINFALL
        ASNOWNEW = 0.                      ! UPDATED ACCUMULATED SNOWFALL
        QI       = QI_COL(L)               ! ICE MIXING RATIO
        QINEW    = 0.                      ! UPDATED ICE MIXING RATIO
        QR       = QR_COL(L)               ! RAIN MIXING RATIO
        QRNEW    = 0.                      ! UPDATED RAIN RATIO
        QW       = QW_COL(L)               ! CLOUD WATER MIXING RATIO
        QWNEW    = 0.                      ! UPDATED CLOUD WATER RATIO
!    
        PCOND  = 0.            ! CONDENSATION (>0) OR EVAPORATION (<0) OF CLOUD WATER (KG/KG)
        PIDEP  = 0.            ! DEPOSITION (>0) OR SUBLIMATION (<0) OF ICE CRYSTALS (KG/KG)
        PIACW  = 0.            ! CLOUD WATER COLLECTION (RIMING) BY PRECIPITATION ICE (KG/KG; >0)
        PIACWI = 0.            ! GROWTH OF PRECIP ICE BY RIMING (KG/KG; >0)
        PIACWR = 0.            ! SHEDDING OF ACCRETED CLOUD WATER TO FORM RAIN (KG/KG; >0)
        PIACR  = 0.            ! FREEZING OF RAIN ONTO LARGE ICE AT SUPERCOOLED TEMPS (KG/KG; >0)
        PICND  = 0.            ! CONDENSATION (>0) ONTO WET, MELTING ICE (KG/KG)
        PIEVP  = 0.            ! EVAPORATION (<0) FROM WET, MELTING ICE (KG/KG)
        PIMLT  = 0.            ! MELTING ICE (KG/KG; >0)
        PRAUT  = 0.            ! CLOUD WATER AUTOCONVERSION TO RAIN (KG/KG; >0)
        PRACW  = 0.            ! CLOUD WATER COLLECTION (ACCRETION) BY RAIN (KG/KG; >0)
        PREVP  = 0.            ! RAIN EVAPORATION (KG/KG; <0)
!---------------------------------------------    
! DOUBLE CHECK INPUT HYDROMETEOR MIXING RATIOS
!--------------------------------------------- 
!    
!------------------------------------------ 
! MAIN MICROPHYSICS CALCULATIONS NOW FOLLOW 
!------------------------------------------ 
!
!--------------------------------------------------------------- 
! CALCULATE A FEW VARIABLES, WHICH ARE USED MORE THAN ONCE BELOW
!--------------------------------------------------------------- 
!
!--------------------------------------------------------------------------------
! LATENT HEAT OF VAPORIZATION AS A FUNCTION OF TEMPERATURE FROM BOLTON (1980,JAS)
!--------------------------------------------------------------------------------    
        XLV    = 3.148E6 - 2370 * TK        ! LATENT HEAT OF VAPORIZATION (LV)
        XLF    = XLS - XLV                  ! LATENT HEAT OF FUSION (LF)
        XLV1   = XLV * RCP                  ! LV/CP
        XLF1   = XLF * RCP                  ! LF/CP
        TK2    = 1. / (TK*TK)               ! 1./TK**2
        XLV2   = XLV * XLV * QSW * TK2 / RV ! LV**2*QSW/(RV*TK**2)
        DENOMW = 1. + XLV2 * RCP            ! DENOMINATOR TERM, CLAUSIUS-CLAPEYRON CORRECTION
!-------------------------------------------------     
! BASIC THERMODYNAMIC QUANTITIES
! DYNVIS - DYNAMIC VISCOSITY  [ KG/(M*S) ]
! THERM_COND - THERMAL CONDUCTIVITY  [ J/(M*S*K) ]
! DIFFUS - DIFFUSIVITY OF WATER VAPOR  [ M**2/S ]
!-------------------------------------------------     
        TFACTOR    = TK ** 1.5 / (TK + 120.)
        DYNVIS     = 1.496E-6 * TFACTOR
        THERM_COND = 2.116E-3 * TFACTOR
        DIFFUS     = 8.794E-5 * TK ** 1.81 / PP
!--------------------------------------------------------------------------------------------------    
! AIR RESISTANCE TERM FOR THE FALL SPEED OF ICE FOLLOWING THE BASIC RESEARCH BY HEYMSFIELD, 
! KAJIKAWA, OTHERS
!--------------------------------------------------------------------------------------------------   
        GAMMAS = (1.E5 / PP) ** C1
!--------------------------------------------------------------------------------------------------   
! AIR RESISTANCE FOR RAIN FALL SPEED (BEARD, 1985, JAS, P.470)
!-------------------------------------------------------------------------------------------------- 
        GAMMAR = (RHO0 / RHO) ** .4
!-------------------------------------  
! IMPORTANT MICROPHYSICS DECISION TREE 
!-------------------------------------   
!
!---------------------------------------------------  
! DETERMINE IF CONDITIONS SUPPORTING ICE ARE PRESENT
!--------------------------------------------------- 
        IF (TC < 0. .OR. QI > CLIMIT .OR. ASNOW > CLIMIT) THEN
            ICE_LOGICAL = .TRUE. 
        ELSE
            ICE_LOGICAL = .FALSE. 
            QLICE = 0.
            QTICE = 0.
        END IF
!-----------------------------  
! DETERMINE IF RAIN IS PRESENT
!----------------------------- 
        RAIN_LOGICAL = .FALSE. 
!
        IF (ARAIN > CLIMIT .OR. QR > CLIMIT) RAIN_LOGICAL = .TRUE. 
!    
        IF (ICE_LOGICAL) THEN
!--------------------------------------------------------------------------------------------------     
! IMPORTANT:  ESTIMATE TIME-AVERAGED PROPERTIES.
!
! FLARGE  - RATIO OF NUMBER OF LARGE ICE TO TOTAL (LARGE & SMALL) ICE
! FSMALL  - RATIO OF NUMBER OF SMALL ICE CRYSTALS TO LARGE ICE PARTICLES
!           SMALL ICE PARTICLES ARE ASSUMED TO HAVE A MEAN DIAMETER OF 50 MICRONS.
! XSIMASS - USED FOR CALCULATING SMALL ICE MIXING RATIO
!
! TOT_ICE - TOTAL MASS (SMALL & LARGE) ICE BEFORE MICROPHYSICS, WHICH IS THE SUM OF THE TOTAL MASS 
!           OF LARGE ICE IN THE CURRENT LAYER AND THE INPUT FLUX OF ICE FROM ABOVE
! PILOSS  - GREATEST LOSS (<0) OF TOTAL (SMALL & LARGE) ICE BY SUBLIMATION, REMOVING ALL OF THE ICE 
!           FALLING FROM ABOVE AND THE ICE WITHIN THE LAYER
! RIMEF1  - RIME FACTOR, WHICH IS THE MASS RATIO OF TOTAL (UNRIMED & RIMED) ICE MASS TO THE UNRIMED
!           ICE MASS (>=1)
! VRIMEF  - THE VELOCITY INCREASE DUE TO RIME FACTOR OR MELTING (RATIO, >=1)
! VSNOW   - FALL SPEED OF RIMED SNOW W/ AIR RESISTANCE CORRECTION
! EMAIRI  - EQUIVALENT MASS OF AIR ASSOCIATED LAYER AND WITH FALL OF SNOW INTO LAYER
! XLIMASS - USED FOR CALCULATING LARGE ICE MIXING RATIO
! FLIMASS - MASS FRACTION OF LARGE ICE
! QTICE   - TIME-AVERAGED MIXING RATIO OF TOTAL ICE
! QLICE   - TIME-AVERAGED MIXING RATIO OF LARGE ICE
! NLICE   - TIME-AVERAGED NUMBER CONCENTRATION OF LARGE ICE
! NSMICE  - NUMBER CONCENTRATION OF SMALL ICE CRYSTALS AT CURRENT LEVEL
! 
! ASSUMED NUMBER FRACTION OF LARGE ICE PARTICLES TO TOTAL (LARGE & SMALL) ICE PARTICLES, WHICH IS 
! BASED ON A GENERAL IMPRESSION OF THE LITERATURE.
!--------------------------------------------------------------------------------------------------   
            WVQW = WV + QW                ! WATER VAPOR AND CLOUD WATER
!        
            IF (TC >= 0. .OR. WVQW < QSIGRD) THEN
!-----------------------------------------------------------------------              
! ELIMINATE SMALL ICE PARTICLE CONTRIBUTIONS FOR MELTING AND SUBLIMATION
!-----------------------------------------------------------------------   
                FLARGE = 1.
            ELSE
!--------------------------------------------------------------------------------------------------            
! ENHANCED NUMBER OF SMALL ICE PARTICLES DURING DEPOSITIONAL GROWTH 
! (EFFECTIVE ONLY WHEN 0C > T >= T_ICE [-10C])
!--------------------------------------------------------------------------------------------------   
                FLARGE = .2
!-------------------------------------------------------------             
! LARGER NUMBER OF SMALL ICE PARTICLES DUE TO RIME SPLINTERING
!------------------------------------------------------------- 
                IF (TC >= -8. .AND. TC <= -3.) FLARGE = .5 * FLARGE
            
            END IF  ! END IF (TC >= 0. .OR. WVQW < QSIGRD)
!
            FSMALL  = (1. - FLARGE) / FLARGE
            XSIMASS = RRHO * MASSI(MDIMIN) * FSMALL
!
            IF (QI <= CLIMIT .AND. ASNOW <= CLIMIT) THEN
                INDEXS  = MDIMIN
                TOT_ICE = 0.
                PILOSS  = 0.
                RIMEF1  = 1.
                VRIMEF  = 1.
                VEL_INC = GAMMAS
                VSNOW   = 0.
                EMAIRI  = THICK
                XLIMASS = RRHO * RIMEF1 * MASSI(INDEXS)
                FLIMASS = XLIMASS / (XLIMASS + XSIMASS)
                QLICE   = 0.
                QTICE   = 0.
                NLICE   = 0.
                NSMICE  = 0.
            ELSE
!--------------------------------------------------------------------------------------------------        
! FOR T<0C MEAN PARTICLE SIZE FOLLOWS HOUZE ET AL. (JAS, 1979, P. 160), CONVERTED FROM FIG. 5 PLOT 
! OF LAMDAS. SIMILAR SET OF RELATIONSHIPS ALSO SHOWN IN FIG. 8 OF RYAN (BAMS, 1996, P. 66).
!--------------------------------------------------------------------------------------------------  
                DUM     = XMIMAX * EXP(.0536 * TC)
                INDEXS  = MIN(MDIMAX, MAX(MDIMIN, INT(DUM)))
                TOT_ICE = THICK * QI + BLEND * ASNOW
                PILOSS  = -TOT_ICE / THICK
                LBEF    = MAX(1,L-1)
                DUM1    = RIMEF_COL(LBEF)
                DUM2    = RIMEF_COL(L)
                RIMEF1  = (DUM2 * THICK * QI + DUM1 * BLEND * ASNOW) / TOT_ICE
                RIMEF1  = MIN(RIMEF1, RFMAX)
!
                DO IPASS=0,1
                    IF (RIMEF1 <= 1.) THEN
                        RIMEF1 = 1.
                        VRIMEF = 1.
                    ELSE
                        IXS  = MAX(2, MIN(INDEXS/100, 9))
                        XRF  = 10.492 * ALOG (RIMEF1)
                        IXRF = MAX(0, MIN(INT(XRF), NRIME))
!
                        IF (IXRF >= NRIME) THEN
                            VRIMEF = VEL_RF(IXS,NRIME)
                        ELSE
                            VRIMEF =  VEL_RF(IXS,IXRF)   + (XRF-FLOAT(IXRF))                      &
    &                              * (VEL_RF(IXS,IXRF+1) - VEL_RF(IXS,IXRF))
                        END IF
!
                    END IF             ! END IF (RIMEF1 <= 1.)
!
                    VEL_INC = GAMMAS  * VRIMEF
                    VSNOW   = VEL_INC * VSNOWI(INDEXS)
                    EMAIRI  = THICK   + BLDTRH * VSNOW
                    XLIMASS = RRHO    * RIMEF1 * MASSI(INDEXS)
                    FLIMASS = XLIMASS / (XLIMASS + XSIMASS)
                    QTICE   = TOT_ICE / EMAIRI
                    QLICE   = FLIMASS * QTICE
                    NLICE   = QLICE   / XLIMASS
                    NSMICE  = FSMALL  * NLICE
                
                    IF ( (NLICE >= NLIMIN .AND. NLICE <= NLIMAX) .OR. IPASS == 1) THEN
                        EXIT
                    ELSE
!-------------------------------------------------------------------------------------------------   
! REDUCE EXCESSIVE ACCUMULATION OF ICE AT UPPER LEVELS ASSOCIATED WITH STRONG GRID-RESOLVED ASCENT
!                    
! FORCE NLICE TO BE BETWEEN NLIMIN AND NLIMAX
!------------------------------------------------------------------------------------------------- 
                    
                        DUM = MAX(NLIMIN, MIN(NLIMAX, NLICE))
                        XLI = RHO * (QTICE / DUM - XSIMASS) / RIMEF1
!
                        IF (XLI <= MASSI(MDIMIN)) THEN
                            INDEXS = MDIMIN
                        ELSE IF (XLI <= MASSI(450) ) THEN
                            DLI = 9.5885E5 * XLI ** .42066 ! DLI IN MICRONS
                            INDEXS = MIN(MDIMAX, MAX(MDIMIN, INT(DLI)))
                        ELSE IF (XLI <= MASSI(MDIMAX)) THEN
                            DLI = 3.9751E6 * XLI ** .49870 ! DLI IN MICRONS
                            INDEXS = MIN(MDIMAX, MAX(MDIMIN, INT(DLI)))
                        ELSE
                            INDEXS = MDIMAX
!--------------------------------------------------------------------------------------------------                        
! 8/22/01: INCREASE DENSITY OF LARGE ICE IF MAXIMUM LIMITS ARE REACHED FOR NUMBER CONCENTRATION 
! (NLIMAX) AND MEAN SIZE (MDIMAX).
! DONE TO INCREASE FALL OUT OF ICE.
!--------------------------------------------------------------------------------------------------  
                            IF (DUM >= NLIMAX)                                                    &
    &                         RIMEF1 = RHO * (QTICE / NLIMAX - XSIMASS) / MASSI(INDEXS)
!
                        END IF ! END IF (XLI <= MASSI(MDIMIN) )
                    END IF ! END IF ( (NLICE >= NLIMIN .AND. NLICE <= NLIMAX) 
                END DO     ! END DO IPASS=0,1
            END IF         ! END IF (QI <= CLIMIT .AND. ASNOW <= CLIMIT)
        END IF             ! END IF (ICE_LOGICAL)
!------------------------------- 
! CALCULATE INDIVIDUAL PROCESSES 
!-------------------------------  
!
!---------------------------------------------------------- 
! CLOUD WATER AUTOCONVERSION TO RAIN AND COLLECTION BY RAIN
!----------------------------------------------------------  	
        IF (QW > CLIMIT .AND. TC >= T_ICE) THEN
!--------------------------------------------------------------------------------------------------        
! QW0 COULD BE MODIFIED BASED ON LAND/SEA PROPERTIES, PRESENCE OF CONVECTION, ETC.  
! THIS IS WHY QAUT0 AND CRAUT ARE PASSED INTO THE SUBROUTINE AS EXTERNALLY DETERMINED PARAMETERS. 
! CAN BE CHANGED IN THE FUTURE IF DESIRED.
!--------------------------------------------------------------------------------------------------    
            QW0 = QAUT0 * RRHO
            PRAUT = MAX(0., QW - QW0) * CRAUT
!
            IF (QLICE > CLIMIT) THEN
!----------------------------------------------------------             
! COLLECTION OF CLOUD WATER BY LARGE ICE PARTICLES ("SNOW")
! PIACWI=PIACW FOR RIMING, PIACWI=0 FOR SHEDDING
!---------------------------------------------------------- 
                FWS   = MIN(1., CIACW * VEL_INC * NLICE * ACCRI(INDEXS) / PP ** C1)
                PIACW = FWS * QW
!
                IF (TC < 0.) PIACWI = PIACW   ! LARGE ICE RIMING
            END IF                            ! END IF (QLICE > CLIMIT)
        END IF                                ! END IF (QW > CLIMIT .AND. TC >= T_ICE)
!------------------------------------------------------------------------
! LOOP AROUND SOME OF THE ICE-PHASE PROCESSES IF NO ICE SHOULD BE PRESENT
!------------------------------------------------------------------------
    
        IF (ICE_logical .EQV. .FALSE. ) GOTO 20
!----------------------------------------------------     
! NOW THE PRETZEL LOGIC OF CALCULATING ICE DEPOSITION
!----------------------------------------------------   
        IF (TC < T_ICE .AND. (WV > QSIgrd .OR. QW > CLIMIT)) THEN
!--------------------------------------------------------------------------------------------------         
! ADJUST TO ICE SATURATION AT T<T_ICE (-10C) IF SUPERSATURATED.
! SOURCES OF ICE DUE TO NUCLEATION AND CONVECTIVE DETRAINMENT ARE EITHER POORLY UNDERSTOOD, POORLY
! RESOLVED AT TYPICAL NWP RESOLUTIONS, OR ARE NOT REPRESENTED 
! (E.G., NO DETRAINED CONDENSATE IN BMJ CU SCHEME).
!--------------------------------------------------------------------------------------------------    
            PCOND = -QW
            DUM1  =  TK + XLV1 * PCOND                 ! UPDATED (DUMMY) TEMPERATURE (DEG K)
            DUM2  =  WV + QW                           ! UPDATED (DUMMY) WATER VAPOR MIXING RATIO
            DUM   =  1000. * FPVS(DUM1)                ! UPDATED (DUMMY) SATURATION VAPOR PRESSURE W/R/T ICE
            DUM   =  RHgrd * EPS * DUM / (PP - DUM)    ! UPDATED (DUMMY) SATURATION MIXING RATIO W/R/T ICE
!
            IF (DUM2 > DUM) PIDEP = DEPOSIT(PP, RHGRD, DUM1, DUM2)
            DWVI = 0.                                  ! USED ONLY FOR DEBUGGING
!        
        ELSE IF (TC < 0.) THEN
!--------------------------------------------------------------------         
! THESE QUANTITIES ARE HANDY FOR ICE DEPOSITION/SUBLIMATION
! PIDEP_MAX - MAX DEPOSITION OR MINIMUM SUBLIMATION TO ICE SATURATION
!--------------------------------------------------------------------  
            DENOMI    = 1. + XLS2 * QSI * TK2
            DWVI      = MIN(WVQW, QSW) - QSI
            PIDEP_MAX = MAX(PILOSS, DWVI / DENOMI)
!
            IF (QTICE > 0.) THEN
!--------------------------------------------------------------------------------------------------            
! CALCULATE ICE DEPOSITION/SUBLIMATION
! SFACTOR - [VEL_INC**.5]*[SCHMIDT**(1./3.)]*[(RHO/DYNVIS)**.5], WHERE SCHMIDT 
! (SCHMIDT NUMBER) = DYNVIS/(RHO*DIFFUS)
! UNITS: SFACTOR - S**.5/M; ABI - M**2/S; NLICE - M**-3; VENTIL, VENTIS - M**-2; VENTI1 - M;
! VENTI2 - M**2/S**.5; DIDEP - UNITLESS
!-------------------------------------------------------------------------------------------------- 
                SFACTOR = VEL_INC ** .5 * (RHO / (DIFFUS * DIFFUS * DYNVIS)) ** C2
                ABI     = 1. / (RHO * XLS3 * QSI * TK2 / THERM_COND + 1. / DIFFUS)
!--------------------------------------------------------------------------------------------------              
! VENTIL - NUMBER CONCENTRATION * VENTILATION FACTORS FOR LARGE ICE
! VENTIS - NUMBER CONCENTRATION * VENTILATION FACTORS FOR SMALL ICE
!            
! VARIATION IN THE NUMBER CONCENTRATION OF ICE WITH TIME IS NOT ACCOUNTED FOR IN THESE CALCULATIONS
! (COULD BE IN THE FUTURE).
!-------------------------------------------------------------------------------------------------- 
                VENTIL = (VENTI1(INDEXS) + SFACTOR * VENTI2(INDEXS)) * NLICE
                VENTIS = (VENTI1(MDIMIN) + SFACTOR * VENTI2(MDIMIN)) * NSMICE
!
                DIDEP  = ABI * (VENTIL + VENTIS) * DTPH
!------------------------------------------------             
! ACCOUNT FOR CHANGE IN WATER VAPOR SUPPLY W/TIME
!------------------------------------------------             
                IF (DIDEP >= XRATIO) DIDEP = (1. - EXP(-DIDEP * DENOMI)) / DENOMI
!
                IF (DWVI > 0.) THEN
                    PIDEP = MIN(DWVI * DIDEP, PIDEP_MAX)
                ELSE IF (DWVI < 0.) THEN
                    PIDEP = MAX(DWVI * DIDEP, PIDEP_MAX)
                END IF
!            
            ELSE IF (WVQW > QSI .AND. TC <= T_ICE_INIT) THEN
!--------------------------------------------------------------------------------------------------             
! ICE NUCLEATION IN NEAR WATER-SATURATED CONDITIONS. ICE CRYSTAL GROWTH DURING TIME STEP CALCULATED
! USING MILLER & YOUNG (1979, JAS)
! THESE DEPOSITION RATES COULD DRIVE CONDITIONS BELOW WATER SATURATION, WHICH IS THE BASIS OF THESE
! CALCULATIONS.
! INTENDED TO APPROXIMATE MORE COMPLEX AND COMPUTATIONALLY INTENSIVE CALCULATIONS.
!-------------------------------------------------------------------------------------------------- 
                INDEX_MY = MAX(MY_T1, MIN(INT(.5 - TC), MY_T2))
!--------------------------------------------------------------------------------------------------             
! DUM1 IS THE SUPERSATURATION W/R/T ICE AT WATER-SATURATED CONDITIONS
! DUM2 IS THE NUMBER OF ICE CRYSTALS NUCLEATED AT WATER-SATURATED CONDITIONS 
! BASED ON MEYERS ET AL. (JAM, 1992).
!            
! PREVENT UNREALISTICALLY LARGE ICE INITIATION (LIMITED BY PIDEP_MAX) IF DUM2 VALUES ARE INCREASED 
! IN FUTURE EXPERIMENTS
!-------------------------------------------------------------------------------------------------- 
                DUM1  = QSW / QSI - 1.
                DUM2  = 1.E3 * EXP(12.96 * DUM1 - .039)
                PIDEP = MIN(PIDEP_MAX, DUM2 * MY_GROWTH(INDEX_MY) * RRHO)
            END IF       ! END IF (QTICE > 0.)
!        
        END IF           ! END IF (TC < T_ICE .AND. (WV > QSIGRD .OR. QW > CLIMIT))
!
     20 CONTINUE      ! JUMP HERE IF CONDITIONS FOR ICE ARE NOT PRESENT
!-------------------------     
! CLOUD WATER CONDENSATION
!-------------------------    
        IF (TC >= T_ICE .AND. (QW > CLIMIT .OR. WV > QSWGRD)) THEN
            IF (PIACWI == 0. .AND. PIDEP == 0.) THEN
                PCOND = CONDENSE(PP, QW, RHGRD, TK, WV)
            ELSE
!-------------------------------------------------------            
! MODIFY CLOUD CONDENSATION IN RESPONSE TO ICE PROCESSES
!-------------------------------------------------------   
                DUM     = XLV * QSWGRD * RCPRV * TK2
                DENOMWI = 1. + XLS * DUM
                DENOMF  = XLF * DUM
                DUM     = MAX(0., PIDEP)
                PCOND   = (WV - QSWGRD - DENOMWI * DUM - DENOMF * PIACWI) / DENOMW
                DUM1    = -QW
                DUM2    = PCOND - PIACW
!
                IF (DUM2 < DUM1) THEN
!------------------------                   
! LIMIT CLOUD WATER SINKS
!------------------------ 
                    DUM    = DUM1 / DUM2
                    PCOND  = DUM  * PCOND
                    PIACW  = DUM  * PIACW
                    PIACWI = DUM  * PIACWI
                END IF ! END IF (DUM2 < DUM1)
            END IF     ! END IF (PIACWI == 0. .AND. PIDEP == 0.)
        END IF         ! END IF (TC >= T_ICE .AND. (QW > CLIMIT .OR. WV > QSWGRD))
!--------------------------------------------------------------------------------------------------     
! LIMIT FREEZING OF ACCRETED RIME TO PREVENT TEMPERATURE OSCILLATIONS, A CRUDE SCHUMANN-LUDLAM 
! LIMIT (P. 209 OF YOUNG, 1993).
!--------------------------------------------------------------------------------------------------   
        TCC = TC + XLV1 * PCOND + XLS1 * PIDEP + XLF1 * PIACWI
!
        IF (TCC > 0.) THEN
            PIACWI = 0.
            TCC = TC + XLV1 * PCOND + XLS1 * PIDEP
        END IF
 !   
        IF (TC > 0. .AND. TCC > 0. .AND. ICE_LOGICAL) THEN
!--------------------------------------------------------------------------------------------------         
! CALCULATE MELTING AND EVAPORATION/CONDENSATION
! UNITS: SFACTOR - S**.5/M; ABI    - M**2/S    ; NLICE - M**-3; VENTIL - M**-2;
!        VENTI1  - M      ; VENTI2 - M**2/S**.5; CIEVP - /S
!-------------------------------------------------------------------------------------------------- 
            SFACTOR = VEL_INC ** .5 * (RHO / (DIFFUS * DIFFUS * DYNVIS)) ** C2
            VENTIL  = NLICE * (VENTI1(INDEXS) + SFACTOR * VENTI2(INDEXS))
            AIEVP   = VENTIL * DIFFUS * DTPH
!
            IF (AIEVP < XRATIO) THEN
                DIEVP = AIEVP
            ELSE
                DIEVP = 1. - EXP(-AIEVP)
            END IF
!
            QSW0 = EPS * ESW0 / (PP - ESW0)
            DWV0 = MIN(WV,QSW) - QSW0
            DUM  = QW + PCOND
!
            IF (WV < QSW .AND. DUM <= CLIMIT) THEN
!--------------------------------------------------------------------------------------------------           
! EVAPORATION FROM MELTING SNOW (SINK OF SNOW) OR SHEDDING OF WATER CONDENSED ONTO MELTING SNOW 
! (SOURCE OF RAIN)
!-------------------------------------------------------------------------------------------------- 
                DUM   = DWV0 * DIEVP
                PIEVP = MAX(MIN(0., DUM), PILOSS)
                PICND = MAX(0., DUM)
            END IF            ! END IF (WV < QSW .AND. DUM <= CLIMIT)
!
            PIMLT = THERM_COND * TCC * VENTIL * RRHO * DTPH / XLF
!------------------------------------------------------------          
! LIMIT MELTING TO PREVENT TEMPERATURE OSCILLATIONS ACROSS 0C
!------------------------------------------------------------  
            DUM1  = MAX(0., (TCC + XLV1 * PIEVP) / XLF1)
            PIMLT = MIN(PIMLT, DUM1)
!--------------------------------------------------- 
! LIMIT LOSS OF SNOW BY MELTING (>0) AND EVAPORATION
!---------------------------------------------------     
            DUM = PIEVP - PIMLT
!
            IF (DUM < PILOSS) THEN
                DUM1  = PILOSS / DUM
                PIMLT = PIMLT * DUM1
                PIEVP = PIEVP * DUM1
            END IF           ! END IF (DUM > QTICE)
        END IF               ! END IF (TC > 0. .AND. TCC > 0. .AND. ICE_LOGICAL)
!--------------------------------------------------------------------------------------------------    
! IMPORTANT:  ESTIMATE TIME-AVERAGED PROPERTIES.
!    
!TOT_RAIN - TOTAL MASS OF RAIN BEFORE MICROPHYSICS, WHICH IS THE SUM OF THE TOTAL MASS OF RAIN IN 
!           THE CURRENT LAYER AND THE INPUT FLUX OF RAIN FROM ABOVE
! VRAIN1   - FALL SPEED OF RAIN INTO GRID FROM ABOVE (WITH AIR RESISTANCE CORRECTION)
! QTRAIN   - TIME-AVERAGED MIXING RATIO OF RAIN (KG/KG)
! PRLOSS   - GREATEST LOSS (<0) OF RAIN, REMOVING ALL RAIN FALLING FROM ABOVE AND THE RAIN WITHIN 
!            THE LAYER
! RQR      - RAIN CONTENT (KG/M**3)
! INDEXR   - MEAN SIZE OF RAIN DROPS TO THE NEAREST 1 MICRON IN SIZE
! N0R      - INTERCEPT OF RAIN SIZE DISTRIBUTION (TYPICALLY 10**6 M**-4)
!--------------------------------------------------------------------------------------------------  
        TOT_RAIN = 0.
        VRAIN1   = 0.
        QTRAIN   = 0.
        PRLOSS   = 0.
        RQR      = 0.
        N0R      = 0.
        INDEXR1  = INDEXR    ! FOR DEBUGGING ONLY
        INDEXR   = MDRMIN
!
        IF (RAIN_LOGICAL) THEN
            IF (ARAIN <= 0.) THEN
                INDEXR = MDRMIN
                VRAIN1 = 0.
            ELSE
!--------------------------------------------------------------------------------------------------             
! INDEXR (RELATED TO MEAN DIAMETER) & N0R COULD BE MODIFIED BY LAND/SEA PROPERTIES, PRESENCE OF 
! CONVECTION, ETC.
!            
! RAIN RATE NORMALIZED TO A DENSITY OF 1.194 KG/M**3
!--------------------------------------------------------------------------------------------------          
                RR = ARAIN / (DTPH * GAMMAR)
!            
                IF (RR <= RR_DRMIN) THEN
!------------------------------------------------------------------------------------------------                 
! ASSUME FIXED MEAN DIAMETER OF RAIN (0.2 MM) FOR LOW RAIN RATES, INSTEAD VARY N0R WITH RAIN RATE
!------------------------------------------------------------------------------------------------ 
                
                    INDEXR = MDRMIN
                ELSE IF (RR <= RR_DR1) THEN
!--------------------------------------------------------------------------------------------------                
! BEST FIT TO MASS-WEIGHTED FALL SPEEDS (V) FROM RAIN LOOKUP TABLES FOR MEAN DIAMETERS (DR) BETWEEN
! 0.05 AND 0.10 MM:
! V(DR)=5.6023E4*DR**1.136, V IN M/S AND DR IN M
! RR = PI*1000.*N0R0*5.6023E4*DR**(4+1.136) = 1.408E15*DR**5.136, RR IN KG/(M**2*S)
! DR (M) = 1.123E-3*RR**.1947 -> DR (MICRONS) = 1.123E3*RR**.1947
!-------------------------------------------------------------------------------------------------- 
                    INDEXR = INT(1.123E3 * RR ** .1947 + .5)
                    INDEXR = MAX(MDRMIN, MIN(INDEXR, MDR1))
!
                ELSE IF (RR <= RR_DR2) THEN
!--------------------------------------------------------------------------------------------------                 
! BEST FIT TO MASS-WEIGHTED FALL SPEEDS (V) FROM RAIN LOOKUP TABLES FOR MEAN DIAMETERS (DR) BETWEEN
! 0.10 AND 0.20 MM:
! V(DR)=1.0867E4*DR**.958, V IN M/S AND DR IN M
! RR = PI*1000.*N0R0*1.0867E4*DR**(4+.958) = 2.731E14*DR**4.958, RR IN KG/(M**2*S)
! DR (M) = 1.225E-3*RR**.2017 -> DR (MICRONS) = 1.225E3*RR**.2017
!--------------------------------------------------------------------------------------------------           
                    INDEXR = INT(1.225E3 * RR ** .2017 + .5)
                    INDEXR = MAX(MDR1, MIN(INDEXR, MDR2))
!
                ELSE IF (RR <= RR_DR3) THEN
!--------------------------------------------------------------------------------------------------                 
! BEST FIT TO MASS-WEIGHTED FALL SPEEDS (V) FROM RAIN LOOKUP TABLES FOR MEAN DIAMETERS (DR) BETWEEN
! 0.20 AND 0.32 MM:
! V(DR)=2831.*DR**.80, V IN M/S AND DR IN M
! RR = PI*1000.*N0R0*2831.*DR**(4+.80) = 7.115E13*DR**4.80, RR IN KG/(M**2*S)
! DR (M) = 1.3006E-3*RR**.2083 -> DR (MICRONS) = 1.3006E3*RR**.2083
!--------------------------------------------------------------------------------------------------  
                    INDEXR = INT(1.3006E3 * RR ** .2083 + .5)
                    INDEXR = MAX(MDR2, MIN(INDEXR, MDR3))
!
                ELSE IF (RR <= RR_DRMAX) THEN
!--------------------------------------------------------------------------------------------------                
! BEST FIT TO MASS-WEIGHTED FALL SPEEDS (V) FROM RAIN LOOKUP TABLES FOR MEAN DIAMETERS (DR) BETWEEN
! 0.32 AND 0.45 MM:
! V(DR)=944.8*DR**.6636, V IN M/S AND DR IN M
! RR = PI*1000.*N0R0*944.8*DR**(4+.6636) = 2.3745E13*DR**4.6636, RR IN KG/(M**2*S)
! DR (M) = 1.355E-3*RR**.2144 -> DR (MICRONS) = 1.355E3*RR**.2144
!-------------------------------------------------------------------------------------------------- 
                    INDEXR = INT(1.355E3 * RR ** .2144 + .5)
                    INDEXR = MAX(MDR3, MIN(INDEXR, MDRMAX))
                ELSE
!--------------------------------------------------------------------------------------------------  
! ASSUME FIXED MEAN DIAMETER OF RAIN (0.45 MM) FOR HIGH RAIN RATES, INSTEAD VARY N0R WITH RAIN RATE
!-------------------------------------------------------------------------------------------------- 
                    INDEXR = MDRMAX
                END IF              ! END IF (RR <= RR_DRMIN) ETC.
                VRAIN1 = GAMMAR * VRAIN(INDEXR)
            END IF                  ! END IF (ARAIN <= 0.)
!
            INDEXR1  = INDEXR       ! FOR DEBUGGING ONLY
            TOT_RAIN = THICK * QR + BLEND * ARAIN
            QTRAIN   = TOT_RAIN / (THICK + BLDTRH * VRAIN1)
            PRLOSS   =-TOT_RAIN / THICK
            RQR      = RHO * QTRAIN
!-------------------------------------------         
! RQR - TIME-AVERAGED RAIN CONTENT (KG/M**3)
!------------------------------------------- 
            IF (RQR <= RQR_DRMIN) THEN
                N0R    = MAX(N0RMIN, CN0R_DMRMIN * RQR)
                INDEXR = MDRMIN
            ELSE IF (RQR >= RQR_DRMAX) THEN
                N0R    = CN0R_DMRMAX * RQR
                INDEXR = MDRMAX
            ELSE
                N0R    = N0R0
                INDEXR = MAX(XMRMIN, MIN(CN0R0 * RQR ** .25, XMRMAX))
            END IF
!        
            IF (TC < T_ICE) THEN
                PIACR = -PRLOSS
            ELSE
                DWVR = WV - PCOND - QSW
                DUM  = QW + PCOND
!
                IF (DWVR < 0. .AND. DUM <= CLIMIT) THEN
!--------------------------------------------------------------------------------------------------                
! RAIN EVAPORATION
!                
! RFACTOR - [GAMMAR**.5]*[SCHMIDT**(1./3.)]*[(RHO/DYNVIS)**.5], 
! WHERE SCHMIDT (SCHMIDT NUMBER) = DYNVIS/(RHO*DIFFUS)
!                
! UNITS: RFACTOR - S**.5/M; ABW    - M**2/S    ; VENTR - M**-2   ;N0R - M**-4; 
!         VENTR1 - M**2   ; VENTR2 - M**3/S**.5; CREVP - UNITLESS
!-------------------------------------------------------------------------------------------------- 
                    RFACTOR = GAMMAR ** .5 * (RHO / (DIFFUS * DIFFUS * DYNVIS)) ** C2
                    ABW     = 1. / (RHO * XLV2 / THERM_COND + 1. / DIFFUS)
!-------------------------------------------------------------------------------------------------- 
! NOTE THAT VENTR1, VENTR2 LOOKUP TABLES DO NOT INCLUDE THE 1/DAVG MULTIPLIER AS IN THE ICE TABLES
!--------------------------------------------------------------------------------------------------  
                    VENTR = N0R * (VENTR1(INDEXR) + RFACTOR * VENTR2(INDEXR))
                    CREVP = ABW * VENTR * DTPH
!
                    IF (CREVP < XRATIO) THEN
                        DUM = DWVR * CREVP
                    ELSE
                        DUM=DWVR * (1. - EXP(-CREVP * DENOMW)) / DENOMW
                    END IF
!
                    PREVP = MAX(DUM, PRLOSS)
!
                ELSE IF (QW > CLIMIT) THEN
                    FWR   = CRACW * GAMMAR * N0R * ACCRR(INDEXR)
                    PRACW = MIN(1.,FWR) * QW
                END IF       ! END IF (DWVR < 0. .AND. DUM <= CLIMIT)
!            
                IF (TC < 0. .AND. TCC < 0.) THEN
!--------------------------------------------------------------                 
! BIGGS (1953) HETEOROGENEOUS FREEZING (E.G., LIN ET AL., 1983)
!--------------------------------------------------------------
                    PIACR = CBFR * N0R * RRHO * (EXP(ABFR * TC) - 1.) * (FLOAT(INDEXR)) ** 7
!
                    IF (QLICE > CLIMIT) THEN
!-------------------------------------------                     
! FREEZING OF RAIN BY COLLISIONS W/LARGE ICE
!------------------------------------------- 
                        DUM  = GAMMAR * VRAIN(INDEXR)
                        DUM1 = DUM - VSNOW
!--------------------------------------------------------------------------------------------------                    
! DUM2 - DIFFERENCE IN SPECTRAL FALL SPEEDS OF RAIN AND LARGE ICE, PARAMETERIZED FOLLOWING EQ. (48)
! ON P. 112 OF MURAKAMI (J. METEOR. SOC. JAPAN, 1990)
!--------------------------------------------------------------------------------------------------  
                        DUM2 = (DUM1 * DUM1 + .04 * DUM * VSNOW) ** .5
                        DUM1 = 5.E-12 * INDEXR * INDEXR + 2.E-12 * INDEXR * INDEXS                &
    &                        + .5E-12 * INDEXS * INDEXS
!
                        FIR = MIN(1., CIACR*NLICE*DUM1*DUM2)
!--------------------------------------------------------------                     
! FUTURE ?  SHOULD COLLECTION BY SMALL ICE SHOULD BE INCLUDED ?
!-------------------------------------------------------------- 
                        PIACR = MIN(PIACR+FIR*QTRAIN, QTRAIN)
                    END IF   ! END IF (QLICE > CLIMIT)
!
                    DUM = PREVP - PIACR
!
                    If (DUM < PRLOSS) THEN
                        DUM1  = PRLOSS / DUM
                        PREVP = DUM1   * PREVP
                        PIACR = DUM1   * PIACR
                    END IF ! END IF (DUM < PRLOSS)
                END IF     ! END IF (TC < 0. .AND. TCC < 0.)
            END IF         ! END IF (TC < T_ICE)
        END IF             ! END IF (RAIN_LOGICAL)
!---------------------- 
! MAIN BUDGET EQUATIONS 
!---------------------- 
!
!-------------------------------------------------------------- 
! UPDATE FIELDS, DETERMINE CHARACTERISTICS FOR NEXT LOWER LAYER 
!-------------------------------------------------------------- 
!
!------------------------------------- 
! CAREFULLY LIMIT SINKS OF CLOUD WATER
!-------------------------------------
        DUM1 = PIACW + PRAUT + PRACW - MIN(0.,PCOND)
!
        IF (DUM1 > QW) THEN
            DUM    = QW  / DUM1
            PIACW  = DUM * PIACW
            PIACWI = DUM * PIACWI
            PRAUT  = DUM * PRAUT
            PRACW  = DUM * PRACW
!
            IF (PCOND < 0.) PCOND = DUM * PCOND
        END IF
!
        PIACWR = PIACW - PIACWI          ! TC >= 0C
!-----------------------------------------     
! QWNEW - UPDATED CLOUD WATER MIXING RATIO
!----------------------------------------- 
        DELW  = PCOND - PIACW - PRAUT - PRACW
        QWnew = QW + DELW
!
        IF (QWNEW <= CLIMIT) QWNEW = 0.
!
        IF (QW > 0. .AND. QWNEW /= 0.) THEN
            DUM = QWNEW / QW
!
            IF (DUM < TOLER) QWNEW = 0.
        END IF
!-------------------------------------------------    
! UPDATE TEMPERATURE AND WATER VAPOR MIXING RATIOS
!-------------------------------------------------  
        DELT = XLV1 * (PCOND + PIEVP +  PICND  + PREVP)                                           &
    &        + XLS1 *  PIDEP + XLF1  * (PIACWI + PIACR - PIMLT)
!
        TNEW = TK + DELT
!    
        DELV  = -PCOND - PIDEP - PIEVP - PICND - PREVP
        WVNEW = WV + DELV
!--------------------------------------------------------------------------------------------------   
! UPDATE ICE MIXING RATIOS
!
! TOT_ICENEW - TOTAL MASS (SMALL AND LARGE) ICE AFTER MICROPHYSICS, WHICH IS THE SUM OF THE TOTAL 
!              MASS OF LARGE ICE IN THE CURRENT LAYER AND THE FLUX OF ICE OUT OF THE GRID BOX BELOW
! RIMEF      - RIME FACTOR, WHICH IS THE MASS RATIO OF TOTAL (UNRIMED AND RIMED) ICE MASS TO THE 
!              UNRIMED ICE MASS (>=1)
! QINEW      - UPDATED MIXING RATIO OF TOTAL (LARGE & SMALL) ICE IN LAYER
!              TOT_ICENEW=QINEW*THICK+BLDTRH*QLICENEW*VSNOW 
!              BUT QLICENEW=QINEW*FLIMASS, 
!              SO TOT_ICENEW=QINEW*(THICK+BLDTRH*FLIMASS*VSNOW)
! ASNOWNEW   - UPDATED ACCUMULATION OF SNOW AT BOTTOM OF GRID CELL
!--------------------------------------------------------------------------------------------------
        DELI  = 0.
        RIMEF = 1.
!
        IF (ICE_LOGICAL) THEN
            DELI = PIDEP + PIEVP + PIACWI + PIACR - PIMLT
            TOT_ICENEW = TOT_ICE + THICK * DELI
!
            IF (TOT_ICE > 0. .AND. TOT_ICENEW /= 0.) THEN
                DUM = TOT_ICENEW / TOT_ICE
                IF (DUM < TOLER) TOT_ICENEW = 0.
            END IF
!
            IF (TOT_ICENEW <= CLIMIT) THEN
                TOT_ICENEW = 0.
                RIMEF      = 1.
                QINEW      = 0.
                ASNOWNEW   = 0.
            ELSE
!----------------------------------            
! UPDATE RIME FACTOR IF APPROPRIATE
!----------------------------------  
                DUM = PIACWI + PIACR
                IF (DUM <= CLIMIT .AND. PIDEP <= CLIMIT) THEN
                    RIMEF = RIMEF1
                ELSE
!--------------------------------------------------------------- 
! RIME FACTOR, RIMEF = (TOTAL ICE MASS)/(TOTAL UNRIMED ICE MASS)
! DUM1 - TOTAL ICE MASS, RIMED & UNRIMED
! DUM2 - ESTIMATED MASS OF *UNRIMED* ICE
!---------------------------------------------------------------
                    DUM1 = TOT_ICE + THICK  * (PIDEP + DUM)
                    DUM2 = TOT_ICE / RIMEF1 + THICK  * PIDEP
!
                    IF (DUM2 <= 0.) THEN
                        RIMEF = RFmax
                    ELSE
                        RIMEF = MIN(RFMAX, MAX(1., DUM1 / DUM2))
                    END IF
                END IF   ! END IF (DUM <= CLIMIT .AND. PIDEP <= CLIMIT)
!
                QINEW = TOT_ICENEW / (THICK + BLDTRH * FLIMASS * VSNOW)
!
                IF (QINEW <= CLIMIT) QINEW = 0.
!
                IF (QI > 0. .AND. QINEW /= 0.) THEN
                    DUM = QINEW / QI
                    IF (DUM < TOLER) QINEW = 0.
                END IF
!
                ASNOWNEW = BLDTRH * FLIMASS * VSNOW * QINEW
!
                IF (ASNOW > 0. .AND. ASNOWNEW /= 0.) THEN
                    DUM = ASNOWNEW / ASNOW
                    IF (DUM < TOLER) ASNOWNEW = 0.
                END IF
            END IF     ! END IF (TOT_ICENEW <= CLIMIT)
        END IF         ! END IF (ICE_LOGICAL)
!--------------------------------------------------------------------------------------------------    
! UPDATE RAIN MIXING RATIOS
!
! TOT_RAINNEW - TOTAL MASS OF RAIN AFTER MICROPHYSICS CURRENT LAYER AND THE INPUT FLUX OF ICE FROM 
!               ABOVE
! VRAIN2      - TIME-AVERAGED FALL SPEED OF RAIN IN GRID AND BELOW (WITH AIR RESISTANCE CORRECTION)
! QRNEW       - UPDATED RAIN MIXING RATIO IN LAYER
!               TOT_RAINNEW=QRNEW*(THICK+BLDTRH*VRAIN2)
! ARAINNEW    - UPDATED ACCUMULATION OF RAIN AT BOTTOM OF GRID CELL
!--------------------------------------------------------------------------------------------------
        DELR        = PRAUT    + PRACW + PIACWR - PIACR + PIMLT + PREVP + PICND
        TOT_RAINNEW = TOT_RAIN + THICK * DELR
!
        IF (TOT_RAIN > 0. .AND. TOT_RAINNEW /= 0.) THEN
            DUM = TOT_RAINNEW / TOT_RAIN
            IF (DUM < TOLER) TOT_RAINNEW = 0.
        END IF
!
        IF (TOT_RAINnew <= CLIMIT) THEN
            TOT_RAINNEW = 0.
            VRAIN2      = 0.
            QRNEW       = 0.
            ARAINNEW    = 0.
        ELSE
!--------------------------------------------------------         
! 1ST GUESS TIME-AVERAGED RAIN RATE AT BOTTOM OF GRID BOX
!--------------------------------------------------------   
            RR = TOT_RAINNEW / (DTPH * GAMMAR)
!--------------------------------------------------------------------------------------------------        
! USE SAME ALGORITHM AS ABOVE FOR CALCULATING MEAN DROP DIAMETER (IDR, IN MICRONS), WHICH IS USED 
! TO ESTIMATE THE TIME-AVERAGED FALL SPEED OF RAIN DROPS AT THE BOTTOM OF THE GRID LAYER. 
! THIS ISNT PERFECT, BUT THE ALTERNATIVE IS SOLVING A TRANSCENDENTAL EQUATION THAT IS NUMERICALLY 
! INEFFICIENT AND NASTY TO PROGRAM (CODED IN EARLIER VERSIONS OF GSMCOLUMN PRIOR TO 8-22-01).
!--------------------------------------------------------------------------------------------------
            IF (RR <= RR_DRMIN) THEN
                IDR = MDRMIN
            ELSE IF (RR <= RR_DR1) THEN
                IDR = INT(1.123E3 * RR ** .1947 + .5)
                IDR = MAX(MDRMIN, MIN(IDR, MDR1))
            ELSE IF (RR <= RR_DR2) THEN
                IDR = INT(1.225E3 * RR ** .2017 + .5)
                IDR = MAX(MDR1, MIN(IDR, MDR2))
            ELSE IF (RR <= RR_DR3) THEN
                IDR = INT(1.3006E3 * RR ** .2083 + .5)
                IDR = MAX(MDR2, MIN(IDR, MDR3) )
            ELSE IF (RR <= RR_DRMAX) THEN
                IDR = INT(1.355E3 * RR ** .2144 + .5)
                IDR = MAX(MDR3, MIN(IDR, MDRMAX))
            ELSE
                IDR = MDRMAX
            END IF  ! END IF (RR <= RR_DRMIN)
!
            VRAIN2 = GAMMAR * VRAIN(IDR)
            QRNEW  = TOT_RAINNEW / (THICK + BLDTRH * VRAIN2)
!
            IF (QRNEW <= CLIMIT) QRNEW = 0.
            IF (QR > 0. .AND. QRNEW /= 0.) THEN
                DUM = QRNEW / QR
                IF (DUM < TOLER) QRNEW = 0.
            END IF
!
            ARAINNEW = BLDTRH * VRAIN2 * QRNEW
            IF (ARAIN > 0. .AND. ARAINNEW /= 0.) THEN
                DUM = ARAINNEW / ARAIN
                IF (DUM < TOLER) ARAINNEW = 0.
            END IF
        END IF
!    
        WCNEW = QWNEW + QRNEW + QINEW
!--------------------------------- 
! BEGIN DEBUGGING AND VERIFICATION  
!--------------------------------- 
!
!------------------------------------------------------------------------------------ 
! QT, QTNEW - TOTAL WATER (VAPOR AND CONDENSATE) BEFORE AND AFTER MICROPHYSICS, RESP.
!------------------------------------------------------------------------------------ 
        QT     = THICK * (WV + WC) + ARAIN + ASNOW
        QTNEW  = THICK * (WVNEW + WCNEW) + ARAINNEW + ASNOWNEW
        BUDGET = QT - QTNEW
!---------------------------------------------------------------------------     
! ADDITIONAL CHECK ON BUDGET PRESERVATION, ACCOUNTING FOR TRUNCATION EFFECTS
!--------------------------------------------------------------------------- 
        DBG_LOGICAL = .FALSE. 
!
        IF ((WVNEW < CLIMIT .OR. DBG_LOGICAL) .AND. PRINT_DIAG) THEN
        
            WRITE(6,"(/2(A,I4),2(A,I2))") '{} I=',I_INDEX,' J=',J_INDEX,' L=',L,' LSFC=',LSFC
!        
            ESW    = 1000. * FPVS0(TNEW)
            QSWNEW = EPS * ESW / (PP - ESW)
!
            IF (TC < 0. .OR. TNEW < 0.) THEN
                ESI    = 1000. * FPVS(TNEW)
                QSINEW = EPS * ESI / (PP - ESI)
            ELSE
                QSI    = QSW
                QSINEW = QSWNEW
            END IF
!
            WSNEW = QSINEW
!
            WRITE(6,"(4(A12,G11.4,1x))")                                                          &
    &       '{} TCOLD=',TC,'TCNEW=',TNEW-T0C,'P=',.01*PP,'RHO=',RHO,                              &
    &       '{} THICK=',THICK,'RHOLD=',WV/WS,'RHNEW=',WVNEW/WSNEW,                                &
    &       'RHGRD=',RHGRD,                                                                       &
    &       '{} RHWOLD=',WV/QSW,'RHWNEW=',WVNEW/QSWNEW,'RHIOLD=',WV/QSI,                          &
    &       'RHINEW=',WVNEW/QSINEW,                                                               &
    &       '{} QSWOLD=',QSW,'QSWNEW=',QSWNEW,'QSIOLD=',QSI,'QSINEW=',QSINEW,                     &
    &       '{} WSOLD=',WS,'WSNEW=',WSNEW,'WVOLD=',WV,'WVNEW=',WVNEW,                             &
    &       '{} WCOLD=',WC,'WCNEW=',WCNEW,'QWOLD=',QW,'QWNEW=',QWNEW,                             &
    &       '{} QIOLD=',QI,'QINEW=',QINEW,'QROLD=',QR,'QRNEW=',QRNEW,                             &
    &       '{} ARAINOLD=',ARAIN,'ARAINNEW=',ARAINNEW,'ASNOWOLD=',ASNOW,                          &
    &       'ASNOWNEW=',ASNOWNEW,                                                                 &
    &       '{} TOT_RAIN=',TOT_RAIN,'TOT_RAINNEW=',TOT_RAINNEW,                                   &
    &       'TOT_ICE=',TOT_ICE,'TOT_ICENEW=',TOT_ICENEW,                                          &
    &       '{} BUDGET=',BUDGET,'QTOLD=',QT,'QTNEW=',QTNEW
        
            WRITE(6,"(4(A12,G11.4,1x))")                                                          &
    &       '{} DELT=',DELT,'DELV=',DELV,'DELW=',DELW,'DELI=',DELI,                               &
    &       '{} DELR=',DELR,'PCOND=',PCOND,'PIDEP=',PIDEP,'PIEVP=',PIEVP,                         &
    &       '{} PICND=',PICND,'PREVP=',PREVP,'PRAUT=',PRAUT,'PRACW=',PRACW,                       &
    &       '{} PIACW=',PIACW,'PIACWI=',PIACWI,'PIACWR=',PIACWR,'PIMLT=',PIMLT,                   &
    &       '{} PIACR=',PIACR
!        
            IF (ICE_LOGICAL) WRITE(6,"(4(A12,G11.4,1x))")                                         &
    &       '{} RIMEF1=',RIMEF1,'GAMMAS=',GAMMAS,'VRIMEF=',VRIMEF,                                &
    &       'VSNOW=',VSNOW,                                                                       &
    &       '{} INDEXS=',FLOAT(INDEXS),'FLARGE=',FLARGE,'FSMALL=',FSMALL,                         &
    &       'FLIMASS=',FLIMASS,                                                                   &
    &       '{} XSIMASS=',XSIMASS,'XLIMASS=',XLIMASS,'QLICE=',QLICE,                              &
    &       'QTICE=',QTICE,                                                                       &
    &       '{} NLICE=',NLICE,'NSMICE=',NSMICE,'PILOSS=',PILOSS,                                  &
    &       'EMAIRI=',EMAIRI,                                                                     &
    &       '{} RIMEF=',RIMEF
!        
            IF (TOT_RAIN > 0. .OR. TOT_RAINNEW > 0.)                                              &
    &       WRITE(6,"(4(A12,G11.4,1x))")                                                          &
    &       '{} INDEXR1=',FLOAT(INDEXR1),'INDEXR=',FLOAT(INDEXR),                                 &
    &       'GAMMAR=',GAMMAR,'N0R=',N0R,                                                          &
    &       '{} VRAIN1=',VRAIN1,'VRAIN2=',VRAIN2,'QTRAIN=',QTRAIN,'RQR=',RQR,                     &
    &       '{} PRLOSS=',PRLOSS,'VOLR1=',THICK+BLDTRH*VRAIN1,                                     &
    &       'VOLR2=',THICK+BLDTRH*VRAIN2
!        
            IF (PRAUT > 0.) WRITE(6,"(A12,G11.4,1x)") '{} QW0=',QW0
!        
            IF (PRACW > 0.) WRITE(6,"(A12,G11.4,1x)") '{} FWR=',FWR
!        
            IF (PIACR > 0.) WRITE(6,"(A12,G11.4,1x)") '{} FIR=',FIR
!        
            DUM = PIMLT + PICND - PREVP - PIEVP
!
            IF (DUM > 0. .OR. DWVI /= 0.)                                                         &
    &       WRITE(6,"(4(A12,G11.4,1x))")                                                          &
    &       '{} TFACTOR=',TFACTOR,'DYNVIS=',DYNVIS,                                               &
    &       'THERM_CON=',THERM_COND,'DIFFUS=',DIFFUS
!        
            IF (PREVP < 0.) WRITE(6,"(4(A12,G11.4,1x))")                                          &
    &       '{} RFACTOR=',RFACTOR,'ABW=',ABW,'VENTR=',VENTR,'CREVP=',CREVP,                       &
    &       '{} DWVR=',DWVR,'DENOMW=',DENOMW
!        
            IF (PIDEP /= 0. .AND. DWVI /= 0.)                                                     &
    &        WRITE(6,"(4(A12,G11.4,1x))")                                                         &
    &       '{} DWVI=',DWVI,'DENOMI=',DENOMI,'PIDEP_MAX=',PIDEP_MAX,                              &
    &       'SFACTOR=',SFACTOR,                                                                   &
    &       '{} ABI=',ABI,'VENTIL=',VENTIL,'VENTIL1=',VENTI1(INDEXS),                             &
    &       'VENTIL2=',SFACTOR*VENTI2(INDEXS),                                                    &
    &       '{} VENTIS=',VENTIS,'DIDEP=',DIDEP
!        
            IF (PIDEP > 0. .AND. PCOND /= 0.)                                                     &
    &       WRITE(6,"(4(A12,G11.4,1x))")                                                          &
    &       '{} DENOMW=',DENOMW,'DENOMWI=',DENOMWI,'DENOMF=',DENOMF,                              &
    &       'DUM2=',PCOND-PIACW
!        
            IF (FWS > 0.) WRITE(6,"(4(A12,G11.4,1x))")                                            &
    &       '{} FWS=',FWS
!        
            DUM = PIMLT + PICND - PIEVP
!
            IF (DUM > 0.) WRITE(6,"(4(A12,G11.4,1x))")                                            &
    &       '{} SFACTOR=',SFACTOR,'VENTIL=',VENTIL,'VENTIL1=',VENTI1(INDEXS),                     &
    &       'VENTIL2=',SFACTOR*VENTI2(INDEXS),                                                    &
    &       '{} AIEVP=',AIEVP,'DIEVP=',DIEVP,'QSW0=',QSW0,'DWV0=',DWV0
!        
        END IF
!------------------------------------------- 
! WATER BUDGET STATISTICS AND MAXIMUM VALUES 
!-------------------------------------------
        IF (PRINT_DIAG) THEN
            ITdX = MAX(ITLO, MIN( INT(TNEW - T0C), ITHI))
!
            IF (QINEW > CLIMIT)                               NSTATS(ITDX,1) = NSTATS(ITDX,1) + 1
            IF (QINEW > CLIMIT .AND.  QRNEW + QWNEW > CLIMIT) NSTATS(ITDX,2) = NSTATS(ITDX,2) + 1
            IF (QWNEW > CLIMIT)                               NSTATS(ITDX,3) = NSTATS(ITDX,3) + 1
            IF (QRNEW > CLIMIT)                               NSTATS(ITDX,4) = NSTATS(ITDX,4) + 1
!        
            QMAX(ITDX,1) = MAX(QMAX(ITDX,1), QINEW)
            QMAX(ITDX,2) = MAX(QMAX(ITDX,2), QWNEW)
            QMAX(ITDX,3) = MAX(QMAX(ITDX,3), QRNEW)
            QMAX(ITDX,4) = MAX(QMAX(ITDX,4), ASNOWNEW)
            QMAX(ITDX,5) = MAX(QMAX(ITDX,5), ARAINNEW)
!
            QTOT(ITDX,1)  = QTOT(ITDX,1) + QINEW * THICK
            QTOT(ITDX,2)  = QTOT(ITDX,2) + QWNEW * THICK
            QTOT(ITDX,3)  = QTOT(ITDX,3) + QRNEW * THICK
!        
            QTOT(ITDX,4)  = QTOT(ITDX,4)  + PCOND  * THICK
            QTOT(ITDX,5)  = QTOT(ITDX,5)  + PICND  * THICK
            QTOT(ITDX,6)  = QTOT(ITDX,6)  + PIEVP  * THICK
            QTOT(ITDX,7)  = QTOT(ITDX,7)  + PIDEP  * THICK
            QTOT(ITDX,8)  = QTOT(ITDX,8)  + PREVP  * THICK
            QTOT(ITDX,9)  = QTOT(ITDX,9)  + PRAUT  * THICK
            QTOT(ITDX,10) = QTOT(ITDX,10) + PRACW  * THICK
            QTOT(ITDX,11) = QTOT(ITDX,11) + PIMLT  * THICK
            QTOT(ITDX,12) = QTOT(ITDX,12) + PIACW  * THICK
            QTOT(ITDX,13) = QTOT(ITDX,13) + PIACWI * THICK
            QTOT(ITDX,14) = QTOT(ITDX,14) + PIACWR * THICK
            QTOT(ITDX,15) = QTOT(ITDX,15) + PIACR  * THICK
!        
            QTOT(ITDX,16) = QTOT(ITDX,16) + (WVNEW - WV) * THICK
            QTOT(ITDX,17) = QTOT(ITDX,17) + (QWNEW - QW) * THICK
            QTOT(ITDX,18) = QTOT(ITDX,18) + (QINEW - QI) * THICK
            QTOT(ITDX,19) = QTOT(ITDX,19) + (QRNEW - QR) * THICK
!
            QTOT(ITDX,20) = QTOT(ITDX,20) + (ARAINNEW - ARAIN)
            QTOT(ITDX,21) = QTOT(ITDX,21) + (ASNOWNEW - ASNOW)
!
            IF (QINEW > 0.) QTOT(ITdX,22) = QTOT(ITDX,22) + QINEW * THICK / RIMEF
!        
        END IF
!-------------- 
! UPDATE ARRAYS 
!-------------- 
            T_COL(L) = TNEW                             ! UPDATED TEMPERATURE   
           QV_COL(L) = MAX(EPSQ, WVNEW  / (1. + WVNEW)) ! UPDATED SPECIFIC HUMIDITY
           WC_COL(L) = MAX(EPSQ, WCNEW)                 ! UPDATED TOTAL CONDENSATE MIXING RATIO
           QI_COL(L) = MAX(EPSQ, QINEW)                 ! UPDATED ICE MIXING RATIO
           QR_COL(L) = MAX(EPSQ, QRNEW)                 ! UPDATED RAIN MIXING RATIO
           QW_COL(L) = MAX(EPSQ, QWNEW)                 ! UPDATED CLOUD WATER MIXING RATIO
        RIMEF_COL(L) = RIMEF                            ! UPDATED RIME FACTOR
!
        ASNOW        = ASNOWNEW                         ! UPDATED ACCUMULATED SNOW
        ARAIN        = ARAINNEW                         ! UPDATED ACCUMULATED RAIN
!
        CONTINUE   ! END "L" LOOP THROUGH MODEL LEVELS
!
 10 END DO

!------------------- 
! RETURN TO GSMDRIVE 
!------------------- 
    RETURN
!
    END SUBROUTINE GSMCOLUMN
!
!
!
!---------------------------------------------------- 
! PRODUCES ACCURATE CALCULATION OF CLOUD CONDENSATION 
!---------------------------------------------------- 
    REAL FUNCTION CONDENSE(PP, QW, RHGRD, TK, WV)
!
    USE F77KINDS
!--------------------------------------------------------------------------------------------------
! THE ASAI (1965) ALGORITHM TAKES INTO CONSIDERATION THE RELEASE OF LATENT HEAT IN INCREASING THE 
! TEMPERATURE AND IN INCREASING THE SATURATION MIXING RATIO (FOLLOWING THE CLAUSIUS-CLAPEYRON EQN.)
!--------------------------------------------------------------------------------------------------
    INTEGER, PARAMETER :: HIGH_PRES = SELECTED_REAL_KIND(15)
!
    REAL   (KIND=HIGH_PRES), PARAMETER :: CLIMIT   = 1.E-20
    REAL   (KIND=HIGH_PRES), PARAMETER :: RHLIMIT  =  .001
    REAL   (KIND=HIGH_PRES), PARAMETER :: RHLIMIT1 = -RHLIMIT
!
    REAL   (KIND=R4KIND), PARAMETER :: CP    = 1004.6
    REAL   (KIND=R4KIND), PARAMETER :: RD    =  287.04
    REAL   (KIND=R4KIND), PARAMETER :: RV    =  461.5
    REAL   (KIND=R4KIND), PARAMETER :: EPS   =  RD / RV
    REAL   (KIND=R4KIND), PARAMETER :: RCP   =  1. / CP
    REAL   (KIND=R4KIND), PARAMETER :: RCPRV = RCP / RV
!
    REAL   (KIND=HIGH_PRES)                                                                     ::&
    & COND    , SSAT    , WCDUM
!---------------------------------- 
! LV (T) IS FROM BOLTON (JAS, 1980)
!----------------------------------
    XLV      = 3.148E6 - 2370. * TK
    XLV1     = XLV * RCP
    XLV2     = XLV * XLV * RCPRV
    TDUM     = TK
    WVDUM    = WV
    WCDUM    = QW
    ESW      = 1000. * FPVS0(TDUM)                 ! SATURATION VAPOR PRESS W/R/T ICE
    WS       = RHGRD * EPS * ESW / (PP - ESW)      ! SATURATION MIXING RATIO
    DWV      = WVDUM - WS                          ! DEFICIT GRID-SCALE WATER VAPOR MIXING RATIO
    SSAT     = DWV / WS                            ! SUPERSATURATION RATIO
    CONDENSE = 0.
!
    DO WHILE ((SSAT < RHLIMIT1 .AND. WCDUM > CLIMIT) .OR. SSAT > RHLIMIT)
        COND     = DWV / (1. + XLV2 * WS / (TDUM*TDUM))       ! ASAI (1965, J. JAPAN)
        COND     = MAX(COND, -WCDUM)                          ! LIMIT CLOUD WATER EVAPORATION
        TDUM     = TDUM + XLV1 * COND                         ! UPDATED TEMPERATURE
        WVDUM    = WVDUM    - COND                            ! UPDATED WATER VAPOR MIXING RATIO
        WCDUM    = WCDUM    + COND                            ! UPDATED CLOUD WATER MIXING RATIO
        CONDENSE = CONDENSE + COND                            ! TOTAL CLOUD WATER CONDENSATION
        ESW      = 1000. * FPVS0(TDUM)                        ! UPDATED SATURATION VAPOR PRESS W/R/T WATER
        WS       = RHGRD * EPS * ESW / (PP - ESW)             ! UPDATED SATURATION MIXING RATIO W/R/T WATER
        DWV      = WVDUM - WS                                 ! DEFICIT GRID-SCALE WATER VAPOR MIXING RATIO
        SSAT     = DWV / WS                                   ! GRID-SCALE SUPERSATURATION RATIO
    END DO
!
    RETURN
!
    END FUNCTION CONDENSE
!
!
!
!------------------------------------
! CALCULATE ICE DEPOSITION AT T<T_ICE 
!------------------------------------
    REAL FUNCTION DEPOSIT(PP, RHGRD, TDUM, WVDUM)
!
    USE F77KINDS
!--------------------------------------------------------------------------------------------------
! ALSO USES THE ASAI (1965) ALGORITHM BUT USES A DIFFERENT TARGET VAPOR PRESSURE FOR THE ADJUSTMENT
!--------------------------------------------------------------------------------------------------
    INTEGER                , PARAMETER :: HIGH_PRES=SELECTED_REAL_KIND(15)
!
    REAL   (KIND=HIGH_PRES), PARAMETER :: RHLIMIT  =      .001
    REAL   (KIND=HIGH_PRES), PARAMETER :: RHLIMIT1 = -RHLIMIT
!
    REAL   (KIND=R4KIND)   , PARAMETER :: CP       =  1004.6 
    REAL   (KIND=R4KIND)   , PARAMETER :: RD       =   287.04
    REAL   (KIND=R4KIND)   , PARAMETER :: RV       =   461.5
    REAL   (KIND=R4KIND)   , PARAMETER :: XLS      =     2.834E6 
    REAL   (KIND=R4KIND)   , PARAMETER :: EPS      = RD / RV
    REAL   (KIND=R4KIND)   , PARAMETER :: RCP      =     1. / CP
    REAL   (KIND=R4KIND)   , PARAMETER :: RCPRV    = RCP / RV
    REAL   (KIND=R4KIND)   , PARAMETER :: XLS1     = XLS * RCP
    REAL   (KIND=R4KIND)   , PARAMETER :: XLS2     = XLS * XLS * RCPRV
!
    REAL   (KIND=HIGH_PRES)                                                                     ::&
    & DEP     , SSAT
!
    ESI     = 1000. * FPVS(TDUM)             ! SATURATION VAPOR PRESS W/R/T ICE
    WS      = RHGRD * EPS * ESI / (PP - ESI) ! SATURATION MIXING RATIO
    DWV     = WVDUM - WS                     ! DEFICIT GRID-SCALE WATER VAPOR MIXING RATIO
    SSAT    = DWV   / WS                     ! SUPERSATURATION RATIO
    DEPOSIT = 0.
!
    DO WHILE (SSAT > RHLIMIT .OR. SSAT < RHLIMIT1)
!--------------------------------------------------------------------------------------------------    
! NOTE THAT XLVS2=LS*LV/(CP*RV)=LV*WS/(RV*T*T)*(LS/CP*DEP1), WHERE WS IS THE SATURATION MIXING 
! RATIO FOLLOWING CLAUSIUS-CLAPEYRON (SEE ASAI,1965; YOUNG,1993,P.405)
!--------------------------------------------------------------------------------------------------   
        DEP     = DWV / (1. + XLS2 * WS / (TDUM * TDUM)) ! ASAI (1965, J. JAPAN)
!
        TDUM    = TDUM    + XLS1 * DEP               ! UPDATED TEMPERATURE
        WVDUM   = WVDUM   - DEP                      ! UPDATED ICE MIXING RATIO
        DEPOSIT = DEPOSIT + DEP                      ! TOTAL ICE DEPOSITION
        ESI     = 1000.   * FPVS(TDUM)               ! UPDATED SATURATION VAPOR PRESS W/R/T ICE
        WS      = RHGRD   * EPS  * ESI / (PP - ESI)  ! UPDATED SATURATION MIXING RATIO W/R/T ICE
        DWV     = WVDUM   - WS                       ! DEFICIT GRID-SCALE WATER VAPOR MIXING RATIO
        SSAT    = DWV     / WS                       ! GRID-SCALE SUPERSATURATION RATIO
    END DO
!
    RETURN
!
    END FUNCTION DEPOSIT
