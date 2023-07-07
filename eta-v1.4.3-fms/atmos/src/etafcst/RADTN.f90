    SUBROUTINE RADTN
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE RADTN
!>
!> SUBROUTINE: RADTN - THE OUTER RADIATION DRIVER
!> PROGRAMMER: BLACK
!> ORG: W/NP22
!> DATE: 93-12-??
!>
!> ABSTRACT:
!> RADTN PRIMARILY SERVES TO SET UP THE ARRAYS NEEDED AS INPUT FOR RADFS (THE INNER RADIATION 
!> DRIVER). GROUPS OF MODEL COLUMNS ARE SENT TO RADFS BUT FIRST THEY ARE "LIFTED" SO THAT THE LOWEST
!> LAYER ABOVE THE GROUND HAS A VERTICAL INDEX VALUE OF LM NOT LMH.
!> THIS ROUTINE IS CALLED AS OFTEN AS DESIRED (EVERY 1 TO 2 HOURS) FOR BOTH THE SHORT AND LONGWAVE 
!> EFFECTS. THE RESULTING TEMPERATURE TENDENCIES, TOTAL DOWNWARD AND SHORTWAVE UPWARD FLUXES ARE 
!> COLLECTED.
!> THE INITIAL GROUND POTENTIAL TEMPERATURE IS ALSO COMPUTED HERE AND IS SIMPLY AN ADIABATIC 
!> EXTRAPOLATION FROM THE LOWEST MID- LAYER VALUE ABOVE THE GROUND.
!>
!> PROGRAM HISTORY LOG:
!> 87-09-??  BLACK      - ORIGINATOR
!> 92-10-??  BALDWIN    - VARIOUS CLOUD EFFECTS WERE INCLUDED WHICH WERE ALREADY IN THE MRF
!> 93-11-??  ZHAO       - TIED TO UPDATED GFDL RADIATION SCHEME USING MODEL-PREDICTED CLOUD
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 95-04-13  BLACK      - PARALLELIZED THE LARGE LOOP STEPPING THROUGH THE DOMAIN THAT CALLS RADFS
!> 95-10-10  ZHAO       - I) THE CALCULATION OF CLOUD FRACTION WAS CHANGED TO USE BOTH CLOUD 
!>                        WATER/ICE MIXING RATIO AND RELATIVE HUMIDITY (RANDALL, 1994);
!>                        II) THE CLOUD INPUTS WERE CHANGED TO USE CLOUD FRACTION IN EACH MODEL 
!>                        LAYER AFTER Y.T. HOU (1995).
!> 96-06-03  ZHAO       - SNOW ALBEDO IS CHANGED ACCORDING TO SUGGESTIONS FROM KEN MITCHELL AND FEI
!>                        CHEN
!> 96-07-23  ZHAO       - ADD CALL TO SOLARD TO CALCULATE THE NON-DIMENSIONAL SUN-EARTH DISTANCE RAD1 
!>                        WHICH WILL BE USED IN RADFS TO COMPUTE SOLAR CONSTANT SOLC ON EACH DAY
!> 96-07-26  BLACK      - ADDED OZONE COMPUTATIONS
!> 97-05-19  ZHAO       - DIAGNOSTIC CLOUDS (LOW, MIDDLE, AND HIGH) ARE MODIFIED TO USE THE MAXIMUM
!>                        OF CONVECTIVE AND STRATIFORM. THIS WILL REPLACE THE PREVIOUS SCHEME WHICH
!>                        USES ONLY CONVECTIVE CLOUDS AT CONVECTIVE POINTS. THIS WILL AFFECT
!>                        CFRACL, CFRACM, CFRACH, AND WILL AFFECT THE TOTAL CLOUD FRACTION 
!>                        CALCULATION IN THE POST PROCESSORS.
!> 98-??-??  TUCCILLO   - ADDED PARALLELISM FOR CLASS VIII
!> 98-10-27  BLACK      - PARALLELISM INTO NEWEST VERSION
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> NOTE: CONVECTIVE CLOUDS ARE ADDED IN THIS SUBROUTINE FOR USE IN THE ETA MODEL IN WHICH 
!>       MODEL-PREDICTED CLOUDS ARE NOT USED IN THE CONVECTIVE PRECIPITATION PROCESSES. 
!>       FOR USE WITH THE VERSION OF THE ETA MODEL IN WHICH THE MODEL-PREDICTED CLOUDS ARE LINKED
!>       INTO THE MODEL'S CONVECTIVE PRECIPITATION PROCESSES, JUST SET: 
!>       CNCLD = .FALSE. 
!>       QINGYUN ZHAO  94-9-12
!>
!> INPUT ARGUMENT LIST:
!> NONE
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE 
!>
!> USE MODULES: ACMCLD
!>              ACMRDL
!>              ACMRDS
!>              CLDWTR
!>              CNVCLD
!>              CTLBLK
!>              CUINIT
!>              DYNAM
!>              F77KINDS
!>              GLB_TABLE
!>              INDX
!>              LOOPS
!>              MAPPINGS
!>              MASKS
!>              MPPCOM 
!>              PARMETA
!>              PARMSOIL
!>              PARM_TBL
!>              PHYS     
!>              PVRBLS     
!>              RD1TIM
!>              SOIL
!>              SWRSAV
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>
!> DRIVER     : EBU
!>              NEWFLT
!>
!> CALLS      : OZON2D
!>              RADFS
!>              ZENITH
!>              ZERO2
!>--------------------------------------------------------------------------------------------------
    USE ACMCLD
    USE ACMRDL
    USE ACMRDS
    USE CLDWTR
    USE CNVCLD
    USE CTLBLK
    USE CUINIT
    USE DYNAM
    USE F77KINDS
    USE GLB_TABLE
    USE INDX
    USE LOOPS
    USE MAPOT
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE PARMETA
    USE PARMSOIL
    USE PARM_TBL
    USE PHYS
    USE PVRBLS
    USE RD1TIM
    USE SOIL
    USE SWRSAV
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
#include "sp.h"
!
    REAL   (KIND=R4KIND), PARAMETER :: CAPA  =  0.28589641
    REAL   (KIND=R4KIND), PARAMETER :: RTD   = 57.2957795
    REAL   (KIND=R4KIND), PARAMETER :: WA    =   .10
    REAL   (KIND=R4KIND), PARAMETER :: WG    =  1. - WA
!
    INTEGER(KIND=I4KIND), PARAMETER :: KSMUD = 0
!------ 
! CLOUD 
!------
    REAL   (KIND=R4KIND), PARAMETER :: A1    = 610.78
    REAL   (KIND=R4KIND), PARAMETER :: A2    =  17.2693882
    REAL   (KIND=R4KIND), PARAMETER :: A3    = 273.16
    REAL   (KIND=R4KIND), PARAMETER :: A4    =  35.86
    REAL   (KIND=R4KIND), PARAMETER :: PQ0   = 379.90516
!------ 
! CLOUD 
!------
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM  = IM * JM - JM / 2

    REAL   (KIND=R4KIND), PARAMETER :: SLPM   =     1.01325E5
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ1  =     1.E-5 
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ   =     2.E-12
    REAL   (KIND=R4KIND), PARAMETER :: EPSO3  =     1.E-10
    REAL   (KIND=R4KIND), PARAMETER :: HPINC  =     1.E1
    REAL   (KIND=R4KIND), PARAMETER :: CLDRH0 =     0.80 
    REAL   (KIND=R4KIND), PARAMETER :: TRESH  =     1.00
    REAL   (KIND=R4KIND), PARAMETER :: RNRM   =     1. / (TRESH - CLDRH0)
    REAL   (KIND=R4KIND), PARAMETER :: CLDRH2 =     0.90
    REAL   (KIND=R4KIND), PARAMETER :: TRESH2 =     1.00
    REAL   (KIND=R4KIND), PARAMETER :: RNRM2  =     1. / (TRESH2 - CLDRH2)
    REAL   (KIND=R4KIND), PARAMETER :: CLAPSE =    -0.0005  
    REAL   (KIND=R4KIND), PARAMETER :: CLPSE  =    -0.0006 
    REAL   (KIND=R4KIND), PARAMETER :: DCLPS  =    -0.0001
    REAL   (KIND=R4KIND), PARAMETER :: CM1    =  2937.4
    REAL   (KIND=R4KIND), PARAMETER :: CM2    =     4.9283
    REAL   (KIND=R4KIND), PARAMETER :: CM3    =    23.5518
    REAL   (KIND=R4KIND), PARAMETER :: EPS    =     0.622
    REAL   (KIND=R4KIND), PARAMETER :: PBOT   = 10000.0
    REAL   (KIND=R4KIND), PARAMETER :: STBOL  =     5.67E-8
    REAL   (KIND=R4KIND), PARAMETER :: PI2    =     2. * 3.14159265
    REAL   (KIND=R4KIND), PARAMETER :: RLAG   =    14.8125
    REAL   (KIND=R4KIND), PARAMETER :: US     =     1.
    REAL   (KIND=R4KIND), PARAMETER :: CLIMIT =     1.0E-20
!
    INTEGER(KIND=I4KIND), PARAMETER :: K15 = SELECTED_REAL_KIND(15)
    REAL(K15)                                                                                   ::&
    & PROD    , DDX     , EEX
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & CALL1   , SHORT   , LONG    , BITX    , BITY    , BITZ    , BITW    , BIT1    , BIT2    ,   &
    & BITC    , BITS    , BITCP1  , BITSP1
!
    LOGICAL(KIND=L4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & BCLD    ,  BTEMP1
!----------- 
! CONVECTION
!----------- 
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & CNCLD
!----------- 
! CONVECTION
!-----------
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & ICVB    , ICVT     
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & PSFC    , TSKN    , ALBDO   , XLAT    , COSZ    , SLMSK   ,  CV     , SV      , FLWUP   ,   &
    & FSWDN   , FSWUP   , FSWDNS  , FSWUPS  , FLWDNS  , FLWUPS
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                                         ::&
    & TENDK
!
    REAL   (KIND=R4KIND), DIMENSION(0:LM)                                                       ::&
    & CLDAMT
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, 3)                                             ::&
    & CLDCFR
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, 3)                                             ::&
    & MBOT    , MTOP
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & CLDF
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LM)                                            ::&
    & TENDS   , TENDL   , PMID    , TMID    , QMID    , THMID   , OZN     , POZN    
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & PDSL    , FNE     , FSE     , TL      , PBOTL   , PTOPL   , PBOTM   , PTOPM   , PBOTH   ,   &
    & PTOPH   , TOT 
!
    REAL   (KIND=R4KIND), DIMENSION(9)                                                          ::&
    & CC      , PPT
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & CSTR    , TAUC    , CVB     , CVT     , TAUDAR
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & PINT    , EMIS
!
    REAL   (KIND=R4KIND), DIMENSION(LP1)                                                        ::&
    & PHALF
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & CAMT 
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & NCLDS   , KCLD
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & ITYP    , KTOP    , KBTM
! 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, NB, LP1)                                       ::&
    & RRCL    , TTCL
!------ 
! CLOUD 
!------
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, LM)                                            ::&
    & IW      
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LM)                                            ::&
    & CCR     , CSMID   , WMID    , HMID    ,CCMID
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & BMID    , UMID 
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NFILE   , I       , J       , NTSPH   , NRADPP  , ITIMSW  , ITIMLW  , JD      , II      ,   &
    & K       , IR      , N       , LML     , LVLIJ   , KNTLYR  , LL2     , IWKL    , LBASE   ,   &
    & L400    , LTROP   , NMOD    , NC      , NLVL    , MALVL   , LLTOP   , LLBOT   , KBT1    ,   &
    & KBT2    , KTH1    , KTH2    , KTOP1   , NBAND   , NCLD    , NKTP    , NBTM    , KS 
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & PLOMD   , PMDHI   , PHITP   , P400    , PLBTM   , UTIM    , CLSTP   , TIME    , DAYI    ,   &
    & HOUR    , ADDL    , RANG    , RCOS1   , RSIN1   , RCOS2   , EXNER   , TIMES   , APES    ,   &
    & HH      , TKL     , QKL     , CWMKL   , TMT0    , TMT15   , AI      , BI      , PP      ,   &
    & QW      , QI      , QINT    , U00KL   , FIQ     , FIW     , QC      , RQKL    , ARG     ,   &
    & CCR_DUM , DDP     , DTHDP   , CLPFIL  , PMOD    , CC1     , CC2     , P1      , P2      ,   &
    & CLDMAX  , CL1     , CL2     , CR1     , DCPL    , QSUM    , PRS1    , PRS2    , DELP    ,   &
    & TCLD    , DD      , EE      , FF      , AA      , BB      , GG      , DENOM   , FCTRA   ,   &
    & FCTRB   , PDSLIJ  , CFRAVG  , DPCL
!------ 
! CLOUD 
!------
    DATA PLOMD /64200./, PMDHI /35000./, PHITP /15000./, P400 /40000./, PLBTM /105000./
    DATA NFILE /14/
    DATA CC /0.  , 0.1, 0.2 , 0.3, 0.4, 0.5,  0.6,  0.7,  0.8/
    DATA PPT /.14,  .31, .70, 1.6, 3.4, 7.7, 17. , 38. , 85. /
!
    UTIM  = 1.
    CNCLD = .TRUE. 
!-------------------------------------------------
! ASSIGN THE PRESSURES FOR CLOUD DOMAIN BOUNDARIES
!-------------------------------------------------
    PTOPC(1) = PLBTM
    PTOPC(2) = PLOMD
    PTOPC(3) = PMDHI
    PTOPC(4) = PHITP
!------------------------------- 
! FIND THE 'SEA LEVEL PRESSURE'.
!------------------------------- 
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            PDSL(I,J) = RES(I,J) * PD(I,J)
        END DO
    END DO
!------------------------------------------------------------------
! THE FOLLOWING CODE IS EXECUTED EACH TIME THE RADIATION IS CALLED.
!------------------------------------------------------------------
!
!-----------  
! CONVECTION 
!-----------  
!--------------------------------------------------------------------------------------------------  
! NRADPP IS THE NUMBER OF TIME STEPS TO ACCUMULATE CONVECTIVE PRECIP FOR RADIATION
! NOTE: THIS WILL NOT WORK IF NRADS AND NRADL ARE DIFFERENT UNLESS THEY ARE INTEGER MULTIPLES OF 
!       EACH OTHER CLSTP IS THE NUMBER OF HOURS OF THE ACCUMULATION PERIOD
!-------------------------------------------------------------------------------------------------- 
    NTSPH  = NINT(3600. / DT)
    NRADPP = MIN(NRADS,NRADL)
    CLSTP  = 1.0 * NRADPP / NTSPH
!-----------  
! CONVECTION 
!-----------  
!
!-----------------------------------------------------------------     
! STATE WHETHER THE SHORT OR LONGWAVE COMPUTATIONS ARE TO BE DONE.
!-----------------------------------------------------------------
    SHORT = .FALSE. 
    LONG  = .FALSE. 
!
    IF (MOD(NTSD,NRADS) == 1) SHORT = .TRUE. 
    IF (MOD(NTSD,NRADL) == 1) LONG  = .TRUE. 
!
    IF (MYPE == 0) THEN
        IF (SHORT) THEN
            WRITE(0,*) 'RADTN: CALCULATE SHORTWAVE, NTSD', NTSD
            WRITE(6,*) 'RADTN: CALCULATE SHORTWAVE, NTSD', NTSD
        END IF
        IF (LONG) THEN
            WRITE(0,*) 'RADTN: CALCULATE LONGWAVE, NTSD', NTSD
            WRITE(6,*) 'RADTN: CALCULATE LONGWAVE, NTSD', NTSD
        END IF
    END IF
!
    ITIMSW = 0
    ITIMLW = 0
!
    IF (SHORT) ITIMSW=1
    IF (LONG)  ITIMLW=1
!--------------------------------------------- 
! FLAG FOR RESETTING CUPPT,HTOP,HBOT IN CHKOUT
!--------------------------------------------- 
    IF (MOD(NTSD,NRADPP) == 1) CURAD = .TRUE. 
!-------------------------------------------------------------------------------------------------- 
! FIND THE MEAN COSINE OF THE SOLAR ZENITH ANGLE BETWEEN THE CURRENT TIME AND THE NEXT TIME 
! RADIATION IS CALLED. 
! ONLY AVERAGE IF THE SUN IS ABOVE THE HORIZON.
!-------------------------------------------------------------------------------------------------- 
    TIME = (NTSD-1) * DT
!
    CALL ZENITH(TIME, DAYI, HOUR)
!
    JD = INT(DAYI + 0.50)
    ADDL = 0.
!
    IF (MOD(IDAT(3),4) == 0) ADDL = 1.
!
    RANG  = PI2 * (DAYI - RLAG) / (365.25 + ADDL)
    RSIN1 = SIN(RANG)
    RCOS1 = COS(RANG)
    RCOS2 = COS(2. * RANG)
!
    IF (SHORT) THEN
!
!$omp parallel do private (I      , J)
!
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                CZMEAN(I,J) = 0.
                TOT   (I,J) = 0.
            END DO
        END DO
!    
        DO II=0,NRADS,NPHS
            TIMES = (NTSD-1) * DT + II * DT
!
            CALL ZENITH(TIMES,DAYI,HOUR)
!
!$omp parallel do private (I      , J)
!
            DO J=MYJS,MYJE
                DO I=MYIS,MYIE
                    IF (CZEN(I,J) > 0.) THEN
                        CZMEAN(I,J) = CZMEAN(I,J) + CZEN(I,J)
                        TOT   (I,J) = TOT(I,J) + 1.
                    END IF
                END DO
            END DO
        END DO
!
!$omp parallel do private (I      , J)
!
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                IF(TOT(I,J) > 0.) CZMEAN(I,J) = CZMEAN(I,J) / TOT(I,J)
            END DO
        END DO
    END IF
!
!$omp parallel do
!$omp private (AA       , ALBDO   , APES    , BB      , BCLD    , BIT1    , BIT2    , BITC    ,   &
!$omp          BITCP1   , BITS    , BITSP1  , BITW    , BITX    , BITY    , BITZ    , BMID    ,   &
!$omp          BTEMP1   , CAMT    , CC1     , CC2     , CCMID   , CCR     , CFRAVG  , CL1     ,   &
!$omp          CL2      , CLDAMT  , CLDCFR  , CLDMAX  , CLPFIL  , COSZ    , CR1     , CSMID   ,   &
!$omp          CSTR     , CV      , CWMKL   , DD      , DDP     , DELP    , DENOM   , DPCL    ,   &
!$omp          DTHDP    , EE      , EMIS    , EXNER   , FCTRA   , FCTRB   , FF      , FIQ     ,   &
!$omp          FIW      , FLWDNS  , FLWUP   , FLWUPS  , FSWDN   , FSWDNS  , FSWUP   , FSWUPS  ,   &
!$omp          GG       , HH      , HMID    , I       , ICVB    , ICVT    , IR      , ITYP    ,   &
!$omp          IW       , IWKL    , J       , KBT1    , KBT2    , KBTM    , KCLD    , KNTLYR  ,   &
!$omp          KTH1     , KTH2    , KTOP    , KTOP1   , K       , L400    , LBASE   , LIN     ,   &
!$omp          LL2      , LLBOT   , LLTOP   , LML     , LTROP   , LVLIJ   , MALVL   , MBOT    ,   &
!$omp          MTOP     , N       , NBAND   , NBTM    , NC      , NCLD    , NCLDS   , NKTP    ,   &
!$omp          NLVL     , NMOD    , OZN     , P1      , P2      , PDSLIJ  , PINT    , PMID    ,   &
!$omp          PMOD     , POZN    , PP      , PRS1    , PRS2    , PSFC    , QC      , QI      ,   &
!$omp          QINT     , QKL     , QMID    , QSUM    , QW      , RQKL    , RRCL    , SLMSK   ,   &
!$omp          SNOFAC   , SV      , TAUC    , TAUDAR  , TCLD    , TENDL   , TENDS   , THMID   ,   &
!$omp          TKL      , TMID    , TMT0    , TMT15   , TSKN    , TTCL    , U00KL   , UMID    ,   &
!$omp          WMID     , XLAT    )
!
!-------------------------------------------------------------
! THIS IS THE BEGINNING OF THE PRIMARY LOOP THROUGH THE DOMAIN
!-------------------------------------------------------------
!
    DO 700 J=MYJS,MYJE
!
        DO 125 K=1,LM
            DO I=MYIS,MYIE
                IR = IRAD(I)
                TMID (I,K) = T(I,J,1)
                QMID (I,K) = EPSQ
                CSMID(I,K) = 0.
                WMID (I,K) = 0.
                CCMID(I,K) = 0.
                IW   (I,K) = 0.
                CCR  (I,K) = 0.
                HMID (I,K) = 0.
                OZN  (I,K) = EPSO3
                TENDS(I,K) = 0.
                TENDL(I,K) = 0.
            END DO
    125 END DO
!    
        DO 140 N=1,3
            DO I=MYIS,MYIE
                CLDCFR(I,N) = 0.
                MTOP  (I,N) = 0
                MBOT  (I,N) = 0
            END DO
    140 END DO
!--------------------------------------------------------------------------------------------
! FILL IN WORKING ARRAYS WHERE VALUES AT L=LM ARE THOSE THAT ARE ACTUALLY AT ETA LEVEL L=LMH.
!--------------------------------------------------------------------------------------------
        DO 200 I=MYIS,MYIE
            IR      = IRAD(I)
            LML     = LMH (I,J)
            LVLIJ   = LVL (I,J)
            BMID(I) = HBM2(I,J)
            UMID(I) = U00 (I,J)
!        
            DO K=1,LML
                 PMID(I,K+LVLIJ)   = AETA(K)   * PDSL(I,J) + PT
                 PINT(I,K+LVLIJ+1) = ETA(K+1)  * PDSL(I,J) + PT
                EXNER              = (1.E5 / PMID(I,K+LVLIJ)) ** CAPA
                 TMID(I,K+LVLIJ)   =   T(I,J,K)
                THMID(I,K+LVLIJ)   =   T(I,J,K) * EXNER
                 QMID(I,K+LVLIJ)   =   Q(I,J,K)
                 WMID(I,K+LVLIJ)   = CWM(I,J,K)
                 HMID(I,K+LVLIJ)   = HTM(I,J,K)
            END DO
!--------------------------------------------------------------- 
! FILL IN ARTIFICIAL VALUES ABOVE THE TOP OF THE DOMAIN.
! PRESSURE DEPTHS OF THESE LAYERS IS 1 HPA.
! TEMPERATURES ABOVE ARE ALREADY ISOTHERMAL WITH (TRUE) LAYER 1.
!---------------------------------------------------------------
            IF (LVLIJ > 0) THEN
                KNTLYR = 0
!            
                DO K=LVLIJ,1,-1
                    KNTLYR       = KNTLYR + 1
                    PMID (I,K  ) = PT - REAL(2 * KNTLYR - 1) * 0.5 * HPINC
                    PINT (I,K+1) = PMID(I,K) + 0.5 * HPINC
                    EXNER        = (1.E5 / PMID(I,K)) ** CAPA
                    THMID(I,K  ) = TMID(I,K) * EXNER
                END DO
            END IF
!        
            IF (LVLIJ == 0) THEN
                PINT(I,1)=PT
            ELSE
                PINT(I,1) = PMID(I,1) - 0.5 * HPINC
            END IF
    200 END DO
!--------------------------------------------------------------------------------------------------
! FILL IN THE SURFACE PRESSURE, SKIN TEMPERATURE, GEODETIC LATITUDE, ZENITH ANGLE, SEA MASK, AND 
! ALBEDO. THE SKIN TEMPERATURE IS NEGATIVE OVER WATER.
!--------------------------------------------------------------------------------------------------
        DO 250 I=MYIS,MYIE
             PSFC(I)  =  PD(I,J) + PT
            APES      = (PSFC(I) * 1.E-5) ** CAPA
             TSKN(I)  = THS(I,J) * APES * (1. - 2. * SM(I,J))
            SLMSK(I)  =  SM(I,J)
! -------------------------------------------------------------------- 
! TURN OFF SNOW ALBEDO CALCULATION SINCE IT IS NOW CALCULATED IN SFLX.
! -------------------------------------------------------------------- 
            ALBDO(I) = ALBEDO(I,J)
!        
             XLAT(I) =   GLAT(I,J) * RTD
             COSZ(I) = CZMEAN(I,J)
    250 END DO
! ------------------------ 
! STRATIFORM CLOUD SECTION 
! ------------------------ 
!
! ------------------------------------------------------------------------------------------------- 
! CALCULATE STRATIFORM CLOUD COVERAGE AT EACH MODEL GRID POINT WHICH WILL BE USED IN THE MODEL 
! RADIATION PARAMETERIZATION SCHEME.
! ------------------------------------------------------------------------------------------------- 
!
! ---------------  
! QW, QI AND QINT 
! --------------- 
        DO 280 I=MYIS,MYIE
            LML   = LMH(I,J)
            LVLIJ = LVL(I,J)
!        
            DO 275 K=1,LML
                LL2   = K + LVLIJ
                HH    = HMID(I,LL2) * BMID(I)
                TKL   = TMID(I,LL2)
                QKL   = QMID(I,LL2)
                CWMKL = WMID(I,LL2)
                TMT0  = (TKL - 273.16)   * HH
                TMT15 = AMIN1(TMT0,-15.) * HH
                AI    = 0.008855
                BI    = 1.
!            
                IF (TMT0 < -20.) THEN
                    AI = 0.007225
                    BI = 0.9674
                END IF
!            
                PP   = PMID(I,LL2)
                QW   = HH * PQ0 / PP * EXP(HH * A2 * (TKL - A3) / (TKL - A4))
                QI   = QW * (BI + AI * AMIN1(TMT0, 0.))
                QINT = QW * (1. - 0.00032 * TMT15 * (TMT15 + 15.))
!
                IF (TMT0 <= -40.) QINT = QI
! ----------------------  
! ICE-WATER ID NUMBER IW 
! ---------------------- 
                U00KL = UMID(I) + UL(K) * (0.95 - UMID(I)) * UTIM
                IF (TMT0 < -15.0) THEN
                    FIQ = QKL - U00KL * QI
                    IF (FIQ > 0. .OR. CWMKL > CLIMIT) THEN
                        IW(I,LL2) = 1
                    ELSE
                        IW(I,LL2) = 0
                    END IF
                END IF
!            
                IF (TMT0 >= 0.) THEN
                    IW(I,LL2) = 0
                END IF
!           
                IF (TMT0 < 0.0 .AND. TMT0 >= -15.0) THEN
                    IW(I,LL2) = 0
                    IF (IW(I,LL2-1) == 1 .AND. CWMKL > CLIMIT) IW(I,LL2) = 1
                END IF
!           
                IWKL = IW(I,LL2)
! --------------------------------             
! THE SATURATION SPECIFIC HUMIDITY 
! --------------------------------             
                FIW = FLOAT(IWKL)
                QC  = (1. - FIW) * QINT + FIW * QI
! ---------------------              
! THE RELATIVE HUMIDITY 
! ---------------------             
                IF (QC <= EPSQ1 .OR. QKL <= EPSQ1) THEN
                    RQKL = 0.
                ELSE
                    RQKL = QKL / QC
                END IF
! ---------------------             
! CLOUD COVER RATIO CCR 
! ---------------------             
                IF (RQKL >= 0.9999) THEN
                    CCR(I,LL2) = AMIN1(US,RQKL) 
! ------------------------------------------------------------------------------------------------- 
! FERRIER, 6/17/02: EMERGENCY CHANGE TO ELIMINATE VERY SMALL CLOUD AMOUNTS (SOON TO BE REPLACED BY 
!                   CORRECTIONS FROM ETAZ PARALLEL)
! -------------------------------------------------------------------------------------------------               
                ELSE IF (CWMKL > 1.E-7) THEN
                    ARG     = -1000. * CWMKL / (US - RQKL)
                    ARG     = AMAX1(ARG,-25.)
                    CCR_DUM = RQKL * (1. - EXP(ARG))
                    IF (CCR_DUM > .01) CCR(I,LL2) = CCR_DUM
                END IF
!
                CSMID(I,LL2) = AMIN1(US, CCR(I,LL2))
!
        275 END DO
    280 END DO
! ------------------------------------------------------------------------------------------------
! NOW CHECK THE CLOUDS PRODUCED ABOVE TO MAKE SURE THEY ARE GOOD ENOUGH FOR RADIATION CALCULATIONS
! ------------------------------------------------------------------------------------------------
!
! ---------------------------------- 
! NO STRATIFORM CLOUDS FOR THIS TYPE
! ---------------------------------- 
        DO 350 I=MYIS,MYIE
!        
            LML   = LMH(I,J)
            LVLIJ = LVL(I,J)
! --------------------------------------------------- 
! ZERO OUT CLDAMT IF LAND AND BELOW PBOT ABOVE GROUND
! --------------------------------------------------- 
            IF (SM(I,J) < 0.5) THEN
                DO K=1,LML
                    LL2 = LML - K + 1 + LVLIJ
                    DDP = PSFC(I) - PMID(I,LL2)
!
                    IF (DDP >= PBOT) GOTO 290
                    CSMID(I,LL2) = 0.
                END DO
!
                290 CONTINUE
            END IF
! -------------------------------------------------------------------------------------------------  
! CHECK FOR OCEAN STRATUS (LOW CLOUD) LOOK ONLY OVER OCEAN AND ONLY IF AN INVERSION 
! (DTHDP.LE.-0.05) IS PRESENT WITH AT LEAST 2 CLOUD FREE LAYERS ABOVE IT
! ------------------------------------------------------------------------------------------------- 
            IF (SM(I,J) > 0.5) THEN
! ----------------------            
! FIND BASE OF INVERSION
! ----------------------              
                LBASE = LM
!
                DO K=1,LML-1
                    LL2   = LML - K + 1 + LVLIJ
                    DTHDP = (THMID(I,LL2-1) - THMID(I,LL2)) / (PMID(I,LL2-1) - PMID(I,LL2))
!
                    IF (DTHDP <= CLAPSE) THEN
                        LBASE = LL2
                        GOTO 300
                    END IF
                END DO
!
            300 CONTINUE
! --------------------------------------             
! CHECK 2 LAYERS ABOVE LBASE FOR DRYNESS
! --------------------------------------             
                IF (CSMID(I,LBASE-1) <= 0. .AND. CSMID(I,LBASE-2) <= 0. .AND. LBASE < LM) THEN
                    IF (DTHDP > CLPSE) THEN
                        CLPFIL = 1. - ((CLPSE - DTHDP) / DCLPS)
                    ELSE
                        CLPFIL = 1.
                    END IF
!                
                    DO K=1,LML
                        LL2 = LML - K + 1 + LVLIJ
                        DDP = PSFC(I) - PMID(I,LL2)
!
                        IF (DDP >= PBOT) GO TO 310
                        CSMID(I,LL2) = CSMID(I,LL2) * CLPFIL
                    END DO
!
            310 CONTINUE
! -------------------------------------------------------------------------------------------------                  
! IF NO INVERSION OR IF CLDS EXIST IN EITHER OF THE 2 LAYERS ABOVE INVERSION, ZERO OUT CLOUD BELOW 
! PBOT
! -------------------------------------------------------------------------------------------------                
                ELSE
                    DO K=1,LML
                        LL2  = LML - K + 1 + LVLIJ
                        DDP = PSFC(I) - PMID(I,LL2)
!
                        IF (DDP >= PBOT) GOTO 320
!
                        CSMID(I,LL2) = 0.
                    END DO
!
                320 CONTINUE
!                
                END IF
!
            END IF
! ---------------------------------------  
! REMOVE HIGH CLOUDS ABOVE THE TROPOPAUSE
! --------------------------------------- 
            L400 = LM
            DO K=1,LML
                LL2 = LML - K + 1 + LVLIJ
                IF (PMID(I,LL2) <= 40000.0) THEN
                    L400 = LL2
                    GOTO 330
                END IF
            END DO
!
        330 CONTINUE
!        
            LTROP = LM
            DO LL2=L400,2,-1
                DTHDP = (THMID(I,LL2-1) - THMID(I,LL2)) / (PMID(I,LL2-1) - PMID(I,LL2))
!            
                IF (DTHDP < -0.0025 .OR. QMID(I,LL2) <= EPSQ1) THEN
                    LTROP = LL2
                    GOTO 340
                END IF
            END DO
!
        340 IF (LTROP < LM) THEN
                DO LL2=LTROP,1,-1
                    CSMID(I,LL2) = 0.
                END DO
            END IF
350 END DO
!--------------------------------  
! END OF STRATIFORM CLOUD SECTION
!--------------------------------
!
!----------- 
! CONVECTION 
!-----------
!
!--------------------------------------------------------------------------------------------------
! CONVECTIVE CLOUD SECTION
!
! THIS PART WAS MODIFIED TO COMPUTE CONVECTIVE CLOUDS AT EACH MODEL LAYER BASED ON CONVECTIVE 
! PRECIPITATION RATES. 
! CURRENTLY, CLOUDS ARE SET TO 0.75*CV(I) BELOW 400MB AND 0.90*CV(I) ABOVE 400MB TO ACCOUNT FOR 
! CIRRUS CAP Q.ZHAO   95-3-22
!
! NON-PRECIPITATING CLOUD FRACTION OF 20 PERCENT IS ADDED AT AT POINTS WHERE THE SHALLOW AND DEEP 
! CONVECTIONS ACCUR.
! Q. ZHAO  97-5-2
!   
!COMPUTE THE CONVECTIVE CLOUD COVER FOR RADIATION
!--------------------------------------------------------------------------------------------------    
        IF (CNCLD) THEN
!
            DO 375 I=MYIS,MYIE
                IF (HBOT(I,J) - HTOP(I,J) > 1.0) THEN
                    SV(I) = 0.0
                ELSE
                    SV(I) = 0.0
                END IF
!            
                PMOD = CUPPT(I,J) * 24.0 * 1000.0 / CLSTP
                NMOD = 0
!            
                DO NC=1,9
                    IF (PMOD > PPT(NC)) NMOD = NC
                END DO
!---------------------------------------------------            
! CLOUD TOPS AND BOTTOMS COME FROM CUCNVC
! ADD LVL TO BE CONSISTENT WITH OTHER WORKING ARRAYS
!---------------------------------------------------              
                IF (NMOD == 0) THEN
                    CV(I) = 0.
                ELSE IF (NMOD == 9) THEN
                    CV(I) = CC(9)
                ELSE
                    CC1 = CC(NMOD  )
                    CC2 = CC(NMOD+1)
!
                    P1 = PPT(NMOD  )
                    P2 = PPT(NMOD+1)
!
                    CV(I) = CC1 + (CC2 - CC1) * (PMOD - P1) / (P2 - P1)
                END IF
!            
                CV(I) = AMAX1(SV(I), CV(I))
                CV(I) = AMIN1(1.0,   CV(I))
!            
                IF (CV(I) == 0.0) THEN
                    ICVT(I) = 0
                    ICVB(I) = 0
                ELSE
                    ICVT(I) = INT(HTOP(I,J) + 0.50) + LVL(I,J)
                    ICVB(I) = INT(HBOT(I,J) + 0.50) + LVL(I,J)
                END IF
        375 END DO
!---------------------------------   
! MAKE SURE CLOUDS ARE DEEP ENOUGH
!--------------------------------- 
            DO I=MYIS,MYIE
                BCLD  (I) = CV(I) > 0. .AND. (ICVB(I) - ICVT(I)) >= 1
                BTEMP1(I) = BCLD(I)
            END DO
!---------------------------------- 
! COMPUTE CONVECTIVE CLOUD FRACTION
!---------------------------------- 
            DO 390 I=MYIS,MYIE
                IF (BCLD(I)) THEN
                    LML   = LMH(I,J)
                    LVLIJ = LVL(I,J)
!                
                    DO K=1,LML
                        LL2 = K + LVLIJ
                        IF (LL2 > ICVB(I) .OR. LL2 < ICVT(I)) THEN
                            CCMID(I,LL2) = 0.
                        ELSE
                            CCMID(I,LL2) = CV(I)
                        END IF
                        CCMID(I,LL2) = AMIN1(1.0, CCMID(I,LL2))
                    END DO
                END IF
        390 END DO
!---------------------------------------- 
! REMOVE HIGH CLOUDS ABOVE THE TROPOPAUSE
!---------------------------------------- 
            L400 = LM
!
            DO 425 I=MYIS,MYIE
                LML   = LMH(I,J)
                LVLIJ = LVL(I,J)
!            
                DO K = 1, LML
                    LL2 = LML - K + 1 + LVLIJ
                    IF (PMID(I,LL2) <= 40000.0) THEN
                        L400 = LL2
                        GOTO 400
                    END IF
                END DO
!
            400 CONTINUE
!            
                LTROP = LM
!
                DO LL2=L400,2,-1
                    DTHDP = (THMID(I,LL2-1) - THMID(I,LL2)) / (PMID(I,LL2-1) - PMID(I,LL2))
                    IF (DTHDP < -0.0025 .OR. QMID(I,LL2) <= EPSQ1) THEN
                        LTROP = LL2
                        GOTO 410
                    END IF
                END DO
!
            410 IF (LTROP < LM) THEN
                    DO LL2=LTROP,1,-1
                        CCMID(I,LL2) = 0.
                    END DO
                END IF
        425 END DO
!
    END IF
!-----------  
! CONVECTION 
!-----------    
!
!--------------------------------
! END OF CONVECTIVE CLOUD SECTION 
!--------------------------------
!
!--------------------------------------------------------------------------------------------------
! DETERMINE THE FRACTIONAL CLOUD COVERAGE FOR HIGH, MID AND LOW OF CLOUDS FROM THE CLOUD COVERAGE 
! AT EACH LEVEL
!
! NOTE: THIS IS FOR DIAGNOSTICS ONLY 
!--------------------------------------------------------------------------------------------------
        DO 500 I=MYIS,MYIE
!        
            CSTR(I) = 0.0
!        
            DO K=0,LM
                CLDAMT(K) = 0.
            END DO
!---------------------------         
! NOW GOES LOW, MIDDLE, HIGH
!---------------------------          
            DO 480 NLVL=1,3
                CLDMAX = 0.
                MALVL  = LM
                LLTOP  = LTOP(NLVL) + LVL(I,J)
!--------------------------------------------------------------------------------------------------  
! GO TO THE NEXT CLOUD LAYER IF THE TOP OF THE CLOUD-TYPE IN QUESTION IS BELOW GROUND OR IS IN THE 
! LOWEST LAYER ABOVE GROUND.
!-------------------------------------------------------------------------------------------------- 
                IF (LLTOP >= LM) GOTO 480
!            
                IF (NLVL > 1) THEN
                    LLBOT = LTOP(NLVL-1) - 1 + LVL(I,J)
                    LLBOT = MIN(LLBOT, LM1)
                ELSE
                    LLBOT = LM1
                END IF
!            
                DO 435 K=LLTOP,LLBOT
                    CLDAMT(K) = AMAX1(CSMID(I,K), CCMID(I,K))
                    IF (CLDAMT(K) > CLDMAX) THEN
                        MALVL  = K
                        CLDMAX = CLDAMT(K)
                    END IF
            435 END DO
!-------------------------------------------------------------------------------------------------- 
! NOW, CALCULATE THE TOTAL CLOUD FRACTION IN THIS PRESSURE DOMAIN USING THE METHOD DEVELOPED BY 
! Y.H., K.A.C. AND A.K. (NOV., 1992).
! IN THIS METHOD, IT IS ASSUMED THAT SEPERATED CLOUD LAYERS ARE RADOMLY OVERLAPPED AND ADJACENT 
! CLOUD LAYERS ARE MAXIMUM OVERLAPPED.
! VERTICAL LOCATION OF EACH TYPE OF CLOUD IS DETERMINED BY THE THICKEST CONTINUING CLOUD LAYERS 
! IN THE DOMAIN.
!--------------------------------------------------------------------------------------------------
                CL1  = 0.0
                CL2  = 0.0
                KBT1 = LLBOT
                KBT2 = LLBOT
                KTH1 = 0
                KTH2 = 0
!            
                DO 450 LL2=LLTOP,LLBOT
                    K = LLBOT - LL2 + LLTOP
                    BIT1 = .FALSE. 
                    CR1 = CLDAMT(K)
                    BITX = (PINT(I,K) >= PTOPC(NLVL+1)) .AND. (PINT(I,K) < PTOPC(NLVL)) .AND.     &
    &                      (CLDAMT(K) > 0.0)
                    BIT1 = BIT1 .OR. BITX
!
                    IF (.NOT. BIT1) GOTO 450
!-------------------------------------------------------------------------------------------------- 
! BITY=T: FIRST CLOUD LAYER; BITZ=T:CONSECUTIVE CLOUD LAYER
! NOTE:  WE ASSUME THAT THE THICKNESS OF EACH CLOUD LAYER IN THE DOMAIN IS LESS THAN 200 MB TO 
! AVOID TOO MUCH COOLING OR HEATING. 
! SO WE SET CTHK(NLVL)=200*E2. BUT THIS LIMIT MAY WORK WELL FOR CONVECTIVE CLOUDS. MODIFICATION 
! MAY BE NEEDED IN THE FUTURE.
!--------------------------------------------------------------------------------------------------
                    BITY = BITX .AND. (KTH2 <= 0)
                    BITZ = BITX .AND. (KTH2 >  0)
!                
                    IF (BITY) THEN
                        KBT2 = K
                        KTH2 = 1
                    END IF
!                
                    IF (BITZ) THEN
                        KTOP1 = KBT2 - KTH2 + 1
                        DPCL  = PMID(I,KBT2) - PMID(I,KTOP1)
                        IF (DPCL < CTHK(NLVL)) THEN
                            KTH2 = KTH2 + 1
                        ELSE
                            KBT2 = KBT2 - 1
                        END IF
                    END IF
!
                    IF (BITX) CL2 = AMAX1(CL2, CR1)
!--------------------------------------------------------------------------------- 
! AT THE DOMAIN BOUNDARY OR SEPARATED CLD LAYERS, RANDOM OVERLAP.
! CHOOSE THE THICKEST OR THE LARGEST FRACTION AMT AS THE CLD LAYER IN THAT DOMAIN.
!--------------------------------------------------------------------------------- 
                    BIT2 = .FALSE. 
                    BITY = BITX .AND. (CLDAMT(K-1) <= 0.0 .OR. PINT(I,K-1) < PTOPC(NLVL+1))
                    BITZ = BITY .AND. CL1 > 0.0
                    BITW = BITY .AND. CL1 <= 0.0
                    BIT2 = BIT2 .OR.  BITY
!
                    IF (.NOT. BIT2) GOTO 450
!                
                    IF (BITZ) THEN
                        KBT1 = INT((CL1 * KBT1 + CL2 * KBT2) / (CL1 + CL2))
                        KTH1 = INT((CL1 * KTH1 + CL2 * KTH2) / (CL1 + CL2)) + 1
!
                        CL1 = CL1 + CL2 - CL1 * CL2
                    END IF
!                
                    IF (BITW) THEN
                        KBT1 = KBT2
                        KTH1 = KTH2
                        CL1  = CL2
                    END IF
!                
                    IF (BITY) THEN
                        KBT2 = LLBOT
                        KTH2 = 0
                        CL2  = 0.0
                    END IF
!
            450 END DO
!
                CLDCFR(I,NLVL) = AMIN1(1.0, CL1)
                MTOP  (I,NLVL) = MIN(KBT1, KBT1 - KTH1 + 1)
                MBOT  (I,NLVL) = KBT1
!
        480 END DO
!
    500 END DO
!-------------------------------- 
! SET THE UN-NEEDED TAUDAR TO ONE
!-------------------------------- 
        DO I=MYIS,MYIE
            TAUDAR(I) = 1.0
        END DO
!--------------------------------------------------------------------------------------------------  
! NOW, CALCULATE THE CLOUD RADIATIVE PROPERTIES AFTER DAVIS (1982), HARSHVARDHAN ET AL (1987)
! AND Y.H., K.A.C. AND A.K. (1993).
!    
! UPDATE: THE FOLLOWING PARTS ARE MODIFIED, AFTER Y.T.H. (1994), TO CALCULATE THE RADIATIVE 
!         PROPERTIES OF CLOUDS ON EACH MODEL LAYER BOTH CONVECTIVE AND STRATIFORM CLOUDS ARE 
!         USED IN THIS CALCULATIONS.
!   
! QINGYUN ZHAO   95-3-22
!--------------------------------------------------------------------------------------------------
!   
!--------------------------------- 
! INITIALIZE ARRAYS FOR USES LATER
!--------------------------------- 
        DO 600 I=MYIS,MYIE
            LML   = LMH(I,J)
            LVLIJ = LVL(I,J)
!------------------------------------------------------------------------------------------------
! NOTE: LAYER=1 IS THE SURFACE, AND LAYER=2 IS THE FIRST CLOUD LAYER ABOVE THE SURFACE AND SO ON.
!------------------------------------------------------------------------------------------------
            EMIS(I,1) = 1.0
            KTOP(I,1) = LP1
            KBTM(I,1) = LP1
            CAMT(I,1) = 1.0
            ITYP(I,1) = 0
            KCLD(I)   = 2
!        
            DO NBAND=1,NB
                RRCL(I,NBAND,1) = 0.0
                TTCL(I,NBAND,1) = 1.0
            END DO
!        
            DO 510 K=2,LP1
                ITYP(I,K) = 0
                CAMT(I,K) = 0.0
                KTOP(I,K) = 1
                KBTM(I,K) = 1
                EMIS(I,K) = 0.0
!            
                DO NBAND=1,NB
                    RRCL(I,NBAND,K) = 0.0
                    TTCL(I,NBAND,K) = 1.0
                END DO
!
        510 END DO
!--------------------------------------------------------------------------------------------------
! NOW CALCULATE THE AMOUNT, TOP, BOTTOM AND TYPE OF EACH CLOUD LAYER
! CLOUD TYPE=1: STRATIFORM CLOUD
!       TYPE=2: CONVECTIVE CLOUD
! WHEN BOTH CONVECTIVE AND STRATIFORM CLOUDS EXIST AT THE SAME POINT, SELECT CONVECTIVE CLOUD 
! (TYPE=2), IN OTHER WORDS, CONVECTIVE CLOUDS HAVE THE HIGHER PRIORITY THAN STRATIFORM CLOUDS.
! CLOUD LAYERS ARE SEPARATED BY:
! 1. NO-CLOUD LAYER
! 2. DIFFERENT CLOUD TYPE
! NOTE: THERE IS ONLY ONE CONVECTIVE CLOUD LAYER IN ONE COLUMN.
! KTOP AND KBTM ARE THE TOP AND BOTTOM OF EACH CLOUD LAYER IN TERMS O ETA MODEL LEVEL.
!--------------------------------------------------------------------------------------------------
            DO 540 K=2,LML
                LL2 = LML - K+1 + LVLIJ
                BITC   = CCMID(I,LL2) > 0.1
                BITS   = CSMID(I,LL2) > 0.1
                BITCP1 = CCMID(I,LL2+1) > 0.1
                BITSP1 = CSMID(I,LL2+1) > 0.1
                BIT1   = BITS .OR. BITC
!
                IF (BIT1) THEN
!
                    IF (ITYP(I,KCLD(I)) == 0) THEN
                        CAMT(I,KCLD(I)) = CSMID(I,LL2)
                        ITYP(I,KCLD(I)) = 1
                        KBTM(I,KCLD(I)) = LL2
!                    
                        IF (BITC) THEN
                            CAMT(I,KCLD(I)) = CCMID(I,LL2)
                            ITYP(I,KCLD(I)) = 2
                        END IF
                    ELSE
                        IF (BITC) THEN
                            IF (BITCP1) THEN
                                CAMT(I,KCLD(I)) = AMAX1(CAMT(I, KCLD(I)), CCMID(I,LL2))
                            ELSE
                                KCLD(I)           = KCLD(I) + 1
                                CAMT(I,KCLD(I))   = CCMID(I,LL2)
                                ITYP(I,KCLD(I))   = 2
                                KTOP(I,KCLD(I)-1) = LL2 + 1
                                KBTM(I,KCLD(I))   = LL2
                            END IF
                        ELSE
                            IF (BITCP1) THEN
                                KCLD(I)           = KCLD(I) + 1
                                CAMT(I,KCLD(I))   = CSMID(I,LL2)
                                ITYP(I,KCLD(I))   = 1
                                KTOP(I,KCLD(I)-1) = LL2 + 1
                                KBTM(I,KCLD(I))   = LL2
                            ELSE
                                CAMT(I,KCLD(I)) = AMAX1(CAMT(I, KCLD(I)), CSMID(I,LL2))
                            END IF
                        END IF
                    END IF
!
                ELSE
!
                    IF (BITCP1 .OR. BITSP1) THEN
                        KCLD(I)           = KCLD(I) + 1
                        KTOP(I,KCLD(I)-1) = LL2 + 1
                        ITYP(I,KCLD(I)  ) = 0
                        CAMT(I,KCLD(I)  ) = 0.0
                    END IF
!
                END IF
!
        540 END DO
!-----------------------------------------------------------------------------------
! THE REAL NUMBER OF CLOUD LAYERS IS (THE FIRST IS THE GROUNG; THE LAST IS THE SKY):
!-----------------------------------------------------------------------------------
            NCLDS(I) =  KCLD(I) - 2
            NCLD     = NCLDS(I)
!----------------------------------------- 
! NOW CALCULATE CLOUD RADIATIVE PROPERTIES
!----------------------------------------- 
            IF (NCLD >= 1) THEN
!-------------------------------------------------------------- 
! NOTE: THE FOLLOWING CALCULATIONS, THE UNIT FOR PRESSURE IS MB 
!-------------------------------------------------------------- 
                DO 580 NC=2,NCLD+1
!                
                    TAUC(I) = 0.0
                    QSUM    = 0.0
                    NKTP    = LP1
                    NBTM    = 0
                    BITX    = CAMT(I,NC) > 0.1
                    NKTP    = MIN(NKTP, KTOP(I,NC))
                    NBTM    = MAX(NBTM, KBTM(I,NC))
!                
                    DO 560 LL2=NKTP,NBTM
                        IF (LL2 >= KTOP(I,NC) .AND. LL2 <= KBTM(I,NC) .AND. BITX) THEN
                            PRS1 = PINT(I,LL2)   * 0.01
                            PRS2 = PINT(I,LL2+1) * 0.01
                            DELP = PRS2 - PRS1
                            TCLD = TMID(I,LL2) - 273.16
                            QSUM = QSUM  + QMID(I,LL2) * DELP                                     &
    &                            * (PRS1 + PRS2) / (120.1612 * SQRT(TMID(I,LL2)))
!-------------------------------------------------------------- 
! FOR CONVECTIVE CLOUD OR STARTIFORM CLOUD WITH TOP ABOVE 500MB
!-------------------------------------------------------------- 
                            IF (ITYP(I,NC) == 2 .OR. PINT(I,KTOP(I,NC)) <= PTOPC(3)) THEN
                                IF (TCLD <= -10.0) THEN
                                    TAUC(I) = TAUC(I) + DELP                                      &
    &                                       * AMAX1(0.1E-3, 2.0E-6 * (TCLD + 82.5) ** 2)
                                ELSE
                                    TAUC(I) = TAUC(I) + DELP * AMIN1(0.08, 6.949E-3 * TCLD + 0.1)
                                END IF
                            ELSE
!----------------------------------  
! FOR LOW AND MID STRATIFORM CLOUDS
!----------------------------------
                                IF (TCLD <= -20.0) THEN
                                    TAUC(I) = TAUC(I) + DELP                                      &
    &                                       * AMAX1(0.1E-3, 2.56E-5 * (TCLD + 82.5) ** 2)
                                ELSE
                                    TAUC(I) = TAUC(I) + DELP * 0.1
                                END IF
                            END IF
                        END IF
!
                560 END DO
!                
                    IF (BITX) EMIS(I,NC) = 1.0 - EXP(-0.75 * TAUC(I))
!
                    IF (QSUM >= EPSQ1) THEN
                        DO 570 NBAND=1,NB
                            IF (BITX) THEN
                                PROD = ABCFF(NBAND) * QSUM
                                DDX  = TAUC(I) / (TAUC(I) + PROD)
                                EEX  = 1.0 - DDX
                                IF (ABS(EEX) >= 1.E-8) THEN
                                    DD = DDX
                                    EE = EEX
                                    FF = 1.0 - DD * 0.85
                                    AA = MIN(50.0,SQRT(3.0 * EE * FF) * TAUC(I))
                                    AA = EXP(-AA)
                                    BB = FF / EE
                                    GG = SQRT(BB)
                                    DD = (GG + 1.0) * (GG + 1.0) - (GG - 1.0)                     &
    &                                  * (GG - 1.0) * AA * AA
                                    RRCL(I,NBAND,NC) =   MAX(0.1E-5, (BB-1.0) * (1.0-AA*AA)/DD)
                                    TTCL(I,NBAND,NC) = AMAX1(0.1E-5, 4.0 * GG * AA/DD)
                                END IF
                            END IF
                    570 END DO
                    END IF
!
            580 END DO
!            
            END IF
!        
    600 END DO
!--------------------------- 
! COMPUTE OZONE AT MIDLAYERS 
!--------------------------- 
!
!-------------------------------------------------------------------------------------------------- 
! MODIFY PRESSURES SO THAT THE ENTIRE COLUMN OF OZONE (TO 0 MB) IS INCLUDED IN THE MODEL COLUMN 
! EVEN WHEN PT > 0 MB
!-------------------------------------------------------------------------------------------------- 
        DO K=1,LM
            DO I=MYIS,MYIE
                DENOM     = 1. / (PINT(I,LP1) - PINT(I,1))
                FCTRA     =              PINT(I,LP1) * DENOM
                FCTRB     = -PINT(I,1) * PINT(I,LP1) * DENOM
                POZN(I,K) =  PMID(I,K) * FCTRA + FCTRB
            END DO
        END DO
!    
        CALL OZON2D(LM, POZN, XLAT, RSIN1, RCOS1, RCOS2, OZN)
!----------------------------------------------------------  
! NOW THE VARIABLES REQUIRED BY RADFS HAVE BEEN CALCULATED.
!----------------------------------------------------------
!
!------------------------------- 
! CALL THE GFDL RADIATION DRIVER
!------------------------------- 
        CALL RADFS (PSFC , PMID , PINT , QMID , TMID , OZN   , TSKN  , SLMSK , ALBDO, XLAT  ,     &
    &               CAMT , ITYP , KTOP , KBTM , NCLDS, EMIS  , RRCL  , TTCL  , COSZ , TAUDAR,     &
    &               1    , 1    , 0    , ETA  , AETA , ITIMSW, ITIMLW, JD    , RAD1 , HOUR  ,     &
    &               TENDS, TENDL, FLWUP, FSWUP, FSWDN, FSWDNS, FSWUPS, FLWDNS, FLWUPS)
!
        DO 650 I=MYIS,MYIE
            PDSLIJ      = PDSL  (I,J)
            PMOD        = CUPPT (I,J) * 24.0 * 1000.0 / CLSTP
            CFRACL(I,J) = CLDCFR(I,1)
            CFRACM(I,J) = CLDCFR(I,2)
            CFRACH(I,J) = CLDCFR(I,3)
!--------------------------------------------------------------------------------------------------        
! ARRAYS ACFRST AND ACFRCV ACCUMULATE AVERAGE STRATIFORM AND CONVECTIVE CLOUD FRACTIONS, 
! RESPECTIVELY.  
! THIS INFORMATION IS PASSED TO THE POST PROCESSOR VIA COMMON BLOCK ACMCLD.
!--------------------------------------------------------------------------------------------------         
            CFRAVG = AMAX1(CFRACL(I,J),AMAX1(CFRACM(I,J),CFRACH(I,J)))
!---------------------------------------------------------
! PREVENT DOUBLE-COUNTING STATISTICS WHEN RESTARTING MODEL
!---------------------------------------------------------        
            IF (NTSD == 1 .OR. NTSD > NSTART+1) THEN
                IF (CNCLD) THEN
                    IF (PMOD <= PPT(1)) THEN
                        ACFRST(I,J) = ACFRST(I,J) + CFRAVG
                        NCFRST(I,J) = NCFRST(I,J) + 1
                    ELSE
                        ACFRCV(I,J) = ACFRCV(I,J) + CFRAVG
                        NCFRCV(I,J) = NCFRCV(I,J) + 1
                    END IF
                ELSE
                    ACFRST(I,J) = ACFRST(I,J) + CFRAVG
                    NCFRST(I,J) = NCFRST(I,J) + 1
                END IF
            END IF
    650 END DO
!--------------------------------------------------------------------------------------------------
! COLLECT ATMOSPHERIC TEMPERATURE TENDENCIES DUE TO RADIATION.
! ALSO COLLECT THE TOTAL SW AND INCOMING LW RADIATION (W/M**2) AND CONVERT TO FORM NEEDED FOR 
! PREDICTION OF THS IN SURFCE.
!--------------------------------------------------------------------------------------------------
        DO 660 I=MYIS,MYIE
            DO K=1,LM
                LL2 = LVL(I,J) + K
                IF (SHORT) RSWTT(I,J,K) = TENDS(I,LL2)
                IF (LONG)  RLWTT(I,J,K) = TENDL(I,LL2)
                IF (LL2 == LM) GOTO 660
            END DO
    660 END DO
!---------------------------------------------------------
! SUM THE LW INCOMING AND SW RADIATION (W/M**2) FOR RADIN.
!---------------------------------------------------------
        DO 675 I=MYIS,MYIE
            IF (LONG) THEN
                SIGT4(I,J) = STBOL * TMID(I,LM) * TMID(I,LM) * TMID(I,LM) * TMID(I,LM)
            END IF
!--------------------------------------------------------------------------------------------------        
! ACCUMULATE VARIOUS LW AND SW RADIATIVE FLUXES FOR POST PROCESSOR. PASSED VIA COMMON ACMRDL AND
! ACMRDS.
!--------------------------------------------------------------------------------------------------        
            IF (LONG) THEN
                RLWIN (I,J) = FLWDNS(I)
                RLWOUT(I,J) = FLWUPS(I)
                RLWTOA(I,J) = FLWUP (I)
            END IF
!
            IF (SHORT) THEN
                RSWIN (I,J) = FSWDNS(I)
                RSWOUT(I,J) = FSWUPS(I)
                RSWTOA(I,J) = FSWUP (I)
            END IF
!
    675 END DO
!--------------------------------- 
! THIS ROW IS FINISHED. GO TO NEXT
!---------------------------------
700 END DO
!------------------------------------------------
! CALLS TO RADIATION THIS TIME STEP ARE COMPLETE.
!------------------------------------------------
!
!----------------------------------------------- 
! HORIZONTAL SMOOTHING OF TEMPERATURE TENDENCIES
!----------------------------------------------- 
    IF (SHORT) THEN
        DO 800 K=1,LM
            CALL ZERO2(TL )
            CALL ZERO2(FNE)
            CALL ZERO2(FSE)
!        
            IF (KSMUD >= 1) THEN
                DO 750 KS=1,KSMUD
!                
                    DO J=MYJS,MYJE
                        DO I=MYIS,MYIE
                            TL(I,J) = RSWTT(I,J,K) * HTM(I,J,K)
                        END DO
                    END DO
!
                    DO J=MYJS,MYJE
                        DO I=MYIS,MYIE
                            FNE(I,J) = (TL(I+IHE(J),J+1) - TL(I,J)) * HTM(I       ,J  ,K)         &
    &                                *                                HTM(I+IHE(J),J+1,K)
                        END DO
                    END DO
!                
                    DO J=MYJS1,MYJE
                        DO I=MYIS,MYIE
                            FSE(I,J) = (TL(I+IHE(J),J-1) - TL(I,J)) * HTM(I+IHE(J),J-1,K)         &
    &                                *                                HTM(I       ,J  ,K)
                        END DO
                    END DO
!                
                    DO J=MYJS2,MYJE2
                        DO I=MYIS,MYIE
                            TL(I,J) = (FNE(I,J)  - FNE(I+IHW(J),J-1)                              &
    &                               +  FSE(I,J)  - FSE(I+IHW(J),J+1))                             &
    &                               * HBM2(I,J)  * 0.125 + TL(I,J)
                        END DO
                    END DO
!                
                    DO J=MYJS,MYJE
                        DO I=MYIS,MYIE
                            RSWTT(I,J,K) = TL(I,J) 
                        END DO
                    END DO
!                
            750 END DO
            END IF
!        
    800 END DO
    END IF
!
    IF (LONG) THEN
!    
        DO 900 K=1,LM
            CALL ZERO2(TL )
            CALL ZERO2(FNE)
            CALL ZERO2(FSE)
!
            IF (KSMUD >= 1) THEN
                DO 850 KS=1,KSMUD
!                
                    DO J=MYJS,MYJE
                        DO I=MYIS,MYIE
                            TL(I,J) = RLWTT(I,J,K) * HTM(I,J,K)
                        END DO
                    END DO
!                
                    DO J=MYJS,MYJE1
                        DO I=MYIS,MYIE
                            FNE(I,J) = (TL(I+IHE(J),J+1)   -  TL(I       ,J  ))                   &
    &                                * HTM(I       ,J  ,K) * HTM(I+IHE(J),J+1,K)
                        END DO
                    END DO
!                
                    DO J=MYJS1,MYJE
                        DO I=MYIS,MYIE
                            FSE(I,J) = (TL(I+IHE(J),J-1)   -  TL(I,J))                            &
    &                                * HTM(I+IHE(J),J-1,K) * HTM(I,J,K)
                        END DO
                    END DO
!                
                    DO J=MYJS2,MYJE2
                        DO I=MYIS,MYIE
                            TL(I,J) = (FNE(I,J) - FNE(I+IHW(J),J-1)                               &
    &                               +  FSE(I,J) - FSE(I+IHW(J),J+1))                              &
    &                               * HBM2(I,J) * 0.125 + TL(I,J)
                        END DO
                    END DO
!                
                    DO J=MYJS,MYJE
                        DO I=MYIS,MYIE
                            RLWTT(I,J,K) = TL(I,J)
                        END DO
                    END DO
!                
            850 END DO
            END IF
!
    900 END DO
    END IF
!
    RETURN
!
    END SUBROUTINE RADTN
