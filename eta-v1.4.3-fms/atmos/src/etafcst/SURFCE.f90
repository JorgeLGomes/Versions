    SUBROUTINE SURFCE(APE, ZINT, CKLQ)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SURFCE
!>
!> SUBPROGRAM: SURFCE - CALCULATE SURFACE CONDITIONS
!> PROGRAMMER: JANJIC   
!> ORG: W/NP22
!> DATE: 95-03-23
!>
!> ABSTRACT:
!> THIS ROUTINE IS THE DRIVER FOR COMPUTATION OF GROUND CONDITIONS.  
!> FOR GCIP, ACCUMULATOR AND OTHER INSTANTANEOUS HOLDING ARRAYS ARE INCLUDED.
!>
!> PROGRAM HISTORY LOG:
!> 95-03-23  JANJIC - ORIGINATOR
!> 95-03-28  BLACK  - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 96-03-29  BLACK  - REMOVED SCRCH COMMON
!> 18-01-15  LUCCI  - MODERNIZATION OF THE CODE, INCLUDING:
!>                    * F77 TO F90/F95
!>                    * INDENTATION & UNIFORMIZATION CODE
!>                    * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                    * DOCUMENTATION WITH DOXYGEN
!>                    * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> APE  - EXNER FUNCTION
!> ZINT - INTERFACE HEIGHTS
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> CKLQ - MASK VALUE
!>
!> USE MODULES: ACMRDL
!>              ACMRDS
!>              ACMPRE
!>              ACMSFC
!>              CTLBLK
!>              DYNAM
!>              F77KINDS
!>              GLB_TABLE
!>              LOOPS
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              OPTIONS
!>              PARMETA
!>              PARMSOIL
!>              PARM_TBL
!>              PHYS
!>              PVRBLS
!>              SOIL
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>         
!> DRIVER     : TURBL
!>
!> CALLS      : SFLX
!>              ZERO2
!>--------------------------------------------------------------------------------------------------
    USE ACMRDL
    USE ACMRDS
    USE ACMPRE
    USE ACMSFC
    USE CTLBLK
    USE DYNAM    , ONLY: AETA, PT
    USE F77KINDS
    USE GLB_TABLE
    USE LOOPS
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE OPTIONS
    USE PARMETA
    USE PARMSOIL
    USE PARM_TBL
    USE PHYS
    USE PVRBLS
    USE SOIL
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND), PARAMETER :: EPSWET =     .001
    REAL   (KIND=R4KIND), PARAMETER :: PQ0    =  379.90516
    REAL   (KIND=R4KIND), PARAMETER :: SEAFC  =     .98
    REAL   (KIND=R4KIND), PARAMETER :: TRESH  =     .95 
    REAL   (KIND=R4KIND), PARAMETER :: A2     =   17.2693882
    REAL   (KIND=R4KIND), PARAMETER :: A3     =  273.16 
    REAL   (KIND=R4KIND), PARAMETER :: A4     =   35.86
    REAL   (KIND=R4KIND), PARAMETER :: T0     =  273.16
    REAL   (KIND=R4KIND), PARAMETER :: T1     =  274.16
    REAL   (KIND=R4KIND), PARAMETER :: CAPA   =    0.28589641
    REAL   (KIND=R4KIND), PARAMETER :: CP     = 1004.6 
    REAL   (KIND=R4KIND), PARAMETER :: STBOL  =    5.67E-8
    REAL   (KIND=R4KIND), PARAMETER :: R      =  287.04
    REAL   (KIND=R4KIND), PARAMETER :: ROW    =    1.E3
    REAL   (KIND=R4KIND), PARAMETER :: G      =    9.81E0
    REAL   (KIND=R4KIND), PARAMETER :: ELWV   =    2.50E6 
    REAL   (KIND=R4KIND), PARAMETER :: ELIV   =    2.834E6
    REAL   (KIND=R4KIND), PARAMETER :: ELIW   =     .334E6
    REAL   (KIND=R4KIND), PARAMETER :: A23M4  = A2   * (A3 - A4)
    REAL   (KIND=R4KIND), PARAMETER :: PQ0SEA = PQ0  * SEAFC
    REAL   (KIND=R4KIND), PARAMETER :: PQ0C   = PQ0  * TRESH
    REAL   (KIND=R4KIND), PARAMETER :: GROW   = G    * ROW  
    REAL   (KIND=R4KIND), PARAMETER :: RLIVWV = ELIV / ELWV
    REAL   (KIND=R4KIND), PARAMETER :: ROWLIW = ROW  * ELIW 
    REAL   (KIND=R4KIND), PARAMETER :: ROWLIV = ROW  * ELIV
!
#include "sp.h"
!----------------------------------------------------- 
! SET LOCAL PARAMETERS DEPENDENT ON GLOBAL PARAMETERS.
!----------------------------------------------------- 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & ZLM     , PS      , ETALM   , PLM     , APELM   , RDSIN   , TLM     , THLM    , QLM     ,   &
    & QLMS    , FFS     , QFC1    , APES    , DQSDT
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(INOUT)       ::&
    & CKLQ
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)         , INTENT(IN)          ::&
    & APE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LP1)        , INTENT(IN)          ::&
    & ZINT
!--------------------------------- 
! EK 18 JAN 2000 - ADD SH2OK ARRAY
!--------------------------------- 
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                                      ::&
    & SMCK    , STCK    , SH2OK
!----------------------------------
! EK 10 FEB 2000 - ADD DECLARATIONS
! DECLARATIONS
!----------------------------------
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & LFIRST  , LFIRSTA
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ICE     , ISLTPK  , IVGTPK  , ISPTPK
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ALB     , ALB2D   , DQSDTK  , DTK     , CHK     , CMCK    , ELFLX   , GFLX    , HFLX    ,   &
    & LWDN    , PLFLX   , PRCP    , PTU     , Q1K     , Q2K     , Q2SAT   , RNOF1K  , RNOF2K  ,   &
    & SFCPRS  , SFCSPD  , SFCTH2  , SFCTMP  , SMELTK  , SNDENS  , SNOALB  , z       , VGFRCK  ,   &
    & SNODPK  , SNOWH   , SOILQM  , SOILQW  , SOLDN   , T1K     , TBOT    
!
    CHARACTER(LEN=80)                                                                             &
    & WORD
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , NS      , LLMH    , INDEX   
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & FDTLIW  , FDTLIV  , FDTW    , FFSK    , FACTRS  , TLMH    , FACTRL  , CHKFF   , CMK     ,   &
    & DUM1    , DUM2    , DUM3    , DUM4    , SCHECK  , TSFC    , TSFC2
!
    DATA    LFIRST  / .TRUE. /
    DATA    LFIRSTA / .TRUE. /
!------------------ 
! START SURFCE HERE
!------------------ 
!
!------------------------------- 
! INITIALIZE SOME WORKING ARRAYS
!------------------------------- 
    CALL ZERO2(QLM)
    CALL ZERO2(QLM)
    CALL ZERO2(QLMS)
!-------------------------------------------
! SET CONSTANTS CALCULATED HERE FOR CLARITY.
!-------------------------------------------
    FDTLIW = DTQ2 / ROWLIW
    FDTLIV = DTQ2 / ROWLIV
      FDTW = DTQ2 / 2.5E9
!------------------------------------------------------
! SET NOAH LSM CONSTANTS AND TIME INDEPENDENT VARIABLES
!  INITIALIZE NOAH LSM HISTORICAL VARIABLES
!------------------------------------------------------
    IF (NTSD < NPHS) THEN
!
!$omp parallel do private (I      , J)
!
        DO 50 J=MYJS,MYJE
            DO 50 I=MYIS,MYIE
                    PS(I,J) = PD(I,J) + PT
                  APES(I,J) = (1.E5 / PS(I,J)) ** CAPA
                PCTSNO(I,J) = -999.0
!------------------------------------------------------
! SET DEFAULT VALUES FOR SEA-ICE OR OCEAN STATES
! OPEN OCEAN, SM=1
! SEA-ICE, SM=0, SICE=1
! LAND, SM=0, SICE=0
! SET ZERO-VALUE FOR SOME OUTPUT DIAGNOSTIC ARRAYS
!------------------------------------------------------
                IF (SM(I,J) < 0.5) THEN
                    IF (SICE(I,J) > 0.5) THEN
!------------- 
! SEA-ICE CASE
!------------- 
                        SMSTAV(I,J) = 1.0
                        SMSTOT(I,J) = 1.0
                        SSROFF(I,J) = 0.0
                        BGROFF(I,J) = 0.0
                           CMC(I,J) = 0.0
!
                        DO NS=1,NSOIL
                             SMC(I,J,NS) = 1.0
                            SH2O(I,J,NS) = 1.0
                        END DO
                    END IF
                ELSE
!----------- 
! WATER CASE
!----------- 
                    SMSTAV(I,J) =   1.0
                    SMSTOT(I,J) =   1.0
                    SSROFF(I,J) =   0.0
                    BGROFF(I,J) =   0.0
                    SOILTB(I,J) = 273.16
                    GRNFLX(I,J) =   0.0
                    SUBSHX(I,J) =   0.0
                    ACSNOW(I,J) =   0.0
                    ACSNOM(I,J) =   0.0
                    SNOPCX(I,J) =   0.0
                       CMC(I,J) =   0.0
                       SNO(I,J) =   0.0
!---------------------------------------------- 
!ADD SI (SNOW DEPTH), SNO=SNOW WATER EQUIVALENT
!---------------------------------------------- 
                    SI(I,J) = 0.0
!
                    DO NS=1,NSOIL
                         SMC(I,J,NS) =   1.0
                        SH2O(I,J,NS) =   1.0
                         STC(I,J,NS) = 273.16
                    END DO
                END IF
!            
     50 END DO
!
    END IF
!---------------------------------- 
! SET LOWEST MODEL LAYER VARIABLES.
!---------------------------------- 
!
!$omp parallel do private (I      , J       , LLMH)
!
    DO 100 J=MYJS2,MYJE2
        DO 100 I=MYIS,MYIE
            LLMH = LMH(I,J)
!
            ETALM(I,J) =  AETA(    LLMH)
            APELM(I,J) =   APE(I,J,LLMH)
              TLM(I,J) =     T(I,J,LLMH)
              QLM(I,J) =     Q(I,J,LLMH)
              ZLM(I,J) = (ZINT(I,J,LLMH) - ZINT(I,J,LLMH+1)) * 0.50
100 END DO
!
!$omp parallel do private (I      , J)
!
    DO 110 J=MYJS2,MYJE2
        DO 110 I=MYIS,MYIE
               PS(I,J) = PD(I,J) + PT
             APES(I,J) = (1.E5 / PS(I,J)) ** CAPA
              PLM(I,J) = ETALM(I,J) * PD(I,J) * RES(I,J) + PT
             QLMS(I,J) = ((1. - SM(I,J)) * PQ0 + SM(I,J) * PQ0SEA)                                &
    &                  /  PLM(I,J) * EXP(A2 * (TLM(I,J) - A3) / (TLM(I,J) - A4))
!
            DQSDT(I,J) = QLMS(I,J) * A23M4 / (TLM(I,J) - A4) ** 2
                  FFSK = AKHS(I,J) * PLM(I,J) * HBM2(I,J) / ((QLM(I,J) * .608 + 1.) * TLM(I,J) * R)
             QFC1(I,J) = APES(I,J) * FFSK * ELWV
              FFS(I,J) = FFSK * CP
110 END DO
!
!$omp parallel do private (I      , J       , FACTRS  , FACTRL  , TLMH)
!
    DO 120 J=MYJS2,MYJE2
        DO 120 I=MYIS,MYIE
!---------------------------------------------------------- 
! COMPUTE RADIN AND RDSIN FOR THIS TIMESTEP
! CZEN IS IN PHYS COMMON AND IS CURRENT FROM CALL TO RDTEMP
!---------------------------------------------------------- 
            IF (CZMEAN(I,J) > 0.) THEN
                FACTRS = CZEN(I,J) / CZMEAN(I,J)
            ELSE
                FACTRS = 0.
            END IF
!        
            IF (SIGT4(I,J) > 0.) THEN
                TLMH   = TLM(I,J)
                FACTRL = STBOL * TLMH * TLMH * TLMH * TLMH / SIGT4(I,J)
            ELSE
                FACTRL = 0.
            END IF
!------------------------------------------------------------         
! EK 10 FEB 2000 - RADIN NO LONGER NEEDED IN SFLX (VIA FK)
! NOW USE RADIN ARRAY FOR INCOMING LONGWAVE 
! PERHAPS CHANGE THE NAME LATER TO E.G. RDLIN FOR CONSISTENCY
!------------------------------------------------------------ 
            RADIN(I,J) = RLWIN(I,J) * FACTRL * HBM2(I,J)
            RDSIN(I,J) = RSWIN(I,J) * FACTRS * HBM2(I,J)
!---------------------------------- 
! DIAGNOSTIC RADIATION ACCUMULATION
!---------------------------------- 
            ASWIN (I,J) = ASWIN (I,J) + RSWIN (I,J) * HBM2(I,J) * FACTRS
            ASWOUT(I,J) = ASWOUT(I,J) - RSWOUT(I,J) * HBM2(I,J) * FACTRS
            ASWTOA(I,J) = ASWTOA(I,J) + RSWTOA(I,J) * HBM2(I,J) * FACTRS
            ALWIN (I,J) = ALWIN (I,J) + RLWIN (I,J) * HBM2(I,J) * FACTRL
            ALWOUT(I,J) = ALWOUT(I,J) - RADOT (I,J) * HBM2(I,J)
            ALWTOA(I,J) = ALWTOA(I,J) + RLWTOA(I,J) * HBM2(I,J)
!----------------------------------------------- 
! CHECK FOR SATURATION AT THE LOWEST MODEL LEVEL
!----------------------------------------------- 
            IF ((QLM(I,J) >= QLMS(I,J)*TRESH) .AND. (QLM(I,J) < QZ0(I,J))) THEN
                CKLQ(I,J) = 0.
            ELSE
                CKLQ(I,J) = HBM2(I,J)
            END IF
120 END DO
!---------------------------------- 
! THS, THLM, CHEATING WET FOR PROFS
!---------------------------------- 
!
!$omp parallel do private (I      , J)    
!
    DO 130 J=MYJS2,MYJE2
        DO 130 I=MYIS,MYIE
!CHOU            THLM(I,J) = TLM (I,J) * APELM(I,J)
!CHOU using Tv
            THLM(I,J) = (TLM (I,J) * (1. + (0.608 * QLM(I,J)))) * APELM(I,J)
!CHOU 
            QFC1(I,J) = QFC1(I,J) * CKLQ (I,J)
130 END DO
!--------------------------------------------------
! EK 10 FEB 2000 - UPDATE THESE PRIVATE STATEMENTS
! ADD NEW VARIABLES: ALB2D, SNOALB, ALB, ISPTPK
!--------------------------------------------------
! EK 10 FEB 2000 - PRIVATE STATEMENTS
! ADD NEW VARIABLES
!--------------------------------------------------
!
!-----------------------------------------------------------------
! EK 18 JAN 2000 - TEMPORARILY SET ISPTPK=1 (2-D FIXED FIELD: X,Y)
! COMES FROM ISLSCP DATA SET 2-D FIXED FIELD
!-----------------------------------------------------------------
    ISPTPK = 1
! ---------------------------------------------------- 
! BEGIN MAIN 'WORKHORSE' LOOP OVER ENTIRE MODEL DOMAIN
! ---------------------------------------------------- 
    DO 160 J=MYJS2,MYJE2
        DO 155 I=MYIS,MYIE
            IF (HBM2(I,J) < 0.5) GOTO 155
            IF (SM(I,J) > 0.5) THEN
                THS(I,J) = SST(I,J) * APES(I,J)
                 QS(I,J) = HBM2(I,J) * PQ0SEA / PS(I,J)                                           &
    &                    * EXP(A2 * (THS(I,J) - A3 * APES(I,J))                                   &
    &                    /          (THS(I,J) - A4 * APES(I,J)))
            END IF
! --------------- 
! LAND OR SEA-ICE
! --------------- 
!
! -------------------------------------------------- 
! LOADING AND UNLOADING NOAH LSM LAND SOIL VARIABLES
! -------------------------------------------------- 
            IF (SM(I,J) < 0.5) THEN
                ICE = NINT(SICE(I,J))
!            
                DTK = DTQ2
                Z = ZLM(I,J)
! -----------------------------------------------------------------
! EK 18 JAN 2000 - ADD LONGWAVE RADIATION CALC NEEDED FOR CALL SFLX
! -----------------------------------------------------------------
                  LWDN = RADIN(I,J)
                 SOLDN = RDSIN(I,J)
                SFCPRS =   PLM(I,J)
                  PRCP =  PREC(I,J) * ROW / DTQ2
                   Q2K =   QLM(I,J)
                 Q2SAT =  QLMS(I,J)
!--------------------------------------------------------------------------------------------------
! Q2K MAY SLIGHTLY EXCEED Q2SAT IN SOME CASES DUE TO ATMOSPHERIC PHYSICS PARAMETERIZATIONS 
! PREVIOUSLY CALLED
!--------------------------------------------------------------------------------------------------
                IF (Q2K > Q2SAT) Q2K = Q2SAT
!
                DQSDTK =  DQSDT(I,J)
                  TBOT =     TG(I,J)
                   CHK =   AKHS(I,J)
                 CHKFF =    FFS(I,J)
                IVGTPK = IVGTYP(I,J)
                ISLTPK = ISLTYP(I,J)
!------------------------------------------------------- 
! MEB  PREVENT ROUTINES IN SFLX FROM GOING OUT OF BOUNDS
!------------------------------------------------------- 
                IF (IVGTPK == 0) IVGTPK = 13
                IF (ISLTPK == 0) ISLTPK =  9
!------------------------------------------------------- 
! MEB  PREVENT ROUTINES IN SFLX FROM GOING OUT OF BOUNDS
!------------------------------------------------------- 
                VGFRCK = VEGFRC(I,J)
                   Q1K =     QS(I,J)
                SFCTMP =   THLM(I,J) / APELM(I,J)
                SFCTH2 =   THLM(I,J) /  APES(I,J)
                   T1K =    THS(I,J) /  APES(I,J)
                  CMCK =    CMC(I,J)
                SNODPK =    SNO(I,J)
!------------------------------------------------------------------------------------ 
! USE 2-D PROGNOSTIC FIELD OF SNOWDEPTH, SI(X,Y) FOR SNOWH (LOCAL SNOWDEPTH VARIABLE)
! EK 17 JAN 2001
!------------------------------------------------------------------------------------ 
                SNOWH = SI(I,J)
!            
                DO 140 NS=1,NSOIL
                    SMCK(NS) = SMC(I,J,NS)
!----------------------------------------------------------------------------------------------- 
! USE 3-D PROGNOSTIC FIELD OF LIQUID SOIL MOISTURE, SH2O(X,Y,4) FOR SH2OK(NS) (LOCAL LIQUID SOIL 
! MOISTURE VARIABLE)
! EK 11 JAN 2001
!----------------------------------------------------------------------------------------------- 
                    SH2OK(NS) = SH2O(I,J,NS)
!                
                     STCK(NS) =  STC(I,J,NS)
            140 END DO
!-----------------------------------------------------------------------------------------------           
! EK 18 JAN 2000 - TEMPORARILY SET ISPTPK=1 (2-D FIXED FIELD: X,Y) COMES FROM ISLSCP DATA SET 
! 2-D FIXED FIELD
! ISPTPK=1
!  SNOALB (FIXED VALUE, MAX SNOW ALBEDO) FROM MXSNAL VIA 2-D FIXED FIELD FROM DAVID ROBINSON
!-----------------------------------------------------------------------------------------------
                SNOALB = MXSNAL(I,J)
!---------------------------------------------------------------------------------- 
! ALB (FIXED VALUE, SNOW-FREE ALBEDO) FROM ALBASE VIA 2-D FIXED FIELD FROM MATTHEWS
!---------------------------------------------------------------------------------- 
                ALB = ALBASE(I,J)
!-----------------------------------------------------------------------------------------------
! SET DYNAMIC ALBEDO FROM THE DYNAMIC ALBEDO 2-D ARRAY, WHICH IS UPDATED ONLY FOR THE LAND IN 
! SFLX, NOT FOR SEA-ICE, SO WE MUST 'PASS THROUGH' ALB2D=0.60 FOR SEA-ICE.
! TURN THIS OFF, AND INSTEAD, DO IT WITHIN SFLX
! ALB2D=ALBEDO(I,J)
!-----------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------
! INITIAL RANGE CHECK FOR VARIABLES/PARAMETERS FOR ENTIRE I,J DOMAIN FOR FIRST TIMESTEP 
! (WHEN LFIRST=TRUE).  SET LFIRST=FALSE AFTER END OF THE 155/160 LOOP.
!-----------------------------------------------------------------------------------------------
                IF (LFIRST) THEN
!----------------------- 
! LAND OR SEA-ICE CHECKS
! LAND CHECKS FIRST
!-----------------------
                    IF (ICE < 0.5) THEN
                        DO INDEX=1,NSOIL
!--------------------------------
! TOTAL SOIL MOISTURE (SMC) CHECK
!--------------------------------
                            IF ((SMCK(INDEX) < 0.02) .OR. (SMCK(INDEX) > 0.468)) THEN
                                WRITE(6,*) 'SMC:INDEX,I,J,MYPE,STCK,SMC,SH2O=',                   &
    &                                       INDEX, I, J, MYPE,                                    &
    &                                       STCK(INDEX), SMCK(INDEX), SH2OK(INDEX)
                            END IF
!-------------------------------------------------------------------- 
! LIQUID SOIL MOISTURE<=TOTAL SOIL MOISTURE (SH2O<=SMC) MAXIMUM CHECK
!-------------------------------------------------------------------- 
                            IF ((SH2OK(INDEX) < 0.02) .OR. (SH2OK(INDEX) > SMCK(INDEX))) THEN
                                WRITE(6,*) 'SH2OL1:INDEX,I,J,MYPE,ICE,STCK,SMC,SH2O=',            &
    &                                       INDEX, I, J, MYPE,                                    &
    &                                       STCK(INDEX), SMCK(INDEX), SH2OK(INDEX)
                            END IF
                        END DO
!------------------------------------- 
! SOIL COLUMN BOTTOM TEMP (TBOT) CHECK
!------------------------------------- 
                        IF ((TBOT < 200.00) .OR. (TBOT > 323.15)) THEN
                            WRITE(6,*) 'TBOTL:INDEX,I,J,MYPE,TBOT=', INDEX, I, J, MYPE, TBOT
                        END IF
!--------------------------------- 
! CANOPY WATER CONTENT (CMC) CHECK
!--------------------------------- 
                        IF ((CMCK < 0.) .OR. (CMCK > 0.5E-3)) THEN
                            WRITE(6,*) 'CMC:INDEX,I,J,MYPE,CMC=', INDEX, I, J, MYPE, CMCK
                        END IF
!---------------- 
! END LAND CHECKS
!---------------- 
                    ELSE
!-------------------- 
! SEA-ICE CHECKS NEXT
!-------------------- 
!
!--------------------------------- 
! SEA-ICE BOTTOM TEMP (TBOT) CHECK
!--------------------------------- 
                        IF ((TBOT < 271.159) .OR. (TBOT > 271.161)) THEN
                            WRITE(6,*) 'TBOTI:INDEX,I,J,MYPE,TBOT=', INDEX, I, J, MYPE, TBOT
                        END IF
!------------------------------------------ 
! SEA-ICE TEMP WITH DEPTH (STC) RANGE CHECK
!------------------------------------------ 
                        DO INDEX=1,4
                            IF ((STCK(INDEX) < 200.00) .OR. (STCK(INDEX) > 274.15)) THEN
                                WRITE(6,*) 'STCI:INDEX,I,J,MYPE,STC=',                            &
    &                                       INDEX, I, J, MYPE,                                    &
    &                                       STCK(INDEX)
                            END IF
                        END DO
!--------------------------------------------- 
! CHECK TO SEE THAT WHEN SEA-ICE, SH2O=SMC=1.0
!--------------------------------------------- 
                        DO INDEX=1,4
                            IF  ((SMCK(INDEX) /=          1.0) .OR.                               &
    &                            (SMCK(INDEX) /= SH2OK(INDEX))) THEN
                                WRITE(6,*) 'SMCI:INDEX,I,J,MYPE,ICE,STCK,SMC,SH2O=',              &
    &                            INDEX, I, J, MYPE, ICE,                                          &
    &                            STCK(INDEX), SMCK(INDEX), SH2OK(INDEX)
                            END IF
                        END DO
!
                    END IF
!---------------------------------------- 
! BOTH LAND AND SEA-ICE CHECKS
! SNOW WATER EQUIVALENT, SNOW DEPTH CHECK
!---------------------------------------- 
                    IF (((SNODPK >  0.) .AND. (SNOWH <= 0.)) .OR.                                 &
    &                   ((SNODPK <= 0.) .AND. (SNOWH >  0.)) .OR.                                 &
    &                    (SNODPK > SNOWH)) THEN
                        WRITE(6,*) 'SNOW:I,J,MYPE,SNODPK,SNOWH=', I, J, MYPE, SNODPK, SNOWH
                    END IF
!------------------- 
! SNOW DENSITY CHECK
!------------------- 
                    IF (SNODPK > 0.) THEN
                        SNDENS = SNODPK / SNOWH
                        IF (SNDENS < 0.05) THEN
                            WRITE(6,*) 'SNDENS<5%:I,J,MYPE,SNODPK,SNOWH,SNDENS=',                 &
    &                                   I, J, MYPE, SNODPK, SNOWH, SNDENS
                        END IF
!
                        IF (SNDENS > 0.40) THEN
                            WRITE(6,*)'SNDENS>40%:I,J,MYPE,SNODPK,SNOWH,SNDENS=',                 &
    &                        I, J, MYPE, SNODPK, SNOWH, SNDENS
                        END IF
                    END IF
!-------------------------- 
! SFC/SKIN TEMP (T1K) CHECK
!-------------------------- 
                    IF ((T1K < 200.00) .OR. (T1K > 323.15)) THEN
                        WRITE(6,*) 'T1:INDEX,I,J,MYPE,T1=', INDEX, I, J, MYPE, T1K
                    END IF
!-------------------------------------------------------------------- 
! CHECK TO SEE THAT 223.15K (-50C) =< SFCTMP,SFCTH2 <= 323.15K (+50C)
! SFCTMP = LOWEST MODEL LEVEL TEMP
! SFCTH2 = LOWEST MODEL LEVEL POT TEMP
!-------------------------------------------------------------------- 
                    IF (((SFCTMP < 223.15) .OR. (SFCTMP > 323.15)) .OR.                           &
    &                   ((SFCTH2 < 223.15) .OR. (SFCTH2 > 323.15))) THEN
                        WRITE(6,*) 'SFCTMP:I,J,MYPE,SFCTMP,SFCTH2=',                              &
    &                               I, J, MYPE, SFCTMP, SFCTH2 
                    END IF
!---------------------------------------------- 
! CHECK TO SEE THAT 0W/M2 =< LWDN  <=  500W/M2
! CHECK TO SEE THAT 0W/M2 =< SOLDN <= 1200W/M2
! LWDN  = DOWNWARD LONGWAVE RADIATION
! SOLDN = DOWNWARD SOLAR RADIATION
!--------------------------------------------- 
                    IF (((LWDN  < 0.) .OR. (LWDN  >  500.)) .OR.                                  &
    &                   ((SOLDN < 0.) .OR. (SOLDN > 1200.))) THEN
                        WRITE(6,*)'LWSOLDN:I,J,MYPE,LWDN,SOLDN=', I, J, MYPE, LWDN, SOLDN
                    END IF
!---------------------------------------------- 
! CHECK TO SEE THAT 0G/KG < Q2K,Q2SAT <= 40G/KG
! CHECK TO SEE THAT Q2K <= Q2SAT
! Q2K   = LOWEST MODEL LEVEL SPEC HUM
! Q2SAT = LOWEST MODEL LEVEL SAT SPEC HUM
!---------------------------------------------- 
                    IF (((Q2K   <= 0.) .OR. (Q2K   > 0.04)) .OR.                                  &
    &                   ((Q2SAT <= 0.) .OR. (Q2SAT > 0.04)) .OR. (Q2K > Q2SAT)) THEN
                        WRITE(6,*) 'Q2:I,J,MYPE,Q2K,Q2SAT=', I, J, MYPE, Q2K, Q2SAT
                    END IF
!-------------------------------------------- 
! CHECK TO SEE THAT 600MB =< SFCPRS <= 1050MB
! SFCPRS = SURFACE PRESSURE (PA)
!-------------------------------------------- 
                    IF ((SFCPRS <= 60000.) .OR. (SFCPRS > 105000.)) THEN
                        WRITE(6,*) 'SFCPRS:I,J,MYPE,SFCPRS=', I, J, MYPE, SFCPRS
                    END IF
!-----------------------------------------------------------------------
! CHECK TO SEE THAT 0 =< PRCP <= 0.04 KG M-2 S-1 (=0.04MM/S = 144MM/HR)
! PRCP = PRECIP RATE (KG M-2 S-1)
!-----------------------------------------------------------------------
                    IF ((PRCP < 0.) .OR. (PRCP > 0.04)) THEN
                        WRITE(6,*) 'PRCP:I,J,MYPE,PRCP=', I, J, MYPE, PRCP
                    END IF
!-----------------------------------------------------------------------
! CHECK TO SEE THAT 0M/S < CHK <= 0.1M/S
! CHK = SFC HEAT EXCHANGE COEFF (M/S)
!-----------------------------------------------------------------------
                    IF ((CHK <= 0.) .OR. (CHK > 0.1)) THEN
                        WRITE(6,*) 'CH:I,J,MYPE,CHK=', I, J, MYPE, CHK
                    END IF
!
                END IF
!-------------------------------- 
! EK 18 JAN 2000 - NEW CALL SFLX
!--------------------------------
!                IF (NTSD.EQ.4361 .AND. I.LE.3 .AND. J.LE.3 .AND. MYPE.EQ.600) THEN
!                IF (NTSD.EQ.8681 .AND. I.LE.3 .AND. J.LE.3 .AND. MYPE.EQ.600) THEN
!		WRITE(2000,*) "Before SFLX", ICE   , DTK   , Z     , NSOIL , SLDPTH, LWDN  , SOLDN , SFCPRS, PRCP  , &
!    &                     SFCTMP, SFCTH2, Q2K   ,         Q2SAT , DQSDTK, IVGTPK, ISLTPK, ISPTPK, &
!    &                     VGFRCK,         TBOT  , ALB   , SNOALB, CMCK  , T1K   , STCK  , SMCK  , &
!    &                     SH2OK , SNOWH , SNODPK, ALB2D , CHK   , CMK   , PLFLX , ELFLX , HFLX  , &
!    &                     GFLX  , RNOF1K, RNOF2K, Q1K   , SMELTK, SOILQW, SOILQM, DUM1  , DUM2  , &
!    &                     DUM3  , DUM4  , I     , J
!                ENDIF
!		
                CALL SFLX(ICE   , DTK   , Z     , NSOIL , SLDPTH, LWDN  , SOLDN , SFCPRS, PRCP  , &
    &                     SFCTMP, SFCTH2, Q2K   , SFCSPD, Q2SAT , DQSDTK, IVGTPK, ISLTPK, ISPTPK, &
    &                     VGFRCK, PTU   , TBOT  , ALB   , SNOALB, CMCK  , T1K   , STCK  , SMCK  , &
    &                     SH2OK , SNOWH , SNODPK, ALB2D , CHK   , CMK   , PLFLX , ELFLX , HFLX  , &
    &                     GFLX  , RNOF1K, RNOF2K, Q1K   , SMELTK, SOILQW, SOILQM, DUM1  , DUM2  , &
    &                     DUM3  , DUM4)
!
!                IF (NTSD.EQ.4361 .AND. I.LE.3 .AND. J.LE.3 .AND. MYPE.EQ.600) THEN
!                IF (NTSD.EQ.8681 .AND. I.LE.3 .AND. J.LE.3 .AND. MYPE.EQ.600) THEN
!		WRITE(2000,*) "After SFLX", ICE   , DTK   , Z     , NSOIL , SLDPTH, LWDN  , SOLDN , SFCPRS, PRCP  , &
!    &                     SFCTMP, SFCTH2, Q2K   ,         Q2SAT , DQSDTK, IVGTPK, ISLTPK, ISPTPK, &
!    &                     VGFRCK,         TBOT  , ALB   , SNOALB, CMCK  , T1K   , STCK  , SMCK  , &
!    &                     SH2OK , SNOWH , SNODPK, ALB2D , CHK   , CMK   , PLFLX , ELFLX , HFLX  , &
!    &                     GFLX  , RNOF1K, RNOF2K, Q1K   , SMELTK, SOILQW, SOILQM, DUM1  , DUM2  , &
!    &                     DUM3  , DUM4  , I     , J
!                ENDIF
!
                IF (LFIRSTA) THEN
!----------------------- 
! LAND OR SEA-ICE CHECKS
! LAND CHECKS FIRST
!----------------------- 
                    IF (ICE < 0.5) THEN
!------------------------------------------------------------------------------------
! ALBEDO CHECKS
! ALB = ALBASE(I,J)  = SNOW FREE ALBEDO
! MIN = 0.11 (MATTHEWS DATA BASE)
! MAX = 0.75 (MATTHEWS DATA BASE)
! ALBEDO(I,J) = DYNAMIC ALBEDO (=ALBASE WHEN SNODPK=0), (=ALB2D ON RETURN FROM SFLX)
! SNOALB = MAXSNAL(I,J) = MAXIMUM SNOW ALBEDO
! MIN = 0.21 (ROBINSON DATA BASE)
! MAX = 0.80 (ROBINSON DATA BASE)
!------------------------------------------------------------------------------------
                        IF ((ALB         > SNOALB     )  .OR.                                     &
    &                       (ALB         > ALBEDO(I,J))  .OR.                                     &
    &                       (ALBEDO(I,J) > SNOALB     )) THEN
                             WRITE(6,*) 'ALBL1:I,J,MYPE,ICE,SNODPK,ALB,ALB2D,SNOALB=',            &
    &                                    I, J, MYPE, ICE, SNODPK, ALB, ALBEDO(I,J), SNOALB
                        END IF
!
                        IF ((ALB    < 0.10)  .OR.                                                 &
    &                       (ALB    > 0.76)  .OR.                                                 &
    &                       (SNOALB < 0.20)  .OR.                                                 &
    &                       (SNOALB > 0.81)) THEN
                             WRITE(6,*) 'ALBL2:I,J,MYPE,ICE,SNODPK,ALB,ALB2D,SNOALB=',            &
    &                                    I, J, MYPE, ICE, SNODPK, ALB, ALBEDO(I,J), SNOALB
                        END IF
!---------------------------------------------------------- 
! VEG,SOIL,SLOPE TYPE, VEG FRACTION, NO. SOIL LAYERS CHECKS
!---------------------------------------------------------- 
                        IF (((IVGTPK <  1) .OR.  (IVGTPK > 16)) .OR.                              &
    &                       ((ISLTPK <  1) .OR.  (ISLTPK > 15)) .OR.                              &
    &                        (ISPTPK /= 1) .OR.                                                   &
    &                       ((VGFRCK < 0.) .OR.  (VGFRCK > 1.)) .OR. (NSOIL  /= 8)) THEN
                            WRITE(6,*) 'LANDSFC:I,J,MYPE,IVGTPK,ISLTPK,ISPTPK,VGFRCK=',           &
    &                                   I, J, MYPE, IVGTPK, ISLTPK, ISPTPK, VGFRCK
                        END IF
!
                        DO INDEX=1,NSOIL
!---------------------------- 
! SOIL TEMP (STC) RANGE CHECK
!----------------------------
                            IF ((STCK(INDEX) < 223.15) .OR. (STCK(INDEX) > 323.15)) THEN
                                WRITE(6,*) 'STCL:INDEX,I,J,MYPE,STC=',                            &
    &                                       INDEX, I, J, MYPE,                                    &
    &                                       STCK(INDEX)
                            END IF
!--------------------------------
! TOTAL SOIL MOISTURE (SMC) CHECK
!--------------------------------
                            IF ((SMCK(INDEX) < 0.02) .OR. (SMCK(INDEX) > 0.468)) THEN
                                WRITE(6,*) 'SMC:INDEX,I,J,MYPE,STCK,SMC,SH2O=',                   &
    &                                       INDEX, I, J, MYPE,                                    &
    &                                       STCK(INDEX), SMCK(INDEX), SH2OK(INDEX)
                            END IF
!--------------------------------------------------------------------
! LIQUID SOIL MOISTURE<=TOTAL SOIL MOISTURE (SH2O<=SMC) MAXIMUM CHECK
!--------------------------------------------------------------------
                            IF ((SH2OK(INDEX) < 0.02) .OR.                                        &
    &                           (SH2OK(INDEX) > SMCK(INDEX))) THEN
                                WRITE(6,*) 'SH2OL1:INDEX,I,J,MYPE,ICE,STCK,SMC,SH2O=',            &
    &                                       INDEX, I, J, MYPE, ICE,                               &
    &                                       STCK(INDEX), SMCK(INDEX), SH2OK(INDEX)                        
                            END IF
!
                        END DO
!------------------------------------- 
! SOIL COLUMN BOTTOM TEMP (TBOT) CHECK
!------------------------------------- 
                        IF ((TBOT < 200.00) .OR. (TBOT > 323.15)) THEN
                            WRITE(6,*) 'TBOTL:INDEX,I,J,MYPE,TBOT=',                              &
    &                                   INDEX, I, J, MYPE, TBOT
                        END IF
!--------------------------------- 
! CANOPY WATER CONTENT (CMC) CHECK
!--------------------------------- 
                        IF ((CMCK < 0.) .OR. (CMCK > 0.5E-3)) THEN
                            WRITE(6,*) 'CMC:INDEX,I,J,MYPE,CMC=',                                 &
    &                                   INDEX, I, J, MYPE, CMCK
                        END IF
!
                    ELSE
!--------------------
! SEA-ICE CHECKS NEXT
!-------------------- 
!
!--------------------------------- 
! SEA-ICE BOTTOM TEMP (TBOT) CHECK
!--------------------------------- 
                        IF ((TBOT < 271.159) .OR. (TBOT > 271.161)) THEN
                            WRITE(6,*) 'TBOTI:INDEX,I,J,MYPE,TBOT=', INDEX, I, J, MYPE, TBOT
                        END IF
!------------------------------------------ 
! SEA-ICE TEMP WITH DEPTH (STC) RANGE CHECK
!------------------------------------------ 
                        DO INDEX=1,4
                            IF ((STCK(INDEX) < 200.00) .OR. (STCK(INDEX) > 274.15)) THEN
                                WRITE(6,*) 'STCI:INDEX,I,J,MYPE,STC=',                            &
    &                                       INDEX, I, J, MYPE, STCK(INDEX)
                            END IF
                        END DO
!------------------------------------------------------------ 
! CHECK TO SEE THAT WHEN SEA-ICE, ALBASE=ALBEDO=0.6, MXSNAL=0
!------------------------------------------------------------ 
                        IF  ((ALB         < 0.59)   .OR.                                          &
    &                        (ALB         > 0.61)   .OR.                                          &
    &                        (ALBEDO(I,J) < 0.59)   .OR.                                          &
    &                        (ALBEDO(I,J) > 0.61)   .OR.                                          &
    &                        (SNOALB      > 1.E-9)) THEN
                            WRITE(6,*) 'ALBI:I,J,MYPE,ICE,SNODPK,ALB,ALB2D,SNOALB=',              &
    &                                   I, J, MYPE, ICE, SNODPK, ALB, ALBEDO(I,J), SNOALB
                        END IF
!--------------------------------------------- 
! CHECK TO SEE THAT WHEN SEA-ICE, SH2O=SMC=1.0
!--------------------------------------------- 
                        DO INDEX=1,4
                            IF  ((SMCK(INDEX) /=          1.0)  .OR.                              &
    &                            (SMCK(INDEX) /= SH2OK(INDEX))) THEN
                                WRITE(6,*) 'SMCI:INDEX,I,J,MYPE,ICE,STCK,SMC,SH2O=',              &
    &                                       INDEX, I, J, MYPE, ICE,                               &
    &                                       STCK(INDEX), SMCK(INDEX), SH2OK(INDEX)
                            END IF
                        END DO
!
                    END IF
!----------------------------------------  
! BOTH LAND AND SEA-ICE CHECKS
! SNOW WATER EQUIVALENT, SNOW DEPTH CHECK
!----------------------------------------  
                    IF (((SNODPK >  0.) .AND. (SNOWH <= 0.)) .OR.                                 &
    &                   ((SNODPK <= 0.) .AND. (SNOWH > 0.))  .OR.  (SNODPK > SNOWH)) THEN
                        WRITE(6,*) 'SNOW:I,J,MYPE,SNODPK,SNOWH=', I, J, MYPE, SNODPK, SNOWH
                    END IF
!------------------- 
! SNOW DENSITY CHECK
!------------------- 
                    IF (SNODPK > 0.) THEN
                        SNDENS = SNODPK / SNOWH
                        IF (SNDENS < 0.05) THEN
                            WRITE(6,*) 'SNDENS<5%:I,J,MYPE,SNODPK,SNOWH,SNDENS=',                 &
    &                                   I, J, MYPE, SNODPK, SNOWH, SNDENS
                        END IF
!
                        IF (SNDENS > 0.40) THEN
                            WRITE(6,*) 'SNDENS>40%:I,J,MYPE,SNODPK,SNOWH,SNDENS=',                 &
    &                                   I, J, MYPE, SNODPK, SNOWH, SNDENS
                        END IF
                    END IF
!--------------------------  
! SFC/SKIN TEMP (T1K) CHECK
!-------------------------- 
                    IF ((T1K < 200.00) .OR. (T1K > 323.15)) THEN
                        WRITE(6,*) 'T1:INDEX,I,J,MYPE,T1=', INDEX, I, J, MYPE, T1K
                    END IF
!--------------------------------------------------------------------- 
! CHECK TO SEE THAT 223.15K (-50C) =< SFCTMP, SFCTH2 <= 323.15K (+50C)
! SFCTMP = LOWEST MODEL LEVEL TEMP
! SFCTH2 = LOWEST MODEL LEVEL POT TEMP
!---------------------------------------------------------------------
                    IF (((SFCTMP < 223.15) .OR. (SFCTMP > 323.15) ) .OR.                           &
    &                   ((SFCTH2 < 223.15) .OR. (SFCTH2 > 323.15) )) THEN
                        WRITE(6,*) 'SFCTMP:I,J,MYPE,SFCTMP,SFCTH2=',                               &
    &                               I, J, MYPE, SFCTMP, SFCTH2
                    END IF
!--------------------------------------------- 
! CHECK TO SEE THAT 0W/M2 =< LWDN  <=  500W/M2
! CHECK TO SEE THAT 0W/M2 =< SOLDN <= 1200W/M2
! LWDN  = DOWNWARD LONGWAVE RADIATION
! SOLDN = DOWNWARD SOLAR RADIATION
!--------------------------------------------- 
                    IF (((LWDN  < 0.) .OR. (LWDN  >  500.))  .OR.                                 &
    &                   ((SOLDN < 0.) .OR. (SOLDN > 1200.))) THEN
                        WRITE(6,*) 'LWSOLDN:I,J,MYPE,LWDN,SOLDN=',                                &
    &                               I, J, MYPE, LWDN, SOLDN
                    END IF
!---------------------------------------------- 
! CHECK TO SEE THAT 0G/KG < Q2K,Q2SAT <= 40G/KG
! CHECK TO SEE THAT Q2K <= Q2SAT
! Q2K   = LOWEST MODEL LEVEL SPEC HUM
! Q2SAT = LOWEST MODEL LEVEL SAT SPEC HUM
!---------------------------------------------- 
                    IF (((Q2K   <= 0.) .OR. (Q2K   > 0.04)) .OR.                                  &
    &                   ((Q2SAT <= 0.) .OR. (Q2SAT > 0.04)) .OR. (Q2K > Q2SAT)) THEN
                        WRITE(6,*) 'Q2:I,J,MYPE,Q2K,Q2SAT=', I, J, MYPE, Q2K, Q2SAT
                    END IF
!---------------------------------------------- 
! CHECK TO SEE THAT 600MB =< SFCPRS <= 1050MB
! SFCPRS = SURFACE PRESSURE (PA)
!---------------------------------------------- 
                    IF ((SFCPRS <= 60000.) .OR. (SFCPRS > 105000.)) THEN
                        WRITE(6,*) 'SFCPRS:I,J,MYPE,SFCPRS=', I, J, MYPE, SFCPRS
                    END IF
!---------------------------------------------------------------------- 
! CHECK TO SEE THAT 0 =< PRCP <= 0.04 KG M-2 S-1 (=0.04MM/S = 144MM/HR)
! PRCP = PRECIP RATE (KG M-2 S-1)
!----------------------------------------------------------------------
                    IF ((PRCP < 0.) .OR. (PRCP > 0.04)) THEN
                        WRITE(6,*) 'PRCP:I,J,MYPE,PRCP=', I, J, MYPE, PRCP
                    END IF
!--------------------------------------- 
! CHECK TO SEE THAT 0M/S < CHK <= 0.1M/S
! CHK = SFC HEAT EXCHANGE COEFF (M/S)
!--------------------------------------- 
                    IF ((CHK <= 0.) .OR. (CHK > 0.1)) THEN
                        WRITE(6,*) 'CH:I,J,MYPE,CHK=', I, J, MYPE, CHK
                    END IF
!
                END IF
!
               SCHECK = Z * CHK
!
                IF (SCHECK <= 1.3E-3) THEN
                    PLFLX = 0.
                    ELFLX = 0.
                END IF
!-------------------------------------------------- 
! GCIP DIAGNOSTICS & MODIFICATION OF QFC1 OVER SNOW
!-------------------------------------------------- 
                SSROFF(I,J) = SSROFF(I,J) + RNOF1K * DTQ2
                BGROFF(I,J) = BGROFF(I,J) + RNOF2K * DTQ2
                SMSTAV(I,J) = SOILQW
                SOILTB(I,J) = TBOT
                SFCEXC(I,J) = CHK
                GRNFLX(I,J) = GFLX
!
                IF (SNO (I,J) > 0. .OR. SICE(I,J) > 0.5) THEN
                   QFC1(I,J) = QFC1(I,J) * RLIVWV
                END IF
!
                IF (SNO(I,J) > 0.) THEN
                    ACSNOM(I,J) = ACSNOM(I,J) + SMELTK
                    SNOPCX(I,J) = SNOPCX(I,J) - SMELTK / FDTLIW
                END IF
!
                POTEVP(I,J) = POTEVP(I,J) + PLFLX * FDTW
                POTFLX(I,J) = POTFLX(I,J) - PLFLX
                SUBSHX(I,J) = SUBSHX(I,J) + GFLX
!------------------------------------ 
! ETA MODEL LOWER BOUNDARY CONDITIONS
!------------------------------------ 
                THS(I,J) = T1K * APES(I,J)
!
                IF (QFC1(I,J) > 0.) QS(I,J) = QLM(I,J) + ELFLX * APES(I,J) / QFC1(I,J)
!--------------------- 
! HISTORICAL VARIABLES
!--------------------- 
!
!-------------------------------------------------- 
! DYNAMIC ALBEDO, ALBEDO, TO BE PASSED TO RADTN.F
!-------------------------------------------------- 
                ALBEDO(I,J) = ALB2D
                   SNO(I,J) = SNODPK
!--------------- 
! SNOW DEPTH, SI
!--------------- 
                    SI(I,J) = SNOWH
                   CMC(I,J) = CMCK
                SMSTOT(I,J) = SOILQM
!
                DO 150 NS=1,NSOIL
                    SMC(I,J,NS) = SMCK(NS)
!---------------------------------- 
! SH2O ARRAY (LIQUID SOIL MOISTURE)
!---------------------------------- 
                    SH2O(I,J,NS) = SH2OK(NS)
                     STC(I,J,NS) =  STCK(NS)
            150 END DO
            END IF
        
    155 END DO
!
160 END DO
!-------------------------------------------------------------------------------------------------
! SET LFIRST=FALSE SO THAT THERE ARE NOT VARIABLE/PARAMETER RANGE CHECKS FOR THE NEXT 155/160 LOOP
!-------------------------------------------------------------------------------------------------
    LFIRSTA = .FALSE. 
    LFIRST  = .FALSE. 
!---------------------------------------------------------------------- 
! VARIABLES TWBS AND QWBS COMPUTED HERE FOR GCIP.
! ACCUMULATE SURFACE HEAT FLUXES HERE.
! FOR GCIP ACCUMULATE ACTUAL AND POTENTIAL EVAPORATION.
! FOR GCIP ACCUMULATE TOTAL SNOW MELT AND THE ASSOCIATED NET HEAT FLUX.
!---------------------------------------------------------------------- 
!
!$omp parallel do private (I      , J)
!
    DO 200 J=MYJS2,MYJE2
        DO 200 I=MYIS,MYIE
            TWBS(I,J) = (THLM(I,J) -  THS(I,J) * (1.-SM(I,J)) - THZ0(I,J) * SM(I,J))              &
    &                 *  FFS (I,J) / APES(I,J)
            QWBS(I,J) = (QLM (I,J) -  QS (I,J) * (1.-SM(I,J)) - QZ0 (I,J) * SM(I,J))              &
    &                 *  QFC1(I,J) / APES(I,J)
!
            SFCSHX(I,J) = SFCSHX(I,J) + TWBS(I,J)
            SFCLHX(I,J) = SFCLHX(I,J) + QWBS(I,J)
            SFCEVP(I,J) = SFCEVP(I,J) - QWBS(I,J) * FDTW
            POTEVP(I,J) = POTEVP(I,J) - QWBS(I,J) * SM(I,J) * FDTW
            POTFLX(I,J) = POTFLX(I,J) + QWBS(I,J) * SM(I,J)
!------------------------------------------       
! IF COLD ENOUGH, IT SNOWS (IN NOAH LSM)...
! FOR GCIP ACCUMULATE TOTAL SNOWFALL.
!------------------------------------------         
            IF (THLM(I,J)/APELM(I,J) <= T0 .AND. SICE(I,J)+SM(I,J) < 0.5) THEN
                ACSNOW(I,J) = ACSNOW(I,J) + PREC(I,J)
!--------------------   
! OTHERWISE IT RAINS.
!-------------------- 
            ELSE
                ACCLIQ(I,J) = ACCLIQ(I,J) + PREC(I,J)
            END IF
!        
!
              PD(I,J) = PD(I,J) - PREC(I,J) * GROW
!
            PREC(I,J) = 0.
200 END DO
!---------------------------- 
! LONGWAVE OUTGOING RADIATION
!---------------------------- 
!
!$omp parallel do private (I      , J       , TSFC    , TSFC2)
!
    DO 210 J=MYJS2,MYJE2
        DO 210 I=MYIS,MYIE
            TSFC = THS(I,J) / APES(I,J)
            TSFC2 = TSFC * TSFC
            RADOT(I,J) = HBM2(I,J) * EPSR(I,J) * STBOL * TSFC2 * TSFC2
210 END DO
!------------------------------------------------------------------  
! INCREMENT TIME STEP COUNTERS FOR USE IN COMPUTING TIME AVE VALUES
!------------------------------------------------------------------ 
    APHTIM = APHTIM + 1.
    ARDSW  = ARDSW  + 1.
    ARDLW  = ARDLW  + 1.
    ASRFC  = ASRFC  + 1.
!
    RETURN
!
    END SUBROUTINE SURFCE
