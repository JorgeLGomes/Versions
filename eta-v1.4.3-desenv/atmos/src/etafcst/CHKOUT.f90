    SUBROUTINE CHKOUT
!>-------------------------------------------------------------------------------------------------- 
!> SUBROUTINE CHKOUT
!>
!> SUBPROGRAM: CHKOUT - POSTS PROFILES AND OUTPUT POST DATA
!> PROGRAMMER: TREADON
!> ORG: W/NP2
!> DATE: 93-02-26
!
!> ABSTRACT: 
!> THIS ROUTINE POSTS PROFILE DATA AND WRITES COMMON BLOCKS TO TEMPORARY FILE FOR USE BY THE POST
!> PROCESSOR.
!> OPTIONALLY, IF RUN UNDER PSHELL THIS ROUTINE WILL SUBMIT POST JOBS AS THE MODEL RUNS.
!> THIS ROUTINE REPLACES ETA MODEL SUBROUTINE OUTMAP.
!>
!> PROGRAM HISTORY LOG:
!> 93-02-26  RUSS TREADON
!> 93-08-30  RUSS TREADON - ADDED DOCBLOC AND DIAGNOSTIC PROFILES
!> 95-03-31  T BLACK      - CONVERTED FROM 1-D TO 2-D IN HORIZONTAL
!> 95-07-31  MIKE BALDWIN - REMOVED SOUNDING DIAGNOSTICS AND BUFR
!> 96-03-13  F MESINGER   - IMPROVED REDUCTION TO SEA LEVEL (TO ACHIEVE EXACT CONSISTENCY WITH THE
!>                          MODELS HYDROSTATIC EQUATION NEXT TO MOUNTAIN SIDES)
!> 96-04-12  MIKE BALDWIN - MODIFIED SOUNDING OUTPUT
!> 96-10-31  T BLACK      - MODIFICATIONS FOR GENERATIONS OF NESTS BCS
!> 98-11-17  T BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
!> 99-05-03  T BLACK      - SLP REDUCTION, BCEX, AND PROFILES REMOVED. EACH PE WRITES ITS OWN 
!>                          MINI-RESTRT FILE
!> 00-08-01  JIM TUCCILLO - QUILT SERVER CAPABILITY ADDED
!> 00-10-11  T BLACK      - MODIFICATIONS FOR RESTART CAPABILITY
!> 18-01-15  LUCCI        - MODERNIZATION OF THE CODE, INCLUDING:
!>                          * F77 TO F90/F95
!>                          * INDENTATION & UNIFORMIZATION CODE
!>                          * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                          * DOCUMENTATION WITH DOXYGEN
!>                          * OPENMP FUNCTIONALITY
!>
!> INPUT  ARGUMENT LIST:
!> NONE
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: ACMCLD
!>              ACMCLH
!>              ACMPRE
!>              ACMRDL
!>              ACMRDS
!>              ACMSFC
!>              BOCO
!>              BUFFER
!>              CLDWTR
!>              CNVCLD
!>              CONTIN
!>              CTLBLK
!>              CUINIT
!>              DYNAM
!>              F77KINDS
!>              GLB_TABLE
!>              INDX
!>              LOOPS
!>              MAPOT
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              NHYDRO
!>              OUTFIL
!>              PARMETA
!>              PARM_TBL
!>              PHYS
!>              PRFHLD
!>              PVRBLS
!>              SOIL
!>              TEMPCOM
!>              TEMPV
!>              TIMCHK
!>              TIMMING
!>              VRBLS
!>
!> DRIVER     : EBU
!>
!> CALLS      : MPI_BARRIER
!>              MPI_COAL
!>              MPI_REDUCE
!>--------------------------------------------------------------------------------------------------
    USE ACMCLD
    USE ACMCLH
    USE ACMPRE
    USE ACMRDL
    USE ACMRDS
    USE ACMSFC
    USE BOCO
    USE BUFFER
    USE CLDWTR
    USE CNVCLD
    USE CONTIN
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
    USE NHYDRO
    USE OUTFIL
    USE PARMETA
    USE PARM_TBL
    USE PHYS
    USE PRFHLD
    USE PVRBLS
    USE SOIL
    USE TEMPCOM
    USE TEMPV
    USE TIMCHK
    USE TIMMING
    USE VRBLS
    USE TOPO
!
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER  :: IMJM  = IM * JM - JM / 2 
    INTEGER(KIND=I4KIND), PARAMETER  :: IMT   = 2  * IM - 1
    INTEGER(KIND=I4KIND), PARAMETER  :: JMT   = JM / 2  + 1
    INTEGER(KIND=I4KIND), PARAMETER  :: NRLX1 = 250
    INTEGER(KIND=I4KIND), PARAMETER  :: NRLX2 = 100
!
    REAL   (KIND=R4KIND), PARAMETER  :: CAPA = 0.285896  
!------------------
! DECLARE VARIABLES
!------------------
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & OUTJ    , STDRD   , MESO    , ONHOUR  , EXBC    , MULTIWRITE
!
    CHARACTER(LEN=8)                                                                            ::&
    & FHR
!
    CHARACTER(LEN=24)                                                                           ::&
    & OUTJOB
!
    CHARACTER(LEN=16)                                                                           ::&
    & ASSIGN
!
    CHARACTER(LEN=8)                                                                            ::&
    & CITAG
!
    CHARACTER(LEN=8)                                                                            ::&
    & ASTMRK  , TMYY
!
    CHARACTER(LEN=200)                                                                          ::&
    & SUBMIT ! CHOU
!
    CHARACTER(LEN=200)                                                                          ::&
    & CHMO
!
    CHARACTER(LEN=32)                                                                           ::&
    & LABEL
!
    INTEGER(KIND=I4KIND), DIMENSION(4)                                                          ::&
    & LABINT
!
    CHARACTER(LEN=88)                                                                           ::&
    & LINE
!
    CHARACTER(LEN=8)    , DIMENSION(85)                                                         ::&
    & LINE1
!
    CHARACTER(LEN=8)                                                                            ::&
    & RESTHR
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & CHX
!
    EQUIVALENCE (LABEL, LABINT)
    EQUIVALENCE (LINE , LINE1 )
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & PSLP    , PDS     , FACTR 
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & SWTTC   , TTND
!
    INTEGER(KIND=I4KIND), DIMENSION(0:INPES*JNPES-1)                                            ::&
    & IKNTS   , IDISP
!
    REAL   (KIND=R4KIND), ALLOCATABLE, DIMENSION(:,:,:)                                         ::&
    & TEMPSOIL
!
    CHARACTER(LEN=48)                                                                           ::&
    & FINFIL
!
    CHARACTER(LEN=8)                                                                            ::&
    & DONE
!---------------------
! DECLARE EQUIVALENCES
!---------------------
    EQUIVALENCE (TTND(1,1,1), SWTTC(1,1,1))
!
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & JSTAT
!
    REAL   (KIND=R8KIND), DIMENSION(LM)                                                         ::&
    & SUMT    , SUMT_0  , SUMT2   , SUMT2_0
!
    REAL   (KIND=R8KIND)                                                                        ::&
    & STDEV   , RMS     , TMEAN
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                                         ::&
    & TMAX    , TMAX_0  ,  TMIN   ,  TMIN_0
!
    REAL   (KIND=R8KIND)                                                                        ::&
    & STRWAIT , ENDWAIT , RTC
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IHS
!
    DATA IHS /MPI_REQUEST_NULL/
!
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & STATUS
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ISERVE
!
    DATA ISERVE / 1 /
!------------------------------------------------
! THE FOLLOWING ARE USED FOR TIMIMG PURPOSES ONLY
!------------------------------------------------
    REAL   (KIND=R8KIND)                                                                        ::&
    & TIMEF
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & MPP_TIM , INIT_TIM
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , NTSPH   , LLMH    , LMHK    , IHR     , JL      , JJ      , IIMAX   ,   &
    & IRTN    , ISTAT   , K       , IER     , N       , IPMAX   , IERR    , IBLOCK  , IRECLEN ,   &
    & ISTART  , IEND    , MYPE_ROW, LUNIN   , LUNOT   , IDX 
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & TLL1    , TERM1   , BTIM0   , BTIMW   , DUMMY   , DIF_TIM , WRT_TIM_0 , WIND
!
    REAL   (KIND=R8KIND)                                                                        ::&
    & TIMES 
!
    REAL   (KIND=R8KIND), PARAMETER :: R3600 =  3600.
!
    NAMELIST /POSTLIST/ CHMO 
    NAMELIST /POSTLIST/ SUBMIT    !CHOU
!------------------
! START CHKOUT HERE
!------------------
!--------------------------------------------------------------------------------------------------
! ON FIRST ENTRY INITIALIZE THE OUTPUT FILE TAG TO ZERO AND DO PRELIMINARY PROFILE DATA ASSIGNMENTS
!--------------------------------------------------------------------------------------------------
    IF (NTSD == 1 .OR. RESTRT) THEN
        ITAG = 0
!    
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                LMHK = LMH(I,J)
                TLL1 = T(I,J,LMHK)
                TLMIN(I,J) = TLL1
                TLMAX(I,J) = TLL1
!Lyra GSM Max Wind
                MAXWU(I,J)   = U10(I,J)
                MAXWV(I,J)   = V10(I,J)
                MAXWIND(I,J) = SQRT((U10(I,J)**2)+(V10(I,J)**2))
!Lyra GSM Max Wind
            END DO
        END DO
    END IF
!--------------------------------------
! UPDATE MAX AND MIN LOWEST LAYER TEMPS
!--------------------------------------
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            LMHK = LMH(I,J)
            TLL1 = T(I,J,LMHK)
            IF (TLL1 < TLMIN(I,J)) TLMIN(I,J) = TLL1
            IF (TLL1 > TLMAX(I,J)) TLMAX(I,J) = TLL1
        END DO
    END DO
!
!Lyra GSM Max Wind
!***  UPDATE MAX WIND
!***
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
           WIND=SQRT((U10(I,J)**2)+(V10(I,J)**2))
           IF ( WIND > MAXWIND(I,J) ) THEN
               MAXWU(I,J)   = U10(I,J)
               MAXWV(I,J)   = V10(I,J)
               MAXWIND(I,J) = WIND
           ENDIF
        ENDDO
    ENDDO
!Lyra GSM Max Wind
!
!---------------------------------------------
! FIGURE OUT JUST WHERE IN THE FORECAST WE ARE
!---------------------------------------------
    NTSPH  = INT(3600. / DT + 0.50)
    TIMES  = (NTSD - 1) * DT
    ONHOUR = .FALSE. 
!
    IF ((MOD(TIMES,R3600) == 0.) .OR. (MOD(TIMES,R3600) > 3600.-DT)) ONHOUR = .TRUE. 
!    IF ( MOD(NTSD,3600/DT) == 1) ONHOUR = .TRUE. 
!--------------------------------------------------------------------------------------------------
! IF THE CURRENT FORECAST TIME IS A FULL HOUR OR EQUALS A FULL BLOWN POST TIME THEN WRITE THE FIELD
! IF NOT, EXIT THIS ROUTINE
!--------------------------------------------------------------------------------------------------
    IF ((NTSD == NSHDE) .OR. ONHOUR) GOTO 100
    IF (NSTART > 0 .AND. NSTART == NSHDE .AND.NTSD == NSHDE) GOTO 100
!--------------------------------------------------------------------------------------------------
! BEGIN: INITIALIZE CONVECTIVE CLOUD FIELDS FOR RADIATION BEFORE RETURNING TO EBU (FERR. 23 JAN 02)
!--------------------------------------------------------------------------------------------------
    IF (CURAD) THEN
        IF (MYPE == 0) THEN
            WRITE(0,"(a)") 'CHKOUT: INITIALIZE CUPPT,HTOP,HBOT'
            WRITE(6,"(a)") 'CHKOUT: INITIALIZE CUPPT,HTOP,HBOT'
	    WRITE(6,*) 'CHKOUT: CURAD, NTSD', CURAD, NTSD
        END IF
!
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                CUPPT(I,J) =   0.
                 HTOP(I,J) = 100.
                 HBOT(I,J) =   0.
            END DO
        END DO
!
        CURAD = .FALSE. 
!
    END IF
!
    RETURN
!-------------------------------------------------------------------------------------------
! IT IS TIME TO WRITE TO THE PROFILE FILE AND/OR WRITE TEMPORARY FILES FOR A FULL BLOWN POST
!-------------------------------------------------------------------------------------------
100 CONTINUE
    IF (MYPE == 0) WRITE(6,*) 'CHKOUT: ONHOUR, NTSD', ONHOUR, NTSD

!--------------------------------------------------
! (CHOU 21 JUN 09)
! PDS IS SURFACE PRESSURE
! TSHLTR HOLDS THE 2M THETA, CONVERT TO TEMPERATURE
! TERM1 IS 2m*G/(Rd*T)
!--------------------------------------------------
!
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            LLMH = LMH(I,J)
            PDS(I,J) = PD(I,J) + PT
            TERM1 = -0.068283 / T(I,J,LLMH)
            PSHLTR(I,J) = PDS(I,J) * EXP(TERM1)
            TSHLTR1(I,J) = TSHLTR(I,J) * (PSHLTR(I,J) * 1.E-5) ** CAPA
                  
            IF (CZMEAN(I,J) > 0.) THEN
                FACTR(I,J) = CZEN(I,J) / CZMEAN(I,J)
            ELSE
                FACTR(I,J) = 0.
            END IF
!
        END DO
    END DO
!-------------------------------------------------------------------------------         
! MAKE SURE POST DOES NOT BLOW UP WHEN COMPUTING RH ON THE GLOBAL N/S BOUNDARIES
!-------------------------------------------------------------------------------        
    IF (MYPE < INPES) THEN
        DO J=1,2
            DO I=MYIS,MYIE
                TSHLTR1(I,J) = TSHLTR1(I,3)
                 QSHLTR(I,J) =  QSHLTR(I,3)
            END DO
        END DO
    END IF
!       
    IF (MYPE >= NPES-INPES) THEN
        DO J=MYJE-1,MYJE
            DO I=MYIS,MYIE
                TSHLTR1(I,J) = TSHLTR1(I,MYJE-2)
                 QSHLTR(I,J) =  QSHLTR(I,MYJE-2)
            END DO
        END DO
    END IF
!         
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            TLMIN(I,J) = AMIN1(TSHLTR1(I,J), TLMIN(I,J))
            TLMAX(I,J) = AMAX1(TSHLTR1(I,J), TLMAX(I,J))
!Lyra GSM Max Wind
            WIND=SQRT((U10(I,J)**2)+(V10(I,J)**2))
            IF ( WIND > MAXWIND(I,J) ) THEN
                 MAXWU(I,J)   = U10(I,J)
                 MAXWV(I,J)   = V10(I,J)
                 MAXWIND(I,J) = AMAX1 ( WIND, MAXWIND(I,J) )
            ENDIF
!Lyra GSM Max Wind
!
        END DO
    END DO
!------------------
! SET FORECAST HOUR
!------------------
    IHR = NTSD / TSPH + 0.5
!----------------------------------------------------------------------------------------------------------------------------------- 
! IF THIS IS NOT A FULL BLOWN OUTPUT TIME, SKIP THE RESTART FILE AND POST JOB WRITES AND GO TO SECTION WHERE ACCUMULATION ARRAYS
! ARE ZEROED OUT IF NECESSARY
!-----------------------------------------------------------------------------------------------------------------------------------
    IF (NTSD /= NSHDE .AND. NSTART+1 /= NSHDE) GOTO 1310
!    IF ( .NOT. ONHOUR) GOTO 1310
!-------------------------------
! COMPUTE TEMPERATURE STATISTICS
!-------------------------------
    BTIM0 = TIMEF()
!
    DO 900 JL=1,LM
!    
         TMAX(JL) = -1.E6
         TMIN(JL) =  1.E6
         SUMT(JL) =  0.
        SUMT2(JL) =  0.
!    
        JJ=0
        DO J=MY_JS_GLB,MY_JE_GLB
            JJ = JJ + 1
!
            IF (MOD(J+1,2) /= 0 .AND. MY_IE_GLB == IM) THEN
                IIMAX = MY_IE_LOC - 1
            ELSE
                IIMAX = MY_IE_LOC
            END IF
!
            DO I=MYIS,IIMAX
                 SUMT(JL) =  SUMT(JL) + T(I,JJ,JL)
                SUMT2(JL) = SUMT2(JL) + T(I,JJ,JL) ** 2
                 TMAX(JL) = AMAX1(TMAX(JL), T(I,JJ,JL))
                 TMIN(JL) = AMIN1(TMIN(JL), T(I,JJ,JL))
            END DO
        END DO
!
900 END DO
!-------------
! GLOBAL STATS
!-------------
    CALL MPI_REDUCE( SUMT,  SUMT_0, LM, MPI_REAL8, MPI_SUM, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_REDUCE(SUMT2, SUMT2_0, LM, MPI_REAL8, MPI_SUM, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_REDUCE( TMAX,  TMAX_0, LM, MPI_REAL , MPI_MAX, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_REDUCE( TMIN,  TMIN_0, LM, MPI_REAL , MPI_MIN, 0, MPI_COMM_COMP, IRTN)
!
    IF (MYPE == 0) THEN
        DO JL=1,LM
            TMEAN = SUMT_0(JL) / DBLE(IMJM)
            STDEV = DSQRT((DBLE(IMJM)*SUMT2_0(JL)-SUMT_0(JL)**2) / DBLE(DBLE(IMJM)*(DBLE(IMJM-1))))
              RMS = DSQRT(SUMT2_0(JL) / DBLE(IMJM))
        END DO
    END IF
!
    STAT_TIM = STAT_TIM + TIMEF() - BTIM0
!----------------------------------------------------------------
! WE REACH THE CODE BELOW ONLY IF IT IS A FULL BLOWN POSTING TIME
! WRITE DATA REQUIRED TO RESTART THE MODEL/INITIALIZE THE POST
!----------------------------------------------------------------
    CALL MPI_BARRIER(MPI_COMM_COMP, ISTAT)
!--------------------------------------------------
! PDS IS SURFACE PRESSURE
! TSHLTR HOLDS THE 2M THETA, CONVERT TO TEMPERATURE
! TERM1 IS 2M*G/(RD*T)
!--------------------------------------------------
!
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            LLMH = LMH(I,J)
            PDS(I,J) = PD(I,J) + PT
            TERM1 = -0.068283 / T(I,J,LLMH)
            PSHLTR(I,J) = PDS(I,J) * EXP(TERM1)
            TSHLTR(I,J) = TSHLTR(I,J) * (PSHLTR(I,J) * 1.E-5) ** CAPA
!
            IF (CZMEAN(I,J) > 0.) THEN
                FACTR(I,J) = CZEN(I,J) / CZMEAN(I,J)
            ELSE
                FACTR(I,J) = 0.
            END IF
!        
        END DO
    END DO
!-------------------------------------------------------------------------------
! MAKE SURE POST DOES NOT BLOW UP WHEN COMPUTING RH ON THE GLOBAL N/S BOUNDARIES
!-------------------------------------------------------------------------------
    IF (MYPE < INPES) THEN
        DO J=1,2
            DO I=MYIS,MYIE
                TSHLTR(I,J) = TSHLTR(I,3)
                QSHLTR(I,J) = QSHLTR(I,3)
            END DO
        END DO
    END IF
!
    IF (MYPE >= NPES-INPES) THEN
        DO J=MYJE-1,MYJE
            DO I=MYIS,MYIE
                TSHLTR(I,J) = TSHLTR(I,MYJE-2)
                QSHLTR(I,J) = QSHLTR(I,MYJE-2)
            END DO
        END DO
    END IF
!----------------------------------------
! SWTTC IS THE CURRENT SW TEMP TENDENCIES
!----------------------------------------
!
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
    DO K=1,LM
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                SWTTC(I,J,K) = RSWTT(I,J,K) * FACTR(I,J)
            END DO
        END DO
    END DO
!----------------------------------------
! TTND IS THE CURRENT RAD TEMP TENDENCIES
!----------------------------------------
    DO K=1,LM
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                TTND(I,J,K) = RLWTT(I,J,K) + SWTTC(I,J,K)
            END DO
        END DO
    END DO
!-----------------------------
! CREATE NAME FOR RESTART FILE
!-----------------------------
    ITAG = NTSD / TSPH + 0.5
!
    CALL GETENV("TMMARK",RESTHR)
!
    IF (RESTHR == '    ') THEN
        WRITE(RSTFIL,1150) ITAG
        1150 FORMAT('RESTRT',I6.6,'.quilt')
    ELSE IF (RESTHR == 'TM00' .AND. IQUILT_GROUP > 0) THEN
        WRITE(RSTFIL,1152) ITAG, MYPE
        1152 FORMAT('RESTRT',I6.6,'.',I3.3)
        MULTIWRITE = .FALSE. 
        IF (NTSD == NTSTM) MULTIWRITE = .TRUE. 
    ELSE
        MULTIWRITE = .FALSE. 
        WRITE(RSTFIL,1155) ITAG, RESTHR
        1155 FORMAT('RESTRT',I6.6,'.quilt.',A4)
    END IF
!--------------------------
! OPEN UNIT TO RESTART FILE
!--------------------------
    LRSTRT  = 8
!
    WRT_TIM = 0.
    BTIMW   = TIMEF()
    BTIM0   = TIMEF()
!
    CLOSE(LRSTRT)
!
    IF (MULTIWRITE) THEN
        OPEN (UNIT=LRSTRT, FILE=RSTFIL, FORM='UNFORMATTED', IOSTAT=IER)
        IF (IER /= 0) WRITE(LIST,*)' LRSTRT OPEN UNIT ERROR IER=',IER
    END IF
!-------------------------------------
! BE SURE THAT THE BUFFER IF AVAILABLE
!-------------------------------------
    CALL COAL(DUMMY,-1)
!----------------------------------------------------
! WRITE DATE AND TIMESTEP INFORMATION TO RESTART FILE
!----------------------------------------------------
    LABEL = 'OMEGA-ALPHA*DT/CP'
!
    IF (MULTIWRITE) WRITE(LRSTRT) RUN, IDAT, IHRST, NTSD, LABEL, ITAG
!
    CALL COAL(RUN  ,1)
    CALL COAL(IDAT ,3)
    CALL COAL(IHRST,1)
    CALL COAL(NTSD ,1)
    CALL COAL(LABEL,8)
!------------------------------
! BEGIN WRITING THE RESTRT FILE
!------------------------------
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((PD(I,J),I=1,MYIE),J=1,MYJE), ((RES(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL( PD(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(RES(1:MYIE,1:MYJE), MYIE*MYJE)
!
    DO K=1,LM
        IF (MULTIWRITE) THEN
            WRITE(LRSTRT) ((OMGALF(I,J,K),I=1,MYIE),J=1,MYJE)
        END IF
!    
        CALL COAL(OMGALF(1:MYIE,1:MYJE,K), MYIE*MYJE)
    END DO
!
    LABEL = 'BND, PD, RES, T, Q, U, V, Q2, TTND, CWM, TRAIN, TCUCN'
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) RUN, IDAT, IHRST, NTSD, LABEL, ITAG, FIRST, IOUT, NSHDE
    END IF
!
    CALL COAL(RUN  ,1)
    CALL COAL(IDAT ,3)
    CALL COAL(IHRST,1)
    CALL COAL(NTSD ,1)
    CALL COAL(LABEL,8)
    CALL COAL(FIRST,1)
    CALL COAL(IOUT ,1)
    CALL COAL(NSHDE,1)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((PD(I,J),I=1,MYIE),J=1,MYJE),                                              &
    &                ((RES(I,J),I=1,MYIE),J=1,MYJE),                                              &
    &                ((FIS(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL(PD (1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(RES(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(FIS(1:MYIE,1:MYJE), MYIE*MYJE)
!-------------------------------------------------
! BOUNDARY CONDITION WRITE CHANGED TO BLANK RECORD
!-------------------------------------------------
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) PDB, TB, QB, UB, VB, Q2B, CWMB
    END IF
!
    CALL COAL(PDB , LB*2)
    CALL COAL(TB  , LB*LM*2)
    CALL COAL(QB  , LB*LM*2)
    CALL COAL(UB  , LB*LM*2)
    CALL COAL(VB  , LB*LM*2)
    CALL COAL(Q2B , LB*LM*2)
    CALL COAL(CWMB, LB*LM*2)
!
    DO K=1,LM
        IF (MULTIWRITE) THEN
            WRITE(LRSTRT) ((T(I,J,K),I=1,MYIE),J=1,MYJE)
        END IF
!
        CALL COAL(T(1:MYIE,1:MYJE,K), MYIE*MYJE)
!  
        IF (MULTIWRITE) THEN
            WRITE(LRSTRT) ((Q(I,J,K),I=1,MYIE),J=1,MYJE)
        END IF
!
        CALL COAL(Q(1:MYIE,1:MYJE,K), MYIE*MYJE)
!    
        IF (MULTIWRITE) THEN
            WRITE(LRSTRT) ((U(I,J,K),I=1,MYIE),J=1,MYJE)
        END IF
!
        CALL COAL(U(1:MYIE,1:MYJE,K), MYIE*MYJE)
    
        IF (MULTIWRITE) THEN
            WRITE(LRSTRT) ((V(I,J,K),I=1,MYIE),J=1,MYJE)
        END IF
!
        CALL COAL(V(1:MYIE,1:MYJE,K), MYIE*MYJE)
!    
        IF (MULTIWRITE) THEN
            WRITE(LRSTRT) ((Q2(I,J,K),I=1,MYIE),J=1,MYJE)
        END IF
!
        CALL COAL(Q2(1:MYIE,1:MYJE,K), MYIE*MYJE)
!    
        IF (MULTIWRITE) THEN
            WRITE(LRSTRT) ((TTND(I,J,K),I=1,MYIE),J=1,MYJE)
        END IF
!
        CALL COAL(TTND(1:MYIE,1:MYJE,K), MYIE*MYJE)
!    
        IF (MULTIWRITE) THEN
            WRITE(LRSTRT) ((CWM(I,J,K),I=1,MYIE),J=1,MYJE)
        END IF
!
        CALL COAL(CWM(1:MYIE,1:MYJE,K), MYIE*MYJE)
!    
        IF (MULTIWRITE) THEN
            WRITE(LRSTRT) ((TRAIN(I,J,K),I=1,MYIE),J=1,MYJE)
        END IF
!
        CALL COAL(TRAIN(1:MYIE,1:MYJE,K), MYIE*MYJE)
!    
        IF (MULTIWRITE) THEN
            WRITE(LRSTRT) ((TCUCN(I,J,K),I=1,MYIE),J=1,MYJE)
        END IF
!
        CALL COAL(TCUCN(1:MYIE,1:MYJE,K), MYIE*MYJE)
!
    END DO
!
    LABEL = 'MISC VARIABLES'
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) RUN, IDAT, IHRST, NTSD, LABEL, ITAG                                ,        &
    &                 ((RSWIN(I,J),I=1,MYIE),J=1,MYJE), ((RSWOUT(I,J),I=1,MYIE),J=1,MYJE),        &
    &                    ((TG(I,J),I=1,MYIE),J=1,MYJE),     ((Z0(I,J),I=1,MYIE),J=1,MYJE),        &
    &                  ((AKMS(I,J),I=1,MYIE),J=1,MYJE),   ((CZEN(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL(RUN,  1)
    CALL COAL(IDAT, 3)
    CALL COAL(IHRST,1)
    CALL COAL(NTSD, 1)
    CALL COAL(LABEL,8)
    CALL COAL( RSWIN(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(RSWOUT(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(    TG(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(    Z0(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(  AKMS(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(  CZEN(1:MYIE,1:MYJE), MYIE*MYJE)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((AKHS(I,J),I=1,MYIE),J=1,MYJE),  ((THS(I,J),I=1,MYIE),J=1,MYJE),           &
    &                   ((QS(I,J),I=1,MYIE),J=1,MYJE), ((TWBS(I,J),I=1,MYIE),J=1,MYJE),           &
    &                 ((QWBS(I,J),I=1,MYIE),J=1,MYJE), ((HBOT(I,J),I=1,MYIE),J=1,MYJE),           &
    &               ((CFRACL(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL(  AKHS(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(   THS(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(    QS(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(  TWBS(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(  QWBS(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(  HBOT(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(CFRACL(1:MYIE,1:MYJE), MYIE*MYJE)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((THZ0(I,J),I=1,MYIE),J=1,MYJE),  ((QZ0(I,J),I=1,MYIE),J=1,MYJE),           &
    &                  ((UZ0(I,J),I=1,MYIE),J=1,MYJE),  ((VZ0(I,J),I=1,MYIE),J=1,MYJE),           &
    &                ((USTAR(I,J),I=1,MYIE),J=1,MYJE), ((HTOP(I,J),I=1,MYIE),J=1,MYJE),           &
    &               ((CFRACM(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL(  THZ0(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(   QZ0(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(   UZ0(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(   VZ0(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL( USTAR(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(  HTOP(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(CFRACM(1:MYIE,1:MYJE), MYIE*MYJE)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((SNO(I,J),I=1,MYIE),J=1,MYJE), ((   SI(I,J),I=1,MYIE),J=1,MYJE),           &
    &              ((CLDEFI(I,J),I=1,MYIE),J=1,MYJE), ((   RF(I,J),I=1,MYIE),J=1,MYJE),           &
    &                ((PSLP(I,J),I=1,MYIE),J=1,MYJE), ((CUPPT(I,J),I=1,MYIE),J=1,MYJE),           &
    &              ((CFRACH(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL(   SNO(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(    SI(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(CLDEFI(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(    RF(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(  PSLP(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL( CUPPT(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(CFRACH(1:MYIE,1:MYJE), MYIE*MYJE)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((SOILTB(I,J),I=1,MYIE),J=1,MYJE), ((SFCEXC(I,J),I=1,MYIE),J=1,MYJE),       &
    &                 ((SMSTAV(I,J),I=1,MYIE),J=1,MYJE), ((SMSTOT(I,J),I=1,MYIE),J=1,MYJE),       &
    &                 ((GRNFLX(I,J),I=1,MYIE),J=1,MYJE), ((PCTSNO(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL(SOILTB(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(SFCEXC(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(SMSTAV(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(SMSTOT(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(GRNFLX(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(PCTSNO(1:MYIE,1:MYJE), MYIE*MYJE)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((RLWIN(I,J),I=1,MYIE),J=1,MYJE), ((RADOT(I,J),I=1,MYIE),J=1,MYJE),         & 
    &                ((CZMEAN(I,J),I=1,MYIE),J=1,MYJE), ((SIGT4(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL( RLWIN(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL( RADOT(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(CZMEAN(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL( SIGT4(1:MYIE,1:MYJE), MYIE*MYJE)
!--------------------------------------------------------------------------------------------------
! U00, UL, AND LC ARE NO LONGER USED AND WILL BE REMOVED WHEN THE NEW CLOUD FIELDS ARE PUT INTO THE 
! RESTART FILES AND THE POST (FERRIER 23 JAN 02)
!--------------------------------------------------------------------------------------------------  
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((U00(I,J),I=1,MYIE),J=1,MYJE), UL, ((LC(I,J),I=1,MYIE),J=1,MYJE),          &
    &                  ((SR(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL(U00(1:MYIE,1:MYJE),MYIE*MYJE)
    CALL COAL(UL,2*LM)
    CALL COAL( LC(1:MYIE,1:MYJE),MYIE*MYJE)
    CALL COAL( SR(1:MYIE,1:MYJE),MYIE*MYJE)
!
    LABEL = 'ACCUMULATED VARIABLES'
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) RUN, IDAT, IHRST, NTSD, LABEL, ITAG,                                        &
    &                 ((PREC(I,J),I=1,MYIE),J=1,MYJE), ((ACPREC(I,J),I=1,MYIE),J=1,MYJE),         &
    &               ((ACCLIQ(I,J),I=1,MYIE),J=1,MYJE), ((CUPREC(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL(RUN  ,1)
    CALL COAL(IDAT ,3)
    CALL COAL(IHRST,1)
    CALL COAL(NTSD ,1)
    CALL COAL(LABEL,8)
!
    CALL COAL(  PREC(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(ACPREC(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(ACCLIQ(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(CUPREC(1:MYIE,1:MYJE), MYIE*MYJE)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((ACFRCV(I,J),I=1,MYIE),J=1,MYJE), ((NCFRCV(I,J),I=1,MYIE),J=1,MYJE),       & 
    &                 ((ACFRST(I,J),I=1,MYIE),J=1,MYJE), ((NCFRST(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL(ACFRCV(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(NCFRCV(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(ACFRST(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(NCFRST(1:MYIE,1:MYJE), MYIE*MYJE)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((ACSNOW(I,J),I=1,MYIE),J=1,MYJE), ((ACSNOM(I,J),I=1,MYIE),J=1,MYJE),       &
    &                 ((SSROFF(I,J),I=1,MYIE),J=1,MYJE), ((BGROFF(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL(ACSNOW(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(ACSNOM(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(SSROFF(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(BGROFF(1:MYIE,1:MYJE), MYIE*MYJE)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((SFCSHX(I,J),I=1,MYIE),J=1,MYJE), ((SFCLHX(I,J),I=1,MYIE),J=1,MYJE),       &
    &                 ((SUBSHX(I,J),I=1,MYIE),J=1,MYJE), ((SNOPCX(I,J),I=1,MYIE),J=1,MYJE),       &
    &                 ((SFCUVX(I,J),I=1,MYIE),J=1,MYJE), ((SFCEVP(I,J),I=1,MYIE),J=1,MYJE),       &
    &                 ((POTEVP(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL(SFCSHX(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(SFCLHX(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(SUBSHX(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(SNOPCX(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(SFCUVX(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(SFCEVP(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(POTEVP(1:MYIE,1:MYJE), MYIE*MYJE)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((ASWIN(I,J),I=1,MYIE),J=1,MYJE), ((ASWOUT(I,J),I=1,MYIE),J=1,MYJE),         &
    &                ((ASWTOA(I,J),I=1,MYIE),J=1,MYJE),  ((ALWIN(I,J),I=1,MYIE),J=1,MYJE),         &
    &                ((ALWOUT(I,J),I=1,MYIE),J=1,MYJE), ((ALWTOA(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL( ASWIN(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(ASWOUT(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(ASWTOA(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL( ALWIN(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(ALWOUT(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(ALWTOA(1:MYIE,1:MYJE), MYIE*MYJE)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ARDSW, ARDLW, ASRFC, AVRAIN, AVCNVC
    END IF
!
    CALL COAL( ARDSW,1)
    CALL COAL( ARDLW,1)
    CALL COAL( ASRFC,1)
    CALL COAL(AVRAIN,1)
    CALL COAL(AVCNVC,1)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((TH10(I,J),I=1,MYIE),J=1,MYJE),    ((Q10(I,J),I=1,MYIE),J=1,MYJE),         &
    &                  ((U10(I,J),I=1,MYIE),J=1,MYJE),    ((V10(I,J),I=1,MYIE),J=1,MYJE),         &
    &               ((TSHLTR(I,J),I=1,MYIE),J=1,MYJE), ((QSHLTR(I,J),I=1,MYIE),J=1,MYJE),         &
    &               ((PSHLTR(I,J),I=1,MYIE),J=1,MYJE),                                            &
!---------
! SM V100M
!---------
    &                ((TH100(I,J),I=1,MYIE),J=1,MYJE),   ((Q100(I,J),I=1,MYIE),J=1,MYJE),         &
    &                 ((U100(I,J),I=1,MYIE),J=1,MYJE),   ((V100(I,J),I=1,MYIE),J=1,MYJE),         &
! Lyra GSM Wind stress
    &             ((XMOMFLUX(I,J),I=1,MYIE),J=1,MYJE),   ((YMOMFLUX(I,J),I=1,MYIE),J=1,MYJE)
! Lyra GSM Wind stress
    END IF
!
    CALL COAL(  TH10(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(   Q10(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(   U10(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(   V10(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(TSHLTR(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(QSHLTR(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(PSHLTR(1:MYIE,1:MYJE), MYIE*MYJE)
!---------
! SM V100M
!---------
    CALL COAL(TH100(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL( Q100(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL( U100(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL( V100(1:MYIE,1:MYJE), MYIE*MYJE)
!
! Lyra GSM Wind stress
    CALL COAL(XMOMFLUX(1:MYIE,1:MYJE),MYIE*MYJE)
    CALL COAL(YMOMFLUX(1:MYIE,1:MYJE),MYIE*MYJE)
! Lyra GSM Wind stress
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) (((SMC(I,J,N),I=1,MYIE),J=1,MYJE),N=1,NSOIL)
    END IF
!
    CALL COAL(SMC(1:MYIE,1:MYJE,1:NSOIL), MYIE*MYJE*NSOIL)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((CMC(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL(CMC(1:MYIE,1:MYJE), MYIE*MYJE)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) (((STC(I,J,N),I=1,MYIE),J=1,MYJE),N=1,NSOIL)
    END IF
    CALL COAL(STC(1:MYIE,1:MYJE,1:NSOIL), MYIE*MYJE*NSOIL)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) (((SH2O(I,J,N),I=1,MYIE),J=1,MYJE),N=1,NSOIL)
    END IF
!
    CALL COAL(SH2O(1:MYIE,1:MYJE,1:NSOIL), MYIE*MYJE*NSOIL)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((ALBEDO(I,J),I=1,MYIE),J=1,MYJE)
    END IF
!
    CALL COAL(ALBEDO(1:MYIE,1:MYJE), MYIE*MYJE)
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((POTFLX(I,J),I=1,MYIE),J=1,MYJE), ((TLMIN(I,J),I=1,MYIE),J=1,MYJE),        &
    &                  ((TLMAX(I,J),I=1,MYIE),J=1,MYJE),                                          &
!Lyra GSM Max Wind
    &                  ((MAXWU(I,J),I=1,MYIE),J=1,MYJE),                                          &
    &                  ((MAXWV(I,J),I=1,MYIE),J=1,MYJE),                                          &
!Lyra GSM Max Wind
    &                 ACUTIM, ARATIM, APHTIM, NHEAT, NPHS, NCNVC, NPREC, NRDSW, NRDLW, NSRFC,     &
    &                 TPH0D, TLM0D, RESTRT
    END IF
!
    CALL COAL(POTFLX(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL( TLMIN(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL( TLMAX(1:MYIE,1:MYJE), MYIE*MYJE)
!Lyra GSM Max Wind
    CALL COAL(MAXWU(1:MYIE,1:MYJE),MYIE*MYJE)
    CALL COAL(MAXWV(1:MYIE,1:MYJE),MYIE*MYJE)
!Lyra GSM Max Wind
    CALL COAL(ACUTIM,1)
    CALL COAL(ARATIM,1)
    CALL COAL(APHTIM,1)
    CALL COAL(NHEAT ,1)
    CALL COAL(NPHS  ,1)
    CALL COAL(NCNVC ,1)
    CALL COAL(NPREC ,1)
    CALL COAL(NRDSW ,1)
    CALL COAL(NRDLW ,1)
    CALL COAL(NSRFC ,1)
    CALL COAL(TPH0D ,1)
    CALL COAL(TLM0D ,1)
    CALL COAL(RESTRT,1)
!
    DO K=1,LM
        IF (MULTIWRITE) THEN
            WRITE(LRSTRT) ((RSWTT(I,J,K),I=1,MYIE),J=1,MYJE)
            WRITE(LRSTRT) ((RLWTT(I,J,K),I=1,MYIE),J=1,MYJE)
        END IF
!
        CALL COAL(RSWTT(1:MYIE,1:MYJE,K), MYIE*MYJE)
        CALL COAL(RLWTT(1:MYIE,1:MYJE,K), MYIE*MYJE)
    END DO
!
    IF (MULTIWRITE) THEN
        WRITE(LRSTRT) ((CNVBOT(I,J),I=1,MYIE),J=1,MYJE)
        WRITE(LRSTRT) ((CNVTOP(I,J),I=1,MYIE),J=1,MYJE)
        WRITE(LRSTRT) ((RSWTOA(I,J),I=1,MYIE),J=1,MYJE)
        WRITE(LRSTRT) ((RLWTOA(I,J),I=1,MYIE),J=1,MYJE)
        CLOSE(LRSTRT)
    END IF
!
    CALL COAL(CNVBOT(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(CNVTOP(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(RSWTOA(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(RLWTOA(1:MYIE,1:MYJE), MYIE*MYJE)
!
    CALL COAL(  HBM2(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(    SM(1:MYIE,1:MYJE), MYIE*MYJE)
    CALL COAL(   SPL(1:LSL),LSL)
    CALL COAL(  DETA(1:LM),LM)
    CALL COAL(PT,1)
    CALL COAL(SPLINE,1)
!------------------------------------------------------------------------------------------------
! AT THIS POINT WE HAVE ACCUMULATED ALL OF THE DATA INTO BUF
! WE WANT TO KNOW THE MAXIMUM AMOUNT ACROSS ALL MPI TASKS
! THIS IS USEFUL IN CASE WE DECIDE TO WRITE A FILE INSTEAD OF SENDING THE DATA TO THE I/O SERVERS
!------------------------------------------------------------------------------------------------
    CALL MPI_ALLREDUCE(IP, IPMAX, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_COMP, IERR)
!-----------------------------------------------------------------
! IPMAX IS THE MAXIMUM NUMBER OF 4 BYTE REALS ACROSS THE MPI TASKS
! LETS COMPUTE A RECLEN THAT IS A MULTIPLE OF 2**18 BYTES
! WE WILL USE THIS WHEN OPENING THE DIRECT ACCESS FILE
!-----------------------------------------------------------------
    IBLOCK = ((IPMAX * 4) / (2 ** 18)) + 1
    IRECLEN = IBLOCK * (2 ** 18)
!------------------------------------------------------
! WE WILL PLACE THE RECLEN IN THE BEGINNING OF THE FILE
! THIS IS HANDY
!------------------------------------------------------
    CALL REPLACE(IRECLEN, 1, 1)
!------------------------------------------------------------------------
! IF WE HAVE ANY I/O SERVERS WE WILL SEND THE DATA TO THEM FOR PROCESSING
!------------------------------------------------------------------------
    IF (IQUILT_GROUP > 0) THEN
    
        IF (MYPE == 0) THEN
            CALL MPI_SEND(ITAG, 1, MPI_INTEGER, 0, 0, MPI_COMM_INTER_ARRAY(ISERVE), IERR)
        END IF
!    
        DO I=0,INUMQ(ISERVE) -1
            CALL PARA_RANGE(0, JNPES-1, INUMQ(ISERVE), I, ISTART, IEND)
            MYPE_ROW = MYPE / INPES
!        
            IF (MYPE_ROW >= ISTART .AND. MYPE_ROW <= IEND ) THEN
               IF (MYPE == 0) WRITE(0,*)  'CALL MPI_ISEND.... ', IP, ITAG
                CALL MPI_ISEND(BUF, IP, MPI_REAL, I, ITAG, MPI_COMM_INTER_ARRAY(ISERVE), IHS, IERR)
            END IF
!        
        END DO
!--------------------------------------------------------------------------------------------------    
! IN CASE WE HAVE MULTIPLE GROUPS OF I/O SERVERS, INCREMENT TO THE NEXT SERVER FOR THE NEXT OUTPUT
! TIME
!--------------------------------------------------------------------------------------------------    
        ISERVE = ISERVE + 1
        IF (ISERVE > IQUILT_GROUP) ISERVE = 1
!---------------------------------------------------------
! APPARENTLY, WE HAVE CHOSEN NOT TO SUPPLY ANY I/O SERVERS
! WE WILL WRITE A DIRECT ACCESS FILE INSTEAD
!---------------------------------------------------------   
    ELSE  
        OPEN(UNIT=LRSTRT,FILE=RSTFIL,FORM='UNFORMATTED',IOSTAT=IER,ACCESS='DIRECT',RECL=IRECLEN)
!
        IF (IER /= 0) WRITE(LIST,*)' LRSTRT OPEN UNIT ERROR IER=', IER
!    
        WRITE(LRSTRT,REC=MYPE+1) (BUF(I),I=1,IP)
        CLOSE(LRSTRT)
!    
    END IF
!
    DIF_TIM = TIMEF() - BTIM0
    WRT_TIM = WRT_TIM + DIF_TIM
!
    CALL MPI_REDUCE(WRT_TIM, WRT_TIM_0, 1, MPI_REAL, MPI_MAX, 0, MPI_COMM_COMP, IERR)
!
    IF (MYPE == 0) THEN
        WRITE(6,*)' SHIPPED OR WROTE DATA, TIME = ', WRT_TIM_0*1.E-03
    END IF
!
    CALL MPI_BARRIER(MPI_COMM_COMP, ISTAT)
!-------------------------------------------------
! SEND SIGNAL THAT ALL TASKS HAVE FINISHED WRITING
!-------------------------------------------------
!
!----------------------------------------------------------------------------
! SM - ALTERACAO FEITA PARA O IO - PARA NAO MAIS ESCREVER O FCSTDONE E OUTJOB
!----------------------------------------------------------------------------
    OUTJ = .FALSE. 
    IF (OUTJ) THEN
        IF (IQUILT_GROUP == 0) THEN
            IF (MYPE == 0) THEN
                DONE = 'DONE'
!
                REWIND(11)
                READ(11,POSTLIST)                        !CHOU
                WRITE(LIST,POSTLIST)                     !CHOU
!                
                WRITE(OUTJOB,1240) ITAG
                1240 FORMAT('OUTJOB_SPECIAL',I6.6,'.ksh')
!
                LUNIN = 67
                LUNOT = 68
                OPEN(LUNIN, FILE='OUTJOB_SPECIAL.ksh')
                OPEN(LUNOT, FILE=OUTJOB)
                1260 FORMAT(A85)
!
                REWIND(LUNIN)
                DO I=1,11
                    READ(LUNIN,1260) LINE
                    IDX = INDEX(LINE,' ') - 1
                    WRITE(LUNOT,1260) LINE
                END DO
                1261 FORMAT('HFCT=',I6.6)
                WRITE(LUNOT,1261) ITAG
!
                1250 READ(LUNIN,1260,END=1290) LINE
                WRITE(LUNOT,1260) LINE
                GO TO 1250
                1290 CONTINUE
!
                CLOSE(LUNOT)
                CLOSE(LUNIN)
!                        
                IDX    = INDEX(SUBMIT,'#') - 1
                SUBMIT = SUBMIT(1:IDX) // OUTJOB
!                         
                CHX  = INDEX(CHMO,'#') - 1
                CHMO = CHMO(1:CHX) // OUTJOB
!
                CALL SYSTEM(CHMO)
!                         
                WRITE(LIST,*) 'CHKOUT:  SUBMIT POST JOB ', SUBMIT
                CALL SYSTEM(SUBMIT)    ! NEC EQUIVALENT
!
                IF (IER /= 0) WRITE(LIST,*)' SIGNAL SENT TO FINFIL:  DONE'
            END IF
        END IF
    END IF
!------------------------------------
! RESET ACCUMULATION COUNTERS TO ZERO
!------------------------------------
    APHTIM = 0.
    ACUTIM = 0.
    ARATIM = 0.
!---------------------------------------
! POST-POSTING UPDATING AND INITIALIZING
!---------------------------------------
!
!----------------------------------------------------------------------------------------------------------------------------
!IF (NTSD.EQ.NSHDE), THEN THIS WAS ALSO A FORECAST OUTPUT TIME.  WE NEED TO INCREMENT NSHDE FOR THE NEXT FORECAST OUTPUT TIME
!----------------------------------------------------------------------------------------------------------------------------
    IF (NTSD == NSHDE .OR. NSTART+1 == NSHDE) THEN
        IOUT = IOUT + 1
        IF (.NOT. RESTRT)   GOTO 1300
        IF (NTSD == NSHDE .OR. NSTART+1 == NSHDE) GOTO 1300
        IOUT = IOUT - 1
        1300 NSHDE = ISHDE(IOUT)
    END IF
!---------------------------
! ZERO ACCUMULATOR ARRAYS
! AVERAGE CLOUD AMOUNT ARRAY
!---------------------------
!    1310 CONTINUE moved to after CHKOUT: ZERO
!
    IF (MOD(NTSD,NCLOD) < NPHS) THEN
        IF (MYPE == 0) WRITE(LIST,*)'CHKOUT: ZERO AVG CLD AMT ARRAY'
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                ACFRCV(I,J) = 0.
                NCFRCV(I,J) = 0
                ACFRST(I,J) = 0.
                NCFRST(I,J) = 0
            END DO
        END DO
    END IF
!---------------------------------------------
! TOTAL AND CONVECTIVE PRECIPITATION ARRAYS
! TOTAL SNOW AND SNOW MELT ARRAYS
! STORM SURFACE AND BASE GROUND RUN OFF ARRAYS
! PRECIPITATION TYPE ARRAY
!---------------------------------------------
    IF (MOD(NTSD,NPREC) < NCNVC) THEN
        IF(MYPE == 0) WRITE(LIST,*) 'CHKOUT: ZERO ACCUM PRECIP ARRAYS'
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                ACPREC(I,J) = 0.
                CUPREC(I,J) = 0.
                ACSNOW(I,J) = 0.
                ACSNOM(I,J) = 0.
                SSROFF(I,J) = 0.
                BGROFF(I,J) = 0.
                SFCEVP(I,J) = 0.
                POTEVP(I,J) = 0.
            END DO
        END DO
    END IF
!--------------------------------------------------
! GRID-SCALE AND CONVECTIVE (LATENT) HEATING ARRAYS
!--------------------------------------------------
    IF (MOD(NTSD,NHEAT) < NCNVC) THEN
        IF (MYPE == 0) WRITE(LIST,*) 'CHKOUT: ZERO ACCUM LATENT HEATING ARRAYS'
        AVRAIN = 0.
        AVCNVC = 0.
        DO K=1,LM
            DO J=MYJS,MYJE
                DO I=MYIS,MYIE
                    TRAIN(I,J,K) = 0.
                    TCUCN(I,J,K) = 0.
                END DO
            END DO
        END DO
    END IF
!----------------------------------------------------------------------------------------------------
! BEGIN: INITIALIZE CONVECTIVE CLOUD FIELDS FOR RADIATION BEFORE RETURNING TO EBU (FERRIER 23 JAN 02)
!----------------------------------------------------------------------------------------------------
    IF (CURAD) THEN
        IF (MYPE == 0) THEN
            WRITE(0,"(a)") 'CHKOUT: INITIALIZE CUPPT,HTOP,HBOT'
            WRITE(6,"(a)") 'CHKOUT: INITIALIZE CUPPT,HTOP,HBOT'
        END IF
!
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                CUPPT(I,J) =   0.
                HTOP(I,J)  = 100.
                HBOT(I,J)  =   0.
            END DO
        END DO
!
        CURAD = .FALSE. 
!
    END IF
!-----------------------------------------------------------------------------------------------------------------------
! RESET CONVECTIVE CLOUD TOP AND BOTTOM ARRAYS (DIAGNOSTIC ONLY; THESE FIELDS ARE NOT CYCLED AND NOT READ WHEN TSTART=0)
!-----------------------------------------------------------------------------------------------------------------------
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            CNVTOP(I,J) = 100.
            CNVBOT(I,J) =   0.
        END DO
    END DO
!----
! END
!----
!
!---------------------------
! LONG WAVE RADIATION ARRAYS
!---------------------------
    IF (MOD(NTSD,NRDLW) < NPHS) THEN
        IF (MYPE == 0) WRITE(LIST,*) 'CHKOUT: ZERO ACCUM LW RADTN ARRAYS'
        ARDLW = 0.
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                 ALWIN(I,J) = 0.
                ALWOUT(I,J) = 0.
                ALWTOA(I,J) = 0.
            END DO
        END DO
    END IF
!----------------------------
! SHORT WAVE RADIATION ARRAYS
!----------------------------
    IF (MOD(NTSD,NRDSW) < NPHS) THEN
        IF (MYPE == 0) WRITE(LIST,*) 'CHKOUT:  ZERO ACCUM SW RADTN ARRAYS'
        ARDSW = 0.
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                 ASWIN(I,J) = 0.
                ASWOUT(I,J) = 0.
                ASWTOA(I,J) = 0.
            END DO
        END DO
    END IF
!---------------------------------------------
! SURFACE SENSIBLE AND LATENT HEAT FLUX ARRAYS
!---------------------------------------------
    IF (MOD(NTSD,NSRFC) < NPHS) THEN
        IF (MYPE == 0) WRITE(LIST,*) 'CHKOUT:  ZERO ACCUM SFC FLUX ARRAYS'
        ASRFC = 0.
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                SFCSHX(I,J) = 0.
                SFCLHX(I,J) = 0.
                SUBSHX(I,J) = 0.
                SNOPCX(I,J) = 0.
                SFCUVX(I,J) = 0.
                POTFLX(I,J) = 0.
            END DO
        END DO
    END IF
!    
    1310 CONTINUE
!    
!-------------------------------------
! RESET THE MAX/MIN TEMPERATURE ARRAYS
!-------------------------------------
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            TLMIN(I,J) =  999.
            TLMAX(I,J) = -999.
!***  RESET THE MAX WIND ARRAYS
!Lyra GSM Max Wind
            MAXWIND(I,J)=0
            WIND=0
!Lyra GSM Max Wind
        END DO
    END DO
!----------------------------------------------------------------
! WRITING IN A SEPARATE FILE FOR EACH MYPE AND EACH OUTPUT PERIOD
!----------------------------------------------------------------
!    IOUT = IOUT -1
!    NSHDE = ISHDE(IOUT)
!    IF (NTSD == NSHDE .OR. NSTART+1 == NSHDE) THEN
!        CALL OUT2RESTRT(ITAG)
!    ENDIF
!    IOUT = IOUT +1
!    NSHDE = ISHDE(IOUT)    
!
    IF (MYPE == 0) THEN
        WRITE(0,"(a)") 'FINISHED CHKOUT'
        WRITE(6,"(a)") 'FINISHED CHKOUT'
    END IF
!---------------
! END OF ROUTINE
!---------------
    RETURN
!
    END SUBROUTINE CHKOUT
!
!
!
    SUBROUTINE COAL(A, LEN)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE COAL
!> 
!> SUBPROGRAM: COAL - ?????
!> PROGRAMMER: ?????          
!> ORG: W/NP22 
!> DATE: ??-??-??
!>     
!> ABSTRACT:
!> ?????
!>     
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????      - ORIGINATOR
!> 06-07-16  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY 
!>
!> INPUT  ARGUMENT LIST:
!> A   -
!> LEN -
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>   
!> CALLS : BUFFER
!>         F77KINDS
!>
!> DRIVER: CHKOUT
!>--------------------------------------------------------------------------------------------------
    USE BUFFER
    USE F77KINDS
!
    IMPLICIT NONE
!
    INCLUDE 'mpif.h'
!
    REAL   (KIND=R4KIND), DIMENSION(*)                                    , INTENT(IN)          ::&
    & A
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LEN       
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IERR    , I
!
    IF (LEN < 0) THEN
        IP = 0
    END IF
!
    IF (IP + LEN > IBUFMAX) THEN
        PRINT *, ' IBUFMAX IN MODULE_BUFFER.f90 IS TOO SMALL, STOPPING'
        PRINT *, ' CHANGE IBUFMAX IN PARMBUF.f90 AND RECOMPILE'
        PRINT *, 'IP + LEN > IBUFMAX',IP, LEN, IBUFMAX
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    END IF
!
    DO I=1,ABS(LEN)
        IP = IP + 1
        BUF(IP) = A(I)
    END DO
!
    END SUBROUTINE COAL
!
!
!
    SUBROUTINE REPLACE(A, LEN, IW)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE REPLACE
!> 
!> SUBPROGRAM: REPLACE - ?????
!> PROGRAMMER: ?????          
!> ORG: W/NP22 
!> DATE: ??-??-??
!>     
!> ABSTRACT:
!> ?????
!>     
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????      - ORIGINATOR
!> 06-07-16  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY 
!>
!> INPUT  ARGUMENT LIST:
!> A   -
!> LEN -
!> IW  -
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>   
!> CALLS : BUFFER
!>         F77KINDS
!>
!> DRIVER: CHKOUT
!>--------------------------------------------------------------------------------------------------
    USE BUFFER
    USE F77KINDS
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND), DIMENSION(*)                                    , INTENT(IN)          ::&
    & A
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LEN     , IW        
!
    INTEGER(KIND=I4KIND)                                                                        ::&  
    & I       , IPP 
!
    IPP = IW
!
    DO I=1,LEN
        BUF(IPP) = A(I)
        IPP = IPP + 1
    END DO
!
    END SUBROUTINE REPLACE          
