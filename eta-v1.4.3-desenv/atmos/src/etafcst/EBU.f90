    PROGRAM EBU
!>--------------------------------------------------------------------------------------------------
!> MAIN PROGRAM EBU
!>
!> MAIN PROGRAM: ETAFCST - EARLY ETA MODEL FORECAST DRIVER
!> PROGRAMMER: JANJIC
!> ORG: NP22
!> DATE: 99-01-20
!>
!> ABSTRACT: 
!> EBU CONTAINS THE PRIMARY RUNSTREAM FOR THE EARLY ETA FORECAST MODEL.  
!> AFTER AN INITIAL CALL TO SUBROUTINE INIT, CALLS ARE MADE TO SUBROUTINES WHICH COMPUTE THE VARIOUS
!> DYNAMICAL AND PHYSICAL PROCESSES IN THE MODEL. 
!> THE VARIABLE 'NTSD' IS THE FUNDAMENTAL TIMESTEP COUNTER AND THUS ITS VALUE DETERMINES WHEN THE 
!> SUBROUTINES ARE CALLED.
!> INFORMATION PERTAINING TO THE SCHEMES USED IN THE MODEL AS WELL AS ADDITIONAL REFERENCES MAY BE
!> FOUND IN "THE STEP-MOUNTAIN ETA  COORDINATE REGIONAL MODEL:  A DOCUMENTATION" (BLACK 1988; 
!> DEVELOPMENT DIVISION) AND "THE NEW NMC MESO-SCALE ETA MODEL: DESCRIPTION AND FORECAST EXAMPLES 
!> (BLACK 1994; WEATHER AND FORECASTING).
!>
!> PROGRAM HISTORY LOG:
!> 87-08-??  JANJIC, BLACK - ORIGINATOR          
!> 93-05-12  TREADON       - DOCBLOCK INSERTED
!> 93-10-25  BLACK         - DOCBLOCK UPDATED
!> 97-03-15  MESINGER      - SPLITTING MODIFIED, TO SEPARATE THE ADJUSTMENT AND THE ADVECTION STEP                        
!> 97-11-19  BLACK         - MODIFIED FOR DISTRIBUTED MEMORY
!> 98-10-20  BLACK         - DISTRIBUTED MEMORY FORM FOR CURRENT OPERATIONAL CODE
!> 00-02-25  TUCCILLO      - INCORPORATED ASYNCHRONOUS I/O SERVERS
!> 00-11-14  BLACK         - INCORPORATED JANJIC NONHYDROSTATIC OPTION
!> 00-11-27  BLACK         - INCORPORATED RESTART CAPABILITY
!> 18-01-15  LUCCI         - MODERNIZATION OF THE CODE, INCLUDING:
!>                           * F77 TO F90/F95
!>                           * INDENTATION & UNIFORMIZATION CODE
!>                           * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                           * DOCUMENTATION WITH DOXYGEN
!>                           * OPENMP FUNCTIONALITY
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
!> USE MODULES: BOCOH
!>              BOCOV
!>              CHKOUT
!>              CLTEND
!>              CUCNVC
!>              DDAMP
!>              DIVHOA
!>              DIVHOAST
!>              EXCH
!>              GOSSIP
!>              GSCOND
!>              GSMDRIVE
!>              HDIFF
!>              HDIFFS
!>              HZADV2
!>              INIT
!>              INITS
!>              MPPINIT
!>              PDNEW
!>              PDTEDT
!>              PGCOR
!>              PRECPD
!>              QUILT
!>              RADTN
!>              RDTEMP
!>              TURBL
!>              VADZ
!>              VTADV
!>              W3TAGB
!>              W3TAGE
!> 
!> DRIVER     : -----
!>
!> CALLS      : BOCOH
!>              BOCOV
!>              CHKOUT
!>              CLTEND
!>              CUCNVC
!>              DDAMP
!>              DIVHOA
!>              DIVHOAST
!>              EPS
!>              EXCH
!>              GETENV
!>              GOSSIP
!>              GSCOND
!>              GSMDRIVE
!>              HDIFF
!>              HDIFFS
!>              HZADV
!>              HZADVS
!>              HZADV2
!>              INIT
!>              INITS
!>              MPI_BARRIER
!>              MPI_FINALIZE
!>              MPI_SEND
!>              MPPINIT
!>              PDNEW 
!>              PDTEDT
!>              PGCOR 
!>              PRECPD
!>              QUILT
!>              RADTN
!>              RDTEMP
!>              SETUP_SERVERS
!>              TURBL
!>              VADZ
!>              VTADV
!>              W3TAGB
!>              W3TAGE
!>            
!>
!> UNIQUE:
!> INIT     - INITIALIZE VARIABLES AT START OF INTEGRATION
!> DIVHOA   - DIVERGENCE, AND HORIZONTAL PART OF THE OMEGA-ALPHA TERM                
!> PGCOR    - PRESSURE GRADIENT AND CORIOLIS FORCE
!> PDTE     - UPDATE SURFACE PRESSURE TENDENCY AND ETADOT
!> VTADV    - VERTICAL ADVECTION
!> HZADV    - HORIZONTAL ADVECTION OF T, U, V, AND TKE
!> HZADV2   - HORIZONTAL ADVECTION OF Q AND CLOUD WATER
!> DDAMP    - APPLY DIVERGENCE DAMPING
!> PDNEW    - UPDATE SURFACE PRESSURE
!> HDIFF    - LATERAL DIFFUSION
!> BOCOH    - UPDATE H POINTS ON THE BOUNDARIES
!> BOCOV    - UPDATE V POINTS ON THE BOUNDARIES
!> RADTN    - RADIATION DRIVER
!> RDTEMP   - APPLY TEMPERATURE TENDENCY DUE TO RADIATION
!> TURBL    - PERFORM THE VERTICAL TURBULENT EXCHANGE
!> SURFACE  - UPDATE SURFACE TEMPERATURE, MOISTURE, AND OTHER GROUND HYDROLOGY            
!> GSCOND   - CLOUD WATER/ICE PHYSICS PARAMETERIZATION (EDAS ONLY)
!> CUCNVC   - CONVECTIVE ADJUSTMENT FOR DEEP OR SHALLOW CONVECTION
!> PRECPD   - GRID SCALE PRECIPITATION (EDAS ONLY)
!> GSMDRIVE - GRID SCALE MICROPHYSICS DRIVER (FREE FORECAST ONLY)
!> ADJPPT   - ADJUST MODEL PRECIPITATION TO OBSERVATIONS
!> CLTEND   - UPDATE TEMPERATURE FROM (GRID- AND SUBGRID-SCALE) CLOUD PROCESSES
!> CHKOUT   - POST PROFILE DATA.  FOR INTERNAL POST, POSTS MODEL OUTPUT.  
!>            FOR EXTERNAL POST, WRITES TEMPORARY FILE CONTAINING ALL MODEL ARRAYS.
!>
!> EXIT STATES:
!> COND =   1 - NORMAL EXIT
!>--------------------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
!         LIMITED AREA ETA MODEL WITH STEP-MOUNTAIN TOPOGRAPHY         
!
! NOAA / NATIONAL CENTERS FOR ENVIRONMENTAL PREDICTION, CAMP SPRINGS, MD
!
!      GEOPHYSICAL FLUID DYNAMICS LABORATORY / NOAA, PRINCETON, NJ     
!
!     UNIVERSITY CORPORATION FOR ATMOSPHERIC RESEARCH, BOULDER, CO.
!                                 AND                                        
!     DEPARTMENT OF METEOROLOGY, UNIVERSITY OF BELGRADE, YUGOSLAVIA 
!----------------------------------------------------------------------- 
!
!--------------------------------------------------------------------------------------------------
!
! REFERENCES FOR THE DYNAMICAL PART OF THE MODEL
!
! STEP-MOUNTAIN ETA COORDINATE: 
! F. MESINGER, 1983, IN RES. ACTIVITIES IN ATMOS. AND OCEANIC MODELING, REP. NO. 5, 
! WMO, GENEVA, 4.9-4.10.
!
! HORIZONTAL ADVECTION, CONTINUITY EQUATION: 
! Z.I. JANJIC, 1984, MWR, 112, NO.6, 1234-1245.
! 
! INTERNAL BOUNDARIES, OMEGA-ALPHA TERM, CODING, PERFORMANCE: 
! MESINGER ET AL., 1988, MWR, 116 NO.7, 1493-1518.
!
! N.B. FOR MORE DETAILS ON THESE TOPICS SEE ALSO: 
!
! 1.  MESINGER, F., AND Z.I. JANJIC, 1985: PROBLEMS AND NUMERICAL METHODS OF THE INCORPORATION OF 
!     MOUNTAINS IN ATMOSPHERIC MODELS.
!     LECTURES IN APPLIED MATHEMATICS, VOL 22, AMER. MATH. SOC.; ALSO, NUMERICAL METHODS FOR 
!     WEATHER PREDICTION, SEMINAR 1983, ECMWF, 103-157;
!     ALSO, SHORT AND MEDIUM-RANGE WEATHER PREDICTION RESEARCH PUBL. SER., NO. 8, 
!     WMO, GENEVA, 175-233.
!
! 2.  JANJIC, Z.I., AND F. MESINGER, 1983: FINITE-DIFFERENCE METHODS FOR THE SHALLOW WATER 
!     EQUATIONS ON VARIOUS HORIZONTAL GRIDS.
!     NUMERICAL METHODS FOR WEATHER PREDICTION, SEMINAR 1983, ECMWF,29-101.
!
! SOME  REFERENCES FOR THE PHYSICS PART OF THE MODEL
!
! JANJIC, Z.I., 1990: THE STEP-MOUNTAIN COORDINATE: 
! PHYSICAL PACKAGE.  MONTHLY WEATHER REVIEW, VOL. 118, NO. 7, 1429-1443.
! JANJIC, Z.I., 1994: THE STEP MOUNTAIN ETA COORDINATE:
! FURTHER DEVELOPMENTS OF THER CONVECTION, VISCOUS SUBLAYER, AND TURBULENCE CLOSURE SCHEMES. 
! MONTHLY WEATHER REVIEW,VOL. 122, 927-945.
!
! ALSO SEE REFERENCES IN PHYSICAL SUBROUTINES 
!
!--------------------------------------------------------------------------------------------------
!
! THIS VERSION OF THE PROGRAM IS WRITTEN IN STANDARD ANSI FORTRAN 90 
!
! PRINCIPAL PROGRAMMERS: 
!
! Z. JANJIC, UNIVERSITY OF BELGRADE, T. BLACK, NCEP 
!
! THE MODEL USES THE SEMI-STAGGERED E GRID IN ARAKAWA NOTATION.
! HORIZONTAL INDEXING IS TWO-DIMENSIONAL. 
!
! H(1,JM)  V(1,JM)  H(2,JM)  V(2,JM) ...... V(IM-1,JM)  H(IM,JM) 
!    .        .        .        .               .          .     
!    .        .        .        .               .          .     
!    .        .        .        .               .          .     
!    .        .        .        .               .          .     
!                                                                
! H(1,3)   V(1,3)   H(2,3)   V(2,3) ....... V(IM-1,3)   H(IM,3)
!
! V(1,2)   H(1,2)   V(2,2)   H(2,2) ....... H(IM-1,2)   V(IM,2) 
!                                                                 
! H(1,1)   V(1,1)   H(2,1)   V(2,1) ....... V(IM-1,1)   H(IM,1) 
!
!
!
! ARRAYS ARE DIMENSIONED (IM,JM). NOTE THAT A PHANTOM COLUMN OF POINTS MUST EXIST ALONG THE 
! EASTERN EDGE FOR THE ARRAYS TO BE COMPLETE.                                             
! 
! THE TOTAL NUMBER OF GRID POINTS IN THE HORIZONTAL EXCLUDING 
! THE PHANTOM COLUMN IS IMJM=IM*JM-JM/2.                      
! 
! AUXILIARY ARRAYS ARE USED TO LOCATE NEIGHBORING GRID POINTS WITH RESPECT TO A GIVEN GRID
! POINT. IHE(J) IS THE INCREMENT TO THE I INDEX NEEDED TO REFER TO THE V POINT EAST OF AN
! H POINT THUS IHE(J)=0 ON ODD ROWS AND =1 ON EVEN ROWS. 
! IHW(J)=IHE(J)-1 IS THE INCREMENT TO THE INDEX OF AN H POINT  TO REFER TO THE V POINT TO THE
! WEST OF THAT H POINT. THE ANALOG EXISTS FOR THE ARRAYS IVE(J) AND IVW(J).
!
! BOUNDARY MASKS AND TOPOGRAPHY MASKS ARE DEFINED FOR VECTOR PROCESSING. THE BOUNDARY MASKS 
! HBM2(K) AND VBM2(K) ARE EQUAL TO ONE EVERYWHERE EXCEPT AT THE TWO OUTERMOST ROWS WHERE THEY
! ARE EQUAL TO ZERO. THE BOUNDARY MASK VBM3(K) IS EQUAL TO ONE EVERYWHERE EXCEPT AT THE THREE 
! OUTERMOST ROWS WHERE IT IS EQUAL TO ZERO. THE TOPOGRAPHY MASKS (HTM(K,L) AND VTM(K,L)) ARE 
! SET TO ZERO UNDERNEATH THE TOPOGRAPHY AND TO ONE ELSWHERE. IN ADDITION, FOR TREATMENT OF 
! PHYSICAL PROCESSES, MAXIMUM VALUES OF THE VERTICAL INDEX ARE DEFINED AND STORED (LMH(K) AND 
! LMV(K).  
!
!--------------------------------------------------------------------------------------------------
! THE NUMBER OF QUILT SERVERS MUST AGREE WITH THE FOLLOWING RELATIONSHIP:
! 0 <=  NUMBER_QUILT_SERVERS <= JNPES 
! WHERE THE NUMBER_QUILT_SERVERS = ( NUMBER_OF MPI_TASKS - INPES*JNPES )
!
! PREFERABLY, THE NUMBER OF QUILT SERVERS DIVIDES EVENLY INTO JNPES
!
! JIM TUCCILLO - AUGUST 2000
!--------------------------------------------------------------------------------------------------
    USE CLDWTR
    USE CONTIN
    USE CTLBLK
    USE EXCHM
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE CMICRO_CONS    , ONLY : RHGRDL, RHGRDS,  VSNOWADST
    USE MPPCOM
    USE NHYDRO
    USE NSOILTYPE
    USE OUTFIL
    USE PARMETA
    USE PVRBLS
    USE TEMPCOM
    USE TIMMING
    USE TOPO
    USE VRBLS
    USE UPDT     
    USE UPDATE_FLDS
!
    IMPLICIT NONE
!
    INCLUDE "EXCHM.h"
!
    LOGICAL(KIND=L4KIND)                                                                          ::&
    & CLTEND_TEST
!
    INCLUDE "mpif.h"
!
#include "sp.h"
!
!  LOGICAL NAMELIST

    LOGICAL(KIND=L4KIND)                                                                          ::&
    &  CUCNVCFLG
!    
    LOGICAL(KIND=L4KIND)                                                                          ::&
    & LDDAMPFLG
!    
    LOGICAL(KIND=L4KIND)                                                                          ::&
    & LFCTTIMCHK
!
    LOGICAL(KIND=L4KIND)                                                                          ::&
    &  SHLCNVFLG
!
    LOGICAL(KIND=L4KIND)                                                                          ::&
    & SLOPE    
!
    LOGICAL(KIND=L4KIND)                                                                          ::&
    & WRITEOUT2RESTRT
!
    INTEGER(KIND=I4KIND)                                                                         ::&  
    & CUCNVCSQM
!
    INTEGER(KIND=I4KIND)                                                                         ::&  
    & FREQOUT2RESTRT
!
    INTEGER(KIND=I4KIND)                                                                         ::&  
    & MICROPHYSSQM   
!
    INTEGER(KIND=I4KIND)                                                                         ::&  
    & SHLCNVSQM
!    
   NAMELIST /ETAINFCTNML/ CUCNVCFLG,CUCNVCSQM,FREQOUT2RESTRT,LCO2,LDDAMPFLG,LFCTTIMCHK             &
    &                    ,LSST,LSSTMNTHLY,MICROPHYSSQM,NSOTYP,SHLCNVFLG                            &
    &                    ,SHLCNVSQM,SLOPE,RHGRDL,RHGRDS,WRITEOUT2RESTRT                            &
    &                    ,VSNOWADST                    
!
!------------------------------------------------ 
! THE FOLLOWING ARE USED FOR TIMIMG PURPOSES ONLY
!------------------------------------------------ 
    REAL   (KIND=R8KIND)                                                                        ::&
    & TIMEF
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    REAL   (KIND=R4KIND)                                                                        ::&
    &  BOCOH_TIM        ,  BOCOV_TIM        ,   CHKOUT_TIM      , CLTEND_TIM        ,             &
    & CUCNVC_TIM        ,  DDAMP_TIM        ,   DIVHOA_TIM      ,    EPS_TIM        ,             &
    &   GOSS_TIM        , GSCOND_TIM        , GSMDRIVE_TIM      ,   HADZ_TIM        ,             &
    &  HDIFF_TIM        ,  HZADV_TIM        ,   HZADV2_TIM      ,   INIT_TIM        ,             &
    &    MPP_TIM        ,  PDNEW_TIM        ,   PDTEDT_TIM      ,  PGCOR_TIM        ,             &
    & PPTADJ_TIM        , PRECPD_TIM        ,    RADTN_TIM      , RDTEMP_TIM        ,             &
    &  TURBL_TIM        ,   VADZ_TIM        ,    VTADV_TIM      ,   TOT2_TIM        ,             &
    &    TOT_TIM        , SHLCNV_TIM
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & BTIMX   , BTIM    , BBTIM   , PCT
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & TSTART
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I        , IER     , IERR              , ICLTEND           , K                 ,            & 
    & CALC_NVEG, CALC_SST, NREC_SST_UP       , NOut
!
    CHARACTER(LEN=8)                                                                            ::&
    & ENVAR   , SRFILE
!--------------------------------------------------------------------------------------------------
! INITIALIZE MPI, SETUP I/O SERVER MECHANICS AND CHECK FOR WHETHER A SUFFICIENT NUMBER OF MPI TASKS
! HAVE BEEN INITIATED.
! IF INSUFFICIENT MPI TASK HAVE BEEN INITIATED THE CODE WILL STOP IN SETUP_SERVERS
!--------------------------------------------------------------------------------------------------
    CALL SETUP_SERVERS(INPES*JNPES  , MYPE          , NPES                , IQUILT_GROUP, INUMQ,  &
    &                  MPI_COMM_COMP, MPI_COMM_INTER, MPI_COMM_INTER_ARRAY)
!
    IF (MYPE == 0) THEN
        CALL W3TAGB('ETAFCST ',0097,0365,0060,'NP22   ')
    END IF
!-------------------------------------------------------------------------------- 
! AT THIS POINT NPES IS THE NUMBER OF MPI TASKS WORKING ON THE MODEL INTEGRATION.
! ALL OTHER TASKS ARE I/O SERVERS. AND AWAY WE GO !
!-------------------------------------------------------------------------------- 
    IF (MYPE >= NPES) THEN
!------------------------   
! FIRE UP THE I/O SERVERS
!------------------------   
        CALL QUILT
!    
    ELSE
!-------------------------------------------------- 
! THESE ARE THE TASKS THAT DO THE MODEL INTEGRATION
!-------------------------------------------------- 
           BOCOH_TIM = 0. 
           BOCOV_TIM = 0.
          CHKOUT_TIM = 0.        
          CLTEND_TIM = 0.        
          CUCNVC_TIM = 0.       
           DDAMP_TIM = 0.
          DIVHOA_TIM = 0.       
             EPS_TIM = 0. 
            EXCH_TIM = 0.
            GOSS_TIM = 0.
          GSCOND_TIM = 0.
        GSMDRIVE_TIM = 0.
            HADZ_TIM = 0.
           HDIFF_TIM = 0.
           HZADV_TIM = 0.
          HZADV2_TIM = 0.
            INIT_TIM = 0.
             MPP_TIM = 0.
             NHB_TIM = 0.
           PDNEW_TIM = 0.
          PDTEDT_TIM = 0.
           PGCOR_TIM = 0.
          PPTADJ_TIM = 0.
          PRECPD_TIM = 0.
           RADTN_TIM = 0.
          RDTEMP_TIM = 0.
             RES_TIM = 0.
          SHLCNV_TIM = 0.
          SURFCE_TIM = 0.
           TURBL_TIM = 0.
            VADZ_TIM = 0.
           VTADV_TIM = 0.
!------------------------------------------------------------- 
! INITIALIZE ALL QUANTITIES ASSOCIATED WITH GRID DECOMPOSITION
!------------------------------------------------------------- 
        BTIMX = TIMEF()
        BTIM  = TIMEF()
!
        CALL MPPINIT
!
        MPP_TIM = MPP_TIM + TIMEF() - BTIM
!----------------------------------------------- 
! INITIALIZE CONSTANTS AND VARIABLES 
! DISTRIBUTE THE VALUES TO THE VARIOUS NODES/PES
!----------------------------------------------- 
        IF (MYPE == 0) THEN
            WRITE(6,*) 'READING ETAINFCTNML'
        END IF
!
        OPEN(UNIT=9, FILE="ETAIN_FCT.nml", STATUS="OLD")
        READ(9,ETAINFCTNML)
        CLOSE(9)
!
        IF (MYPE == 0) THEN
            WRITE(6,*) 'NSOTY: ',NSOTYP
!
            WRITE(6,*) 'SLOPE: ',SLOPE
!
            IF (CUCNVCFLG) THEN
                 SELECT CASE (CUCNVCSQM)
                     CASE (1) 
                       WRITE(6,*) 'USING: Betts-Miller Convective Scheme'
                     CASE (2) 
                       WRITE(6,*) 'USING: Kain-Fritsch Convective Scheme'
                END SELECT
            ELSE
                WRITE(6,*) 'USING: No Convective SchemeBetts-Miller'
            END IF
!
            SELECT CASE (MICROPHYSSQM)
                CASE (1) 
                  WRITE(6,*) 'USING: Ferrier Microphysics Scheme'
                CASE (2) 
                  WRITE(6,*) 'USING: Zhao Microphysics Scheme'
            END SELECT
!
            IF (LSST) THEN
	        WRITE(6,*) 'SST UPDATE: ',LSST
                IF (LSSTMNTHLY) THEN
                    WRITE(6,*) 'SST INTERPOLATED FROM MONTHLY DATA'
                ELSE
                    WRITE(6,*) 'SST UPDATED FROM DAYLY DATA ARQUIVES'
                ENDIF
             END IF
!	    
            WRITE(6,*) 'CO2 UPDATE: ',LCO2
        ENDIF
!
        BBTIM = TIMEF()
!
        IF (SLOPE) THEN
            CALL INITS(TSTART)
        ELSE
            CALL INIT
        END IF
!
        INIT_TIM = TIMEF() - BBTIM
!    
        BTIM = TIMEF()
!
        CALC_SST= 86400 / DT          

        IF (LSST) THEN
          IF ((NTSD == 0) .AND. (.NOT.RESTRT)) THEN
          ELSEIF (LSSTMNTHLY) THEN
              CALL FLDS_UPDATE_DRIVER('UPDATE','SSTM2D')
          ELSE
              CALL FLDS_UPDATE_DRIVER('UPDATE','SST',CALC_SST,RECSST_REF=1)
          END IF
        END IF
        CALL GOSSIP
!
        GOSS_TIM = GOSS_TIM + TIMEF() - BTIM
!------------------------------------------- 
! INVOKE THE LYNCH DIGITAL FILTER IF DESIRED 
!------------------------------------------- 
!        DO NFLT=1,3
!           IF (.NOT. NEST .AND .NFLT .GT. 1 .AND. MYPE .EQ. 0) THEN
!              REWIND NBC
!              READ(NBC)
!              READ(NBC)BCHR
!    
!           END IF
!        CALL DIGFLT
!    
!       END DO
!
        CALL GETENV("TMMARK", ENVAR)
!
        IF (MYPE == 0) WRITE(0,*) "EBU FINDS THAT TMMARK =", ENVAR
!
        IF (ENVAR /= 'TM00') THEN
!
         10 FORMAT('SR.',A4)
!
        END IF
!
        CLTEND_TEST = .TRUE. 
!
        IF (NPHS /= NCNVC) THEN
            CLTEND_TEST = .FALSE. 
            IF (ENVAR /= 'TM00') THEN
                WRITE(0,"(A)") 'WARNING: RESULTS COULD BE IN ERROR !'
            END IF
        END IF
!-----------------------------------------
!READ OUT2RESTRT FILES WHEN RESTRT IS TRUE
!-----------------------------------------
        IF (RESTRT) THEN
	    PRINT*, "Entered READ_OUT2RESTRT, MYPE: ", MYPE
            CALL READ_OUT2RESTRT(TSTART)
	ENDIF
!------------------------------------------------------- 
! SPECIAL CONSIDERATION WHEN NTSD=0 AT START OF FORECAST
!------------------------------------------------------- 
        IF (NTSD == 0) NTSD = 1
!-------------------------------------------------------- 
! CALLED AT BEGINNING EVERY TIME IN ORDER TO TEST CYCLING
!-------------------------------------------------------- 
        IF (NTSD .EQ. 1) THEN
            BTIM = TIMEF()
!
            CALL RADTN
!
            RADTN_TIM = RADTN_TIM + TIMEF() - BTIM
        END IF
!------------------------ 
! GENERATE INITIAL OUTPUT 
!------------------------ 
        BTIM = TIMEF()
!
            CALL CHKOUT
!
        CHKOUT_TIM = CHKOUT_TIM   + TIMEF() - BTIM
!
        IF (NTSD == 1) NTSD = 0
!------------------------- 
! ENTRY INTO THE TIME LOOP 
!------------------------- 
    2000 CONTINUE
        NTSD = NTSD + 1
!
        IF ((MYPE == 0) .AND. (LFCTTIMCHK))                                                          &
    &        WRITE(0,2001) NTSD, (NTSD-1) * DT, (NTSD-1) * DT / 3600
   2001 FORMAT('EBU:  TIMESTEP NTSD=',I7,'  FCST TIME=',F12.0,' S',' AND ',F8.3,' H')
! GSM XC
! UPDATE  SST DATA 
!---------------------------------------
!
        IF (LSST) THEN
           IF ((NTSD == 1) .OR. (RESTRT) .OR. (MOD(NTSD,CALC_SST) == 0)) THEN
              CALC_SST = 86400 / DT
              IF (LSSTMNTHLY) THEN
                 CALL FLDS_UPDATE_DRIVER('UPDATE','SSTM2D')
              ELSE
                 CALL FLDS_UPDATE_DRIVER('UPDATE','SST',CALC_SST,RECSST_REF=1)
              END IF
            END IF
        END IF

!----------------------------
! CHOU 02-02-2008
! UPDATE MONTHLY VEG GREENESS
!----------------------------
!
        CALC_NVEG = 86400 / DT
!
        IF ((NTSD == 1) .OR. (RESTRT) .OR. (MOD(NTSD,CALC_NVEG) == 0)) THEN
!------------------------------
! UPDATE EVERY UPDATE EVERY 24H
!------------------------------
                CALL FLDS_UPDATE_DRIVER('UPDATE          ','VGREEN          ')
        END IF
!
!--------------------------------------------------------------------------------------------- 
! START THE ADJUSTMENT STEP: INTEGRATE FORWARD THE CONTINUITY EQUATION (UPDATE THE MASS FIELD)
!--------------------------------------------------------------------------------------------- 
!
!------------------------------------------------------- 
! DIVERGENCE AND HORIZONTAL PART OF THE OMEGA-ALPHA TERM
!------------------------------------------------------- 
        IF (NTSD > 1)     CALL EXCH(T , LM, U , LM, V, LM, Q, LM, CWM, LM, 2, 2) 
!
        IF (.NOT. HYDRO) THEN
            IF (NTSD > 1) CALL EXCH(DWDT, LM, PINT, LM+1, 5, 5)
        END IF
 !   
        BTIM = TIMEF()
!
        IF (SLOPE) THEN
            CALL DIVHOASTQL
        ELSE
            CALL DIVHOA
        END IF
!
        DIVHOA_TIM = DIVHOA_TIM + TIMEF() - BTIM  
!------------------------------------------------------- 
! PRESSURE TENDENCY, ETA/SIGMA DOT, VERTICAL OMEGA-ALPHA 
!------------------------------------------------------- 
        BTIM = TIMEF()
!
        CALL EXCH(PD, 1, DIV, LM, PINT, LM+1, 2, 2)
!
        EXCH_TIM = EXCH_TIM     + TIMEF() - BTIM
!   
        BTIM = TIMEF()
!
        CALL PDTEDT !CONTAINS CALL TO EXCH
!
        PDTEDT_TIM = PDTEDT_TIM + TIMEF() - BTIM
!-------------------------------------------------------
! DO VERTICAL ADVECTION WITHIN THE FIRST ADJUSTMENT STEP 
!-------------------------------------------------------
!        IF (MOD(NTSD-1,IDTAD) == 0) THEN
        IF (MOD(NTSD,IDTAD) == 1) THEN
            BTIM = TIMEF()
!
            CALL EXCH(ETADT,LM-1,1,1)
!
            EXCH_TIM = EXCH_TIM   + TIMEF() - BTIM
!        
            BTIM = TIMEF()
!
             CALL VTADV
!
             VTADV_TIM = VTADV_TIM + TIMEF() - BTIM
!        
            BTIM = TIMEF()
!
            CALL EXCH(T, LM, U, LM, V, LM, Q, LM, Q2, LM, 1, 1)
!
            EXCH_TIM  = EXCH_TIM  + TIMEF() - BTIM
        END IF
!----------------------------- 
! UPDATING PRESSURE DIFFERENCE 
!----------------------------- 
        BTIM = TIMEF()
!
        CALL PDNEW
!
        PDNEW_TIM = PDNEW_TIM + TIMEF() - BTIM
!----------------------------------------- 
! PDATING BOUNDARY VALUES AT HEIGHT POINTS 
!-----------------------------------------  
        BTIM = TIMEF()
!
!        IF (MOD(NTSD,IDTAD) == 0) THEN
        IF (MOD(NTSD,IDTAD) == 1) THEN
!
            CALL EXCH(T, LM, Q, LM, Q2, LM, 1, 1)
!
        END IF
!
        CALL EXCH(PD, 1, CWM, LM, 1, 1)
!
        EXCH_TIM  =  EXCH_TIM + TIMEF() - BTIM
!    
        BTIM = TIMEF()
!
        CALL BOCOH 
!
        BOCOH_TIM = BOCOH_TIM + TIMEF() - BTIM
!-----------------------------------------------------------------
! INTEGRATE BACKWARD THE MOMENTUM EQUATION (UPDATE THE WIND FIELD)
!-----------------------------------------------------------------
!
!------------------------------------------ 
!PRESSURE GRADIENT AND CORIOLIS FORCE TERMS
!------------------------------------------
        BTIM = TIMEF()
!
        CALL EXCH(PD, 1, T, LM, Q, LM, 2, 2)
        CALL EXCH(CWM, LM, 2, 2)
!
        IF (.NOT. HYDRO) THEN
            CALL EXCH(PINT, LM+1, 5, 5)
        END IF
!    
        EXCH_TIM  = EXCH_TIM  + TIMEF() - BTIM
!    
        BTIM = TIMEF()
!
        CALL PGCOR
!
        PGCOR_TIM = PGCOR_TIM + TIMEF() - BTIM
!    
        BTIM = TIMEF()
!
        CALL EXCH(PDSL, 1, 5, 5)
!
        EXCH_TIM =  EXCH_TIM  + TIMEF() - BTIM
!-------------------  
! DIVERGENCE DAMPING 
!-------------------  
        IF ((MOD(NTSD,NTDDMP) == 0).AND.(LDDAMPFLG))THEN
            BTIM = TIMEF()
!
            CALL EXCH(T, LM, U, LM, V, LM, DIV, LM, 1, 1)
!
            EXCH_TIM  = EXCH_TIM  + TIMEF() - BTIM    
            BTIM = TIMEF()
!
            CALL DDAMP
!
            DDAMP_TIM = DDAMP_TIM + TIMEF() - BTIM
        END IF
!-------------------------------------------- 
! UPDATING BOUNDARY VALUES AT VELOCITY POINTS 
!-------------------------------------------- 
        BTIM = TIMEF()
!
        CALL EXCH(U, LM, V, LM, 1, 1)
!
        EXCH_TIM  = EXCH_TIM  + TIMEF() - BTIM    
        BTIM = TIMEF()
!
        CALL BOCOV
!
        BOCOV_TIM = BOCOV_TIM + TIMEF() - BTIM
!-------------------------------------------------------------------------------------------------- 
! THE ADJUSTMENT STEP IS NOW DONE. MAKE THE REMAINING CALLS WHICH TRADITIONALLY (SO FAR) HAVE BEEN 
! DONE EVERY ADJUSTMENT STEP
!--------------------------------------------------------------------------------------------------
!
!-------------------------------------------- 
! APPLY TEMPERATURE TENDENCY DUE TO RADIATION 
!--------------------------------------------
        BTIM = TIMEF()
!
        CALL RDTEMP
!
        RDTEMP_TIM = RDTEMP_TIM + TIMEF() - BTIM
!------------------  
! LATERAL DIFFUSION 
!------------------ 
        BTIM = TIMEF()
!
        CALL EXCH(T , LM, U, LM, V, LM, Q, LM, 2, 2)
!
        CALL EXCH(Q2, LM, 2, 2)
        CALL EXCH(CWM, LM, 2, 2)
!
        EXCH_TIM   = EXCH_TIM   + TIMEF() - BTIM
!    
        BTIM = TIMEF()
!
        IF (SLOPE) THEN
            CALL HDIFFS
        ELSE
            CALL HDIFF
        END IF
!
        HDIFF_TIM  = HDIFF_TIM  + TIMEF() - BTIM
!--------------------- 
! HORIZONTAL ADVECTION  
!---------------------  
!        IF (MOD(NTSD,IDTAD) == 0) THEN
        IF (MOD(NTSD,IDTAD) == 1) THEN
            BTIM = TIMEF()
!
            CALL EXCH(T , LM, U, LM, V, LM, 4, 4)
!
            CALL EXCH(Q2, LM, 5, 5)
!
            EXCH_TIM  = EXCH_TIM + TIMEF() - BTIM
!        
            BTIM = TIMEF()
!
            IF (SLOPE) THEN
                CALL HZADVS
!                CALL SLADVT
            ELSE
                CALL HZADV
!                CALL HZADV_LM1
            END IF
!
            HZADV_TIM = HZADV_TIM + TIMEF() - BTIM
!        
            BTIM = TIMEF()
!
            CALL EXCH(U, LM, V, LM, Q, LM, CWM, LM, 2, 2)
!
            EXCH_TIM = EXCH_TIM + TIMEF() - BTIM
!----------------------------------------        
! HORIZONTAL ADVECTION OF WATER SUBSTANCE
!----------------------------------------         
            BTIM = TIMEF()
!
            CALL HZADV2
!
            HZADV2_TIM = HZADV2_TIM + TIMEF() - BTIM
        END IF
!--------------------------------------------------------------------------------------------------
! IF THE TIME IS RIGHT, NOW DO VARIOUS PHYSICS CALLS
! (WARNING: TO AVOID ENDING THE INTEGRATION WITH PHYSICS CALLS WHICH HAVE NOT BEEN FOLLOWED BY 
! ADJUSTMENT STEPS, PHYSICS CALLS ARE OFFSET BY HALVES OF VARIOUS CALLING INTERVALS.  
! IT IS ASSUMED THAT THE CALLING INTERVALS, NPHS AND NCNVC, ARE DIVISIBLE BY IDTAD. 
! IF NOT, INTEGRATION WILL END WITH AN INCORRECT NUMBER OF CALLS HAVING BEEN MADE.
!--------------------------------------------------------------------------------------------------
!
!-------------------------------------- 
! TURBULENT PROCESSES AND PRECIPITATION 
!--------------------------------------
!        IF (MOD(NTSD-NPHS/2,NPHS) == 0) THEN
        IF (MOD(NTSD,NPHS) == 1) THEN
!
            IF (MYPE == 0) THEN
                WRITE(0,"(A)") 'EBU:  PHYSICS TIME STEP'
            END IF
!
            BTIM = TIMEF()
!
            CALL EXCH(PD, 1, UZ0, 1, VZ0, 1, T, LM, U, LM, V, LM, Q, LM, 1, 1)
            CALL EXCH(CWM, LM, 1, 1)
!
            EXCH_TIM      = EXCH_TIM  + TIMEF() - BTIM
!        
            BTIM = TIMEF()
!
            CALL TURBL !CONTAINS CALLS TO EXCH
!
            TURBL_TIM = TURBL_TIM + TIMEF() - BTIM
!--------------------------------- 
! STORE ORIGINAL TEMPERATURE ARRAY
!--------------------------------- 
            IF (CLTEND_TEST) THEN
                BTIM = TIMEF()
                ICLTEND = -1
!
                CALL CLTEND(ICLTEND)
!
                CLTEND_TIM = CLTEND_TIM + TIMEF() - BTIM
            END IF
!
            IF (MICROPHYSSQM == 2) THEN
                BTIM = TIMEF()
!
                CALL GSCOND
!
                GSCOND_TIM   = GSCOND_TIM + TIMEF() - BTIM
            END IF
!-------------------------
! SHALLOW CONVECTION
!
            IF (SHLCNVFLG) THEN
                 BTIM = TIMEF()
                 SELECT CASE (SHLCNVSQM)
                     CASE (1)
                       CALL CUCNVC_SHALLOW  !> SHLCNVSQM = 1 => USING scheme based on Betts-Miller Convective Scheme'
!
                       SHLCNV_TIM = SHLCNV_TIM + TIMEF() - BTIM
                     CASE (2)
                       CALL CUCNVC_SHALLOW  !> SHLCNVSQM = 2 => USING .... 
!
                       SHLCNV_TIM = SHLCNV_TIM + TIMEF() - BTIM
                END SELECT
            END IF
!
!------------------------- 
! CONVECTIVE PRECIPITATION
!
            IF (CUCNVCFLG) THEN
                 BTIM = TIMEF()
                 SELECT CASE (CUCNVCSQM)
                     CASE (1) 
                       CALL CUCNVC !> CUCNVCSQM = 1 => USING Betts-Miller Convective Scheme'
!
                       CUCNVC_TIM = CUCNVC_TIM + TIMEF() - BTIM
                     CASE (2) 
                       CALL CUCNVC !> CUCNVCSQM = 2 => USING Kain-Fritsch Convective Scheme'
!
                       CUCNVC_TIM = CUCNVC_TIM + TIMEF() - BTIM
                END SELECT
            END IF
!
!-----------------------------------------------------------------------
! GRIDSCALE MICROPHYSICS (CONDENSATION AND PRECIPITATION; FORECAST ONLY)
!-----------------------------------------------------------------------
            BTIM = TIMEF()
            SELECT CASE (MICROPHYSSQM)
                 CASE (1) 
                   CALL GSMDRIVE !> MICROPHYSSQM = 1 => USING Ferrier Microphysics Scheme'
                   GSMDRIVE_TIM = GSMDRIVE_TIM + TIMEF() - BTIM
                 CASE (2) 
                   CALL PRECPD   !> MICROPHYSSQM = 2 => USING Zhao Microphysics Scheme'
                   PRECPD_TIM     = PRECPD_TIM   + TIMEF() - BTIM
             END SELECT
!
!---------------------------- 
! PRECIPIPTATION ASSIMILATION 
!----------------------------
!
!----------------------------------------------------- 
! CALCULATE TEMP TENDENCIES AND RESTORE ORIGINAL TEMPS 
!----------------------------------------------------- 
            IF (CLTEND_TEST) THEN
                BTIM    = TIMEF()
                ICLTEND = 0
!
                CALL CLTEND(ICLTEND)
!
                CLTEND_TIM = CLTEND_TIM + TIMEF() - BTIM
            END IF
!-------------------- 
! END PHYSICS IF LOOP
!-------------------- 
        END IF
!------------------------------------------------------------ 
! UPDATE TEMP TENDENCIES FROM CLOUD PROCESSES EVERT TIME STEP
!------------------------------------------------------------ 
        IF (CLTEND_TEST) THEN
            BTIM    = TIMEF()
            ICLTEND = 1
!
            CALL CLTEND(ICLTEND)
!
            CLTEND_TIM = CLTEND_TIM + TIMEF() - BTIM
        END IF  
!----------------------------- 
! VERTICAL ADVECTION OF HEIGHT 
!----------------------------- 
        BTIM = TIMEF()
!
        CALL VADZ
!
        VADZ_TIM = VADZ_TIM  + TIMEF() - BTIM
!------------------------------- 
! HORIZONTAL ADVECTION OF HEIGHT 
!------------------------------- 
        IF ( .NOT. HYDRO) THEN
            BTIM = TIMEF()
!
            CALL EXCH(U, LM  , V, LM, 1, 1)
!
            CALL EXCH(Z, LM+1, 2, 2)
!
            EXCH_TIM = EXCH_TIM  + TIMEF() - BTIM
        END IF
!    
        BTIM = TIMEF()
!
        CALL HADZ
!
        HADZ_TIM = HADZ_TIM + TIMEF() - BTIM
!--------------- 
! ADVECTION OF W 
!--------------- 
        IF (HYDRO) THEN
            BTIM = TIMEF()
!
            CALL EXCH(PDSL, 1   , 2, 2)
!
            CALL EXCH(PINT, LM+1, 3, 3)
!
            EXCH_TIM  = EXCH_TIM + TIMEF() - BTIM
        ELSE
            BTIM = TIMEF()
!
            CALL EXCH(PDSL, 1 , 2, 2)
!
            CALL EXCH(U   , LM, V, LM, DWDT, LM, PINT, LM+1, W, LM+1, 3, 3)
!
            EXCH_TIM = EXCH_TIM + TIMEF() - BTIM
        END IF
!    
        BTIM = TIMEF()
!
        CALL EPS
!
        EPS_TIM = EPS_TIM  + TIMEF() - BTIM
!
        IF (NTSD > NSTART+1) THEN
!---------- 
! RADIATION 
!---------- 
            IF (MOD(NTSD,NRADS) == 1 .OR. MOD(NTSD,NRADL) == 1) THEN
                BTIM = TIMEF()
!
                CALL RADTN
!
                RADTN_TIM = RADTN_TIM + TIMEF() - BTIM
!
            END IF
!------------------------------------------------------- 
! IS IT TIME FOR A CHECK POINT ON THE MODEL HISTORY FILE 
!------------------------------------------------------- 
            BTIM = TIMEF()
!
            CALL CHKOUT
!
            CHKOUT_TIM = CHKOUT_TIM + TIMEF() - BTIM
        END IF
!----------------------- 
! CLEAN UP AFTER RESTART 
!----------------------- 
        IF (RESTRT) THEN
            RESTRT = .FALSE. 
        END IF
!----------------------------------------------------------------
! WRITING IN A SEPARATE FILE FOR EACH MYPE AND EACH OUTPUT PERIOD
!----------------------------------------------------------------
!        IF ( NTSD == 1 .OR. MOD((NTSD-1)*DT,3*3600) == 0) THEN
!        IF ( MOD(NTSD,12*3600/DT) == 1) THEN
         IF (WRITEOUT2RESTRT) THEN
             NOut = 24*FREQOUT2RESTRT*3600/DT     
             IF ( MOD(NTSD,NOut) == 1) THEN   ! Out2restrt file each 60 days
                 IF (MYPE == 0) PRINT*, 'OUT2RESTRT CALLED, NTSD= ', NTSD, NTSD*DT/3600
	         CALL OUT2RESTRT(ITAG)
             ENDIF	
        ENDIF
! 
        IF (NTSD < NTSTM) GOTO 2000
!------------------------ 
! EXIT FROM THE TIME LOOP     
!------------------------ 
        2005 CONTINUE
!	
        TOT2_TIM = TIMEF() - BTIMX
        TOT_TIM  =    MPP_TIM +   INIT_TIM +   GOSS_TIM + RADTN_TIM    + CHKOUT_TIM + DIVHOA_TIM  &
    &            + PDTEDT_TIM +  VTADV_TIM +  PDNEW_TIM + BOCOH_TIM    +  PGCOR_TIM +  DDAMP_TIM  &
    &            +  BOCOV_TIM + RDTEMP_TIM +  HDIFF_TIM + HZADV_TIM    + HZADV2_TIM +  TURBL_TIM  &
    &            + GSCOND_TIM + CUCNVC_TIM +   EXCH_TIM + GSMDRIVE_TIM + CLTEND_TIM +   VADZ_TIM  &
    &            +   HADZ_TIM +    EPS_TIM + PRECPD_TIM
!    
        IF (MYPE == 0) THEN
            PCT =      MPP_TIM / TOT_TIM*1.E2
            WRITE(6,*)       ' MPP=',MPP_TIM*1.E-3,'    PCT=',PCT
            PCT =     INIT_TIM / TOT_TIM*1.E2
            WRITE(6,*)      ' INIT=',INIT_TIM*1.E-3,'   PCT=',PCT
            PCT =     GOSS_TIM / TOT_TIM*1.E2
            WRITE(6,*)      ' GOSS=',GOSS_TIM*1.E-3,'   PCT=',PCT
            PCT =    RADTN_TIM / TOT_TIM*1.E2
            WRITE(6,*)     ' RADTN=',RADTN_TIM*1.E-3,'  PCT=',PCT
            PCT =   CHKOUT_TIM / TOT_TIM*1.E2
            WRITE(6,*)    ' CHKOUT=',CHKOUT_TIM*1.E-3,' PCT=',PCT
            PCT =   DIVHOA_TIM / TOT_TIM*1.E2
            WRITE(6,*)    ' DIVHOA=',DIVHOA_TIM*1.E-3,' PCT=',PCT
            PCT =   PDTEDT_TIM / TOT_TIM*1.E2
            WRITE(6,*)  ' PDTEDT=',PDTEDT_TIM*1.E-3,'   PCT=',PCT
            PCT =    VTADV_TIM / TOT_TIM*1.E2
            WRITE(6,*)     ' VTADV=',VTADV_TIM*1.E-3,'  PCT=',PCT
            PCT =    PDNEW_TIM / TOT_TIM*1.E2
            WRITE(6,*)     ' PDNEW=',PDNEW_TIM*1.E-3,'  PCT=',PCT
            PCT =    BOCOH_TIM / TOT_TIM*1.E2
            WRITE(6,*)     ' BOCOH=',BOCOH_TIM*1.E-3,'  PCT=',PCT
            PCT =    PGCOR_TIM / TOT_TIM*1.E2
            WRITE(6,*)     ' PGCOR=',PGCOR_TIM*1.E-3,'  PCT=',PCT
            PCT =   PRECPD_TIM / TOT_TIM*1.E2
            WRITE(6,*)   ' PRECPD=',PRECPD_TIM*1.E-3,'  PCT=',PCT
            PCT =    DDAMP_TIM / TOT_TIM*1.E2
            WRITE(6,*)     ' DDAMP=',DDAMP_TIM*1.E-3,'  PCT=',PCT
            PCT =    BOCOV_TIM / TOT_TIM*1.E2
            WRITE(6,*)     ' BOCOV=',BOCOV_TIM*1.E-3,'  PCT=',PCT
            PCT =   RDTEMP_TIM / TOT_TIM*1.E2
            WRITE(6,*)    ' RDTEMP=',RDTEMP_TIM*1.E-3,' PCT=',PCT
            PCT =    HDIFF_TIM / TOT_TIM*1.E2
            WRITE(6,*)     ' HDIFF=',HDIFF_TIM*1.E-3,'  PCT=',PCT
            PCT =    HZADV_TIM / TOT_TIM*1.E2
            WRITE(6,*)     ' HZADV=',HZADV_TIM*1.E-3,'  PCT=',PCT
            PCT =   HZADV2_TIM / TOT_TIM*1.E2
            WRITE(6,*)    ' HZADV2=',HZADV2_TIM*1.E-3,' PCT=',PCT
            PCT =     VADZ_TIM / TOT2_TIM*1.E2
            WRITE(6,*)       ' VADZ=',VADZ_TIM*1.E-3,'  PCT=',PCT
            PCT =     HADZ_TIM / TOT2_TIM*1.E2
            WRITE(6,*)       ' HADZ=',HADZ_TIM*1.E-3,'  PCT=',PCT
            PCT =      EPS_TIM / TOT2_TIM*1.E2
            WRITE(6,*)         ' EPS=',EPS_TIM*1.E-3,'  PCT=',PCT
            PCT =    TURBL_TIM /  TOT_TIM*1.E2
            WRITE(6,*)     ' TURBL=',TURBL_TIM*1.E-3,'  PCT=',PCT
            PCT =   CUCNVC_TIM /  TOT_TIM*1.E2
            WRITE(6,*)    ' CUCNVC=',CUCNVC_TIM*1.E-3,' PCT=',PCT
            PCT = GSMDRIVE_TIM /  TOT_TIM*1.E2
            WRITE(6,*)' GSMDRIVE=',GSMDRIVE_TIM*1.E-3,' PCT=',PCT
            PCT =   CLTEND_TIM /  TOT_TIM*1.E2
            WRITE(6,*)    ' CLTEND=',CLTEND_TIM*1.E-3,' PCT=',PCT
            PCT =     EXCH_TIM /  TOT_TIM*1.E2
            WRITE(6,*)      ' EXCH=',EXCH_TIM*1.E-3,'   PCT=',PCT
!
            WRITE(6,*) ' TOTAL=', TOT_TIM*1.E-3
            WRITE(6,*)' TOTAL2=',TOT2_TIM*1.E-3
        END IF
!------------------------------------------------------------------
! WE MUST NOW SHUT DOWN THE I/O SERVERS
! THIS IS DONE BY SENDING A -999 TO MPI TASK 0 OF EACH SERVER GROUP
!------------------------------------------------------------------    
        IF (MYPE == 0) THEN
            DO I=1,IQUILT_GROUP
                CALL MPI_SEND(-999, 1, MPI_INTEGER, 0, 0, MPI_COMM_INTER_ARRAY(I), IER)
            END DO
        END IF
!---------------------------------------------------- 
!END IF ON TASKS FOR MODEL INTEGRATION VS I/O SERVING
!---------------------------------------------------- 
    END IF 
!
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
!
    IF (MYPE == 0) THEN
        CALL W3TAGE('ETAFCST ')
    END IF
!
    CALL MPI_FINALIZE(IERR)
!
    END PROGRAM EBU
