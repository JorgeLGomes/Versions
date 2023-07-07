    SUBROUTINE INIT
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE INIT
!> 
!> SUBPROGRAM: INIT - INITIALIZE VARIABLE FOR MODEL RUN
!> PROGRAMMER: JANJIC 
!> ORG: W/NP22
!> DATE: ??-??-??
!>
!> ABSTRACT:  
!> INIT READS IN PRIMARY AND AUXILIARY VARIABLES AND CONSTANTS AND SETS INITIAL VALUES FOR OTHERS.
!>
!> PROGRAM HISTORY LOG:
!> 87-06-??  JANJIC  - ORIGINATOR
!> 92-10-27  DEAVEN  - CHANGED READS OF NHB, NFC, AND NBC TO ACCOMODATE SHORTENED RECORD LENGTHS
!> 95-03-27  BLACK   - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 96-10-31  BLACK   - ADDED NAMELIST BCEXDATA FOR THE NESTS
!> 98-06-10  ROGERS  - MADE Y2K COMPLIANT BY REPLACING CALL TO W3FI13 TO W3DOXDAT
!> 98-09-04  PYLE    - CHANGED TO NOT RE-INITIALIZE TSHLTR AND QSHLTR IF RESTART=TRUE
!> 98-10-21  BLACK   - CHANGES FOR DISTRIBUTED MEMORY
!> 98-11-17  BLACK   - ADDED CODE TO LOCATE THE INNER DOMAIN BOUNDARIES ON THE RELEVANT PEs
!> 00-08-??  BLACK   - MODIFIED FOR RESTART CAPABILITY
!> 18-01-15  LUCCI   - MODERNIZATION OF THE CODE, INCLUDING:
!>                     * F77 TO F90/F95
!>                     * INDENTATION & UNIFORMIZATION CODE
!>                     * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                     * DOCUMENTATION WITH DOXYGEN
!>                     * OPENMP FUNCTIONALITY 
!>
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
!> INPUT FILES:
!> NFC - THE INITIAL VALUES OF SFC PRESSURE, T, Q, U, AND V
!> NHB - A LARGE VARIETY OF ARRAY AND SCALAR CONSTANTS
!> NBC - THE BOUNDARY CONDITIONS AND TENDENCIES
!>
!>                         OR
!>
!> RESTRT - A RESTART FILE WITH ALL NECESSARY QUANTITIES
!>
!> OUTPUT FILES:
!> NONE
!>
!> USE MODULES: ACMCLD
!>              ACMCLH 
!>              ACMPRE
!>              ACMRDL
!>              ACMRDS
!>              ACMSFC
!>              BOCO
!>              CLDWTR
!>              CMICRO_START
!>              CNVCLD
!>              CONTIN
!>              CTLBLK
!>              CUINIT
!>              CUPARM
!>              C_TADJ
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
!>              PARMETA
!>              PARMSOIL
!>              PARM_TBL
!>              PHYS
!>              PPTASM
!>              PVRBLS
!>              RD1TIM
!>              SOIL
!>              TEMPCOM
!>              TEMPV
!>              TIMMING
!>              TOPO
!>              VRBLS
!>              Z0EFFT
!>
!> DRIVER     : EBU
!>
!> CALLS      : DSTRB
!>              GRADFS
!>              MPI_BARRIER
!>              MPI_BCAST
!>              MPI_FINALIZE
!>              O3CLIM
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>              SOLARD
!>              ZENITH
!>              ZERO2
!>              ZERO3
!>--------------------------------------------------------------------------------------------------
    USE ACMCLD
    USE ACMCLH 
    USE ACMPRE 
    USE ACMRDL 
    USE ACMRDS 
    USE ACMSFC
    USE BOCO
    USE CLDWTR
    USE CMICRO_START
    USE CNVCLD
    USE CONTIN
    USE CTLBLK
    USE CUINIT
    USE CUPARM
    USE C_TADJ
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
    USE PARMETA
    USE PARMSOIL
    USE PARM_TBL
    USE PHYS
    USE PPTASM
    USE PVRBLS
    USE RD1TIM
    USE SOIL
    USE TEMPCOM
    USE TEMPV
    USE TIMMING
    USE TOPO
    USE VRBLS
    USE Z0EFFT
!
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!
#include "sp.h"
! 
    REAL   (KIND=R4KIND), PARAMETER :: CM1    = 2937.4
    REAL   (KIND=R4KIND), PARAMETER :: CM2    =    4.9283
    REAL   (KIND=R4KIND), PARAMETER :: CM3    =   23.5518
    REAL   (KIND=R4KIND), PARAMETER :: EPS    =    0.622 
    REAL   (KIND=R4KIND), PARAMETER :: PI2    = 2. * 3.14159265
    REAL   (KIND=R4KIND), PARAMETER :: RLAG   =   14.8125
    REAL   (KIND=R4KIND), PARAMETER :: Q2INI  =     .50 
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ2  =    0.2
    REAL   (KIND=R4KIND), PARAMETER :: EPSWET =    0.0 
    REAL   (KIND=R4KIND), PARAMETER :: Z0LAND =     .10
    REAL   (KIND=R4KIND), PARAMETER :: Z0SEA  =     .001
    REAL   (KIND=R4KIND), PARAMETER :: FCM    =     .00001 
    REAL   (KIND=R4KIND), PARAMETER :: DTR    =    0.1745329E-1
    REAL   (KIND=R4KIND), PARAMETER :: A1     =  610.78
    REAL   (KIND=R4KIND), PARAMETER :: WA     =     .10
    REAL   (KIND=R4KIND), PARAMETER :: WG     =    1.0 - WA
!
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM   = IM * JM - JM / 2
    INTEGER(KIND=I4KIND), PARAMETER :: JMP1   = JM + 1
    INTEGER(KIND=I4KIND), PARAMETER :: IMT    = 2  * IM - 1
    INTEGER(KIND=I4KIND), PARAMETER :: NSTAT  = 1000
!------------------  
! DECLARE VARIABLES
!------------------  
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & RUNB    , EXBC    , INSIDEH , INSIDEV
!
    CHARACTER(32)                                                                               ::&
    & LABEL
!
    CHARACTER(40)                                                                               ::&
    & CONTRL  , FILALL  , FILMST  , FILTMP  , FILTKE  , FILUNV  , FILCLD  , FILRAD  , FILSFC
!
    REAL   (KIND=R4KIND), DIMENSION(LP1)                                                        ::&
    & PHALF
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & NPEBND
!--------------------------------------------------------------------------------------------------
! NOTE: THE DIMENSION OF THE FOLLOWING ARRAYS IS ARBITRARILY CHOSEN TO EXCEED ANY NUMBER OF 
! BOUNDARY POINTS THAT MIGHT EXIST IN ANY INNER DOMAIN
!--------------------------------------------------------------------------------------------------
    REAL   (KIND=R4KIND), DIMENSION(1500)                                                       ::&
    &  HLATI  , HLONI   , VLATI   , VLONI   , THLONI  , THLATI  , TVLONI  , TVLATI
!
    REAL   (KIND=R4KIND), DIMENSION(NSTAT)                                                      ::&
    & TSLAT   , TSLON
!
    INTEGER(KIND=I4KIND), DIMENSION(3)                                                          ::&
    & IDATB
!
    INTEGER(KIND=I4KIND), DIMENSION(8)                                                          ::&
    & INIDAT  , IBCDAT
!
    LOGICAL(KIND=L8KIND)                                                                        ::&
    & RUNBX
!
    INTEGER(KIND=I8KIND), DIMENSION(3)                                                          ::&
    & IDATBX
!
    INTEGER(KIND=I8KIND)                                                                        ::&
    & IHRSTBX               
!------------------------------------------ 
! THE FOLLOWING IS FOR TIMIMG PURPOSES ONLY
!------------------------------------------ 
    REAL   (KIND=R8KIND)                                                                        ::&
    & TIMEF
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NFILE   , IUNWGT  , NMAP    , NRADSH  , NRADLH  , J       , KNT     , I       , N       ,   &
    & NHB     , NHIBU   , LFCSTD  , K       , IERR    , KBI     , KBI2    , LRECBC  , IHRSTB  ,   &
    & IRTN    , ISTART  , NREC    , LLMH    , NS      , KCCO2   , LMHK    , LMVK      
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & PLOMD   , PMDHI   , PHITP   , P400    , PLBTM   , TSTART  , TEND    , TCP     , BTIM    ,   &
    & TBOCO   , BCHR    , TERM1   , PM1     , APEM1   , TSFCK   , PSFCK   , CLOGES  , ESE     ,   &
    & PSUM    , SLPM    , PDIF    , PSS     , TIME    , DAYI    , HOUR    , ADDL    , RANG    ,   &
    & RSIN1   , RCOS1   , RCOS2   , ULM     , VLM     , TLM     , QLM     , PLM2    , APELM   ,   &
    & APELMNW , EXNERR  , THLM    , DPLM    , DZLM    , FAC1    , FAC2  
!
    DATA PLOMD/64200./, PMDHI/35000./, PHITP/15000./, P400/40000./, PLBTM/105000./
!
    DATA NFILE/14/, IUNWGT/40/
!------------------------------------------------------ 
! FLAG FOR INITIALIZING CONVECTIVE CLOUDS FOR RADIATION
!------------------------------------------------------ 
!------------------- 
! DECLARE NAMELISTS.
!-------------------
    NAMELIST /FCSTDATA/                                                                           &
    & TSTART, TEND  , TCP  , RESTRT, SINGLRST, SUBPOST, NMAP, TSHDE, SPL , NPHS , NCNVC , NRADSH, &
    & NRADLH, NTDDMP, TPREC, THEAT , TCLOD   , TRDSW  ,TRDLW, TSRFC, NEST, HYDRO, SPLINE
!-----------------
! START INIT HERE.
!-----------------
!
!-------------------------------------------
! CALCULATE THE I-INDEX EAST-WEST INCREMENTS
!-------------------------------------------
    DO J=1,JM
        IHEG(J) = MOD(J+1,2)
        IHWG(J) = IHEG(J) - 1
        IVEG(J) = MOD(J,2)
        IVWG(J) = IVEG(J) - 1
    END DO
!-------------------------------------------
! CALCULATE THE INDIRECT I INDICES FOR RADTN
!-------------------------------------------
    KNT = 0
!
    DO I=1,IM
        KNT = KNT + 1
        IRADG(KNT) = I
    END DO
!
    DO I=1,IM-1
        KNT = KNT + 1
        IRADG(KNT) = IM + 2 + I
    END DO
!-------------------------------- 
! ZERO OUT LOCALLY INDEXED ARRAYS
!-------------------------------- 
    CALL ZERO2(PDSL)
    CALL ZERO3(T     , LM  )
    CALL ZERO3(Q     , LM  )
    CALL ZERO3(U     , LM  )
    CALL ZERO3(V     , LM  )
    CALL ZERO2(RES)
    CALL ZERO3(RTOP  , LM  )
    CALL ZERO3(OMGALF, LM  )
    CALL ZERO3(DIV   , LM  )
    CALL ZERO3(ETADT , LM-1)
    CALL ZERO3(HTM   , LM  )
    CALL ZERO3(VTM   , LM  )
    CALL ZERO2(HBM2)
    CALL ZERO2(AKMS)
    CALL ZERO2(UZ0)
    CALL ZERO2(VZ0)
    CALL ZERO2(FAD)
!------------------ 
! READ Z0 EFFECTIVE
!------------------
    OPEN(UNIT=22, FILE='ZEFF', FORM='UNFORMATTED')
!
    DO N=1,4
        IF (MYPE == 0) THEN
            READ(22) TEMP1
        END IF
!
        CALL DSTRB(TEMP1,ZEFFIJ,1,4,N)
!
    END DO
!------------------------------------------------ 
! READ "CONSTANT" DATA FROM UNIT CONNECTED TO NHB
!------------------------------------------------ 
    NHB  = 12
    LSL  = LSM
    BTIM = TIMEF()
!
    CALL READ_NHB(NHB)
!
    NHB_TIM = TIMEF() - BTIM
!
    NHIBU = 12
    IF (MYPE == 0) WRITE(LIST,*)'INIT:  READ CONSTANTS FILE'
!--------------------------------------------------------------------------------------- 
! READ NAMELIST FCSTDATA WHICH CONTROLS TIMESTEPS, ACCUMULATION PERIODS, STANDARD OUTPUT
!---------------------------------------------------------------------------------------
    RESTRT = .FALSE. 
    LFCSTD = 11
!
    OPEN(LFCSTD, FILE='FCSTDATA.meso', STATUS='OLD')
    REWIND LFCSTD
    READ(LFCSTD, FCSTDATA)
!
    IF  (MYPE == 0) THEN
        WRITE(6,*) 'HYDRO='  , HYDRO
        WRITE(6,*) 'SPLINE=' , SPLINE
    END IF
!
    IF (MYPE == 0) THEN
        WRITE(LIST,*)'INIT:  READ NAMELIST FCSTDATA - LISTED BELOW'
        WRITE(LIST,*)'  TSTART,TEND  :  ',TSTART, TEND
        WRITE(LIST,*)'  TCP          :  ',TCP
        WRITE(LIST,*)'  RESTRT       :  ',RESTRT
        WRITE(LIST,*)'  HYDRO        :  ',HYDRO
        WRITE(LIST,*)'  SINGLRST     :  ',SINGLRST
        WRITE(LIST,*)'  SUBPOST      :  ',SUBPOST
        WRITE(LIST,*)'  NMAP,NPHS    :  ',NMAP  , NPHS
        WRITE(LIST,*)'  NCNVC        :  ',NCNVC
        WRITE(LIST,*)'  NRADSH,NRADLH:  ',NRADSH, NRADLH
        WRITE(LIST,*)'  NTDDMP       :  ',NTDDMP
        WRITE(LIST,*)'  TPREC,THEAT  :  ',TPREC , THEAT
        WRITE(LIST,*)'  TCLOD,TRDSW  :  ',TCLOD , TRDSW
        WRITE(LIST,*)'  TRDLW,TSRFC  :  ',TRDLW , TSRFC
        WRITE(LIST,*)'  TSHDE (POSTED FORECAST HOURS) BELOW:  '
        WRITE(LIST,75) (TSHDE(K),K=1,NMAP)
        WRITE(LIST,*)'  SPL (POSTED PRESSURE LEVELS) BELOW: '
        WRITE(LIST,80) (SPL(K),K=1,LSM)
        75 FORMAT(14(F4.1,1X))
        80 FORMAT(8(F8.1,1X))
    END IF
!------------------------------------- 
! SET TIME STEPPING RELATED CONSTANTS.
!------------------------------------- 
    FIRST  = .TRUE. 
    NSTART = INT(TSTART * TSPH + 0.5)
    NTSTM  = INT(TEND   * TSPH + 0.5) + 1
    NCP    = INT(TCP    * TSPH + 0.5)
    NPREC  = INT(TPREC  * TSPH + 0.5)
    NHEAT  = INT(THEAT  * TSPH + 0.5)
    NCLOD  = INT(TCLOD  * TSPH + 0.5)
    NRDSW  = INT(TRDSW  * TSPH + 0.5)
    NRDLW  = INT(TRDLW  * TSPH + 0.5)
    NSRFC  = INT(TSRFC  * TSPH + 0.5)
!
    IF (MYPE == 0) THEN
        WRITE(0,*)' NTSTM=',NTSTM,' TSPH=',TSPH,' DT=',DT
    END IF
!------------------------------------------------
! SET VARIOUS PHYSICS PACKAGE TIMESTEP VARIABLES.
!------------------------------------------------
    NRADS = NINT(TSPH) * NRADSH
    NRADL = NINT(TSPH) * NRADLH
    DTQ2  = NPHS * DT
    TDTQ2 = DTQ2 + DTQ2
    DTD   = 0.5  * DTQ2
    TDTD  = DTD  + DTD
    KTM   = INT(DTQ2 / DTD + 0.5)
!
    IF (MYPE == 0) THEN
        WRITE(LIST,*)' '
        WRITE(LIST,*)'SET TIME STEPPING CONSTANTS'
        WRITE(LIST,*)' FIRST             :  ',FIRST
        WRITE(LIST,*)' NSTART,NSTSM,NCP  :  ',NSTART, NTSTM, NCP
        WRITE(LIST,*)' NTDDMP,NPREC,NHEAT:  ',NTDDMP, NPREC, NHEAT
        WRITE(LIST,*)' NCLOD,NRDSW,NRDLW :  ',NCLOD , NRDSW, NRDLW
        WRITE(LIST,*)' NSRFC             :  ',NSRFC
        WRITE(LIST,*)' NRADS,NRADL,KTM   :  ',NRADS , NRADL, KTM
        WRITE(LIST,*)' DTQ2,TDTQ2        :  ',DTQ2  , TDTQ2
        WRITE(LIST,*)' DTD,TDTD          :  ',DTD   , TDTD
        WRITE(LIST,*)' '
    END IF
!-------------------------------------- 
! COMPUTE DERIVED MAP OUTPUT CONSTANTS.
!-------------------------------------- 
    DO K=1,LSL
        ALSL(K) = LOG(SPL(K))
    END DO
!
    DO I=1,NMAP
        ISHDE(I) = INT(TSHDE(I) * TSPH + 0.5) + 1
    END DO
!-------------------------------------- 
! SET UP ARRAY IRAD (INDICES FOR RADTN)
!--------------------------------------
    DO I=MYIS,MYIE
        IRAD(I) = IRADG(I+MY_IS_GLB-1) - MY_IS_GLB + 1
    END DO
!----------------------------------------
! READ INITIAL CONDITIONS OR RESTART FILE.
!----------------------------------------
    BTIM = TIMEF()
!
    IF (SINGLRST) THEN
        CALL READ_RESTRT
    ELSE
        CALL READ_RESTRT2
    END IF
!
    RES_TIM = TIMEF() - BTIM
!----------------------------------------------------------------------------------------------------
! IF NOT RUNNING THE MODEL, PRINT DATE OF INITIAL CONDITIONS JUST READ AND STOP. OTHERWISE, CONTINUE.
!----------------------------------------------------------------------------------------------------
    IF (RUN) GOTO 190
!
    IF (MYPE == 0) THEN
        WRITE(LIST,165) IHRST, IDAT
        WRITE(LIST,166)
!
        CALL MPI_FINALIZE(IERR)
!
        STOP 2
        165 FORMAT('0*** ',I2,' GMT ',2(I2,'/'),I4,' ***')
        166 FORMAT('0F*** NO INITIAL CONDITIONS. RUN TERMINATED.')
    END IF
!--------------------------------------------------------------------------------------------------
! IF THE TIMESTEP COUNTER (NTSD) EXCEEDS THE "STOP MODEL" T TIMESTEP,CONTINUE, STOP EXECUTION. 
! OTHERWISE, CONTINUE.
!--------------------------------------------------------------------------------------------------
 190 IF (NTSD >= NTSTM) THEN
        IF (MYPE == 0) THEN
            WRITE(LIST,165) IHRST, IDAT
            WRITE(LIST,195)
            195 FORMAT('0F*** FORECAST ALREADY DONE. RUN TERMINATED.')
            CALL MPI_FINALIZE(IERR)
            STOP 3
        END IF
    END IF
!-------------------------- 
! READ BOUNDARY CONDITIONS.
!-------------------------- 
    IF (MYPE == 0) THEN
        OPEN(UNIT=NBC, FORM='UNFORMATTED', FILE='BNDY.file')
!
        IF (NEST) THEN
            KBI    = 2 * IM + JM - 3
            KBI2   = KBI - 4
            LRECBC = 4 * (1 + (1 + 6 * LM) * KBI * 2 + (KBI + KBI2) * (LM + 1))
            OPEN (UNIT=NBC, ACCESS='DIRECT', RECL=LRECBC)
        END IF
!    
        IF ( .NOT. NEST) REWIND NBC
!    
        IF (NEST) THEN
            READ(NBC,REC=1) RUNBX, IDATBX, IHRSTBX, TBOCO
       ELSE
            READ(NBC)       RUNBX, IDATBX, IHRSTBX, TBOCO
        END IF
!    
        RUNB   = RUNBX
        IDATB  = IDATBX
        IHRSTB = IHRSTBX
!
        IF (NEST) THEN
            READ(NBC,REC=1) RUNB, IDATB, IHRSTB, TBOCO
        ELSE
            WRITE(6,*) 'READING FROM NBC HERE'
            READ(NBC) RUNB, IDATB, IHRSTB, TBOCO
            WRITE(6,*) 'PAST FROM NBC HERE'
            WRITE(6,*) 'IDATB: ' , IDATB
            WRITE(6,*) 'IHRSTB: ', IHRSTB
        END IF
!
    END IF
!
    CALL MPI_BCAST(RUNB  , 1, MPI_LOGICAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST(IDATB , 3, MPI_INTEGER, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST(IHRSTB, 1, MPI_INTEGER, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST(TBOCO , 1, MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BARRIER                         (MPI_COMM_COMP, IRTN)
!
    IF (MYPE == 0 .AND. .NOT. NEST) THEN
        ISTART = NINT(TSTART)
    
        READ(NBC) BCHR
    205 READ(NBC)
        READ(NBC)
        READ(NBC)
        READ(NBC)
        READ(NBC)
        READ(NBC)
        READ(NBC)
!    
        IF (ISTART == NINT(BCHR)) THEN
            IF (ISTART > 0) READ(NBC) BCHR
            GOTO 215
        ELSE
            READ(NBC) BCHR
        END IF
!    
        IF (ISTART >= NINT(BCHR)) GOTO 205
    END IF
!
    IF (MYPE == 0 .AND. NEST) THEN
        ISTART = NINT(TSTART)
        NREC = 1
!    
    210 NREC = NREC + 1
        READ(NBC,REC=NREC) BCHR
    
        IF (ISTART == NINT(BCHR)) THEN
            IF (ISTART > 0) READ(NBC,REC=NREC+1) BCHR
            GOTO 215
        ELSE
            GOTO 210
        END IF
    END IF
!
215 CONTINUE
!
    CALL MPI_BCAST(BCHR, 1, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BARRIER                    (MPI_COMM_COMP, IRTN)
!
    IF (MYPE == 0) WRITE(LIST,*)'  READ UNIT NBC=', NBC
!-------------------------------------------------
! COMPUTE THE 1ST TIME FOR BOUNDARY CONDITION READ
!-------------------------------------------------
    NBOCO = NINT(BCHR * TSPH)
!
    IF (NTSD == 0) THEN
        IF (MYPE == 0 .AND. .NOT. NEST) THEN
            BACKSPACE NBC
            BACKSPACE NBC
            BACKSPACE NBC
            BACKSPACE NBC
            BACKSPACE NBC
            BACKSPACE NBC
            BACKSPACE NBC
            WRITE(LIST,*)'  BACKSPACE UNIT NBC=', NBC
        END IF
    END IF
!---------------------------------------- 
! SET ARRAYS CONTROLLING POST PROCESSING.
!----------------------------------------
    IF (MYPE == 0) THEN
        WRITE(LIST,*)'INIT:  READ IOUT,NSHDE,NTSD=', IOUT, NSHDE, NTSD
    END IF
!
    DO I=1,NMAP
        IOUT = I
        IF (ISHDE(I) >= NTSD) GOTO 220
    END DO
!
220 NSHDE = ISHDE(IOUT)
!
    IF (MYPE == 0) THEN
        WRITE(LIST,*)'INIT:  SET IOUT,NSHDE =', IOUT, NSHDE,' FOR ISHDE,NTSD=', ISHDE(IOUT), NTSD
    END IF
!----------------------------------------------------------------
! INITIALIZE PHYSICS VARIABLES IF STARTING THIS RUN FROM SCRATCH.
!----------------------------------------------------------------
    IF (NEST) THEN
        DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
            DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
!            
                LLMH = LMH(I,J)
!            
                IF (T(I,J,LLMH) == 0.) THEN
                    T(I,J,LLMH) = T(I,J,LLMH-1)
                END IF
!            
                TERM1 = -0.068283 / T(I,J,LLMH)
                PSHLTR(I,J) = (PD(I,J) + PT) * EXP(TERM1)
            END DO
        END DO
    END IF
!
    IF ( .NOT. RESTRT) THEN
        DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
            DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
                LLMH       = LMH(I,J)
                  PDSL(I,J) = PD(I,J) * RES(I,J)
                  PREC(I,J) =   0.
                ACPREC(I,J) =   0.
                CUPREC(I,J) =   0.
                    Z0(I,J) = SM(I,J) * Z0SEA + (1. - SM(I,J)) * (FIS(I,J) * FCM + Z0LAND)
                    QS(I,J) =   0.
                  AKMS(I,J) =   0.
                  AKHS(I,J) =   0.
                  TWBS(I,J) =   0.
                  QWBS(I,J) =   0.
                CLDEFI(I,J) =   1.
                  HTOP(I,J) = 100.
                  HBOT(I,J) =   0.
!--------------------------------------------------------------------------------------------------
! AT THIS POINT, WE MUST CALCULATE THE INITIAL POTENTIAL TEMPERATURE OF THE SURFACE AND OF THE 
! SUBGROUND.
! EXTRAPOLATE DOWN FOR INITIAL SURFACE POTENTIAL TEMPERATURE. ALSO DO THE SHELTER PRESSURE.
!--------------------------------------------------------------------------------------------------
                PM1      = PDSL(I,J) * AETA(LLMH) + PT
                APEM1    = (1.E5 / PM1) ** CAPA
                THS(I,J) = T(I,J,LLMH) * (1. + 0.608 * Q(I,J,LLMH)) * APEM1
                TSFCK    = T(I,J,LLMH) * (1. + 0.608 * Q(I,J,LLMH))
                PSFCK    = PD(I,J) + PT
!            
                IF (SM(I,J) < 0.5) THEN
                    QS(I,J) = PQ0 / PSFCK * EXP(A2 * (TSFCK - A3) / (TSFCK - A4))
                ELSE IF(SM(I,J) > 0.5) THEN
                    THS(I,J) = SST(I,J) * (1.E5 / (PD(I,J) + PT)) ** CAPA
                END IF
!            
                TERM1 = -0.068283 / T(I,J,LLMH)
                PSHLTR(I,J) = (PD(I,J) + PT) * EXP(TERM1)
!            
                USTAR(I,J) = 0.1
                 THZ0(I,J) = THS(I,J)
                  QZ0(I,J) =  QS(I,J)
                  UZ0(I,J) = 0.
                  VZ0(I,J) = 0.
            
            END DO
        END DO
!------------------------     
! INITIALIZE CLOUD FIELDS
!------------------------     
        DO K=1,LM
            DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
                DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
                    CWM(I,J,K) = 0.
                END DO
            END DO
        END DO
!---------------------------------------     
! INITIALIZE ACCUMULATOR ARRAYS TO ZERO.
!---------------------------------------     
        ARDSW = 0.0
        ARDLW = 0.0
        ASRFC = 0.0
        AVRAIN= 0.0
        AVCNVC= 0.0
!    
        DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
            DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
                ACFRCV(I,J) = 0.
                NCFRCV(I,J) = 0
                ACFRST(I,J) = 0.
                NCFRST(I,J) = 0
                ACSNOW(I,J) = 0.
                ACSNOM(I,J) = 0.
                SSROFF(I,J) = 0.
                BGROFF(I,J) = 0.
                 ALWIN(I,J) = 0.
                ALWOUT(I,J) = 0.
                ALWTOA(I,J) = 0.
                 ASWIN(I,J) = 0.
                ASWOUT(I,J) = 0.
                ASWTOA(I,J) = 0.
                SFCSHX(I,J) = 0.
                SFCLHX(I,J) = 0.
                SUBSHX(I,J) = 0.
                SNOPCX(I,J) = 0.
                SFCUVX(I,J) = 0.
                SFCEVP(I,J) = 0.
                POTEVP(I,J) = 0.
                POTFLX(I,J) = 0.
            END DO
        END DO
!--------------------------------------------------------    
! INITIALIZE SATURATION SPECIFIC HUMIDITY OVER THE WATER.
!--------------------------------------------------------     
        DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
            DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
                IF(SM(I,J) > 0.5)THEN
                    CLOGES = -CM1 / SST(I,J) - CM2 * ALOG10(SST(I,J)) + CM3
                    ESE    = 10. ** (CLOGES + 2.)
                    QS(I,J)= SM(I,J) * EPS * ESE / (PD(I,J) + PT - ESE * (1. - EPS))
                END IF
            END DO
        END DO
!--------------------------------------------------------    
! PAD GROUND WETNESS IF IT IS TOO SMALL.
! 
! DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
! DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
!    WET(I,J)=AMAX1(WET(I,J),EPSWET)
! END DO
! END DO
!
!--------------------------------------------------------------------------------------------------
! INITIALIZE TURBULENT KINETIC ENERGY (TKE) TO A SMALL VALUE (EPSQ2) ABOVE GROUND.  SET TKE TO ZERO 
! IN THE THE LOWEST MODEL LAYER.
! IN THE LOWEST TWO ATMOSPHERIC ETA LAYERS SET TKE TO A SMALL VALUE (Q2INI).
!--------------------------------------------------------------------------------------------------   
        DO K=1,LM1
            DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
                DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
                    Q2(I,J,K) = HTM(I,J,K+1) * HBM2(I,J) * EPSQ2
                END DO
            END DO
        END DO
!    
        DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
            DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
                Q2(I,J,LM)     = 0.
                LLMH           = LMH(I,J)
                Q2(I,J,LLMH-2) = HBM2(I,J) * Q2INI
                Q2(I,J,LLMH-1) = HBM2(I,J) * Q2INI
            END DO
        END DO
!-------------------------------------------------------    
! PAD ABOVE GROUND SPECIFIC HUMIDITY IF IT IS TOO SMALL.
! INITIALIZE LATENT HEATING ACCUMULATION ARRAYS.
!-------------------------------------------------------     
        DO K=1,LM
            DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
                DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
                    IF(Q(I,J,K) < EPSQ)Q(I,J,K) = EPSQ*HTM(I,J,K)
                    TRAIN(I,J,K) = 0.
                    TCUCN(I,J,K) = 0.
                END DO
            END DO
        END DO
!-------------------------------------------     
! END OF SCRATCH START INITIALIZATION BLOCK.
 !-------------------------------------------    
        IF (MYPE == 0) THEN
            WRITE(LIST,*)'INIT:  INITIALIZED ARRAYS FOR CLEAN START'
        END IF
    END IF
!--------------------------------------------------------------------------- 
! RESTART INITIALIZING. CHECK TO SEE IF WE NEED TO ZERO ACCUMULATION ARRAYS.
!--------------------------------------------------------------------------- 
    IF (RESTRT) THEN
!---------------------------     
! AVERAGE CLOUD AMOUNT ARRAY
!---------------------------    
        IF (NSTART == 0) THEN
            IF (MYPE == 0) WRITE(LIST,*)'  ZERO AVG CLD AMT ARRAY'
            DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
                DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
                    ACFRCV(I,J) = 0.
                    NCFRCV(I,J) = 0
                    ACFRST(I,J) = 0.
                    NCFRST(I,J) = 0
                END DO
            END DO
        END IF
!-------------------------------------------------
! GRID-SCALE AND CONVECTIVE LATENT HEATING ARRAYS.
!-------------------------------------------------    
        IF (NSTART == 0) THEN
            IF (MYPE == 0) THEN
                WRITE(LIST,*)'  ZERO ACCUM LATENT HEATING ARRAYS'
            END IF
!        
            AVRAIN = 0.
            AVCNVC = 0.
!
            DO K=1,LM
                DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
                    DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
                        TRAIN(I,J,K) = 0.
                        TCUCN(I,J,K) = 0.
                    END DO
                END DO
            END DO
        END IF
!----------------------------------------------
! IF THIS IS NOT A NESTED RUN, INITIALIZE TKE
!
! IF (.NOT.NEST) THEN
! DO L=1,LM
! DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
! DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
! Q2(I,J,L)=AMAX1(Q2(I,J,L)*HBM2(I,J),EPSQ2)
! END DO
! END DO
! END DO
! END IF
!   
! TOTAL AND CONVECTIVE PRECIPITATION ARRAYS.
! TOTAL SNOW AND SNOW MELT ARRAYS.
! STORM SURFACE AND BASE GROUND RUN OFF ARRAYS.
!----------------------------------------------    
        IF (NSTART == 0) THEN
            IF(MYPE == 0)WRITE(LIST,*)'  ZERO ACCUM PRECIP ARRAYS'
            DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
                DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
                    ACPREC(I,J) = 0.
                    CUPREC(I,J) = 0.
                    ACSNOW(I,J) = 0.
                    ACSNOM(I,J) = 0.
                    SSROFF(I,J) = 0.
                    BGROFF(I,J) = 0.
                END DO
            END DO
        END IF
!----------------------------    
! LONG WAVE RADIATION ARRAYS.
!----------------------------     
        IF (NSTART == 0) THEN
            IF (MYPE == 0)WRITE(LIST,*)'  ZERO ACCUM LW RADTN ARRAYS'
            ARDLW = 0.
            DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
                DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
                     ALWIN(I,J) = 0.
                    ALWOUT(I,J) = 0.
                    ALWTOA(I,J) = 0.
                END DO
            END DO
        END IF
!-----------------------------     
! SHORT WAVE RADIATION ARRAYS.
!-----------------------------    
        IF (NSTART == 0) THEN
            IF (MYPE == 0)WRITE(LIST,*)'  ZERO ACCUM SW RADTN ARRAYS'
            ARDSW=0.
            DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
                DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
                    ASWIN(I,J)  = 0.
                    ASWOUT(I,J) = 0.
                    ASWTOA(I,J) = 0.
                END DO
            END DO
        END IF
!----------------------------------------------    
! SURFACE SENSIBLE AND LATENT HEAT FLUX ARRAYS.
!----------------------------------------------     
        IF (NSTART == 0) THEN
            IF(MYPE == 0)WRITE(LIST,*)'  ZERO ACCUM SFC FLUX ARRAYS'
            ASRFC = 0.
            DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
                DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
                    SFCSHX(I,J) = 0.
                    SFCLHX(I,J) = 0.
                    SUBSHX(I,J) = 0.
                    SNOPCX(I,J) = 0.
                    SFCUVX(I,J) = 0.
                    SFCEVP(I,J) = 0.
                    POTEVP(I,J) = 0.
                    POTFLX(I,J) = 0.
                END DO
            END DO
        END IF
!-------------------------------------------------    
! END IF FOR RESTART FILE ACCUMULATION ZERO BLOCK.
!-------------------------------------------------    
        IF (MYPE == 0) THEN
            WRITE(LIST,*)'INIT:  INITIALIZED ARRAYS FOR RESTART START'
        END IF
    END IF
!--------------------------- 
! INITIALIZE CLOUD CONSTANTS
!--------------------------- 
    DO 350 J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
        DO 350 I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
            U00(I,J) = .95
350 END DO
!--------------------------------------------------------------------------------------------- 
! FLAG FOR INITIALIZING ARRAYS, LOOKUP TABLES AND CONSTANTS USED IN MICROPHYSICS AND RADIATION
!--------------------------------------------------------------------------------------------- 
    MICRO_START = .TRUE. 
!------------------------------------------------------------  
! FLAG FOR INITIALIZING CONVECTIVE CLOUD ARRAYS FOR RADIATION
!------------------------------------------------------------
    CURAD = .FALSE. 
!
    DO 355 K=1,2*LM
        IF (K >= LM-10 .AND. K <= LM) THEN
            UL(K) = 0.1 * FLOAT(K-LM+10)
        ELSE
            UL(K) = 0.
        END IF
355 END DO
!---------------------------------------- 
! SET INDEX ARRAYS FOR UPSTREAM ADVECTION
!----------------------------------------
    KNT = 0
!
    DO J=3,5
        KNT       = KNT + 1
        IHLA(KNT) = 2
        IHHA(KNT) = IM  - 1 - MOD(J+1,2)
        IVLA(KNT) = 2
        IVHA(KNT) = IM  - 1 - MOD(J  ,2)
         JRA(KNT) = J
    END DO
!
    DO J=JM-4,JM-2
        KNT       = KNT + 1
        IHLA(KNT) = 2
        IHHA(KNT) = IM  - 1 - MOD(J+1,2)
        IVLA(KNT) = 2
        IVHA(KNT) = IM  - 1 - MOD(J  ,2)
         JRA(KNT) = J
    END DO
!
    DO J=6,JM-5
        KNT       = KNT + 1
        IHLA(KNT) = 2
        IHHA(KNT) = 2   + MOD(J  ,2)
        IVLA(KNT) = 2 
        IVHA(KNT) = 2   + MOD(J+1,2)
         JRA(KNT) = J
    END DO
!
    DO J=6,JM-5
        KNT       = KNT + 1
        IHLA(KNT) = IM  - 2
        IHHA(KNT) = IM  - 2 + MOD(J  ,2)
        IVLA(KNT) = IM  - 2
        IVHA(KNT) = IM  - 2 + MOD(J+1,2)
         JRA(KNT) = J
    END DO
!-------------------------------------------------
! SET ZERO-VALUE FOR SOME OUTPUT DIAGNOSTIC ARRAYS
!-------------------------------------------------
    IF (NSTART == 0) THEN
!    
        DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
            DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
                PCTSNO(I,J) = -999.0
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
                    GRNFLX(I,J) =   0.
                    SUBSHX(I,J) =   0.0
                    ACSNOW(I,J) =   0.0
                    ACSNOM(I,J) =   0.0
                    SNOPCX(I,J) =   0.0
                       CMC(I,J) =   0.0
                       SNO(I,J) =   0.0
!
                    DO NS=1,NSOIL
                         SMC(I,J,NS) =  1.0
                        SH2O(I,J,NS) =  1.0
                         STC(I,J,NS) =273.16
                    END DO
                END IF
!            
            END DO
        END DO
!    
        APHTIM = 0.0
        ARATIM = 0.0
        ACUTIM = 0.0
!    
    END IF
!--------------------------------------------------------------------------------------------------
! INITIALIZE RADTN VARIABLES
! CALCULATE THE NUMBER OF STEPS AT EACH POINT.
! THE ARRAY 'LVL' WILL COORDINATE VERTICAL LOCATIONS BETWEEN THE LIFTED WORKING ARRAYS AND THE 
! FUNDAMENTAL MODEL ARRAYS.
! LVL HOLDS THE HEIGHT (IN MODEL LAYERS) OF THE TOPOGRAPHY AT EACH GRID POINT.
!--------------------------------------------------------------------------------------------------
    DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
        DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
            LVL(I,J) = LM - LMH(I,J)
        END DO
    END DO
!------------------------------------------------------------------------ 
! DETERMINE MODEL LAYER LIMITS FOR HIGH(3), MIDDLE(2), AND LOW(1) CLOUDS.
! ALSO FIND MODEL LAYER THAT IS JUST BELOW (HEIGHT-WISE) 400 MB. (K400)
!------------------------------------------------------------------------ 
    K400 = 0
    PSUM = PT
    SLPM = 101325.
    PDIF = SLPM - PT
!
    DO K=1,LM
        PSUM = PSUM + DETA(L) * PDIF
        IF (LTOP(3) == 0) THEN
             IF (PSUM > PHITP) LTOP(3) = K
        ELSE IF (LTOP(2) == 0) THEN
             IF (PSUM > PMDHI) LTOP(2) = K
        ELSE IF (K400 == 0) THEN
             IF (PSUM > P400)  K400    = K
        ELSE IF (LTOP(1) == 0) THEN
             IF (PSUM > PLOMD) LTOP(1) = K
        END IF
    END DO
!----------------------------------------------------
! CALL GRADFS ONCE TO CALC. CONSTANTS AND GET O3 DATA
!----------------------------------------------------
    KCCO2 = 0
!------------------------------------------------------------
! CALCULATE THE MIDLAYER PRESSURES IN THE STANDARD ATMOSPHERE
!------------------------------------------------------------
    PSS  = 101325.
    PDIF = PSS - PT
!
    DO K=1,LM1
        PHALF(K+1) = AETA(K) * PDIF + PT
    END DO
!
    PHALF(1)   = 0.
    PHALF(LP1) = PSS
!
    CALL GRADFS(PHALF, KCCO2, NFILE)
!------------------------------------------------------------
! CALL SOLARD TO CALCULATE NON-DIMENSIONAL SUN-EARTH DISTANCE
!------------------------------------------------------------
    IF (MYPE == 0) CALL SOLARD(RAD1)
    CALL MPI_BCAST(RAD1, 1, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
!------------------------------------------------------------------------------ 
! CALL ZENITH SIMPLY TO GET THE DAY OF THE YEAR FOR THE SETUP OF THE OZONE DATA
!------------------------------------------------------------------------------
    TIME = (NTSD - 1) * DT
!
    CALL ZENITH(TIME, DAYI, HOUR)
!
    ADDL = 0.
    IF (MOD(IDAT(3),4) == 0) ADDL = 1.
    RANG  = PI2 * (DAYI - RLAG) / (365.25 + ADDL)
    RSIN1 = SIN(RANG)
    RCOS1 = COS(RANG)
    RCOS2 = COS(2. * RANG)
!
    CALL O3CLIM
!------------------------------------------------- 
! SOME INITIAL VALUES RELATED TO TURBULENCE SCHEME
!------------------------------------------------- 
    DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
        DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
!------------------------------------------------         
! TRY A SIMPLE LINEAR INTERP TO GET 2/10 M VALUES
!------------------------------------------------        
            PDSL(I,J) =  PD(I,J) * RES(I,J)
            LMHK      = LMH(I,J)
            LMVK      = LMV(I,J)
!
            ULM = U(I,J,LMVK)
            VLM = V(I,J,LMVK)
            TLM = T(I,J,LMHK)
            QLM = Q(I,J,LMHK)
!
            PLM2   = PDSL(I,J) * AETA(LMHK) + PT
!
            APELM   = (1.0E5 / PLM2)         ** CAPA
            APELMNW = (1.0E5 / PSHLTR(I,J)) ** CAPA
            EXNERR  = (PSHLTR(I,J) * 1.E-5) ** CAPA
!
            THLM = TLM * APELM
            DPLM = PDSL(I,J) * DETA(LMHK) * 0.5
            DZLM = 287.04 * DPLM * TLM * (1. + 0.608 * QLM) / (9.801 * PLM2)
            FAC1 = 10. / DZLM
            FAC2 = (DZLM - 10.) / DZLM
!
            IF (DZLM <= 10.) THEN
                FAC1 = 1.
                FAC2 = 0.
            END IF
!        
            IF ( .NOT. RESTRT) THEN
                TH10(I,J) = FAC2 * THS(I,J) + FAC1 * THLM
                 Q10(I,J) = FAC2 *  QS(I,J) + FAC1 * QLM
                 U10(I,J) = ULM
                 V10(I,J) = VLM
            END IF
!        
            FAC1 = 2. / DZLM
            FAC2 = (DZLM - 2.) / DZLM
!
            IF (DZLM <= 2.) THEN
                FAC1 = 1.
                FAC2 = 0.
            END IF
!        
            IF ( .NOT. RESTRT .OR. NEST) THEN
                TSHLTR(I,J) = 0.1 * THS(I,J) + 0.9 * THLM
                QSHLTR(I,J) = QLM
            END IF
!--------- 
! SM V100M
!---------
            FAC1 = 100. / DZLM
            FAC2 = (DZLM - 100.) / DZLM
!
            IF (DZLM <= 100.) THEN
                FAC1 = 1.
                FAC2 = 0.
            END IF
!             							      
            IF ( .NOT. RESTRT) THEN
                TH100(I,J) = FAC2 * THS(I,J) + FAC1 * THLM
                 Q100(I,J) = FAC2 *  QS(I,J) + FAC1 * QLM
                 U100(I,J) = ULM
                 V100(I,J) = VLM
            END IF
!---------
! SM V100M
!---------
!
!----------------------------------------------------------------------------------------
! NEED TO CONVERT TO THETA IF IS THE RESTART CASE AS CHKOUT.F WILL CONVERT TO TEMPERATURE
!----------------------------------------------------------------------------------------
            IF (RESTRT) THEN
                TSHLTR(I,J) = TSHLTR(I,J) * APELMNW
            END IF
        END DO
    END DO
!------------------------------------- 
! INITIALIZE NONHYDROSTATIC QUANTITIES
!------------------------------------- 
    DO K=1,LM
!    
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                DWDT(I,J,K) = 1.
            END DO
        END DO
    END DO
!
    IF (SIGMA) THEN
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                PDSL(I,J) = PD(I,J)
            END DO
        END DO
    ELSE
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                PDSL(I,J) = RES(I,J) * PD(I,J)
            END DO
        END DO
    END IF
!
    DO K=1,LP1
!    
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                PINT(I,J,K) = PDSL(I,J) * ETA(K) + PT
                   Z(I,J,K) = PINT(I,J,K)
            END DO
        END DO
    END DO
!------------------------ 
! END OF SUBROUTINE INIT.
!------------------------ 
    IF (MYPE == 0) THEN
        WRITE(LIST,*)'INIT:  EXIT INIT AND START MODEL INTEGRATION'
        WRITE(LIST,*)' '
    END IF
!
    RETURN
!
    END SUBROUTINE INIT
!
!
!
    BLOCK DATA CLOUD
!--------------------------------------------------------------------------------------------------
    USE PARMETA
    USE RD1TIM
!
    END BLOCK DATA CLOUD
