    SUBROUTINE CUCNVC_SHALLOW
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE CUCNVC_SHALLOW
!>
!> SUBPROGRAM: CUCNVC_SHALLOW - CONVECTIVE PRECIPITATION PARAMETERIZATION
!> PROGRAMMER: JANJIC
!> ORG: W/NP2 
!> DATE: 93-11-02
!>
!> ABSTRACT:
!> CUCNVC CALCULATES THE SUB-GRID SCALE CONVECTION INCLUDING DEEP AND SHALLOW CONVECTIVE CLOUDS 
!> FOLLOWING THE SCHEME DESCRIBED BY JANJIC (1994) BUT WITH SIGNIFICANT MODIFICATIONS.
!> IN ADDITION, THE LATENT HEAT RELEASE AND MOISTURE CHANGE DUE TO PRECIPITATING AND NON -
!> PRECIPITATING CLOUDS ARE COMPUTED.
!>
!> PROGRAM HISTORY LOG:
!> 87-09-??  JANJIC     - ORIGINATOR
!> 90-11-21  JANJIC     - TWO SETS OF DSP PROFILES (FAST AND SLOW) REPLACE THE ORIGINAL ONE SET
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 96-03-28  BLACK      - ADDED EXTERNAL EDGE
!> 98-11-02  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
!> 20-04-30  GUSTAVO    - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY 
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
!> USE MODULES: ACMCLH
!>              CNVCLD
!>              CTLBLK
!>              CUPARM
!>              DYNAM
!>              F77KINDS
!>              GLB_TABLE
!>              INDX
!>              LOOPS
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              PARMETA
!>              PARM_TBL
!>              PHYS
!>              PPTASM
!>              PVRBLS
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>  
!> DRIVER     : EBU
!>              NEWFLT
!>
!> CALLS      : TTBLEX
!>              ZERO2
!>-------------------------------------------------------------------------------------------------- 
!
!-----------------------------------------------------------------------------------------
! REFERENCES:                                                 
!
! JANJIC, Z.I., 1994:  THE STEP-MOUNTAIN ETA COORDINATE MODEL:
! FURTHER DEVELOPMENTS OF THE CONVECTION, VISCOUS SUBLAYER AND TURBULENCE CLOSURE SCHEMES.  
! MONTHLY WEATHER REVIEW, VOL. 122, 927-945.
!-----------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! WARNING: THIS SUBROUTINE WILL NOT WORK IF (LM .LT. 12);
! MUST BE CALLED IN THE SAME STEP WITH PROFQ2 BECAUSE PROFQ DEFINES APE;
!-----------------------------------------------------------------------
    USE ACMCLH
    USE CNVCLD
    USE CTLBLK
    USE CUPARM
    USE DYNAM
    USE F77KINDS
    USE GLB_TABLE
    USE INDX
    USE LOOPS
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE PARMETA
    USE PARM_TBL
    USE PHYS
    USE PPTASM
    USE PVRBLS
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & UNIS    , UNIL    , OCT90
!
    NAMELIST /CUPARMDATA/                                                                         &
    & STABS   , STABD   , STABFC  , DTTOP   ,                                                     &
    & RHF     , EPSUP   , EPSDN   , EPSTH   ,                                                     &
    & PBM     , PQM     , PNO     , PONE    , PSH     ,                                           &
    & PFRZ    , PSHU    ,                                                                         &
    & UNIS    , UNIL    , OCT90   ,                                                                &
    & FSS     , EFIMN   , EFMNT   , FCC     ,                                                     &
    & DSPBFL  , DSP0FL  , DSPTFL  , FSL     ,                                                     &
    & DSPBFS  , DSP0FS  , DSPTFS  ,                                                               &
    & TREL    , EPSNTP  , EFIFC   ,                                                               &
    & DSPC    , EPSP    ,                                                                         &
    & STEFI
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
!
!------------------------------
! INSTABILITY FOR TOO LARGE LSH
!------------------------------
    INTEGER(KIND=I4KIND), PARAMETER :: KSMUD = 0
    INTEGER(KIND=I4KIND), PARAMETER :: NROW  = 0
!--------------------------------------------------------------
! INSTABILITY FOR TOO LARGE LSH
!
! LP1=LM+1,LM1=LM-1,LNO=1,LSH=LM/3-1,LSHU=LM/2-1,LQM=LM/5,KBM=3
! LP1=LM+1,LM1=LM-1,LNO=1,LSH=LM/3  ,LSHU=LM/2-1,LQM=LM/5,KBM=3
! LP1=LM+1,LM1=LM-1,LNO=3,LSH=LM/3  ,LSHU=LM/2-1,LQM=LM/5,KBM=3
! LP1=LM+1,LM1=LM-1,LNO=2,LSH=LM/3-1,LSHU=LM/2-1,LQM=LM/5,KBM=3
! LP1=LM+1,LM1=LM-1,LNO=2,LSH=LM/3-2,LSHU=LM/2-1,LQM=LM/5,KBM=3
!--------------------------------------------------------------
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM     = IM    * JM    - JM / 2
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM_LOC = IDIM2 * JDIM2
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                                         ::&
    & TREFK   , QREFK   , PK      , APEK    , TK      ,                                           &
    & THSK    , PSK     , APESK   , QK      , THERK   ,                                           &
    & THVREF  , THEVRF  , THVMOD  , DIFT    , DIFQ    ,                                           &
    & QSATK   , FPK
!
    INTEGER(KIND=I4KIND), DIMENSION(LM)                                                         ::&
    & NTOPD   , NBOTD   , NTOPS   , NBOTS   ,                                                     &
    & NDPTHD  , NDPTHS
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & LTOP    , LBOT    ,                                                                         &
    & IPTB    , ITHTB
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & PTOP    , PBOT    ,                                                                         &
    & PDSL    , APEBT   ,                                                                         &
    & TBT     , Q2BT    ,                                                                         &
    & QQ      , PP      ,                                                                         &
    & PSP     , THBT    ,                                                                         &
    & THESP   , P       ,                                                                         &
    & BTH     , STH     ,                                                                         &
    & T00     , T10     ,                                                                         &
    & T01     , T11     ,                                                                         &
    & WF1     , WF2     ,                                                                         &
    & WF3     , WF4     ,                                                                         &
    & PRECOL
!
    INTEGER(KIND=I4KIND), DIMENSION(IMJM_LOC)                                                   ::&
    &  IBUOY  , JBUOY   ,                                                                         &
    &  IDEEP  , JDEEP   ,                                                                         &
    &  ISHAL  , JSHAL   ,                                                                         &
    &  ILRES  , JLRES   ,                                                                         &
    &  IHRES  , JHRES
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & APE     ,                                                                                   &
    & TREF    ,                                                                                   &
    & TMOD    ,                                                                                   &
    & QMOD
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & DSPB    , DSP0    ,                                                                         &
    & DSPT    ,                                                                                   &
    & TL      ,   QL    ,                                                                         &
    & TNE     ,  TSE    ,                                                                         &
    & QNE     ,  QSE
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & STABS   , STABD   , STABFC  , DTTOP   , RHF     , EPSUP   , EPSDN   , EPSTH   , PBM     ,   &
    & PSHU    , FSS     , EFIMN   , EFMNT   , FCC     , DSPBFL  , DSP0FL  , DSPTFL  , FSL     ,   &
    & DSPBFS  , DSP0FS  , DSPTFS  , TREL    , PQM     , PNO     , PONE    , PFRZ    , PSH     ,   &
    & EPSNTP  , EFIFC   , DSPC    , EPSP    ,                                                     &
    & STEFI   , FCB     , DSPBSL  , DSP0SL  ,                                                     &
    & DSPTSL  , DSPBSS  , DSP0SS  , DSPTSS  ,                                                     &
    & AVGEFI  , SLOPBL  , SLOP0L  , SLOPTL  ,                                                     &
    & SLOPBS  , SLOP0S  , SLOPTS  , SLOPE   ,                                                     &
    & A23M4L  , ELOCP   , CPRLG   , RCP     ,                                                     &
    & DTCNVC  , RDTCNVC , TAUK    , APESTS  , PKL     ,                                           &
    & PSFCK   , QBT     , TTHBT   , TTH     , QQ1     , BQS00K  ,                                 &
    & SQS00K  , BQS10K  , SQS10K  , BQ      , SQ      , TQ      , PP1     , P00K    ,             &
    & P10K    , P01K    , P11K    , TPSP    ,                                                     &
    & APESP   , TTHES   , AETAL   , PRESK   , EFI     , SMK     , RSMK    ,                       &
    & PSHNEW  , PSFCIJ  , DEPMIN  , DEPTH   ,                                                     &
    & DSPBK   , DSP0K   , DSPTK   , TKL     , QKL     , APEKL   , PKB     , PKT     ,             &
    & PK0     , TREFKX  , THERKX  , APEKXX  ,                                                     &
    & THERKY  , APEKXY  , STABDL  , RDP0T   ,                                                     &
    & DTHEM   , DEPWL   , DSP     , SUMDP   , SUMDE   , HCORR   , TSKL    , DHDT    ,             &
    & THSKL   , DENTPY  , AVRGT   , PRECK   , DIFTL   , DIFQL   , AVRGTL  , PBTK    ,             &
    & PTPK    , RHH     , RHMAX   , RHL     , DRHDP   , DRHEAT  , FEFI    , PZ0     ,             &
    & THVMKL  , THTPK   , TTHK    , QQK     , BQK     , SQK     , TQK     , PPK     ,             &
    & PART1   , PART2   , PART3   , DPMIX   , SMIX    , PKXXXX  , SUMDT   , RDPSUM  ,             &
    & TCORR   , PKXXXY  , TRFKL   , PSUM    , QSUM    , POTSUM  , QOTSUM  , OTSUM   ,             &
    & DST     , FPTK    , DPKL    , RTBAR   , ROTSUM  , DSTQ    , DEN     , DQREF   ,             &
    & QRFTP   , QRFKL   , QNEW    , DTDETA  , CUTOP   , CUBOT
! 
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K       , KB      , NDSTN   , NDSTP   , LLMH    ,                       &
    & LMHK    , IQTB    , IQ      , IT      , NSHAL   , NEGDS   , NSMUD   ,                       &
    & LMHIJ   , KNUML   , KNUMH   , KROW    , KS      ,                                           &
    & KOUNT   , KHDEEP  , N       , LTPK    , ITTB    , ITTBK   ,                                 &
    & LBTK    , LB      , LTP1    , LBM1    ,                                                     &
    & L0      , L0M1    , ITER    , LCOR    ,                                                     &
    & LQM     , LTSH    , LSHU    , LTP2    ,                                                     &
    & NDEEP   , NDEPTH  , NNEG    , KHSHAL  
!
    OPEN(11,FILE='CUPARMDATA.dat',STATUS='OLD')
    REWIND 11
    READ(11,CUPARMDATA)
    CLOSE(11)
!
    FCB = 1. - FCC
    DSPBSL = DSPBFL * FSL
    DSP0SL = DSP0FL * FSL
    DSPTSL = DSPTFL * FSL
    DSPBSS = DSPBFS * FSS
    DSP0SS = DSP0FS * FSS
    DSPTSS = DSPTFS * FSS
!
    AVGEFI = (EFIMN + 1.) * .5
!
    SLOPBL = (DSPBFL-DSPBSL) / (1.-EFIMN)
    SLOP0L = (DSP0FL-DSP0SL) / (1.-EFIMN)
    SLOPTL = (DSPTFL-DSPTSL) / (1.-EFIMN)
    SLOPBS = (DSPBFS-DSPBSS) / (1.-EFIMN)
    SLOP0S = (DSP0FS-DSP0SS) / (1.-EFIMN)
    SLOPTS = (DSPTFS-DSPTSS) / (1.-EFIMN)
!
    SLOPE  = (1. - EFMNT) / (1. - EFIMN)
    A23M4L = A2 * (A3 - A4) * ELWV
    ELOCP  = ELIVW / CP
    CPRLG  = CP / (ROW * G * ELWV)
    RCP    = 1. / CP
!
    CALL ZERO2(DSP0)
    CALL ZERO2(DSPB)
    CALL ZERO2(DSPT)
    CALL ZERO2(PSP)
!
    AVCNVC  = AVCNVC + 1.
    ACUTIM  = ACUTIM + 1.
    DTCNVC  = NCNVC * DT
    RDTCNVC = 1. / DTCNVC
    TAUK    = DTCNVC / TREL
!----------------------------------------------
! POSSIBLE FUTURE CHANGE FOR SHALLOW CONVECTION
! TAUKSC=DTCNVC/(5.*TREL)
!----------------------------------------------
!
!------------
! DIAGNOSTICS
!------------
    DO K=1,LM
         NTOPD(K) = 0
         NBOTD(K) = 0
         NTOPS(K) = 0
         NBOTS(K) = 0
        NDPTHS(K) = 0
        NDPTHD(K) = 0
    END DO
!------------
!PREPARATIONS
!------------
!
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
    DO 120 J=MYJS,MYJE
        DO 120 I=MYIS,MYIE
              LBOT(I,J) = LMH(I,J)
             THESP(I,J) = 0.
             PDSL (I,J) = RES(I,J) * PD(I,J)
            PRECOL(I,J) = 0.
            TREF(I,J,1) = T(I,J,1)
120 END DO
!---------------------------------------
! PADDING SPECIFIC HUMIDITY IF TOO SMALL
! RESTORE APE TO SCRATCH ARRAY
!---------------------------------------
!
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (APESTS)
!
    DO 130 K=1,LM
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                APESTS     = PDSL(I,J) * AETA(K) + PT
                APE(I,J,K) = (1.E5 / APESTS) ** CAPA
                IF (Q(I,J,K) < EPSQ) Q(I,J,K) = HTM(I,J,K) * EPSQ
            END DO
        END DO
130 END DO
!----------------------------------
! SEARCH FOR MAXIMUM BUOYANCY LEVEL
!----------------------------------
    DO 170 KB=1,LM
!---------------------------------------
! TRIAL MAXIMUM BUOYANCY LEVEL VARIABLES
!---------------------------------------
!
!------- 
! OPENMP
!------- 
!$omp parallel do
!$omp private (APESP    , BQ      , BQS00K  , BQS10K  , IQ      , IT      , ITTB    , ITTBK   ,   &
!$omp          IQTB     , LMHK    , P00K    , P01K    , P10K    , P11K    , PKL     , PP1     ,   &
!$omp          PSFCK    , QBT     , QQ1     , SQ      , SQS00K  , SQS10K  , TPSP    , TQ      ,   &
!$omp          TTH      , TTHBT   , TTHES)
!
        DO 155 J=MYJS2,MYJE2
            DO 150 I=MYIS1,MYIE1
!            
                PKL   = AETA(KB)   * PDSL(I,J) + PT
                LMHK  = LMH(I,J)
                PSFCK = AETA(LMHK) * PDSL(I,J) + PT
!--------------------------------------------------------------------------------------------------
! NOW SEARCHING OVER A SCALED DEPTH IN FINDING THE PARCEL WITH THE MAX THETA-E INSTEAD OF THE OLD 
! 130 MB
!--------------------------------------------------------------------------------------------------
                IF (KB <= LMHK .AND. PKL >= 0.80*PSFCK) THEN
                    QBT   = Q(I,J,KB)
                    TTHBT = T(I,J,KB) * APE(I,J,KB)
                    TTH   = (TTHBT-THL) * RDTH
                    QQ1   = TTH - AINT(TTH)
                    ITTB  = INT(TTH) + 1
!---------------------------------
! KEEPING INDICES WITHIN THE TABLE
!---------------------------------
                    IF (ITTB < 1) THEN
                        ITTB = 1
                        QQ1  = 0.
                    END IF
!
                    IF (ITTB >= JTB) THEN
                        ITTB = JTB - 1
                        QQ1  = 0.
                    END IF
!
                    CONTINUE
!-------------------------------------------
! BASE AND SCALING FACTOR FOR SPEC. HUMIDITY
!-------------------------------------------
                    ITTBK  = ITTB
                    BQS00K = QS0(ITTBK)
                    SQS00K = SQS(ITTBK)
                    BQS10K = QS0(ITTBK+1)
                    SQS10K = SQS(ITTBK+1)
!---------------------------------------
! SCALING SPEC. HUMIDITY AND TABLE INDEX
!---------------------------------------
                    BQ  = (BQS10K-BQS00K) * QQ1 + BQS00K
                    SQ  = (SQS10K-SQS00K) * QQ1 + SQS00K
                    TQ  = (QBT-BQ) / SQ * RDQ
                    PP1 = TQ - AINT(TQ)
                    IQTB = INT(TQ) + 1
!---------------------------------
! KEEPING INDICES WITHIN THE TABLE
!---------------------------------
                    IF (IQTB < 1) THEN
                        IQTB = 1
                        PP1  = 0.
                    END IF
!
                    IF (IQTB >= ITB) THEN
                        IQTB = ITB - 1 
                        PP1  = 0.
                    END IF
!---------------------------------------------------
! SATURATION PRESSURE AT FOUR SURROUNDING TABLE PTS.
!---------------------------------------------------
                    IQ = IQTB
                    IT = ITTB
                    P00K = PTBL(IQ  ,IT  )
                    P10K = PTBL(IQ+1,IT  )
                    P01K = PTBL(IQ  ,IT+1)
                    P11K = PTBL(IQ+1,IT+1)
!-----------------------------------------
! SATURATION POINT VARIABLES AT THE BOTTOM
!-----------------------------------------
                    TPSP  = P00K + (P10K-P00K) * PP1 + (P01K-P00K) * QQ1 + (P00K-P10K-P01K+P11K)  &
    &                     * PP1 * QQ1
                    APESP = (1.E5/TPSP) ** CAPA
                    TTHES = TTHBT * EXP(ELOCP * QBT * APESP / TTHBT)
!---------------------------
! CHECK FOR MAXIMUM BUOYANCY
!---------------------------
                    IF (TTHES > THESP(I,J)) THEN
                          PSP(I,J) = TPSP
                         THBT(I,J) = TTHBT
                        THESP(I,J) = TTHES
                    END IF
!
                END IF
!
        150 END DO
    155 END DO
170 END DO
!------------------------------------------------
! CHOOSE CLOUD BASE AS MODEL LEVEL JUST BELOW PSP
!------------------------------------------------
    DO 240 K=1,LM1
        AETAL = AETA(K)
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO J=MYJS2,MYJE2
            DO I=MYIS,MYIE
                P(I,J) = PDSL(I,J) * AETAL + PT
                IF (P(I,J) < PSP(I,J) .AND. P(I,J) >= PQM) LBOT(I,J) = K + 1
            END DO
        END DO
!  
240 END DO
!--------------------------------------------------------------
! WARNING: LBOT MUST NOT BE GT LMH(I,J)-1 IN SHALLOW CONVECTION
! MAKE SURE CLOUD BASE IS AT LEAST PONE ABOVE THE SURFACE
!--------------------------------------------------------------
!
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (LMHIJ  , PSFCK)
!
    DO 250 J=MYJS2,MYJE2
        DO 250 I=MYIS,MYIE
            LMHIJ = LMH(I,J)
            PBOT(I,J) = AETA(LBOT(I,J)) * PDSL(I,J) + PT
            PSFCK     = AETA(LMHIJ)     * PDSL(I,J) + PT
            IF (PBOT(I,J) >= PSFCK-PONE .OR. LBOT(I,J) >= LMHIJ) THEN
!            
                DO K=1,LMHIJ-1
                    P(I,J) = AETA(K) * PDSL(I,J) + PT
                    IF (P(I,J) < PSFCK-PONE) LBOT(I,J) = K
                END DO
!            
                PBOT(I,J) = AETA(LBOT(I,J)) * PDSL(I,J) + PT
            END IF
250 END DO
!----------------------   
! CLOUD TOP COMPUTATION
!----------------------
!
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            LTOP(I,J) = LBOT(I,J)
            PTOP(I,J) = PBOT(I,J)
        END DO
    END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!$omp private (IHRES    , ILRES   , IPTB    , ITHTB   , JHRES   , JLRES   , KNUMH   , KNUML   ,   &
!$omp          PP       , PRESK   , QQ      )
!
    DO 290 K=LM,1,-1
!------------------------------------    
! SCALING PRESSURE AND TT TABLE INDEX
!------------------------------------
        KNUML = 0
        KNUMH = 0
!    
        DO 270 J=MYJS2,MYJE2
            DO 270 I=MYIS1,MYIE1
                PRESK = PDSL(I,J) * AETA(K) + PT
                IF (PRESK < PLQ) THEN
                    KNUML = KNUML + 1
                    ILRES(KNUML) = I
                    JLRES(KNUML) = J
                ELSE
                    KNUMH = KNUMH + 1
                    IHRES(KNUMH) = I
                    JHRES(KNUMH) = J
                END IF
    270 END DO
!---------------------------------------------------------------
! COMPUTE PARCEL TEMPERATURE ALONG MOIST ADIABAT FOR PRESSURE<PL
!---------------------------------------------------------------
        IF (KNUML > 0) THEN
            CALL TTBLEX(TREF(IDIM1,JDIM1,K)                                                       &
    &        ,          TTBL, ITB, JTB, KNUML, ILRES, JLRES, PDSL                                 &
    &        ,          AETA(K)                                                                   &
    &        ,          HTM(IDIM1,JDIM1,K)                                                        &
    &        ,          PT, PL                                                                    &
    &        ,          QQ(IDIM1,JDIM1), PP(IDIM1,JDIM1)                                          &
    &        ,          RDP, THE0, STHE, RDTHE                                                    &
    &        ,          THESP(IDIM1,JDIM1), IPTB(IDIM1,JDIM1)                                     &
    &        ,          ITHTB(IDIM1,JDIM1))
        END IF
!---------------------------------------------------------------
! COMPUTE PARCEL TEMPERATURE ALONG MOIST ADIABAT FOR PRESSURE>PL
!---------------------------------------------------------------
        IF (KNUMH > 0) THEN
            CALL TTBLEX(TREF(IDIM1,JDIM1,K)                                                       &
    &       ,           TTBLQ             , ITBQ             , JTBQ , KNUMH                       &
    &       ,           IHRES             , JHRES            , PDSL                               & 
    &       ,           AETA(K)                                                                   &
    &       ,           HTM(IDIM1,JDIM1,K)                                                        &
    &       ,           PT                , PLQ                                                   &
    &       ,              QQ(IDIM1,JDIM1),   PP(IDIM1,JDIM1)                                     &
    &       ,           RDPQ              , THE0Q            , STHEQ, RDTHEQ                      &
    &       ,           THESP(IDIM1,JDIM1), IPTB(IDIM1,JDIM1)                                     &
    &       ,           ITHTB(IDIM1,JDIM1))
        END IF
290 END DO
!---------------    
! BUOYANCY CHECK
!---------------
    DO 295 K=LM,1,-1
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO J=MYJS2,MYJE2
            DO I=MYIS1,MYIE1
                IF (TREF(I,J,K) > T(I,J,K)-DTTOP) LTOP(I,J) = K
            END DO
        END DO
!
295 END DO
!-------------------
! CLOUD TOP PRESSURE
!-------------------
!
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
    DO J=MYJS2,MYJE2
        DO I=MYIS1,MYIE1
            PTOP(I,J) = AETA(LTOP(I,J)) * PDSL(I,J) + PT
        END DO
    END DO
!----------------------------------
! DEFINE AND SMOOTH DSPS AND CLDEFI
! UNIFIED OR SEPARATE LAND/SEA CONV 
!----------------------------------
    IF (UNIS) THEN
!
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (EFI)
!
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                EFI = CLDEFI(I,J)
                DSPB(I,J) = (EFI-EFIMN) * SLOPBS + DSPBSS
                DSP0(I,J) = (EFI-EFIMN) * SLOP0S + DSP0SS
                DSPT(I,J) = (EFI-EFIMN) * SLOPTS + DSPTSS
            END DO
        END DO
!
    ELSE IF ( .NOT. UNIL) THEN
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (EFI)
!
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                EFI=CLDEFI(I,J)
                DSPB(I,J) = ((EFI-EFIMN)*SLOPBS+DSPBSS) * SM(I,J) + ((EFI-EFIMN)*SLOPBL+DSPBSL)   &
    &                     * (1.-SM(I,J))
!
                DSP0(I,J) = ((EFI-EFIMN)*SLOP0S+DSP0SS) * SM(I,J) + ((EFI-EFIMN)*SLOP0L+DSP0SL)   &
    &                     * (1.-SM(I,J))
!
                DSPT(I,J) = ((EFI-EFIMN)*SLOPTS+DSPTSS) * SM(I,J) + ((EFI-EFIMN)*SLOPTL+DSPTSL)   &
    &                     * (1.-SM(I,J)) 
           END DO
        END DO
!
    ELSE
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (EFI)
!
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                EFI=CLDEFI(I,J)
                DSPB(I,J) = ((EFI-EFIMN) * SLOPBL + DSPBSL)
                DSP0(I,J) = ((EFI-EFIMN) * SLOP0L + DSP0SL)
                DSPT(I,J) = ((EFI-EFIMN) * SLOPTL + DSPTSL)
            END DO
        END DO
!
    END IF
!------------------------------------------------
! EXTENDING SEA STRUCTURES INLAND ALONG COASTLINE
!------------------------------------------------
    IF (NROW > 0 .AND. .NOT. UNIS .AND. .NOT. UNIL) THEN
!------- 
! OPENMP
!------- 
!
!$omp parallel do
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                WF1(I,J) = 0.
                WF2(I,J) = 0.
                WF3(I,J) = 0.
                WF4(I,J) = 0.
            END DO
        END DO
!    
        KROW = NROW
!    
        DO 350 KOUNT=1,KROW
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
            DO 345 J=MYJS2,MYJE2
                DO 345 I=MYIS1,MYIE1
                    WF1(I,J) =   (DSPB(I+IHW(J),J-1) +   DSPB(I+IHE(J),J-1)                       &
    &                        +    DSPB(I+IHW(J),J+1) +   DSPB(I+IHE(J),J+1)                       &
    &                        +    4. * DSPB(I,J)) * HBM2(I,J)   * 0.125
!
                    WF2(I,J) =   (DSP0(I+IHW(J),J-1) +   DSP0(I+IHE(J),J-1)                       &
    &                        +    DSP0(I+IHW(J),J+1) +   DSP0(I+IHE(J),J+1)                       &
    &                        +    4. * DSP0(I,J)) * HBM2(I,J)   * 0.125
!
                    WF3(I,J) =   (DSPT(I+IHW(J),J-1) +   DSPT(I+IHE(J),J-1)                       &
    &                        +    DSPT(I+IHW(J),J+1) +   DSPT(I+IHE(J),J+1)                       &
    &                        +    4. * DSPT(I,J)) * HBM2(I,J)   * 0.125
!
                    WF4(I,J) = (CLDEFI(I+IHW(J),J-1) + CLDEFI(I+IHE(J),J-1)                       &
    &                        +  CLDEFI(I+IHW(J),J+1) + CLDEFI(I+IHE(J),J+1)                       &
    &                        +    4. * CLDEFI(I,J)) * HBM2(I,J) * 0.125
        345 END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (RSMK   , SMK)
!
            DO J=MYJS,MYJE
                DO I=MYIS,MYIE
                    SMK = SM(I,J)
                    RSMK = 1. - SMK
                      DSPB(I,J) =   DSPB(I,J) * SMK + WF1(I,J) * RSMK
                      DSP0(I,J) =   DSP0(I,J) * SMK + WF2(I,J) * RSMK
                      DSPT(I,J) =   DSPT(I,J) * SMK + WF3(I,J) * RSMK
                    CLDEFI(I,J) = CLDEFI(I,J) * SMK + WF4(I,J) * RSMK
                END DO
            END DO
!        
    350 END DO
!
    END IF
!-----------------------------------------------
! INITIALIZE CHANGES OF T AND Q DUE TO CONVECTION
!-----------------------------------------------
!
!------- 
! OPENMP
!-------  
!
!$omp parallel do
!
    DO 360 K=1,LM
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                TMOD(I,J,K) = 0.
                QMOD(I,J,K) = 0.
            END DO
        END DO
360 END DO
!-------------------------------------------
! CLEAN UP AND GATHER DEEP CONVECTION POINTS 
!-------------------------------------------
!
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
     DO 380 J=MYJS2,MYJE2
        DO 380 I=MYIS1,MYIE1
            IF (LTOP(I,J) >= LBOT(I,J)) THEN
                LBOT(I,J) = 0
                LTOP(I,J) = LBOT(I,J)
                PTOP(I,J) = PBOT(I,J)
            END IF
!
            IF (HBM2(I,J) < 0.90) THEN
                LBOT(I,J) = 0
                LTOP(I,J) = LBOT(I,J)
                PTOP(I,J) = PBOT(I,J)
            END IF
!
            IF (PTOP(I,J) > PBOT(I,J)-PNO .OR. LTOP(I,J) > LBOT(I,J)-2)                           &
    &           CLDEFI(I,J) = AVGEFI * SM(I,J) + STEFI * (1.-SM(I,J))

!
     380 END DO
!
    KHDEEP = 0
    PSHNEW = 20000.
!
    DO J=MYJS2,MYJE2
        DO I=MYIS1,MYIE1
            PSFCIJ = PD(I,J) + PT
!--------------------------------------------------------------------------------------------------
! DEPTH OF CLOUD REQUIRED TO MAKE THE POINT A DEEP CONVECTION POINT IS NOW A SCALED VALUE OF THE 
! PSFC INSTEAD OF 290 MB EVERYWHERE
!--------------------------------------------------------------------------------------------------
            DEPMIN = PSHNEW * PSFCIJ * 1.E-5
            DEPTH  = PBOT(I,J) - PTOP(I,J)
!
            IF (DEPTH >= DEPMIN) THEN
                KHDEEP = KHDEEP + 1
                IDEEP(KHDEEP) = I
                JDEEP(KHDEEP) = J
            END IF
        END DO
    END DO
!------------------------------------
! HORIZONTAL LOOP FOR DEEP CONVECTION
!------------------------------------
!
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!$omp private (APEK     , APEKL   , APEKXX  , APEKXY  , APESK   , AVRGT   , AVTGTL  , DENTPY  ,   &
!$omp          DEPMIN   , DEPTH   , DEPWL   , DHDT    , DIFQ    , DIFQL   , DIFT    , DIFTL   ,   &
!$omp          DRHEAT   , DRHDP   , DSP     , DSP0K   , DSPBK   , DSPTK   , DTHEM   , EFI     ,   &
!$omp          FEFI     , HCORR   , I       , J       , L0      , L0M1    , LB      , LBM1    ,   &
!$omp          LBTK     , LCOR    , LQM     , LSHU    , LTP1    , LTP2    , LTPK    , LTSH    ,   &
!$omp          PBTK     , PK      , PK0     , PKB     , PKL     , PKT     , PRECK   , PSFCIJ  ,   &
!$omp          PSK      , PTHRS   , PTPK    , QK      , QKL     , QREFK   , QSATK   , RDP0T   ,   &
!$omp          RHH      , RHL     , RHMAX   , SUMDE   , SUMDP   , THERK   , THERKX  , THERKY  ,   &
!$omp          THSK     , THSKL   , TK      , TKL     , TREFK   , TREFKX  , TSKL    )
!
!GSM SHALLOW CONVECTION
!
    NDEEP = 0
!
    DO 620 J=MYJS2,MYJE2
        DO 620 I=MYIS,MYIE
            LTPK = LTOP(I,J)
            LBTK = LBOT(I,J)
            LB   =  LMH(I,J) - 1
            PSFCIJ = PD(I,J) + PT
            DEPMIN = PSHNEW * PSFCIJ * 1.E-5
!
            IF (PTOP(I,J) < PBOT(I,J)-DEPMIN) THEN
                NDEEP = NDEEP + 1
                NDEPTH = LB - LTPK
                NTOPD (LTPK  ) = NTOPD (LTPK  ) + 1
                NBOTD (LB    ) = NBOTD (LB    ) + 1
                NDPTHD(NDEPTH) = NDPTHD(NDEPTH) + 1
            END IF
620 END DO
!
    NNEG = KHDEEP - NDEEP
!--------------------------------- 
! GATHER SHALLOW CONVECTION POINTS
!---------------------------------
    KHSHAL = 0
    NDSTN  = 0
    NDSTP  = 0
!
    DO 630 J=MYJS2,MYJE2
        DO 630 I=MYIS,MYIE
            IF (PTOP(I,J) > PBOT(I,J)-PNO .OR. LTOP(I,J) > LBOT(I,J)-2) GOTO 630
                PSFCIJ = PD(I,J) + PT
                KHSHAL = KHSHAL + 1
                ISHAL(KHSHAL) = I
                JSHAL(KHSHAL) = J
!
630 END DO
!----------------------------------------
! HORIZONTAL LOOP FOR SHALLOW CONVECTION
!----------------------------------------
!
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!$omp private (APEK     , APEKL   , APEKXX  , APEKXY  , BQK     , BQS00K  , BQS10K  , DEN     ,   &
!$omp          DENTPY   , DPKL    , DPMIX   , DQREF   , DST     , DSTQ    , DTDETA  , FPK     ,   &
!$omp          FPTK     , I       , IQ      , IT      , J       , LBM1    , LBTK    , LTP1    ,   &
!$omp          LTPK     , OTSUM   , PART1   , PART2   , PART3   , PK      , PKL     , PKXXXX  ,   &
!$omp          PKXXXY   , POTSUM  , PPK     , PSUM    , PTPK    , PZ0     , QK      , QKL     ,   &
!$omp          QNEW     , QOTSUM  , QQK     , QREFK   , QRFKL   , QRFTP   , QSATK   , QSUM    ,   &
!$omp          RDPSUM   , RTBAR   , SMIX    , SQK     , SQS00K  , SQS10K  , SUMDP   , SUMDT   ,   &
!$omp          TCORR    , THVMKL  , THVREF  , TK      , TKL     , TQK     , TREFK   , TREFKX  ,   &
!$omp          TRFKL    , TTHK)
!
    DO 800 N=1,KHSHAL
!
        I = ISHAL(N)
        J = JSHAL(N)
!------------------- 
! SHALLOW CONVECTION
!------------------- 
        PZ0  =  PD(I,J) + PT
        LLMH = LMH(I,J)
!    
        DO 650 K=1,LLMH
            TKL      = T(I,J,K)
            TK   (K) = TKL
            TREFK(K) = TKL
            QKL      = Q(I,J,K)
            QK   (K) = QKL
            QREFK(K) = QKL
            PKL      = AETA(K) * PDSL(I,J) + PT
            PK   (K) = PKL
            QSATK(K) = PQ0 / PK(K) * EXP(A2 * (TK(K) - A3) / (TK(K) - A4))
            APEKL    = APE(I,J,K)
!-----------------------------
! CHOOSE THE PRESSURE FUNCTION
!
! FPK  (K) =ALOG(PKL)
! FPK  (K) =PKL
! FPK  (K) =-1./PKL
!-----------------------------
            APEK (K)  = APEKL
            THVMKL    = TKL * APEKL * (QKL * 0.608 + 1.)
            THVREF(K) = THVMKL
        650 END DO
!------------------ 
! SHALLOW CLOUD TOP
!------------------ 
        LBTK = LBOT(I,J)
        LBM1 = LBTK - 1
        PTPK = PTOP(I,J)
        LTP1 = LTOP(I,J)
        LTPK = LTOP(I,J) - 1
!
        IF (PTOP(I,J) > PBOT(I,J)-PNO .OR. LTOP(I,J) > LBOT(I,J)-2) THEN
            LBOT(I,J) = 0
            LTOP(I,J) = LBOT(I,J)
            PTOP(I,J) = PBOT(I,J)
            GOTO 800
        END IF
!-----------------------------------------------------
! SCALING POTENTIAL TEMPERATURE AND TABLE INDEX AT TOP
!-----------------------------------------------------
        THTPK = T(I,J,LTPK) * APE(I,J,LTPK)
!    
        TTHK = (THTPK - THL) * RDTH
        QQK  = TTHK - AINT(TTHK)
        IT   = INT(TTHK) + 1
!
        IF (IT < 1) THEN
            IT  = 1
            QQK = 0.
        END IF
!
        IF (IT >= JTB) THEN
            IT  = JTB - 1
            QQK = 0.
        END IF
!--------------------------------------------------
! BASE AND SCALING FACTOR FOR SPEC. HUMIDITY AT TOP
!--------------------------------------------------
        BQS00K = QS0(IT)
        SQS00K = SQS(IT)
        BQS10K = QS0(IT+1)
        SQS10K = SQS(IT+1)
!----------------------------------------------
! SCALING SPEC. HUMIDITY AND TABLE INDEX AT TOP
!----------------------------------------------
        BQK = (BQS10K - BQS00K) * QQK + BQS00K
        SQK = (SQS10K - SQS00K) * QQK + SQS00K
!
        TQK = (Q(I,J,LTPK) - BQK) / SQK * RDQ
!    
        PPK = TQK - AINT(TQK)
        IQ  = INT(TQK) + 1
!
        IF (IQ < 1) THEN
            IQ  = 1
            PPK = 0.
        END IF
!
        IF (IQ >= ITB) THEN
            IQ  = ITB - 1
            PPK = 0.
        END IF
!------------------------------------
! CLOUD TOP SATURATION POINT PRESSURE
!------------------------------------
        PART1 = (PTBL(IQ+1,IT) - PTBL(IQ  ,IT)) * PPK
        PART2 = (PTBL(IQ,IT+1) - PTBL(IQ  ,IT)) * QQK
        PART3 = (PTBL(IQ,IT  ) - PTBL(IQ+1,IT) - PTBL(IQ,IT+1) + PTBL(IQ+1,IT+1)) * PPK * QQK
        PTPK  =  PTBL(IQ,IT  ) + PART1 + PART2 + PART3
!
        DPMIX = PTPK - PSP(I,J)
        IF (ABS(DPMIX) < 3000.) DPMIX = -3000.
!-------------------------- 
! TEMPERATURE PROFILE SLOPE
!-------------------------- 
        SMIX = (THTPK - THBT(I,J)) / DPMIX * STABS
!    
        TREFKX = TREFK(LBTK+1)
        PKXXXX =    PK(LBTK+1)
        PKXXXY =    PK(LBTK)
        APEKXX =  APEK(LBTK+1)
        APEKXY =  APEK(LBTK)
!    
        DO 670 K=LBTK,LTP1,-1
            TREFKX   = ((PKXXXY - PKXXXX) * SMIX + TREFKX * APEKXX) / APEKXY
            TREFK(K) = TREFKX
            APEKXX   = APEKXY
            PKXXXX   = PKXXXY
            APEKXY   = APEK(K-1)
            PKXXXY   =   PK(K-1)
    670 END DO
!-----------------------------------------
! TEMPERATURE REFERENCE PROFILE CORRECTION 
!-----------------------------------------
        SUMDT = 0.
        SUMDP = 0.
!    
        DO K=LTP1,LBTK
            SUMDT = (TK(K) - TREFK(K)) * DETA(K) + SUMDT
            SUMDP = SUMDP + DETA(K)
        END DO
!    
        RDPSUM = 1. / SUMDP
        FPK(LBTK) = TREFK(LBTK)
!    
        TCORR=SUMDT*RDPSUM
!    
        DO K=LTP1,LBTK
            TRFKL    = TREFK(K)+TCORR
            TREFK(K) = TRFKL
            FPK  (K) = TRFKL
        END DO
!--------------------------- 
! HUMIDITY PROFILE EQUATIONS
!---------------------------
        PSUM   = 0.
        QSUM   = 0.
        POTSUM = 0.
        QOTSUM = 0.
        OTSUM  = 0.
        DST    = 0.
        FPTK   = FPK(LTP1)
!    
        DO 700 K=LTP1,LBTK
            DPKL   = FPK(K) - FPTK
            PSUM   = DPKL  * DETA(K) + PSUM
            QSUM   = QK(K) * DETA(K) + QSUM
            RTBAR  = 2. / (TREFK(K) + TK(K))
            OTSUM  = DETA(K) * RTBAR + OTSUM
            POTSUM = DPKL   * RTBAR * DETA(K) + POTSUM
            QOTSUM = QK(K)  * RTBAR * DETA(K) + QOTSUM
            DST    = (TREFK(K) - TK(K)) * RTBAR * DETA(K) + DST
    700 END DO
!    
        PSUM   = PSUM * RDPSUM
        QSUM   = QSUM * RDPSUM
        ROTSUM = 1. / OTSUM
        POTSUM = POTSUM * ROTSUM
        QOTSUM = QOTSUM * ROTSUM
        DST    = DST    * ROTSUM * CP / ELWV
!------------------------------- 
! ENSURE POSITIVE ENTROPY CHANGE 
!-------------------------------
        IF (DST > 0.) THEN
            LBOT(I,J) = 0
            LTOP(I,J) = LBOT(I,J)
            PTOP(I,J) = PBOT(I,J)
            GOTO 800       
        ELSE
            DSTQ = DST * EPSDN
        END IF
!-------------------------------- 
! CHECK FOR ISOTHERMAL ATMOSPHERE
!--------------------------------
        DEN = POTSUM - PSUM
!    
        IF (-DEN/PSUM < 5.E-5) THEN
            LBOT(I,J) = 0
            LTOP(I,J) = LBOT(I,J)
            PTOP(I,J) = PBOT(I,J)
            GOTO 800
!----------------------------------------        
! SLOPE OF THE REFERENCE HUMIDITY PROFILE
!----------------------------------------
        ELSE
            DQREF = (QOTSUM-DSTQ-QSUM) / DEN
        END IF
!--------------------------------------- 
! HUMIDITY DOES NOT INCREASE WITH HEIGHT
!---------------------------------------
        IF (DQREF < 0.) THEN
            LBOT(I,J) = 0
            LTOP(I,J) = LBOT(I,J)
            PTOP(I,J) = PBOT(I,J)
            GOTO 800
        END IF
!--------------------------
! HUMIDITY AT THE CLOUD TOP
!--------------------------
        QRFTP = QSUM - DQREF * PSUM
!----------------- 
! HUMIDITY PROFILE
!-----------------
        DO 720 K=LTP1,LBTK
            QRFKL = (FPK(K) - FPTK) * DQREF + QRFTP
!------------------------------------------
! SUPERSATURATION OR NEGATIVE Q NOT ALLOWED
!-------------------------------------------    
            QNEW = (QRFKL - QK(K)) * TAUK + QK(K)
!
            IF (QNEW > QSATK(K)*STRESH .OR. QNEW < 0.) THEN
                LBOT(I,J) = 0
                LTOP(I,J) = LBOT(I,J)
                PTOP(I,J) = PBOT(I,J)
                GOTO 800
            END IF
!
            THVREF(K) = TREFK(K) * APEK(K) * (QRFKL * 0.608 + 1.)
             QREFK(K) = QRFKL
    720 END DO
!-----------------------------------------------
! ELIMINATE IMPOSSIBLE SLOPES (BETTS, DTHETA/DQ)
!-----------------------------------------------
        DO 730 K=LTP1,LBTK
            DTDETA = (THVREF(K-1) - THVREF(K)) / (AETA(K) - AETA(K-1))
            IF (DTDETA < EPSTH) THEN
                LBOT(I,J) = 0
                LTOP(I,J) = LBOT(I,J)
                PTOP(I,J) = PBOT(I,J)
                GOTO 800
            END IF
    730 END DO
!-------------------- 
! DIAGNOSTICS 
!
! IF (DST.GT.0.) THEN
! NDSTP = NDSTP + 1
! ELSE
! NDSTN = NDSTN + 1
! END IF
!--------------------
        DENTPY = 0.
!    
        DO K=LTP1,LBTK
            DENTPY = ((TREFK(K) - TK(K)) * CP + (QREFK(K) - QK(K)) * ELWV)                        &
    &              / (TK(K) + TREFK(K)) * DETA(K) + DENTPY
        END DO
!--------------------------------------    
! RELAXATION TOWARDS REFERENCE PROFILES 
!--------------------------------------

        DO 750 K=LTP1,LBTK
            TMOD(I,J,K) = (TREFK(K) - TK(K)) * TAUK
            QMOD(I,J,K) = (QREFK(K) - QK(K)) * TAUK
    750 END DO       
800 END DO
!-------------------------- 
! END OF SHALLOW CONVECTION
!--------------------------
!
!------------
! DIAGNOSTICS
!------------
    NSHAL = 0
!
    DO 820 J=MYJS2,MYJE2
        DO 820 I=MYIS,MYIE
            LTPK = LTOP(I,J)
            LBTK = LBOT(I,J)
            PTPK = PTOP(I,J)
            PBTK = PBOT(I,J)
!
            IF (PTPK > PBTK-PNO .OR. LTPK > LBTK-2) GOTO 820
            PSFCIJ = PD(I,J) + PT
            DEPMIN = PSHNEW * PSFCIJ * 1.E-5
!
            IF (PTPK >= PBTK-DEPMIN) THEN
                NSHAL          = NSHAL + 1
                NTOPS(LTPK)    = NTOPS(LTPK) + 1
                NBOTS(LBTK)    = NBOTS(LBTK) + 1
                NDEPTH         = LBTK - LTPK
                NDPTHS(NDEPTH) = NDPTHS(NDEPTH) + 1
            END IF
!
820 END DO
!
    NEGDS = KHSHAL - NSHAL
!---------------------------------------------
! SMOOTHING TEMPERATURE & HUMIDITY CORRECTIONS
!---------------------------------------------
    IF (KSMUD == 0) THEN
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!    
        DO 900 K=1,LM
!-----------------------------------------
! UPDATE THE FUNDAMENTAL PROGNOSTIC ARRAYS
!-----------------------------------------
            DO 830 J=MYJS,MYJE
                DO 830 I=MYIS,MYIE
                    T(I,J,K) = T(I,J,K) + TMOD(I,J,K)
                    Q(I,J,K) = Q(I,J,K) + QMOD(I,J,K)
!---------------------------------------------------------------------------------------------------------------
! ACCUMULATE LATENT HEATING DUE TO CONVECTION.
! SCALE BY THE RECIPROCAL OF THE PERIOD AT WHICH THIS ROUTINE IS CALLED. THIS PERIOD IS THE CONVECTION TIMESTEP.
!---------------------------------------------------------------------------------------------------------------
                    TCUCN(I,J,K) = TCUCN(I,J,K) + TMOD(I,J,K) * RDTCNVC

       830 END DO
 
    900 END DO 
  
!
    ELSE
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (QL     , QNE     , QSE     , TL      , TNE     , TSE)
!    
        DO 910 K=1,LM
!        
            CALL ZERO2(QL)
            CALL ZERO2(QNE)
            CALL ZERO2(QSE)
            CALL ZERO2(TL)
            CALL ZERO2(TNE)
            CALL ZERO2(TSE)
!
            DO J=MYJS,MYJE
                DO I=MYIS,MYIE
                    TL(I,J) = TMOD(I,J,K)
                    QL(I,J) = QMOD(I,J,K)
                END DO
            END DO
!
            NSMUD = KSMUD
!
            DO 870 KS=1,NSMUD
!
                DO J=MYJS,MYJE1
                    DO I=MYIS,MYIE
                        TNE(I,J) = (TL(I+IHE(J),J+1) - TL(I,J)) * HTM(I,J,K) * HTM(I+IHE(J),J+1,K)
                        QNE(I,J) = (QL(I+IHE(J),J+1) - QL(I,J)) * HTM(I,J,K) * HTM(I+IHE(J),J+1,K)
                    END DO
                END DO
!
                DO J=MYJS1,MYJE
                    DO I=MYIS,MYIE
                        TSE(I,J) = (TL(I+IHE(J),J-1) - TL(I,J)) * HTM(I+IHE(J),J-1,K) * HTM(I,J,K)
                        QSE(I,J) = (QL(I+IHE(J),J-1) - QL(I,J)) * HTM(I+IHE(J),J-1,K) * HTM(I,J,K)
                    END DO
                END DO
!
                DO J=MYJS2,MYJE2
                    DO I=MYIS,MYIE
                        TL(I,J) = (TNE(I,J) - TNE(I+IHW(J),J-1) + TSE(I,J) - TSE(I+IHW(J),J+1))   &
    &                           * HBM2(I,J) * 0.125 + TL(I,J)
                        QL(I,J) = (QNE(I,J) - QNE(I+IHW(J),J-1) + QSE(I,J) - QSE(I+IHW(J),J+1))   &
    &                           * HBM2(I,J) * 0.125 + QL(I,J)
                    END DO
                END DO
!
        870 END DO
!-----------------------------------------
! UPDATE THE FUNDAMENTAL PROGNOSTIC ARRAYS
!-----------------------------------------
            DO J=MYJS,MYJE
                DO I=MYIS,MYIE
                    T(I,J,K) = T(I,J,K) + TL(I,J)
                    Q(I,J,K) = Q(I,J,K) + QL(I,J)
                END DO
            END DO
!---------------------------------------------------------------------------------------------------------------
! ACCUMULATE LATENT HEATING DUE TO CONVECTION.
! SCALE BY THE RECIPROCAL OF THE PERIOD AT WHICH THIS ROUTINE IS CALLED. THIS PERIOD IS THE CONVECTION TIMESTEP.
!---------------------------------------------------------------------------------------------------------------
            DO J=MYJS,MYJE
                DO I=MYIS,MYIE
                    TCUCN(I,J,K) = TCUCN(I,J,K) + TL(I,J) * RDTCNVC
                END DO
            END DO
!        
    910 END DO
!    
    END IF
!----------------------------------------
! SAVE CLOUD TOP AND BOTTOM FOR RADIATION
!----------------------------------------
!
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            IF (LTOP(I,J) > 0 .AND. LTOP(I,J) < LBOT(I,J)) THEN
!
                CUTOP = FLOAT(LTOP(I,J))
                CUBOT = FLOAT(LBOT(I,J))
!
                  HTOP(I,J) = MIN(CUTOP,  HTOP(I,J))
                  HBOT(I,J) = MAX(CUBOT,  HBOT(I,J))
                CNVTOP(I,J) = MIN(CUTOP,CNVTOP(I,J))
                CNVBOT(I,J) = MAX(CUBOT,CNVBOT(I,J))
            END IF
        END DO
    END DO
!---------------------------------------------------------------------------------------------------------------------
! DIAGNOSTICS
!
! WRITE(LIST,950) NTSD, NSHAL, NDEEP, NNEG, NEGDS, NDSTN, NDSTP
!
! DO 940 L=1,LM
!    WRITE(LIST,952) L
!    WRITE(LIST,954) NBOTS(L), NTOPS(L), NDPTHS(L), NBOTD(L), NTOPD(L), NDPTHD(L)
! 940 CONTINUE
!
! 950 FORMAT(' NTSD=',I3,I8,' SHALLOW, ',I4,' DEEP, ',I4,' NEG., ',I4,' NEG. SHALL.,',I4,' DST.LT.0, ',I4,' DST.GT.0')
! 952 FORMAT(' LAYER (FROM TOP),',I2)
! 954 FORMAT('     NBOTS=',I4,'     NTOPS=',I4,'     NDPTHS=',I4,'     NBOTD=',I4,'     NTOPD=',I4,'     NDPTHD=',I4)
!---------------------------------------------------------------------------------------------------------------------
    RETURN
!
    END SUBROUTINE CUCNVC_SHALLOW
