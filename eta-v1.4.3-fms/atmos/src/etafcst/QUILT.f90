    SUBROUTINE QUILT
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE QUILT
!>
!> SUBROUTINE: QUILT - I/O SERVERS
!> PROGRAMMER: TUCCILLO
!> ORG: IBM
!> DATE: 00-01-20
!>
!> ABSTRACT:
!> I/O SERVERS
!>
!> PROGRAM HISTORY LOG:
!> 00-01-20  TUCCILLO - ORIGINATOR
!> 00-11-02  BLACK    - SLP FOR NEST BOUNDARIES
!> 00-12-12  BLACK    - RESTART CAPABILITY
!> 18-01-15  LUCCI    - MODERNIZATION OF THE CODE, INCLUDING:
!>                      * F77 TO F90/F95
!>                      * INDENTATION & UNIFORMIZATION CODE
!>                      * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                      * DOCUMENTATION WITH DOXYGEN
!>                      * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NONE
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT FILES:  NONE
!>
!> OUTPUT FILES:  NONE
!>
!> USE MODULES: BUFFER
!>              F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MPPCOM
!>              PARA
!>              PARMETA
!>              PARMSOIL
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : EBU
!>
!> CALLS      : COLLECT
!>              DECOAL
!>              GETENV
!>              MPI_ABORT
!>              MPI_BCAST
!>              MPI_FIRST
!>              MPI_RECV
!>              SLP
!>              SLPSIG
!>              SLPSIGSPLINE
!-------------------------------------------------------------------------------------------------- 
!
!--------------------------------------------------------------------------------- 
! THIS CODE ASSUMES THAT NSOIL IS GE TO 4. IF THIS IS NOT TRUE,THE CODE WILL STOP. 
! THE EQUIVALENCING IS THE PROBLEM.
!--------------------------------------------------------------------------------- 
    USE BUFFER
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARA
    USE PARMETA
    USE PARMSOIL
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: LB = 2 * IM + JM - 3
!
    REAL   (KIND=R4KIND), DIMENSION(IM, JM)                                                     ::&
    & DUM1    , DUM2    , DUM3    , DUM4    , DUM5    , DUM6    , DUM7
!---------
! SM v100M
!---------
    REAL   (KIND=R4KIND), DIMENSION(IM, JM)                                                     ::&
    & DUM8    , DUM9    , DUM10   , DUM11
!---------
! SM v100M
!---------
! Lyra GSM Wind stress
    REAL   (KIND=R4KIND), DIMENSION(IM, JM)                                                     ::&
    & DUM12    , DUM13
! Lyra GSM Wind stress
!
    REAL   (KIND=R4KIND), DIMENSION(IM, JM, NSOIL)                                              ::&
    & DUMS
!
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & STATUS
!
    EQUIVALENCE (DUM1(1,1), DUMS(1,1,1))
    EQUIVALENCE (DUM2(1,1), DUMS(1,1,2))
    EQUIVALENCE (DUM3(1,1), DUMS(1,1,3))
    EQUIVALENCE (DUM4(1,1), DUMS(1,1,4))
!
    REAL   (KIND=R4KIND), DIMENSION(:,:), ALLOCATABLE                                           ::&
    & PDOMG   , RESOMG  , PD      , RES     , FIS     ,                                           &
    & RSWIN   , RSWOUT  , TG      , Z0      , AKMS    ,                                           &
    & CZEN    , AKHS    , THS     , QS      , TWBS    ,                                           &
    & QWBS    , HBOT    , CFRACL  , THZ0    , QZ0     ,                                           &
    & UZ0     , VZ0     , USTAR   , HTOP    , CFRACM  ,                                           &
    & SNO     , SI      , CLDEFI  , RF      , PSLP    ,                                           &
    & CUPPT   , CFRACH  , SOILTB  , SFCEXC  ,                                                     &
    & SMSTAV  , SMSTOT  , GRNFLX  , PCTSNO  ,                                                     &
    & RLWIN   , RADOT   , CZMEAN  , SIGT4   ,                                                     &
    & U00     , SR      ,   PREC  , ACPREC  , ACCLIQ  ,                                           &
    & CUPREC  , ACFRCV  , ACFRST  , SFCSHX  ,                                                     &
    & ACSNOW  , ACSNOM  , SSROFF  , BGROFF  ,                                                     &
    & SFCLHX  , SUBSHX  , SNOPCX  , SFCUVX  ,                                                     &
    & SFCEVP  , POTEVP  ,  ASWIN  , ASWOUT  ,                                                     &
    & ASWTOA  , ALWIN   , ALWOUT  , ALWTOA  ,                                                     &
!---------
! SM v100M
!---------
    & TH100   , Q100    , U100    , V100    ,                                                     &
!---------
! SM v100M
!---------
! Lyra GSM Wind stress
    & XMOMFLUX, YMOMFLUX,                                                                         &
! Lyra GSM Wind stress
!
    & TH10    , Q10     , U10     , V10     , TSHLTR  ,                                           &
    & QSHLTR  , PSHLTR  , CMC     , POTFLX  ,                                                     &
    & TLMIN   , TLMAX   , RSWTOA  , RLWTOA  ,                                                     &
    & CNVBOT  , CNVTOP  , ALBEDO  ,                                                               &
!
!Lyra GSM Max wind
    & MAXWU , MAXWV,                                                                              &
!Lyra GSM Max wind
!
    & SM      , HBM2 
    REAL   (KIND=R4KIND), DIMENSION(:), ALLOCATABLE                                             ::&
    & DETA
!
    REAL   (KIND=R4KIND), DIMENSION(2*LM)                                                       ::&
    & UL
!
    REAL   (KIND=R4KIND), DIMENSION(:,:,:), ALLOCATABLE                                         ::&
    & OMGALF  , T       , Q       , U       , V       , Q2      , TTND    , CWM     , TRAIN   ,   &
    & TCUCN   , RSWTT   , RLWTT   , SMC     , STC     , SH2O
!
    REAL   (KIND=R4KIND), DIMENSION(LB, 2)                                                      ::&
    & PDB
!
    REAL   (KIND=R4KIND), DIMENSION(LB, LM, 2)                                                  ::&
    & TB      , QB      , UB      , VB      , Q2B     , CWMB
!
    INTEGER(KIND=I4KIND), DIMENSION(3)                                                          ::&
    & IDAT
!
    INTEGER(KIND=I4KIND), DIMENSION(:,:), ALLOCATABLE                                           ::&
    & LC      , NCFRCV  , NCFRST
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & RUN     , FIRST   , HYDRO   , SIGMA
!
    CHARACTER(LEN= 56) :: RSTFIL1 , RSTFIL2 , ENVAR   , FINFIL
    CHARACTER(LEN=  8) :: RESTHR  , BLANK
    CHARACTER(LEN= 32) :: LABEL
    CHARACTER(LEN= 80) :: FNAME
    CHARACTER(LEN= 16) :: DONE 
!
    CHARACTER(LEN=200) :: SUBMIT !CHOU
    CHARACTER(LEN=200) :: CHMO
    CHARACTER(LEN= 88) :: LINE
    CHARACTER(LEN= 32) :: OUTJOB
!           
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & LME     , OUTJ
!
    DATA LRSTRT1/21/, LRSTRT2/61/, NHB/12/, BLANK/'    '/
!
    REAL   (KIND=R8KIND)                                                                        ::&
    & TIMEF   , IST     , ISP     , RTC     , IST2    , ISP2    , ICUM
!
    REAL   (KIND=R4KIND), DIMENSION(999999)                                                     ::&
    & TSHDE
!
    REAL   (KIND=R4KIND), DIMENSION(LSM)                                                        ::&
    & SPL
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & RESTRT  , SINGLRST, SUBPOST , NEST    , SPLINE
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & LRSTRT1 , LRSTRT2 , NHB     , NMAP    , NPHS    , NCNVC   , NRADSH  , NRADLH  , NTDDMP  ,   &
    & LFCSTD  , IERR    , NFCST   , NBC     , LIST    , IDTAD   , IHOUR   , IER     , IXXX    ,   &
    & IPE     , IS      , IE      , JS      , JE      , LEN_CH  , IDUM    , IHRST   , NTSD    ,   &
    & L       , IOUT    , NSHDE   , NHEAT   , NPREC   , NRDSW   , NSRFC   , KPATH   , KRST    ,   &
    & ITAG    , LFINFIL , NRDLW  
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & TSTART  , TEND    , TCP     , TPREC   , THEAT   , TCLOD   , TRDSW   , TRDLW   , TSRFC   ,   &
    & DT      , ARDSW   , ARDLW   , ASRFC   , AVRAIN  , AVCNVC  , ACUTIM  , ARATIM  , APHTIM  ,   &
    & TPH0D   , TLM0D   , PT       
!-----------------
! DECLARE NAMELIST
!-----------------
    NAMELIST /FCSTDATA/                                                                           &
    & TSTART, TEND  , TCP  , RESTRT, SINGLRST, SUBPOST, NMAP , TSHDE, SPL, NPHS  , NCNVC , NRADSH,&
    & NRADLH, NTDDMP, TPREC, THEAT , TCLOD   , TRDSW  , TRDLW, TSRFC, NEST, HYDRO, SPLINE
!
    NAMELIST /POSTLIST/ CHMO
    NAMELIST /POSTLIST/ SUBMIT !CHOU
!
    CALL MPI_FIRST
!-----------------------------------------------------------
! READ NAMELIST FCSTDATA TO FIND OUT IF THIS IS A NESTED RUN
!-----------------------------------------------------------
    LFCSTD = 11
    OPEN(LFCSTD, FILE='FCSTDATA.meso', STATUS='OLD')
    READ(LFCSTD, FCSTDATA)
!
    IF (NSOIL < 4) THEN
        PRINT*, ' NSOIL IS LESS THAN 4. CHANGE THE EQUIVALENCES'
        PRINT*, ' STOPPING'
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    END IF
!
    IF (ME == 0) THEN
        LME = .TRUE. 
    ELSE
        LME = .FALSE. 
    END IF
!------------------------------ 
! READ NHB FILE TO OBTAIN SIGMA
!------------------------------          
    OPEN(UNIT=NHB, FORM='UNFORMATTED', FILE='CNST.file')
    REWIND NHB
    READ(NHB) NFCST, NBC, LIST, DT, IDTAD, SIGMA
    IF (ME == 0) THEN
        PRINT*,'IN QUILT, SIGMA= ', SIGMA
    END IF
!
    ALLOCATE   (PDOMG(IM,MY_JSD:MY_JED))
    ALLOCATE  (RESOMG(IM,MY_JSD:MY_JED))
    ALLOCATE  (OMGALF(IM,MY_JSD:MY_JED,1:LM))
    ALLOCATE      (PD(IM,MY_JSD:MY_JED))
    ALLOCATE     (RES(IM,MY_JSD:MY_JED))
    ALLOCATE     (FIS(IM,MY_JSD:MY_JED))
    ALLOCATE       (T(IM,MY_JSD:MY_JED,1:LM))
    ALLOCATE       (Q(IM,MY_JSD:MY_JED,1:LM))
    ALLOCATE       (U(IM,MY_JSD:MY_JED,1:LM))
    ALLOCATE       (V(IM,MY_JSD:MY_JED,1:LM))
    ALLOCATE      (Q2(IM,MY_JSD:MY_JED,1:LM))
    ALLOCATE    (TTND(IM,MY_JSD:MY_JED,1:LM))
    ALLOCATE     (CWM(IM,MY_JSD:MY_JED,1:LM))
    ALLOCATE   (TRAIN(IM,MY_JSD:MY_JED,1:LM))
    ALLOCATE   (TCUCN(IM,MY_JSD:MY_JED,1:LM))
    ALLOCATE   (RSWIN(IM,MY_JSD:MY_JED))
    ALLOCATE  (RSWOUT(IM,MY_JSD:MY_JED))
    ALLOCATE      (TG(IM,MY_JSD:MY_JED))
    ALLOCATE      (Z0(IM,MY_JSD:MY_JED))
    ALLOCATE    (AKMS(IM,MY_JSD:MY_JED))
    ALLOCATE    (CZEN(IM,MY_JSD:MY_JED))
    ALLOCATE    (AKHS(IM,MY_JSD:MY_JED))
    ALLOCATE     (THS(IM,MY_JSD:MY_JED))
    ALLOCATE      (QS(IM,MY_JSD:MY_JED))
    ALLOCATE    (TWBS(IM,MY_JSD:MY_JED))
    ALLOCATE    (QWBS(IM,MY_JSD:MY_JED))
    ALLOCATE    (HBOT(IM,MY_JSD:MY_JED))
    ALLOCATE  (CFRACL(IM,MY_JSD:MY_JED))
    ALLOCATE    (THZ0(IM,MY_JSD:MY_JED))
    ALLOCATE     (QZ0(IM,MY_JSD:MY_JED))
    ALLOCATE     (UZ0(IM,MY_JSD:MY_JED))
    ALLOCATE     (VZ0(IM,MY_JSD:MY_JED))
    ALLOCATE   (USTAR(IM,MY_JSD:MY_JED))
    ALLOCATE    (HTOP(IM,MY_JSD:MY_JED))
    ALLOCATE     (SNO(IM,MY_JSD:MY_JED))
    ALLOCATE      (SI(IM,MY_JSD:MY_JED))
    ALLOCATE  (CLDEFI(IM,MY_JSD:MY_JED))
    ALLOCATE      (RF(IM,MY_JSD:MY_JED))
    ALLOCATE    (PSLP(IM,MY_JSD:MY_JED))
    ALLOCATE   (CUPPT(IM,MY_JSD:MY_JED))
    ALLOCATE  (CFRACH(IM,MY_JSD:MY_JED))
    ALLOCATE  (CFRACM(IM,MY_JSD:MY_JED))
    ALLOCATE  (SOILTB(IM,MY_JSD:MY_JED))
    ALLOCATE  (SFCEXC(IM,MY_JSD:MY_JED))
    ALLOCATE  (SMSTAV(IM,MY_JSD:MY_JED))
    ALLOCATE  (SMSTOT(IM,MY_JSD:MY_JED))
    ALLOCATE  (GRNFLX(IM,MY_JSD:MY_JED))
    ALLOCATE  (PCTSNO(IM,MY_JSD:MY_JED))
    ALLOCATE   (RLWIN(IM,MY_JSD:MY_JED))
    ALLOCATE   (RADOT(IM,MY_JSD:MY_JED))
    ALLOCATE  (CZMEAN(IM,MY_JSD:MY_JED))
    ALLOCATE   (SIGT4(IM,MY_JSD:MY_JED))
    ALLOCATE     (U00(IM,MY_JSD:MY_JED))
    ALLOCATE      (LC(IM,MY_JSD:MY_JED))
    ALLOCATE      (SR(IM,MY_JSD:MY_JED))
    ALLOCATE    (PREC(IM,MY_JSD:MY_JED))
    ALLOCATE  (ACPREC(IM,MY_JSD:MY_JED))
    ALLOCATE  (ACCLIQ(IM,MY_JSD:MY_JED))
    ALLOCATE  (CUPREC(IM,MY_JSD:MY_JED))
    ALLOCATE  (ACFRCV(IM,MY_JSD:MY_JED))
    ALLOCATE  (NCFRCV(IM,MY_JSD:MY_JED))
    ALLOCATE  (ACFRST(IM,MY_JSD:MY_JED))
    ALLOCATE  (NCFRST(IM,MY_JSD:MY_JED))
    ALLOCATE  (ACSNOW(IM,MY_JSD:MY_JED))
    ALLOCATE  (ACSNOM(IM,MY_JSD:MY_JED))
    ALLOCATE  (SSROFF(IM,MY_JSD:MY_JED))
    ALLOCATE  (BGROFF(IM,MY_JSD:MY_JED))
    ALLOCATE  (SFCSHX(IM,MY_JSD:MY_JED))
    ALLOCATE  (SFCLHX(IM,MY_JSD:MY_JED))
    ALLOCATE  (SUBSHX(IM,MY_JSD:MY_JED))
    ALLOCATE  (SNOPCX(IM,MY_JSD:MY_JED))
    ALLOCATE  (SFCUVX(IM,MY_JSD:MY_JED))
    ALLOCATE  (SFCEVP(IM,MY_JSD:MY_JED))
    ALLOCATE  (POTEVP(IM,MY_JSD:MY_JED))
    ALLOCATE   (ASWIN(IM,MY_JSD:MY_JED))
    ALLOCATE  (ASWOUT(IM,MY_JSD:MY_JED))
    ALLOCATE  (ASWTOA(IM,MY_JSD:MY_JED))
    ALLOCATE   (ALWIN(IM,MY_JSD:MY_JED))
    ALLOCATE  (ALWOUT(IM,MY_JSD:MY_JED))
    ALLOCATE  (ALWTOA(IM,MY_JSD:MY_JED))
!---------
! SM v100M
!---------
    ALLOCATE   (TH100(IM,MY_JSD:MY_JED))
    ALLOCATE    (Q100(IM,MY_JSD:MY_JED))
    ALLOCATE    (U100(IM,MY_JSD:MY_JED))
    ALLOCATE    (V100(IM,MY_JSD:MY_JED))
!---------
! SM v100M
!---------
! Lyra GSM Wind stress
    ALLOCATE(XMOMFLUX(IM,MY_JSD:MY_JED))
    ALLOCATE(YMOMFLUX(IM,MY_JSD:MY_JED))
! Lyra GSM Wind stress
!
    ALLOCATE    (TH10(IM,MY_JSD:MY_JED))
    ALLOCATE     (Q10(IM,MY_JSD:MY_JED))
    ALLOCATE     (U10(IM,MY_JSD:MY_JED))
    ALLOCATE     (V10(IM,MY_JSD:MY_JED))
    ALLOCATE  (TSHLTR(IM,MY_JSD:MY_JED))
    ALLOCATE  (QSHLTR(IM,MY_JSD:MY_JED))
    ALLOCATE  (PSHLTR(IM,MY_JSD:MY_JED))
    ALLOCATE     (SMC(IM,MY_JSD:MY_JED,1:NSOIL))
    ALLOCATE     (CMC(IM,MY_JSD:MY_JED))
    ALLOCATE     (STC(IM,MY_JSD:MY_JED,1:NSOIL))
    ALLOCATE    (SH2O(IM,MY_JSD:MY_JED,1:NSOIL))
    ALLOCATE  (ALBEDO(IM,MY_JSD:MY_JED))
    ALLOCATE  (POTFLX(IM,MY_JSD:MY_JED))
    ALLOCATE   (TLMIN(IM,MY_JSD:MY_JED))
    ALLOCATE   (TLMAX(IM,MY_JSD:MY_JED))
!Lyra GSM Max wind
    ALLOCATE(MAXWU(IM,MY_JSD:MY_JED))
    ALLOCATE(MAXWV(IM,MY_JSD:MY_JED))
!Lyra GSM Max wind
!
    ALLOCATE   (RSWTT(IM,MY_JSD:MY_JED,1:LM))
    ALLOCATE   (RLWTT(IM,MY_JSD:MY_JED,1:LM))
!
    ALLOCATE (CNVBOT(MY_ISD:MY_IED,MY_JSD:MY_JED))
    ALLOCATE (CNVTOP(MY_ISD:MY_IED,MY_JSD:MY_JED))
    ALLOCATE (RSWTOA(MY_ISD:MY_IED,MY_JSD:MY_JED))
    ALLOCATE (RLWTOA(MY_ISD:MY_IED,MY_JSD:MY_JED))
!
    ALLOCATE (HBM2(IM,MY_JSD:MY_JED))
    ALLOCATE   (SM(IM,MY_JSD:MY_JED))
!
    ALLOCATE (DETA(1:LM))
!------------------------------- 
! LOOP OVER ALL THE OUTPUT TIMES
!-------------------------------
    666 CONTINUE
!
    IF (ME == 0) THEN
        CALL MPI_RECV(IHOUR, 1, MPI_INTEGER, 0, 0, MPI_COMM_INTER, STATUS, IER)
        PRINT*,' IHOUR IN QUILT = ', IHOUR
    END IF
!
    CALL MPI_BCAST(IHOUR, 1, MPI_INTEGER, 0, MPI_COMM_COMP, IER)
!
    IF (IHOUR == -999) GOTO 667
!
    DO 200 IXXX=1,JEND(ME)-JSTA(ME)+1
!--------------------------------------------------------------------- 
! RECEIVE ALL THE DATA FROM CHKOUT FROM THE APPROPRIATE FORECAST TASKS
!---------------------------------------------------------------------
        IST = TIMEF()
!
        CALL MPI_RECV(BUF, IBUFMAX, MPI_REAL, MPI_ANY_SOURCE, IHOUR, MPI_COMM_INTER, STATUS, IER)
!
        IPE = STATUS(MPI_SOURCE)
!    
        IF (IER /= 0) THEN
            PRINT*,' ERROR FROM MPI_REC = ', IER
        END IF
!
        IS = MY_IS_GLB_A(IPE)
        IE = MY_IE_GLB_A(IPE)
        JS = MY_JS_GLB_A(IPE)
        JE = MY_JE_GLB_A(IPE)
!
        LEN_CH = (IE - IS + 1) * (JE - JS + 1)
!-----------------------------------------------------------------------     
! EXTRACT RECORD LENGTH - LETS KEEP THIS BECAUSE IT IS POTENTIALLY HANDY
!----------------------------------------------------------------------- 
        CALL DECOAL(IDUM               , -1    )   
        CALL DECOAL(RUN                ,  1    )
        CALL DECOAL(IDAT               ,  3    )
        CALL DECOAL(IHRST              ,  1    )
        CALL DECOAL(NTSD               ,  1    )
        CALL DECOAL(LABEL              ,  8    )
        CALL DECOAL( PDOMG(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(RESOMG(IS:IE,JS:JE), LEN_CH)
!    
        DO L=1,LM
            CALL DECOAL(OMGALF(IS:IE,JS:JE,L), LEN_CH)
        END DO
!    
        CALL DECOAL(RUN             , 1      )
        CALL DECOAL(IDAT            , 3      )
        CALL DECOAL(IHRST           , 1      )
        CALL DECOAL(NTSD            , 1      )
        CALL DECOAL(LABEL           , 8      )
        CALL DECOAL(FIRST           , 1      )
        CALL DECOAL(IOUT            , 1      )
        CALL DECOAL(NSHDE           , 1      )
        CALL DECOAL(PD(IS:IE,JS:JE) , LEN_CH )
        CALL DECOAL(RES(IS:IE,JS:JE), LEN_CH )
        CALL DECOAL(FIS(IS:IE,JS:JE), LEN_CH )
        CALL DECOAL(PDB             , LB*2   )
        CALL DECOAL(TB              , LB*LM*2)
        CALL DECOAL(QB              , LB*LM*2)
        CALL DECOAL(UB              , LB*LM*2)
        CALL DECOAL(VB              , LB*LM*2)
        CALL DECOAL(Q2B             , LB*LM*2)
        CALL DECOAL(CWMB            , LB*LM*2)
!    
        DO L=1,LM
            CALL DECOAL(    T(IS:IE,JS:JE,L), LEN_CH)
            CALL DECOAL(    Q(IS:IE,JS:JE,L), LEN_CH)
            CALL DECOAL(    U(IS:IE,JS:JE,l), LEN_CH)
            CALL DECOAL(    V(IS:IE,JS:JE,l), LEN_CH)
            CALL DECOAL(   Q2(IS:IE,JS:JE,L), LEN_CH)
            CALL DECOAL( TTND(IS:IE,JS:JE,L), LEN_CH)
            CALL DECOAL(  CWM(IS:IE,JS:JE,L), LEN_CH)
            CALL DECOAL(TRAIN(IS:IE,JS:JE,L), LEN_CH)
            CALL DECOAL(TCUCN(IS:IE,JS:JE,L), LEN_CH)
        END DO
!    
        CALL DECOAL(RUN                , 1     )
        CALL DECOAL(IDAT               , 3     )
        CALL DECOAL(IHRST              , 1     )
        CALL DECOAL(NTSD               , 1     )
        CALL DECOAL(LABEL              , 8     )
        CALL DECOAL( RSWIN(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(RSWOUT(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(    TG(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(    Z0(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(  AKMS(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(  CZEN(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(  AKHS(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(   THS(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(    QS(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(  TWBS(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(  QWBS(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(  HBOT(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(CFRACL(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(  THZ0(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(   QZ0(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(   UZ0(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(   VZ0(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL( USTAR(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(  HTOP(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(CFRACM(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(   SNO(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(    SI(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(CLDEFI(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(    RF(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(  PSLP(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL( CUPPT(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(CFRACH(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(SOILTB(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(SFCEXC(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(SMSTAV(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(SMSTOT(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(GRNFLX(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(PCTSNO(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL( RLWIN(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL( RADOT(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(CZMEAN(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL( SIGT4(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(   U00(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(UL                 , 2*LM  )
        CALL DECOAL(    LC(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(    SR(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(RUN                , 1     )
        CALL DECOAL(IDAT               , 3     )
        CALL DECOAL(IHRST              , 1     )
        CALL DECOAL(NTSD               , 1     )
        CALL DECOAL(LABEL              , 8     )
        CALL DECOAL(  PREC(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(ACPREC(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(ACCLIQ(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(CUPREC(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(ACFRCV(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(NCFRCV(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(ACFRST(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(NCFRST(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(ACSNOW(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(ACSNOM(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(SSROFF(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(bgroff(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(SFCSHX(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(SFCLHX(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(SUBSHX(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(SNOPCX(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(SFCUVX(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(SFCEVP(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(POTEVP(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL( ASWIN(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(ASWOUT(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(ASWTOA(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL( ALWIN(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(ALWOUT(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(ALWTOA(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(ARDSW              , 1     )
        CALL DECOAL(ARDLW              , 1     )
        CALL DECOAL(ASRFC              , 1     )
        CALL DECOAL(AVRAIN             , 1     )
        CALL DECOAL(AVCNVC             , 1     )
        CALL DECOAL(  TH10(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(   Q10(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(   U10(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(   V10(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(TSHLTR(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(QSHLTR(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(PSHLTR(IS:IE,JS:JE), LEN_CH)
!---------
! SM v100 M
!---------
        CALL DECOAL( TH100(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(  Q100(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(  U100(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(  V100(IS:IE,JS:JE), LEN_CH)
!---------
! SM v100M
!---------
! Lyra GSM Wind stress
        CALL DECOAL(XMOMFLUX(IS:IE,JS:JE),LEN_CH)
        CALL DECOAL(YMOMFLUX(IS:IE,JS:JE),LEN_CH)
! Lyra GSM Wind stress
!
        CALL DECOAL(    SMC(IS:IE,JS:JE,1:NSOIL), LEN_CH*NSOIL)
        CALL DECOAL(    CMC(IS:IE,JS:JE)        , LEN_CH)
        CALL DECOAL(    STC(IS:IE,JS:JE,1:NSOIL), LEN_CH*NSOIL)
        CALL DECOAL(   SH2O(IS:IE,JS:JE,1:NSOIL), LEN_CH*NSOIL)
        CALL DECOAL( ALBEDO(IS:IE,JS:JE)        , LEN_CH)
        CALL DECOAL( POTFLX(IS:IE,JS:JE)        , LEN_CH)
        CALL DECOAL(  TLMIN(IS:IE,JS:JE)        , LEN_CH)
        CALL DECOAL(  TLMAX(IS:IE,JS:JE)        , LEN_CH)
!
!Lyra GSM Max wind
        CALL DECOAL(  MAXWU(IS:IE,JS:JE)        , LEN_CH)
        CALL DECOAL(  MAXWV(IS:IE,JS:JE)        , LEN_CH)
!Lyra 
!
        CALL DECOAL(ACUTIM,1)
        CALL DECOAL(ARATIM,1)
        CALL DECOAL(APHTIM,1)
        CALL DECOAL(NHEAT ,1)
        CALL DECOAL(NPHS  ,1)
        CALL DECOAL(NCNVC ,1)
        CALL DECOAL(NPREC ,1)
        CALL DECOAL(NRDSW ,1)
        CALL DECOAL(NRDLW ,1)
        CALL DECOAL(NSRFC ,1)
        CALL DECOAL(TPH0D ,1)
        CALL DECOAL(TLM0D ,1)
        CALL DECOAL(RESTRT,1)
!    
        DO L=1,LM
            CALL DECOAL(RSWTT(IS:IE,JS:JE,L), LEN_CH)
            CALL DECOAL(RLWTT(IS:IE,JS:JE,L), LEN_CH)
        END DO
!    
        CALL DECOAL(CNVBOT(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(CNVTOP(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(RSWTOA(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(RLWTOA(IS:IE,JS:JE), LEN_CH)
!-----------------------
! ADDED FOR SLPSIG STUFF
!-----------------------
        CALL DECOAL(HBM2(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL(  SM(IS:IE,JS:JE), LEN_CH)
        CALL DECOAL( SPL(1:LSM), LSM)
        CALL DECOAL(DETA(1:LM) , LM)
        CALL DECOAL(PT         , 1)
        CALL DECOAL(SPLINE     , 1)
!       WRITE(6,*) 'DECOALED IN QUILT THAT SPLINE= ', SPLINE
!
200 END DO
!-----------------------------------------------------
! BEFORE WRITING OUT THE RESTRT FILE, COMPUTE THE MSLP
!-----------------------------------------------------
    IF (SIGMA) THEN
        IF (SPLINE) THEN
            PRINT*,'CALLING SLPSIGSPLINE'
            CALL SLPSIGSPLINE(PD, FIS, T, Q, SPL, LSM, DETA, PT, PSLP)
        ELSE
            PRINT*,'CALLING SLPSIG'
            CALL SLPSIG(PD, FIS, SM, T, Q, CWM, HBM2, U00, SPL, LSM, UL, DETA, PT, PSLP)
        END IF
    ELSE
        PRINT*,'CALLING SLP'
!
        CALL SLP(NHB, PD, RES, FIS, T, Q, NTSD, NEST, PSLP)
!
    END IF
!---------------------------------- 
! WRITE OUT THE GLOBAL RESTRT FILE.
!---------------------------------- 
!
!--------------------------------------------------- 
! GENERATE THE NAME OF THE GLOBAL OUTPUT RESTRT FILE
!--------------------------------------------------- 
    ENVAR = ' '
    CALL GETENV("RSTFNL", ENVAR)
    CALL GETENV("TMMARK", RESTHR)
    RESTHR='t00s'
    KPATH = INDEX(ENVAR,' ') -1
!
    IF (KPATH <= 0) KPATH = LEN(ENVAR)

    IF (RESTHR == '    ') THEN
        WRITE(RSTFIL2,280) IHOUR
        280 FORMAT('RESTRT', I6.6)
    ELSE
        WRITE(RSTFIL2,285) IHOUR, RESTHR
        285 FORMAT('RESTRT',I6.6,'.',A4)
    END IF
!
    KRST = INDEX(RSTFIL2,' ') -1
    IF (KRST <= 0) KRST = LEN(RSTFIL2)
!------------------------------------- 
! OPEN UNIT TO THE GLOBAL RESTART FILE
!------------------------------------- 
    CLOSE(LRSTRT2)
!
    IF (ENVAR(1:4) == BLANK) THEN
        OPEN(UNIT=LRSTRT2, FILE=RSTFIL2, FORM='UNFORMATTED', IOSTAT=IER)
    ELSE
        FNAME = ENVAR(1:KPATH) // RSTFIL2(1:KRST)
        OPEN(UNIT=LRSTRT2, FILE=FNAME, FORM='UNFORMATTED', IOSTAT=IER)
    END IF
!
    IF (LME) WRITE(LRSTRT2) RUN, IDAT, IHRST, NTSD, LABEL, IHOUR
!
    CALL COLLECT(PDOMG , DUM1)
    CALL COLLECT(RESOMG, DUM2)
!
    IF (LME) WRITE(LRSTRT2) DUM1, DUM2
!
    DO L=1,LM
        CALL COLLECT(OMGALF(:,:,L), DUM1)
        IF (LME) WRITE(LRSTRT2) DUM1
    END DO
!
    IF (LME) WRITE(LRSTRT2) RUN, IDAT, IHRST, NTSD, LABEL, IHOUR, FIRST, IOUT, NSHDE
!
    CALL COLLECT(PD ,DUM1)
    CALL COLLECT(RES,DUM2)
    CALL COLLECT(FIS,DUM3)
!
    IF (LME) WRITE(LRSTRT2) DUM1, DUM2, DUM3
    IF (LME) WRITE(LRSTRT2) PDB , TB  , QB  , UB, VB, Q2B, CWMB
!
    DO L=1,LM
        CALL COLLECT     (T(:,:,L),DUM1)
        IF (LME) WRITE(LRSTRT2) DUM1
        CALL COLLECT     (Q(:,:,L),DUM1)
        IF (LME) WRITE(LRSTRT2) DUM1
        CALL COLLECT     (U(:,:,L),DUM1)
        IF (LME) WRITE(LRSTRT2) DUM1
        CALL COLLECT     (V(:,:,L),DUM1)
        IF (LME) WRITE(LRSTRT2) DUM1
        CALL COLLECT    (Q2(:,:,L),DUM1)
        IF (LME) WRITE(LRSTRT2) DUM1
        CALL COLLECT  (TTND(:,:,L),DUM1)
        IF (LME) WRITE(LRSTRT2) DUM1
        CALL COLLECT   (CWM(:,:,L),DUM1)
        IF (LME) WRITE(LRSTRT2) DUM1
        CALL COLLECT (TRAIN(:,:,L),DUM1)
        IF (LME) WRITE(LRSTRT2) DUM1
        CALL COLLECT (TCUCN(:,:,L),DUM1)
        IF (LME) WRITE(LRSTRT2) DUM1
    END DO
!
    CALL COLLECT(RSWIN , DUM1)
    CALL COLLECT(RSWOUT, DUM2)
    CALL COLLECT(TG    , DUM3)
    CALL COLLECT(Z0    , DUM4)
    CALL COLLECT(AKMS  , DUM5)
    CALL COLLECT(CZEN  , DUM6)
!
    IF (LME) WRITE(LRSTRT2) RUN , IDAT, IHRST, NTSD, LABEL, IHOUR,                                 &
    &                       DUM1, DUM2, DUM3 , DUM4, DUM5 , DUM6
!
    CALL COLLECT(AKHS  ,DUM1)
    CALL COLLECT(THS   ,DUM2)
    CALL COLLECT(QS    ,DUM3)
    CALL COLLECT(TWBS  ,DUM4)
    CALL COLLECT(QWBS  ,DUM5)
    CALL COLLECT(HBOT  ,DUM6)
    CALL COLLECT(CFRACL,DUM7)
!
    IF (LME) WRITE(LRSTRT2) DUM1, DUM2, DUM3, DUM4, DUM5, DUM6, DUM7
!
    CALL COLLECT(THZ0  ,DUM1)
    CALL COLLECT(QZ0   ,DUM2)
    CALL COLLECT(UZ0   ,DUM3)
    CALL COLLECT(VZ0   ,DUM4)
    CALL COLLECT(USTAR ,DUM5)
    CALL COLLECT(HTOP  ,DUM6)
    CALL COLLECT(CFRACM,DUM7)
!
    IF (LME) WRITE(LRSTRT2) DUM1, DUM2, DUM3, DUM4, DUM5, DUM6, DUM7
!
    CALL COLLECT(SNO   ,DUM1)
    CALL COLLECT(SI    ,DUM2)
    CALL COLLECT(CLDEFI,DUM3)
    CALL COLLECT(RF    ,DUM4)
    CALL COLLECT(PSLP  ,DUM5)
    CALL COLLECT(CUPPT ,DUM6)
    CALL COLLECT(CFRACH,DUM7)
!
    IF (LME) WRITE(LRSTRT2) DUM1, DUM2, DUM3, DUM4, DUM5, DUM6, DUM7
!
    CALL COLLECT(SOILTB,DUM1)
    CALL COLLECT(SFCEXC,DUM2)
    CALL COLLECT(SMSTAV,DUM3)
    CALL COLLECT(SMSTOT,DUM4)
    CALL COLLECT(GRNFLX,DUM5)
    CALL COLLECT(PCTSNO,DUM6)
!
    IF (LME) WRITE(LRSTRT2) DUM1, DUM2, DUM3, DUM4, DUM5, DUM6
!
    CALL COLLECT(RLWIN ,DUM1)
    CALL COLLECT(RADOT ,DUM2)
    CALL COLLECT(CZMEAN,DUM3)
    CALL COLLECT(SIGT4 ,DUM4)
!
    IF (LME) WRITE(LRSTRT2) DUM1, DUM2, DUM3, DUM4
!
    CALL COLLECT(U00,DUM1)
    CALL COLLECT(LC ,DUM2)
    CALL COLLECT(SR ,DUM3)
!
    IF (LME) WRITE(LRSTRT2) DUM1, UL, DUM2, DUM3
!
    CALL COLLECT(PREC  ,DUM1)
    CALL COLLECT(ACPREC,DUM2)
    CALL COLLECT(ACCLIQ,DUM3)
    CALL COLLECT(CUPREC,DUM4)
!
    IF (LME) WRITE(LRSTRT2) RUN, IDAT, IHRST, NTSD, LABEL, IHOUR, DUM1, DUM2, DUM3, DUM4
!
    CALL COLLECT(ACFRCV,DUM1)
    CALL COLLECT(NCFRCV,DUM2)
    CALL COLLECT(ACFRST,DUM3)
    CALL COLLECT(NCFRST,DUM4)
!
    IF (LME) WRITE(LRSTRT2) DUM1, DUM2, DUM3, DUM4
!
    CALL COLLECT(ACSNOW,DUM1)
    CALL COLLECT(ACSNOM,DUM2)
    CALL COLLECT(SSROFF,DUM3)
    CALL COLLECT(BGROFF,DUM4)
!
    IF (LME) WRITE(LRSTRT2) DUM1, DUM2, DUM3, DUM4
!
    CALL COLLECT(SFCSHX,DUM1)
    CALL COLLECT(SFCLHX,DUM2)
    CALL COLLECT(SUBSHX,DUM3)
    CALL COLLECT(SNOPCX,DUM4)
    CALL COLLECT(SFCUVX,DUM5)
    CALL COLLECT(SFCEVP,DUM6)
    CALL COLLECT(POTEVP,DUM7)
!
    IF (LME) WRITE(LRSTRT2) DUM1, DUM2, DUM3, DUM4, DUM5, DUM6, DUM7
!
    CALL COLLECT(ASWIN ,DUM1)
    CALL COLLECT(ASWOUT,DUM2)
    CALL COLLECT(ASWTOA,DUM3)
    CALL COLLECT(ALWIN ,DUM4)
    CALL COLLECT(ALWOUT,DUM5)
    CALL COLLECT(ALWTOA,DUM6)
!
    IF (LME) WRITE(LRSTRT2) DUM1, DUM2, DUM3, DUM4, DUM5, DUM6
!
    IF (LME) WRITE(LRSTRT2) ARDSW, ARDLW, ASRFC, AVRAIN, AVCNVC
!
    CALL COLLECT(TH10  ,DUM1)
    CALL COLLECT(Q10   ,DUM2)
    CALL COLLECT(U10   ,DUM3)
    CALL COLLECT(V10   ,DUM4)
    CALL COLLECT(TSHLTR,DUM5)
    CALL COLLECT(QSHLTR,DUM6)
    CALL COLLECT(PSHLTR,DUM7)
!---------
! SM V100M
!---------
    CALL COLLECT(TH100,DUM8 )
    CALL COLLECT(Q100 ,DUM9 )
    CALL COLLECT(U100 ,DUM10)
    CALL COLLECT(V100 ,DUM11)
!---------
! SM V100M
!---------          
! Lyra Wind stress
    CALL COLLECT(XMOMFLUX,DUM12)
    CALL COLLECT(YMOMFLUX,DUM13)
! Lyra Wind stress
!
    IF (LME) WRITE(LRSTRT2) DUM1, DUM2, DUM3, DUM4, DUM5, DUM6, DUM7, DUM8, DUM9, DUM10, DUM11, DUM12, DUM13
!---------
! SM V100M
!---------   
    DO L=1,NSOIL
        CALL COLLECT(SMC(:,:,L), DUMS(:,:,L))
    END DO
!
    IF (LME) WRITE(LRSTRT2) DUMS
!
    CALL COLLECT(CMC, DUM1)
!
    IF (LME) WRITE(LRSTRT2) DUM1
!
    DO L=1,NSOIL
        CALL COLLECT(STC(:,:,L), DUMS(:,:,L))
    END DO
!
    IF (LME) WRITE(LRSTRT2) DUMS
!
    DO L=1,NSOIL
        CALL COLLECT(SH2O(:,:,L), DUMS(:,:,L))
    END DO
!
    IF (LME) WRITE(LRSTRT2) DUMS
!
    CALL COLLECT(ALBEDO, DUM1)
!
    IF (LME) WRITE(LRSTRT2) DUM1
!
    CALL COLLECT(POTFLX,DUM1)
    CALL COLLECT(TLMIN ,DUM2)
    CALL COLLECT(TLMAX ,DUM3)
!
! Lyra GSM Max wind
    CALL COLLECT(MAXWU ,DUM4)
    CALL COLLECT(MAXWV ,DUM5)    
! Lyra GSM Max wind
!
    IF (LME) WRITE(LRSTRT2) DUM1 , DUM2 , DUM3 , DUM4  , DUM5  , ACUTIM, ARATIM, APHTIM, NHEAT , NPHS, NCNVC,     &
    &                       NPREC, NRDSW, NRDLW, NSRFC , TPH0D , TLM0D , RESTRT
!
    DO L=1,LM
        CALL COLLECT(RSWTT(:,:,L), DUM1)
!
        IF (LME) WRITE(LRSTRT2)    DUM1
!
        CALL COLLECT(RLWTT(:,:,L), DUM1)
!
        IF (LME) WRITE(LRSTRT2)    DUM1
    END DO
!
    CALL COLLECT(CNVBOT(:,:), DUM1)
!
    IF (LME) WRITE(LRSTRT2)   DUM1
!
    CALL COLLECT(CNVTOP(:,:), DUM1)
!
    IF (LME) WRITE(LRSTRT2)   DUM1
!
    CALL COLLECT(RSWTOA(:,:), DUM1)
!
    IF (LME) WRITE(LRSTRT2)   DUM1
!
    CALL COLLECT(RLWTOA(:,:), DUM1)
!
    IF (LME) WRITE(LRSTRT2)   DUM1
!
    CLOSE(LRSTRT2)

!----------------------------------------------------------------------------------
! SM - ALTERA\E7\E3O FEITA PARA O IO - PARA N\E3O MAIS ESCREVER O OUTJOB E FCSTDONE
!----------------------------------------------------------------------------------          
    OUTJ = .TRUE. 
!
    IF (OUTJ) THEN
!              
        IF (LME) THEN
            DONE = 'DONE'
            ITAG = IHOUR
!
            WRITE(FINFIL,1190) ITAG, RESTHR
!
            1190 FORMAT('FCSTDONE',I6.6,'.',A4)
            LFINFIL = 91
            CLOSE(LFINFIL)
!
            OPEN(UNIT=LFINFIL, FILE=FINFIL, FORM='UNFORMATTED', IOSTAT=IER)
            WRITE(LFINFIL) DONE
            CLOSE(LFINFIL)
!
            IF (IER /= 0) WRITE(LIST,*)' SIGNAL SENT TO FINFIL:  DONE'
        END IF
    END IF
!
    GOTO 666
    667 CONTINUE
!
    PRINT*,' QUILT I/O SERVER SHUTTING DOWN NOW'
    END SUBROUTINE QUILT
!
!
    SUBROUTINE DECOAL(A, LEN_CH)
!--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE BUFFER
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & A(*)
!
    IF (LEN_CH < 0) THEN
        IP = 0
    END IF
!
    DO I=1,ABS(LEN_CH)
        IP   = IP + 1
        A(I) = BUF(IP)
    END DO
!
    END SUBROUTINE DECOAL

