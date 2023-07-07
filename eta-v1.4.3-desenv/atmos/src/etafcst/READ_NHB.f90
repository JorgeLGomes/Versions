    SUBROUTINE READ_NHB(NHB)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE READ_NHB
!>
!> SUBPROGRAM: READ_NHB - READ AND DISTRIBUTE NHB FILE
!> PROGRAMMER: BLACK
!> ORG: W/NP2
!> DATE: 98-10-22
!>
!> ABSTRACT:
!> READ_NHB READS IN QUANTITIES FROM THE NHB FILE AND THEN DISTRIBUTES THEM TO THE OTHER NODES/PES 
!> AS NEEDED
!>
!> PROGRAM HISTORY LOG:
!> 97-??-??  MEYS       - ORIGINATOR
!> 97-08-??  BLACK      - REWROTE FOR BENCHMARK
!> 98-??-??  TUCCILLO   - MODIFIED FOR SINGLE OR DOUBLE PRECISION
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NHB - FILE NUMBER OF THE NHB FILE
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
!>              CLDWTR
!>              CNVCLD
!>              CONTIN
!>              CTLBLK
!>              DYNAM
!>              F77KINDS
!>              GLB_TABLE
!>              INDX
!>              LOOPS
!>              MAPOT
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              PARMETA
!>              PARMSOIL
!>              PARM_TBL
!>              PHYS
!>              PVRBLS
!>              SLOPES
!>              SOIL
!>              TEMPCOM
!>              TOPO     
!>              VRBLS     
!>
!> DRIVER     : INIT
!>              INITS
!>
!> CALLS      : DSTRB
!>              IDSTRB
!>              MPI_BARRIER
!>              MPI_BCAST
!>--------------------------------------------------------------------------------------------------
    USE ACMCLD
    USE ACMCLH
    USE ACMPRE
    USE ACMRDL
    USE ACMRDS
    USE ACMSFC
    USE BOCO
    USE CLDWTR
    USE CNVCLD
    USE CONTIN
    USE CTLBLK
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
    USE SLOPES
    USE SOIL
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND), PARAMETER :: G      =    9.8
    REAL   (KIND=R4KIND), PARAMETER :: CM1    = 2937.4
    REAL   (KIND=R4KIND), PARAMETER :: CM2    =    4.9283
    REAL   (KIND=R4KIND), PARAMETER :: CM3    =   23.5518
    REAL   (KIND=R4KIND), PARAMETER :: EPS    =    0.622
    REAL   (KIND=R4KIND), PARAMETER :: Q2INI  =     .50
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ2  =    0.12
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ   =    2.E-12
    REAL   (KIND=R4KIND), PARAMETER :: EPSWET =    0.0
    REAL   (KIND=R4KIND), PARAMETER :: Z0LAND =     .10
    REAL   (KIND=R4KIND), PARAMETER :: Z0SEA  =     .001
    REAL   (KIND=R4KIND), PARAMETER :: FCM    =     .00001
    REAL   (KIND=R4KIND), PARAMETER :: DTR    =    0.1745329E-1
    REAL   (KIND=R4KIND), PARAMETER :: H360   =  360.0
    REAL   (KIND=R4KIND), PARAMETER :: H1905  =  190.5
    REAL   (KIND=R4KIND), PARAMETER :: H105   =  105.0
!
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM   = IM * JM - JM / 2
    REAL   (KIND=R4KIND), PARAMETER :: JMP1   = JM + 1
    REAL   (KIND=R4KIND), PARAMETER :: IMT    =  2 * IM - 1
!------------------   
! DECLARE VARIABLES
!------------------ 
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & RUNB
! 
    CHARACTER(32)                                                                               ::&
    & LABEL
!
    CHARACTER(40)                                                                               ::&
    & CONTRL  , FILALL  , FILMST  , FILTMP  , FILTKE  , FILUNV  , FILCLD  , FILRAD  , FILSFC
!
    INTEGER(KIND=I4KIND), DIMENSION(3)                                                          ::&
    & IDATB      
!
    INCLUDE "mpif.h"
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!
    INTEGER(KIND=I4KIND), DIMENSION(IM, JM)                                                     ::&
    & ILMH
!
#ifdef DP_REAL
    INTEGER(KIND=I8KIND), DIMENSION(IM, JM)                                                     ::&
    & ITEMPX
!
    INTEGER(KIND=I8KIND)                                                                        ::&
    & NFCSTX  , NBCX    , LISTX   , IDTADX  , KHLAX   , KHHAX   , KVLAX   , KVHAX   , KHL2X   ,   &
    & KHH2X   , KVL2X   , KVH2X   , IXMX    , IYMX
!
    LOGICAL(KIND=L8KIND)                                                                        ::&
    & SIGMAX
#endif
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & NHB
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NMAP    , NRADSH  , NRADLH  , KVL2    , KHLA    , KHHA    , KVLA    , KVHA    , KHL2    ,   &
    & KVH2    , IRTN    , I       , J       , JG      , IHL     , IHH     , IG      , K       ,   &
    & LMHN    , LMHS    , LMHE    , LMHW    , KHH2    , N       , NSST_UP , ISST_UP     
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & TSTART  , TEND    , TCP     , TDDAMP 
!------------------ 
! DECLARE NAMELIST.
!------------------ 
    NAMELIST /FCSTDATA/                                                                           &
    & TSTART, TEND , TCP  , RESTRT, SUBPOST, NMAP, TSHDE, SPL, NPHS, NCNVC, NRADSH, NRADLH, TPREC,&
    & THEAT , TCLOD, TRDSW, TRDLW , TSRFC
!
    IF (MYPE == 0) THEN
!
        OPEN(UNIT=NHB, FORM='UNFORMATTED', FILE='CNST.file')
!
#ifdef DP_REAL
!
        READ(NHB) NFCSTX, NBCX , LISTX, DT, IDTADX, SIGMAX, KHLAX, KHHAX, KVLAX, KVHAX, KHL2X,    &
    &             KHH2X , KVL2X, KVH2X
!
        NFCST = NFCSTX
        NBC   = NBCX
        LIST  = LISTX
        IDTAD = IDTADX
        SIGMA = SIGMAX
        KHLA  = KHLAX
        KHHA  = KHHAX
        KVLA  = KVLAX
        KVHA  = KVHAX
        KHL2  = KHL2X
        KHH2  = KHH2X
        KVL2  = KVL2X
        KVH2  = KVH2X
!
#else
!
        READ(NHB) NFCST, NBC, LIST, DT, IDTAD, SIGMA, KHLA, KHHA, KVLA, KVHA, KHL2, KHH2, KVL2,   &
    &             KVH2
!
#endif
!
    END IF
!
    CALL MPI_BCAST  (NFCST, 1, MPI_INTEGER, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (NBC  , 1, MPI_INTEGER, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (LIST , 1, MPI_INTEGER, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (DT   , 1, MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (IDTAD, 1, MPI_INTEGER, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (SIGMA, 1, MPI_LOGICAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BARRIER(                          MPI_COMM_COMP, IRTN)
!
    LIST = 6
!--------------- 
! DISTRIBUTE LMH
!--------------- 
    IF (MYPE == 0) THEN
!
#ifdef DP_REAL
!
        READ(NHB) ITEMPX
        ITEMP = ITEMPX
!
#else
!
        READ(NHB) ITEMP
!
#endif
!
        ILMH = ITEMP
    END IF
!
    CALL IDSTRB(ITEMP, LMH)
!---------------  
! DISTRIBUTE LMV
!---------------
    IF (MYPE == 0) THEN
!
#ifdef DP_REAL
!
        READ(NHB) ITEMPX
        ITEMP = ITEMPX
!
#else
!
        READ(NHB) ITEMP
!
#endif
!
    END IF
!
    CALL IDSTRB(ITEMP, LMV)
!----------------
! DISTRIBUTE HBM2
!----------------
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, HBM2, 1, 1, 1)
!---------------------
! FILL HBM3 ON EACH PE
!---------------------
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            HBM3(I,J) = 0.
        END DO
    END DO
!
    DO J=MYJS,MYJE
        JG = J + MY_JS_GLB - 1
!
        IF (JG >= 4 .AND. JG <= JM-3) THEN
            IHL = 2 - IHWG(JG)
            IHH = IM - 2
            DO I=MYIS,MYIE
                IG = I + MY_IS_GLB - 1
                IF (IG >= IHL .AND. IG <= IHH) HBM3(I,J) = 1.
            END DO
        END IF
!
    END DO
!---------------- 
! DISTRIBUTE VBM2
!----------------
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, VBM2, 1, 1, 1)
!----------------
! DISTRIBUTE VBM3
!----------------
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, VBM3, 1, 1, 1)
!-------------- 
! DISTRIBUTE SM
!-------------- 
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, SM, 1, 1, 1)
!---------------- 
! DISTRIBUTE SICE
!----------------
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, SICE, 1, 1, 1)
!--------------- 
! DISTRIBUTE HTM
!---------------
    DO K=1,LM
        IF (MYPE == 0) THEN
            READ(NHB) TEMP1
        END IF
        CALL DSTRB(TEMP1, HTM, 1, LM, K)
    END DO
!---------------
! DISTRIBUTE VTM
!---------------
    DO K=1,LM
        IF (MYPE == 0) THEN
            READ(NHB) TEMP1
        END IF
        CALL DSTRB(TEMP1, VTM, 1, LM, K)
    END DO
!
    IF (MYPE == 0)THEN
        READ(NHB) DY , CPGFV, EN, ENT, R, PT, TDDAMP, F4D, F4Q, EF4T, DETA, RDETA, AETA, F4Q2,    &
    &             ETA, DFL  , EM, EMT
    END IF
!
    CALL MPI_BCAST  (DY      , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (CPGFV   , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (EN      , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (ENT     , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (R       , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (PT      , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (TDDAMP  , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (F4D     , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (F4Q     , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (EF4T    , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  ( DETA(1), LM , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (RDETA(1), LM , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  ( AETA(1), LM , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  ( F4Q2(1), LM , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (  ETA(1), LP1, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (  DFL(1), LP1, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (   EM(1), JAM, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (  EMT(1), JAM, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BARRIER(                            MPI_COMM_COMP, IRTN)
!--------------  
! DISTRIBUTE DX
!--------------  
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, DX, 1, 1, 1)
!----------------- 
! DISTRIBUTE WPDAR
!-----------------
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, WPDAR, 1, 1, 1)
!-----------------
! DISTRIBUTE CPGFU
!-----------------
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, CPGFU, 1, 1, 1)
!----------------
! DISTRIBUTE CURV
!----------------
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, CURV, 1, 1, 1)
!--------------- 
! DISTRIBUTE FCP
!--------------- 
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, FCP, 1, 1, 1)
!---------------- 
! DISTRIBUTE FDIV
!---------------- 
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, FDIV, 1, 1, 1)
!--------------- 
! DISTRIBUTE FAD
!---------------
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, FAD, 1, 1, 1)
!------------- 
! DISTRIBUTE F
!------------- 
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, F, 1, 1, 1)
!----------------- 
! DISTRIBUTE DDMPU
!----------------- 
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, DDMPU, 1, 1, 1)
!----------------- 
! DISTRIBUTE DDMPV
!----------------- 
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, DDMPV, 1, 1, 1)
!--------------------
! DISTRIBUTE PT, GLAT
!--------------------
    IF (MYPE == 0) THEN
        READ(NHB) PT ,TEMP1
    END IF
!
    CALL DSTRB(TEMP1, GLAT, 1, 1, 1)
    CALL MPI_BCAST(PT , 1, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
!---------------- 
! DISTRIBUTE GLON
!----------------
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, GLON, 1, 1, 1)
!
    IF (MYPE == 0) THEN
        READ(NHB)PLQ, RDPQ, RDTHEQ, STHEQ, THE0Q
    END IF
!
    CALL MPI_BCAST  (PLQ     , 1   , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (RDPQ    , 1   , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (RDTHEQ  , 1   , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (STHEQ(1), ITBQ, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (THE0Q(1), ITBQ, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BARRIER(                             MPI_COMM_COMP, IRTN)
!
    IF (MYPE == 0) THEN
        READ(NHB) ROS, CS, DS, ROI, CI, DI, PL, THL, RDQ, RDTH, RDP, RDTHE, DETA, AETA, DFRLG,    &
    &             QS0, SQS, STHE, THE0
    END IF
!
    CALL MPI_BCAST  (ROS     , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (CS      , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (DS      , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (ROI     , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (CI      , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (DI      , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (PL      , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (THL     , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (RDQ     , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (RDTH    , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (RDP     , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (RDTHE   , 1  , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (DETA (1), LM , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (AETA (1), LM , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (DFRLG(1), LP1, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (  QS0(1), JTB, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (  SQS(1), JTB, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  ( STHE(1), ITB, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  ( THE0(1), ITB, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BARRIER(                            MPI_COMM_COMP, IRTN)
!------------------
! DISTRIBUTE MXSNAL
!------------------
    IF (MYPE == 0) THEN
        READ(NHB)TEMP1
    END IF
!
    CALL DSTRB(TEMP1, MXSNAL, 1, 1, 1)
!---------------- 
! DISTRIBUTE EPSR
!---------------- 
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, EPSR, 1, 1, 1)
!-------------- 
! DISTRIBUTE TG
!-------------- 
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, TG, 1, 1, 1)
!---------------- 
! DISTRIBUTE GFFC
!----------------
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, GFFC, 1, 1, 1)
!--------------- 
! DISTRIBUTE SST
!---------------
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, SST, 1, 1, 1)
!------------------
! DISTRIBUTE ALBASE
!------------------
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, ALBASE, 1, 1, 1)
!---------------- 
! DISTRIBUTE HDAC
!----------------
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, HDAC, 1, 1, 1)
!-----------------
! DISTRIBUTE HDACV
!-----------------
    IF (MYPE == 0) THEN
        READ(NHB)TEMP1
    END IF
!
    CALL DSTRB(TEMP1, HDACV, 1, 1, 1)
!-----------------
! DISTRIBUTE TTBLQ
!-----------------
    IF (MYPE == 0) THEN
        READ(NHB) TTBLQ
    END IF
!
    CALL MPI_BCAST  (TTBLQ(1,1), ITBQ * JTBQ, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
!
    CALL MPI_BARRIER(                                      MPI_COMM_COMP, IRTN)
!
    IF (MYPE == 0) THEN
!
#ifdef DP_REAL
!
        READ(NHB) PTBL, TTBL, R  , PT  , TSPH, WBD  , SBD  , TLM0D, TPH0D, DLMD, DPHD, CMLD,      &
    &             DP30, X1P ,Y1P , IXMX, IYMX, DETA , AETA , ETA
        IXM = IXMX
        IYM = IYMX
!
#else
!
        READ(NHB) PTBL, TTBL, R  , PT , TSPH , WBD  , SBD, TLM0D, TPH0D, DLMD, DPHD, CMLD, DP30,  &
    &             X1P , Y1P , IXM, IYM, DETA , AETA , ETA
!
#endif
!
    END IF
!
    CALL MPI_BCAST  (PTBL(1,1), ITB*JTB, MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (TTBL(1,1), JTB*ITB, MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (R        , 1      , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (PT       , 1      , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (TSPH     , 1      , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (WBD      , 1      , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (SBD      , 1      , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (TLM0D    , 1      , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (TPH0D    , 1      , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (DLMD     , 1      , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (DPHD     , 1      , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (CMLD     , 1      , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (DP30     , 1      , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (X1P      , 1      , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (Y1P      , 1      , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (IXM      , 1      , MPI_INTEGER, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (IYM      , 1      , MPI_INTEGER, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (DETA (1) , LM     , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (AETA (1) , LM     , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST  (ETA  (1) , LP1    , MPI_REAL   , 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BARRIER(                                    MPI_COMM_COMP, IRTN)
!------------------
! DISTRIBUTE IVGTYP
!------------------
    IF (MYPE == 0) THEN
!
#ifdef DP_REAL
!
        READ(NHB) ITEMPX
!
        ITEMP = ITEMPX
!
#else
!
        READ(NHB) ITEMP
!
#endif
!
    END IF
!
    CALL IDSTRB(ITEMP, IVGTYP)
!------------------
! DISTRIBUTE ISLTYP
!------------------
    IF (MYPE == 0) THEN
!
#ifdef DP_REAL
!
        READ(NHB) ITEMPX
!
        ITEMP = ITEMPX
!
#else
!
        READ(NHB) ITEMP
!
#endif
!
    END IF
!
    CALL IDSTRB(ITEMP, ISLTYP)
!------------------
! DISTRIBUTE ISLOPE
!------------------
    IF (MYPE == 0) THEN
!
#ifdef DP_REAL
!
        READ(NHB) ITEMPX
!
        ITEMP = ITEMPX
!
#else
!
        READ(NHB) ITEMP
!
#endif
!
    END IF
!
    CALL IDSTRB(ITEMP, ISLOPE)
!------------------
! DISTRIBUTE VEGFRC
!------------------
    IF (MYPE == 0) THEN
        READ(NHB) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, VEGFRC, 1, 1, 1)
!------------------
! DISTRIBUTE SLDPTH
!------------------
    IF (MYPE == 0) THEN
        READ(NHB) SLDPTH
    END IF
!
    CALL MPI_BCAST(SLDPTH(1), NSOIL, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
!------------------
! DISTRIBUTE RTDPTH
!------------------
    IF (MYPE == 0) THEN
        READ(NHB) RTDPTH
    END IF
!
    CALL MPI_BCAST(RTDPTH(1), NSOIL, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
!
    IF (MYPE == 0) THEN
        DO J=1,JM
            DO I=1,IM
                ITEMP(I,J) = 0
            END DO
        END DO
!
        DO J=8,JM-7
            DO I=4+MOD(J+1,2),IM-4
!
                ITEMP(I,J) = 0
                LMHN = ILMH(I,J+1)
                LMHS = ILMH(I,J-1)
                LMHE = ILMH(I+MOD(J,2),J)
                LMHW = ILMH(I-1+MOD(J,2),J)
!
                IF (LMHE < LMHN .AND. LMHE < LMHW .AND. LMHE < LMHS) THEN
                    ITEMP(I,J) = 1
                ELSE IF (LMHE == LMHN .AND. LMHE < LMHW .AND. LMHE < LMHS) THEN
                    ITEMP(I,J) = 2
                ELSE IF (LMHN < LMHW  .AND. LMHN < LMHS .AND. LMHN < LMHE) THEN
                    ITEMP(I,J) = 3
                ELSE IF (LMHN == LMHW .AND. LMHN < LMHS .AND. LMHN < LMHE) THEN
                    ITEMP(I,J) = 4
                ELSE IF (LMHW < LMHS  .AND. LMHW < LMHE .AND. LMHW < LMHN) THEN
                    ITEMP(I,J) = 5
                ELSE IF (LMHW == LMHS .AND. LMHW < LMHE .AND. LMHW < LMHN) THEN
                    ITEMP(I,J) = 6
                ELSE IF (LMHS < LMHE  .AND. LMHS < LMHN .AND. LMHS < LMHW) THEN
                    ITEMP(I,J) = 7
                ELSE IF (LMHS == LMHE .AND. LMHS < LMHN .AND. LMHS < LMHW) THEN
                    ITEMP(I,J) = 8
                ELSE
                    ITEMP(I,J) = 0
                END IF
            END DO
        END DO
    END IF
!
    CALL IDSTRB(ITEMP, ISLD)

    RETURN
!
    END SUBROUTINE READ_NHB

