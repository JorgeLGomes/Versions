    SUBROUTINE PRECPD
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE PRECPD
!>
!> SUBROUTINE: PRECPD - LARGE SCALE PRECIPITATION
!> PROGRAMMER: ZHAO
!> ORG: W/NP22
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> PRECPD COMPUTES THE GRID SCALE PRECIPITATION.
!>
!> PROGRAM HISTORY LOG:
!> 94-??-??  ZHAO       - ORIGINATOR
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 95-11-20  ABELES     - PARALLEL OPTIMIZATION
!> 96-03-29  BLACK      - REMOVED SCRCH COMMON
!> 96-07-18  ZHAO       - NEW WMIN CALCULATION
!> 06-09-25  BALDWIN    - NEW SR CALCULATION
!> 98-11-02  BLACK      - MODIFICATION FOR DISTRIBUTED MEMORY
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
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
!> OUTPUT FILES:
!> NONE
!>
!> USE MODULES: ACMCLH 
!>              CLDWTR
!>              CTLBLK
!>              DYNAM
!>              F77KINDS
!>              GLB_TABLE
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
!> CALLS      : SGETMO           
!>--------------------------------------------------------------------------------------------------
    USE ACMCLH
    USE CLDWTR
    USE CTLBLK
    USE DYNAM    , ONLY: DETA, AETA, PT
    USE F77KINDS
    USE GLB_TABLE
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
    REAL   (KIND=R4KIND), PARAMETER :: A1    =  610.78
    REAL   (KIND=R4KIND), PARAMETER :: A2    =   17.2693882
    REAL   (KIND=R4KIND), PARAMETER :: A3    =  273.16  
    REAL   (KIND=R4KIND), PARAMETER :: A4    =   35.86
    REAL   (KIND=R4KIND), PARAMETER :: PQ0   =  379.90516
    REAL   (KIND=R4KIND), PARAMETER :: TRESH =     .95
    REAL   (KIND=R4KIND), PARAMETER :: R     =  287.04 
    REAL   (KIND=R4KIND), PARAMETER :: C0    =    0.15 
    REAL   (KIND=R4KIND), PARAMETER :: CP    = 1004.6
    REAL   (KIND=R4KIND), PARAMETER :: ELWV  =    2.50E6 
    REAL   (KIND=R4KIND), PARAMETER :: ELIV  =    2.834E6 
    REAL   (KIND=R4KIND), PARAMETER :: ROW   =    1.E3 
    REAL   (KIND=R4KIND), PARAMETER :: G     =    9.8    
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ  =    2.E-12  
    REAL   (KIND=R4KIND), PARAMETER :: DLDT  = 2274.0
    REAL   (KIND=R4KIND), PARAMETER :: ELIW  = ELIV - ELWV 
!
    REAL   (KIND=R4KIND), PARAMETER :: ARCP  = A2  * (A3 - A4) / CP
    REAL   (KIND=R4KIND), PARAMETER :: RCP   = 1.  / CP
    REAL   (KIND=R4KIND), PARAMETER :: PQ0C  = PQ0 * TRESH
    REAL   (KIND=R4KIND), PARAMETER :: RROG  = 1.  / (ROW * G)
    REAL   (KIND=R4KIND), PARAMETER :: RROW  = 1.  / ROW
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM = IM * JM - JM / 2
    INTEGER(KIND=I4KIND), PARAMETER :: LTOP = 1
    INTEGER(KIND=I4KIND), PARAMETER :: LBOT = LM
!
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM_LOC = IDIM2 * JDIM2
    INTEGER(KIND=I4KIND), PARAMETER :: LDA      = (IDIM2 - IDIM1 + 1) * (JDIM2 - JDIM1 + 1)
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & NOZ
!
    INTEGER(KIND=I4KIND), DIMENSION(IMJM_LOC)                                                   ::&
    & IPREC   , JPREC
!
    REAL   (KIND=R4KIND), DIMENSION(LM,IDIM1:IDIM2, JDIM1:JDIM2)                                ::&
    &     T_T ,    Q_T  ,                                                                         &
    & TRAIN_T ,  HTM_T  ,                                                                         &
    &   CWM_T , TLAT_T
!
    REAL   (KIND=R4KIND), DIMENSION(IMJM_LOC)                                                   ::&
    & PRECRL1 , PRECSL1   
!
    INTEGER(KIND=I4KIND), DIMENSION(IMJM_LOC)                                                   ::&
    & IWL1
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & KE      , INIT    , MI0
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K       , N       , NPRE    , LML     , IWL     , IWK 

    REAL   (KIND=R4KIND)                                                                        ::&  
    & DTPH    , RDTPH   , TWODT   , RTWODT  , US      , EPS     , CCLIMIT , CWS     , CSM1    ,   &
    & CRS1    , CRS2    , CR      , AA2     , UTIM    , TTEMP   , WFIX    , WMIN    , HBM2K   ,   &
    & PDSL    , CLIMIT  , CONSTA  , U00IJ   , DETAL   , ULL     , AETAL   , PRECRL  , PRECSL  ,   &
    & PRAUT   , PSAUT   , PRACW   , PSACI   , ERR     , ERS     , PSM     , PSM1    , PSM2    ,   &
    & PPR     , PPS     , CPDR    , HH      , PID     , CONDE   , RCONDE  , TT      , QQ      ,   &
    & WW      , HTMK    , TTLAT   , U00KL   , WMINK   , PRECRK  , PRECSK  , TK      , QK      ,   &
    & TMT0    , TMT15   , AI      , BI      , QW      , QI      , QINT    , FI      , FIW     ,   &
    & QC      , RQ      , CCR     , RQKLL   , CWMK    , EXPF    , AA1     , AMAXCM  , TMT0K   ,   &
    & U00KLT  , AMAXRQ  , HHT     , QTEMP   , TMT0T   , TMT15T  , QWT     , QIT     , QINTT   ,   &
    & QCT     , RQT     , RQTT    , ERK     , RPRS    , ERRT    , ERST    , AMAXPS  , PRECRS  ,   &
    & PRECSS  , TOTPPT

!------------------------- 
! PREPARATORY CALCULATIONS 
!------------------------- 
    DTPH    = NPHS * DT
    RDTPH   = 1. / DTPH
    TWODT   = DTPH
    RTWODT  = 1. / TWODT
    KE      = 2.0E-5
    US      = 1.
    EPS     = 0.622E0
    CCLIMIT = 1.0E-3
    CLIMIT  = 1.0E-20
    CWS     = 0.025
    CSM1    = 5.00000E-8
    CRS1    = 5.00000E-6
    CRS2    = 6.66600E-10
    CR      = 5.0E-4
    MI0     = 5.0E-4
    AA2     = 1.25E-3
!
    AVRAIN = AVRAIN + 1.
    ARATIM = ARATIM + 1
!---------------------------------------- 
! PADDING CLOUD MIXING RATIO IF TOO SMALL
!----------------------------------------
!
!$omp parallel do
!
    DO 20 K=1,LM
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                CWM(I,J,K) = CWM(I,J,K) * HTM(I,J,K) * HBM2(I,J)
!
                IF (CWM(I,J,K) < 0.) CWM(I,J,K) = 0.
!---------------------------------------                  
! PADDING SPECIFIC HUMIDITY IF TOO SMALL
!---------------------------------------
                IF (Q(I,J,K) < EPSQ) Q(I,J,K) = EPSQ * HTM(I,J,K)
            END DO
        END DO
 20 END DO
!
    UTIM = 1.
!
!$omp parallel do 
!
    DO N=1,IMJM_LOC
           IWL1(N) = 0
        PRECRL1(N) = 0.
        PRECSL1(N) = 0.
    END DO
!----------------------------------------------
! CHOOSE THE COLUMNS WHERE PREC CAN BE PRODUCED
!----------------------------------------------
    NPRE = 0
!
    DO 35 J=MYJS2,MYJE2
        DO 35 I=MYIS,MYIE
!        
            DO K=2,LM
                TTEMP = 0.025  * (T(I,J,K) - 273.16)
                WFIX  = 0.9814 * EXP(0.01873 * K)
                WMIN  = 0.1E-3 * EXP(TTEMP) * WFIX
!
                IF (CWM(I,J,K) > WMIN) GOTO 33
!
            END DO
!        
            GOTO 35
!
         33 NPRE = NPRE + 1
!
            IPREC(NPRE) = I
            JPREC(NPRE) = J
 35 END DO
!----------------- 
! TRANSPOSE ARRAYS
!-----------------
!
!$omp parallel sections
!
!$omp section
!
    CALL SGETMO(T     , LDA, LDA, LM, T_T    ,LM)
!
!$omp section
!
    CALL SGETMO(Q     , LDA, LDA, LM, Q_T    ,LM)
!
!$omp section
!
    CALL SGETMO(CWM   , LDA, LDA, LM, CWM_T  ,LM)
!
!$omp section
!
    CALL SGETMO(HTM   , LDA, LDA, LM, HTM_T  ,LM)
!
!$omp section
!
    CALL SGETMO(TRAIN , LDA, LDA, LM, TRAIN_T,LM)
!
!$omp section
!
    CALL SGETMO(TLATGS, LDA, LDA, LM, TLAT_T ,LM)
!
!$omp end parallel sections
!
!-------------------------------------- 
! BEGINING OF PRECIPITATION CALCULATION 
!-------------------------------------- 
!
!-------------------------------------------- 
! LOOP OVER ALL POSSIBLE PRECIPITATION POINTS
!-------------------------------------------- 
!
!$omp parallel do
!$omp private (AAI      , AETAL   , AI      , AMAXCM  , AMAXPS  , AMAXRQ  , BI      , CCR     ,   &
!$omp          CONDE    , CONSTA  , CPDR    , CS      , CWMK    , DETAL   , ERK     , ERR     ,   &
!$omp          ERRT     , ERS     , ERST    , EXPF    , FI      , FIW     , HBM2K   , HH      ,   &
!$omp          HTMK     , I       , IWL     , J       , LML     , MI0     , PDSL    , PID     ,   &
!$omp          PPR      , PPS     , PRACW   , PRAUT   , PRECRK  , PRECRL  , PRECSK  , PRECSL  ,   &
!$omp          PRECSS   , PSACI   , PSAUT   , PSM     , PSM1    , PSM2    , QC      , QCT     ,   &
!$omp          QI       , QINT    , QINTT   , QIT     , QK      , QQ      , QTEMP   , QW      ,   &
!$omp          QWT      , RCONDE  , RPRS    , RQ      , RQKLL   , RQT     , RQTT    , TK      ,   &
!$omp          TMT0     , TMT0K   , TMT0T   , TMT15   , TMT15T  , TOTPPT  , TT      , TTEMP   ,   &
!$omp          TTLAT    , U00IJ   , U00KL   , U00KLT  , ULL     , WFIX    , WMINK   , WW      )
!
    DO 300 N=1,NPRE
!    
        I = IPREC(N)
        J = JPREC(N)
!
        HBM2K  = HBM2(I,J)
        PDSL   =  RES(I,J) * PD(I,J)
!
        CONSTA = PDSL / G * TWODT
        LML    = LM - LMH(I,J)
        U00IJ  = U00(I,J)
!    
        DO 180 K=2,LM
!        
            DETAL = DETA(K)
            ULL   =   UL(K)
            AETAL = AETA(K)
!
            WFIX = 0.9814 * EXP(0.01873 * K)
!        
            PRECRL = 0.
            PRECSL = 0.
            PRAUT  = 0.
            PSAUT  = 0.
            PRACW  = 0.
            PSACI  = 0.
            ERR    = 0.
            ERS    = 0.
            PSM    = 0.
            PSM1   = 0.
            PSM2   = 0.
            PPR    = 0.
            PPS    = 0.
            CPDR   = 0.
            HH     = 0.
            PID    = 0.
            IWL    = 0.
            CONDE  = 0.
            RCONDE = 0.
!
            TT    =    T_T(K,I,J)
            QQ    =    Q_T(K,I,J)
            WW    =  CWM_T(K,I,J)
            HTMK  =  HTM_T(K,I,J)
            TTLAT = TLAT_T(K,I,J)
!
            U00KL = U00IJ  + UL(K + LML) * (0.95 - U00IJ) * UTIM
            TTEMP = 0.025  * (TT - 273.16)
            WMINK = 0.1E-3 * EXP(TTEMP)  * WFIX
!------------------------------------------------------ 
! CHOOSE THE POINTS WHERE PRECIPITATION CAN BE PRODUCED 
!------------------------------------------------------ 
            PRECRK = AMAX1(0., PRECRL1(N))
            PRECSK = AMAX1(0., PRECSL1(N))
!
            HH = HTMK * HBM2K
!
            IF (WW < WMINK .AND. (PRECRK+PRECSK) == 0.) THEN
                PID = 0.
            ELSE
                PID = HH
            END IF
!---------------- 
! QW, QI AND QINT
!----------------
            IF (PID == 1.) THEN
                CONDE  = CONSTA * DETAL
                RCONDE = 1. / CONDE
                TK     = TT
                QK     = QQ
                TMT0   = (TK - 273.16)     * HH
                TMT15  = AMIN1(TMT0, -15.) * HH
                AI     = 0.008855
                BI     = 1.
!
                IF (TMT0 < -20.) THEN
                    AI = 0.007225
                    BI = 0.9674
                END IF
!            
                QW   = HH *  PQ0 / (PDSL * AETAL + PT) * EXP(HH * A2 * (TK - A3) / (TK - A4))
                QI   = QW * (BI  + AI * AMIN1(TMT0, 0.))
                QINT = QW * (1. - 0.00032 * TMT15 * (TMT15 + 15.))
!
                IF (TMT0 <= -40.) QINT = QI
!-----------------------
! ICE-WATER ID NUMBER IW 
!-----------------------
                IF (TMT0 < -15.) THEN
                    FI = QK - U00KL * QI
                    IF (FI > 0. .OR. WW > CLIMIT) THEN
                        IWL = 1
                    ELSE
                        IWL = 0
                    END IF
                END IF
!            
                IF (TMT0 < 0.0 .AND. TMT0 >= -15.0) THEN
                    IWL = 0
                    IF (IWL1(N) == 1 .AND. WW > CLIMIT) IWL = 1
                END IF
!            
                IF (TMT0 >= 0.) THEN
                    IWL = 0
                END IF
!--------------------------------
! THE SATUATION SPECIFIC HUMIDITY 
!--------------------------------
                FIW = FLOAT(IWL)
                QC  = (1. - FIW) * QINT + FIW * QI
!---------------------- 
! THE RELATIVE HUMIDITY 
!---------------------- 
                IF (QC <= 0.) THEN
                    RQ = 1.E-10
                ELSE
                    RQ = QK / QC
                END IF
!---------------------- 
! CLOUD COVER RATIO CCR
!---------------------- 
                IF (RQ <= U00KL) THEN
                    CCR = 0.
                ELSE
                    RQKLL = AMIN1(US, RQ)
                    CCR   = 1. - SQRT((US - RQKLL) / (US - U00KL))
                END IF
!---------------------------------------------------- 
! CORRECT CCR IF IT IS TOO SMALL IN LARGE CWM REGIONS
!---------------------------------------------------- 
                IF (CCR >= 0.01 .AND. CCR <= 0.2 .AND. WW >= 0.2E-3) THEN
                    CCR = AMIN1(1., WW * 1.0E3)
                END IF
            END IF
!
         60 CONTINUE
!------------------------------- 
! PRECIPITATION PRODUCTION RATES 
!       AUTO-CONVERT RATES 
!------------------------------- 
            IF (PID == 1.) THEN
                IWK  = IWL
                CWMK = AMAX1(0., WW - CLIMIT)
                MI0  = WMINK
!            
                IF (IWK == 1) THEN
                    EXPF  = EXP(0.025 * TMT0)
                    AA1   = 1.E-3 * EXPF
                    PSAUT = AA1 * AMAX1(0., CWMK - MI0)
                    CPDR  = -PSAUT * TWODT
!
                    IF (-CPDR >= CWMK) THEN
                        CPDR  = -CWMK
                        PSAUT = -CPDR * RTWODT
                    END IF
                ELSE
                    AMAXCM = AMAX1(0., CWMK - MI0)
                    PRAUT  = C0 * AMAXCM * AMAXCM
                    CPDR   = -PRAUT * TWODT
!
                    IF (-CPDR >= CWMK) THEN
                        CPDR  = -CWMK
                        PRAUT = -CPDR * RTWODT
                    END IF
                END IF
                PPR = PRAUT * CONDE
                PPS = PSAUT * CONDE
            END IF
!        
            IF (PID == 1.) THEN
                WW = CPDR * HH + WW
!
                PRECRL = PRECRL1(N) + PPR * HH
                PRECSL = PRECSL1(N) + PPS * HH
            END IF
!----------- 
! ACCRETIONS 
!-----------  
            IF (PID == 1.) THEN
                IWK  = IWL
                CWMK = WW
!
                PRECRK = AMAX1(0., PRECRL1(N))
                PRECSK = AMAX1(0., PRECSL1(N))
!
                IF (IWK == 1) THEN
                    EXPF  = EXP(0.025 * TMT0)
                    CS    = AA2    * EXPF
                    PSACI = CS     * AMAX1(0., CWMK) * PRECSK
                    CPDR  = -PSACI * TWODT
!
                    IF (-CPDR >= CWMK) THEN
                        CPDR  = -CWMK
                        PSACI = -CPDR * RTWODT
                    END IF
                ELSE
                    PRACW = CR     * AMAX1(0., CWMK) * (PRECRK + PRECSK)
                    CPDR  = -PRACW * TWODT
!
                    IF (-CPDR >= CWMK) THEN
                        CPDR  = -CWMK
                        PRACW = -CPDR * RTWODT
                    END IF
                END IF
                PPR = PRACW * CONDE
                PPS = PSACI * CONDE
            END IF
!        
            IF (PID == 1.) THEN
                WW = CPDR * HH + WW
                PRECRL = PRECRL + PPR * HH
                PRECSL = PRECSL + PPS * HH
            END IF
!------------------------------------------  
! EVAPORATION/CONDENSATION OF PRECIPITATION 
! ERR AND ERS POSITIVE--EVAPORATION
! ERR AND ERS NEGTIVE---CONDENSATION
!------------------------------------------
            IF (PID == 1.0) THEN
                QK    = QQ
                TMT0K = TMT0
!
                IF (TMT0K < -30.) TMT0K = -30.
                PRECRK = AMAX1(0., PRECRL)
                PRECSK = AMAX1(0., PRECSL)
!------------------------------------------------------------ 
! INCREASE THE EVAPORATION/CONDENSATION FOR STRONG/LIGHT PREC
!------------------------------------------------------------ 
                U00KLT = U00KL
                AMAXRQ = AMAX1(0., U00KL - RQ)
                ERR    = KE * AMAXRQ * PRECRK ** 0.5
!            
                IF (TMT0 >= 0.) THEN
                    ERS = 0.
                ELSE
                    ERS = (CRS1 + CRS2 * TMT0K) * AMAXRQ * PRECSK / U00KLT
                END IF
!            
                IF (ERR + ERS <= 1.E-20) GO TO 125
!------------------------------------  
! CORRECT IF OVER-EVAPO./COND. OCCURS
!------------------------------------ 
                HHT   = HH * TWODT
                TTEMP = TT - RCP * (ELWV * ERR + ELIV * ERS) * HHT
                QTEMP = QQ + HHT * (       ERR        + ERS)
!
                TMT0T = (TTEMP - 273.16) * HH
!
                IF (TMT0T < -30.) TMT0T = -30.
!
                TMT15T = AMIN1(TMT0T, -15.) * HH
                AI = 0.008855
                BI = 1.
!            
                IF (TMT0T < -20.) THEN
                    AI = 0.007225
                    BI = 0.9674
                END IF
!            
                QWT = HH  * PQ0 / (PDSL * AETAL + PT) * EXP(HH * A2 * (TTEMP - A3) / (TTEMP - A4))
                QIT = QWT * (BI + AI * AMIN1(TMT0T, 0.))
                QINTT=QWT * (1. - 0.00032 * TMT15T * (TMT15T + 15.))
!
                IF (TMT0T <= -40.) QINTT = QIT
                FIW = FLOAT(IWL)
                QCT = (1. - FIW) * QINTT + FIW * QIT
!            
                IF (QCT <= 1.E-10) THEN
                     RQT = 1.E-10
                    RQTT = 1.E-10
                ELSE
                    RQT  = QTEMP / QCT
                    RQTT = QQ    / QCT
                END IF
!            
                IF (RQT <= U00KL) GO TO 125
!            
                ERK  = (U00KL - RQTT) * QCT * RTWODT
                RPRS = ERK / (PRECRK + PRECSK)
!
                ERRT =PRECRK  * RPRS
                ERST = PRECSK * RPRS
!
                ERR = AMAX1(0., 0.5 * (ERR + ERRT))
                ERS = AMAX1(0., 0.5 * (ERS + ERST))
!            
            125 CONTINUE
!            
                PPR = -ERR * CONDE
                PPS = -ERS * CONDE
!            
                IF (-PPR >= PRECRK) THEN
                    PPR = -PRECRK
                    ERR = -PPR * RCONDE
                END IF
!            
                IF (-PPS >= PRECSK) THEN
                    PPS = -PRECSK
                    ERS = -PPS * RCONDE
                END IF
!            
            END IF
!        
            IF (PID == 1.) THEN
                PRECRL = PRECRL + PPR * HH
                PRECSL = PRECSL + PPS * HH
            END IF
!-------------------- 
! MELTING OF THE SNOW 
!-------------------- 
            IF (PID == 1.) THEN
                CWMK   = WW
                AMAXPS = AMAX1(0., PRECSL)
!            
                IF (TMT0 > 0.) THEN
                    PSM1 = CSM1 * TMT0 * TMT0 * AMAXPS
                    PSM2 = CWS  * CR   * CWMK * AMAXPS
                    PSM  = PSM1 + PSM2
                ELSE
                    PSM1 = 0.
                    PSM2 = 0.
                    PSM  = 0.
                END IF
!            
                PPR =  PSM * CONDE
                PPS = -PSM * CONDE
!            
                IF (-PPS >= AMAXPS) THEN
                    PPS  = -AMAXPS
                    PPR  =  AMAXPS
!
                    PSM1 = -PPS * RCONDE
                    PSM2 = 0.
                    PSM  = PSM1
                END IF
!            
            END IF
!        
            IF (PID == 1.) THEN
                PRECRL = PRECRL + PPR * HH
                PRECSL = PRECSL + PPS * HH
            END IF
!--------------- 
! UPDATE T AND Q 
!---------------
            IF (PID == 1.) THEN
                HHT= HH * TWODT
                TT = -RCP * (ELWV * ERR + ELIV * ERS + ELIW * PSM1) * HHT + TT
                QQ = (ERR + ERS) * HHT + QQ
!
                TTLAT = -RCP  * (ELWV  * ERR + ELIV * ERS + ELIW * PSM1 - ELWV * (PRAUT + PRACW)    &
    &                 -  ELIV * (PSAUT + PSACI))    * HHT + TTLAT
            END IF
!        
            IF (HH == 1.) THEN
                   IWL1(N) = IWL
                PRECRL1(N) = PRECRL
                PRECSL1(N) = PRECSL
            END IF
!--------------------------------------------------------------------------------------------------        
! ACCUMULATE LATENT HEATING DUE TO GRID-SCALE PRECIP/EVAP. SCALE BY THE RECIPROCAL OF THE PERIOD AT
! WHICH THIS ROUTINE IS CALLED.
! THIS PERIOD IS THE PHYSICS TIMESTEP.
!-------------------------------------------------------------------------------------------------- 
            TRAIN_T(K,I,J) = TRAIN_T(K,I,J) + (TT-T_T(K,I,J)) * RDTPH
                T_T(K,I,J) = TT
                Q_T(K,I,J) = QQ
              CWM_T(K,I,J) = WW
             TLAT_T(K,I,J) = TTLAT
    180 END DO
!------------------------- 
! THE PRECIPITATION ON SFC 
!------------------------- 
        PRECRS = PRECRL1(N) * RROW
        PRECSS = PRECSL1(N) * RROW
!    
         APREC(I,J) = PRECRS + PRECSS
          PREC(I,J) =   PREC(I,J) + PRECRS + PRECSS
        ACPREC(I,J) = ACPREC(I,J) + APREC(I,J)
!-------------------------------------------- 
! THE SNOW AND RAIN RATIO OF SFC PREC 
! SR IS THE RATIO OF SNOW TO THE TOTAL PRECIP 
! IF TOTAL PRECIP IS ZERO, SR IS ZERO- 
!-------------------------------------------- 
        TOTPPT = PRECRS + PRECSS
        IF (TOTPPT > 1.E-8) THEN
            SR(I,J) = PRECSS / TOTPPT
        ELSE
            SR(I,J) = 0.
        END IF
!
300 END DO
!--------------- 
! TRANSPOSE BACK
!--------------- 
!
!$omp parallel sections
!$omp section
!
    CALL SGETMO(T_T    , LM, LM, LDA, T     , LDA)
!
!$omp section 
!
    CALL SGETMO(Q_T    , LM, LM, LDA, Q     , LDA)
!
!$omp section
!
    CALL SGETMO(CWM_T  , LM, LM, LDA, CWM   , LDA)
!
!$omp section
!
    CALL SGETMO(TRAIN_T, LM, LM, LDA, TRAIN , LDA)
!
!$omp section
!
    CALL SGETMO(TLAT_T , LM, LM, LDA, TLATGS, LDA)
!
!$omp end parallel sections
!
    RETURN
!
    END SUBROUTINE PRECPD
