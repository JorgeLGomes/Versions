    SUBROUTINE GSCOND
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE GSCOND
!> 
!> SUBPROGRAM: GSCOND - GRID SCALE CONDENSATION AND EVAPORATION
!> PROGRAMMER: ZHAO
!> ORG: W/NP22
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> GSCOND COMPUTES THE GRID SCALE EVAPORATION AND CONDENSATION
!>
!> PROGRAM HISTORY LOG:
!> 94-??-??  ZHAO       - ORIGINATOR
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 95-03-28  BLACK      - ADDED EXTERNAL EDGE
!> 98-11-02  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
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
!>              PVRBLS
!>              TEMPCOM
!>              TEMPV
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
    USE DYNAM    , ONLY: AETA, PT
    USE F77KINDS
    USE GLB_TABLE
    USE LOOPS
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE PARMETA
    USE PARM_TBL
    USE PHYS
    USE PVRBLS
    USE TEMPCOM
    USE TEMPV
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
    REAL   (KIND=R4KIND), PARAMETER :: CP    = 1004.6 
    REAL   (KIND=R4KIND), PARAMETER :: ELWV  =    2.50E6
    REAL   (KIND=R4KIND), PARAMETER :: ELIV  =    2.834E6
    REAL   (KIND=R4KIND), PARAMETER :: ROW   =    1.E3
    REAL   (KIND=R4KIND), PARAMETER :: G     =    9.8 
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ  =    2.E-12
    REAL   (KIND=R4KIND), PARAMETER :: DLDT  = 2274.0
    REAL   (KIND=R4KIND), PARAMETER :: TM10  =  263.16
    REAL   (KIND=R4KIND), PARAMETER :: R     =  287.04
    REAL   (KIND=R4KIND), PARAMETER :: CPR   = CP * R
    REAL   (KIND=R4KIND), PARAMETER :: RCPR  = 1. / CPR
!
    REAL   (KIND=R4KIND), PARAMETER :: ARCP  = A2  * (A3-A4) / CP
    REAL   (KIND=R4KIND), PARAMETER :: RCP   = 1.  / CP
    REAL   (KIND=R4KIND), PARAMETER :: PQ0C  = PQ0 * TRESH
    REAL   (KIND=R4KIND), PARAMETER :: RROG  = 1.  / (ROW*G)
!
    INTEGER(KIND=R4KIND), PARAMETER :: IMJM  = IM  * JM - JM / 2
    INTEGER(KIND=R4KIND), PARAMETER :: LTOP  = 1
    INTEGER(KIND=R4KIND), PARAMETER :: LBOT  = LM
    INTEGER(KIND=R4KIND), PARAMETER :: LDA   = (IDIM2 - IDIM1 + 1) * (JDIM2 - JDIM1 + 1)
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & NOZ
!
    INTEGER(KIND=I4KIND), DIMENSION(LM)                                                         ::&
    & IW 
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & PDSL
!
    REAL   (KIND=R4KIND), DIMENSION(LM, IDIM1:IDIM2, JDIM1:JDIM2)                               ::&
    & T_T     , T0_T    , Q_T     , Q0_T    , TRAIN_T , CWM_T   , HTM_T
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & MR      , KE      , INIT
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    REAL   (KIND=R4KIND)                                                                        ::&
    & DTPH    , RDTPH   , TWODT   , RTWODT  , C0      , C1      , C2      , US      , EPS     ,   &
    & CCLIMIT , UTIM    , HBM2IJ  , U00IJ   , P0IJ    , RESIJ   , PDSLIJ  , TKL     , QKL     ,   &
    & CWMKL   , COND    , CLIMIT  , E0      , HH      , TMT0    , TMT15   , AI      , BI      ,   &
    & QW      , QI      , QINT    , U00KL   , FI      , THH     , PP      , PP0     , AT      ,   &
    & AQ      , AP      , FIW     , ELV     , QC      , RQKL    , CCR     , RQKLL   , CCRKL   ,   &
    & EC      , US00    , CCRKL1  , AA      , AB      , AC      , AD      , AE      , AF      ,   &
    & AG      , CONDK   , QTEMP   , RQTMP   , CONE0
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K       , LMHIJ   , IWKL    , LML        
#include "sp.h"
!------------------------- 
! PREPARATORY CALCULATIONS 
!------------------------- 
    DTPH    = NPHS * DT
    RDTPH   =   1. / DTPH
    TWODT   = DTPH
    RTWODT  =   1. / TWODT
    C0      =   1.5E-4
    C1      = 300.
    C2      =   0.5
    MR      =   3.0E-4
    KE      =   2.0E-5
    US      =   1.
    EPS     =   0.622
    CCLIMIT =   1.0E-3
    CLIMIT  =   1.0E-20
!-----------------------------------------------
! PADDING SPECIFIC HUMIDITY AND CWM IF TOO SMALL 
!-----------------------------------------------
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
    DO 30 K=1,LM
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                IF (Q(I,J,K) < EPSQ)       Q(I,J,K) =   EPSQ*HTM(I,J,K)
                IF (CWM(I,J,K) < CLIMIT) CWM(I,J,K) = CLIMIT*HTM(I,J,K)
            END DO
        END DO
 30 END DO
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            PDSL(I,J) = RES(I,J) * PD(I,J)
        END DO
    END DO
!
    IW(1) = 0
    UTIM  = 1.
!------------------------------------------------
! BEGINNING OF GRID-SCALE CONDENSATION/EVAP. LOOP 
!------------------------------------------------
!
!----------------- 
! TRANSPOSE ARRAYS
!----------------- 
!------- 
! OPENMP
!------- 
!
!$omp parallel sections
!$omp section
!
    CALL SGETMO(T    , LDA, LDA, LM, T_T     , LM)
    CALL SGETMO(Q    , LDA, LDA, LM, Q_T     , LM)
    CALL SGETMO(HTM  , LDA, LDA, LM, HTM_T   , LM)
    CALL SGETMO(CWM  , LDA, LDA, LM, CWM_T   , LM)
!
!$omp section
!
    CALL SGETMO(T0   , LDA, LDA, LM, T0_T    , LM)
    CALL SGETMO(Q0   , LDA, LDA, LM, Q0_T    , LM)
    CALL SGETMO(TRAIN, LDA, LDA, LM, TRAIN_T , LM)
!
!$omp end parallel sections
!
!---------------- 
! QW, QI AND QINT 
!---------------- 
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!$omp private (AA       , AB      , AC      , AD      , AE      , AF      , AG      , AI      ,   &
!$omp          AP       , AQ      , AT      , BI      , CCR     , CCRKL   , CCRKL1  , COND    ,   &
!$omp          CONDK    , CONE0   , CWMKL   , E0      , EC      , ELV     , FI      , FIW     ,   &
!$omp          HBM2IJ   , HH      , IWKL    , LMHIJ   , LML     , P0IJ    , PDSLIJ  , PP      ,   &
!$omp          PP0      , QC      , QI      , QINT    , QKL     , QTEMP   , QW      , RESIJ   ,   &
!$omp          RQKL     , RQKLL   , RQTMP   , THH     , TKL     , TMT0    , TMT15   , U00IJ   ,   &
!$omp          U00KL    , US00    )
!$omp firstprivate (IW)
!
    IF (MYPE == 1) THEN
        PRINT*,'MYIS, MYIE, MYJS2, MYJE2=', MYIS, MYIE, MYJS2, MYJE2
    END IF
!
    DO 100 J=MYJS2,MYJE2
        DO 100 I=MYIS,MYIE        
            LMHIJ  =  LMH(I,J)
            HBM2IJ = HBM2(I,J)
            U00IJ  =  U00(I,J)
            P0IJ   =   P0(I,J)
            RESIJ  =  RES(I,J)
            PDSLIJ = PDSL(I,J)
!
            DO 90 K=2,LM            
                TKL   =   T_T(K,I,J)
                QKL   =   Q_T(K,I,J)
                CWMKL = CWM_T(K,I,J)
!            
                COND  = 0.
                E0    = 0.
                LML   = LM-LMHIJ
                HH    = HTM_T(K,I,J)     * HBM2IJ
                TMT0  = (TKL-273.16)     * HH
                TMT15 = AMIN1(TMT0,-15.) * HH
                AI    = 0.008855
                BI    = 1.
!            
                IF (TMT0 < -20.) THEN
                    AI = 0.007225
                    BI = 0.9674
                END IF
!            
                QW   = HH * PQ0 / (PDSLIJ * AETA(K) + PT) * EXP(HH*A2 * (T(I,J,K) - A3)           &
    &                /                                                  (T(I,J,K) - A4))
!
                QI   = QW * (BI + AI * AMIN1(TMT0,0.))
                QINT = QW * (1. - 0.00032 * TMT15 * (TMT15 + 15.))
!
                IF (TMT0 <= -40.) QINT = QI
!----------------------- 
! ICE-WATER ID NUMBER IW 
!----------------------- 
                IF (TMT0 < -15.) THEN
                    U00KL = U00IJ    + UL(K+LML) * (0.95-U00IJ) * UTIM
                    FI    = Q(I,J,K) - U00KL     * QI
                    IF (FI > 0. .OR. CWMKL > CLIMIT) THEN
                        IW(K) = 1
                    ELSE
                        IW(K) = 0
                    END IF
                END IF
!            
                IF (TMT0 >= 0.) THEN
                    IW(K) = 0
                END IF
!            
                IF (TMT0 < 0.0 .AND. TMT0 >= -15.) THEN
                    IW(K) = 0
                    IF (IW(K-1) == 1 .AND. CWMKL > CLIMIT) IW(K) = 1
                END IF
!-------------------------------------------------------
! CONDENSATION AND EVAPORATION OF CLOUD AT, AQ AND DP/DT 
!-------------------------------------------------------
                THH   = TWODT  * HH
                PP    = PDSLIJ *            AETA(K) + PT
                PP0   = P0IJ   * RESIJ    * AETA(K) + PT
                AT    = (TKL-T0_T(K,I,J)) * RTWODT
                AQ    = (QKL-Q0_T(K,I,J)) * RTWODT
                AP    = (PP-PP0)          * RTWODT
                IWKL  = IW(K)
                U00KL = U00IJ + UL(K+LML) * (0.95-U00IJ) * UTIM
!-------------------------------- 
! THE SATUATION SPECIFIC HUMIDITY 
!-------------------------------- 
                FIW = FLOAT(IWKL)
                ELV = (1.-FIW) * ELWV + FIW * ELIV
                QC  = (1.-FIW) * QINT + FIW * QI
!---------------------- 
! THE RELATIVE HUMIDITY 
!---------------------- 
                IF (QC <= 0.) THEN
                    RQKL = 0.
                ELSE
                    RQKL = QKL / QC
                END IF
!---------------------- 
! CLOUD COVER RATIO CCR 
!---------------------- 
                IF (RQKL <= U00KL) THEN
                    CCR   = 0.
                ELSE
                    RQKLL = AMIN1(US,RQKL)
                    CCR   = 1. -SQRT((US-RQKLL) / (US-U00KL))
                END IF
!----------------------------------------------------
! CORRECT CCR IF IT IS TOO SMALL IN LARGE CWM REGIONS 
!---------------------------------------------------- 
                IF (CCR >= 0.01 .AND. CCR <= 0.2 .AND. CWMKL >= 0.2E-3) THEN
                    CCR   = AMIN1(1.,CWMKL*1.E3)
                END IF
!            
                CCRKL = CCR
!-------------------------------------------------------
! GIVE UP THIS POINT  IF NO CLOUD NOR CONDENSATION EXIST 
!-------------------------------------------------------
                IF (CCRKL <= CCLIMIT .AND. CWMKL <= CLIMIT) GOTO 90
!--------------------------- 
! EVAPORATION OF CLOUD WATER 
!---------------------------
                EC = 0.
!
                IF (CCRKL <= CCLIMIT .AND. CWMKL > CLIMIT) THEN
                    EC = QC * (U00KL-RQKL) * RTWODT
                    E0 = AMAX1(EC,0.0)
                    E0 = AMIN1(CWMKL*RTWODT,E0) * HH
                    E0 = AMAX1(0.,E0)
                END IF
!---------------------- 
! CONDENSATION OF CLOUD 
!---------------------- 
                IF (CCRKL <= 0.20 .OR. QC <= EPSQ) THEN
                    COND = 0.
                    GOTO 80
                END IF
!------------------------------------------------------
! THE EQS. FOR COND. HAS BEEN REORGANIZED TO REDUCE CPU 
!------------------------------------------------------
                US00   = US - U00KL
                CCRKL1 = 1. - CCRKL
!
                AA     = EPS     * ELV    * PP    * QKL
                AB     = CCRKL   * CCRKL1 * QC    * US00
                AC     = AB      + 0.5    * CWMKL
                AD     = AB      * CCRKL1
                AE     = CPR     * TKL    * TKL
                AF     = AE      * PP
                AG     = AA      * ELV
                AI     = CP      * AA
!
                COND   = (AC - AD) * (AF * AQ - AI * AT + AE * QKL * AP) / (AC * (AF + AG))
!---------------------------------------------- 
! CHECK AND CORRECT IF OVER CONDENSATION OCCURS 
!----------------------------------------------
                CONDK = (QKL - U00KL * QC * 0.1) * RTWODT
!
                IF (COND > CONDK) THEN
                    COND = CONDK
                END IF
!------------------------------------------------
! CHECK AND CORRECT IF SUPERSATUATION IS TOO HIGH 
!------------------------------------------------ 
                QTEMP = QKL - AMAX1(0.,(COND-E0)) * THH
                RQTMP = QTEMP / QC
!
                IF (RQTMP >= 1.10) THEN
                    COND = (QKL-1.10*QC) * RTWODT
                END IF
!
                IF (COND < 0.) THEN
                    COND = 0.
                END IF
!----------------------- 
! UPDATE OF T, Q AND CWM 
!----------------------- 
             80 CONTINUE
!
                IF (MYPE == 1 .AND. I == 25 .AND. J == 31 .AND. K == 31) THEN
                    PRINT*,'COND, E0=', COND, E0
                END IF
!
                CONE0 = COND - E0
!
                CWM_T(K,I,J) = CONE0 * THH + CWMKL         
!--------------------------------------------------------------------------------------------------
! ACCUMULATE LATENT HEATING DUE TO GRID-SCALE PRECIP/EVAP. SCALE BY THE RECIPROCAL OF THE PERIOD AT
! WHICH THIS ROUTINE IS CALLED.
! THIS PERIOD IS THE PHYSICS TIMESTEP.
!--------------------------------------------------------------------------------------------------
                    T_T(K,I,J) =  ELV   * RCP * CONE0 * THH + TKL
                TRAIN_T(K,I,J) =  ELV   * RCP * CONE0 * THH * RDTPH + TRAIN_T(K,I,J)
                    Q_T(K,I,J) = -CONE0 * THH + QKL
!
                IF (CWM_T(K,I,J) <= 0.) CWM_T(K,I,J) = 0.
         90 END DO
100 END DO
!------------------------------ 
! SAVE T, Q AND P FOR THIS STEP 
!------------------------------ 
!
!--------------------------------- 
! TRANSPOSE BACK THE NEEDED ARRAYS
!--------------------------------- 
!------- 
! OPENMP
!------- 
!
!$omp parallel sections
!$omp section
!
    CALL SGETMO(T_T    , LM, LM, LDA, T    , LDA)
    CALL SGETMO(Q_T    , LM, LM, LDA, Q    , LDA)
!
!$omp section
!
    CALL SGETMO(TRAIN_T, LM, LM, LDA, TRAIN, LDA)
    CALL SGETMO(CWM_T  , LM, LM, LDA, CWM  , LDA)
!
!$omp end parallel sections
!
!$omp parallel do
!
    DO 125 K=1,LM
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                Q0(I,J,K) = Q(I,J,K)
                T0(I,J,K) = T(I,J,K)
            END DO
        END DO
125 END DO
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            P0(I,J) = PD(I,J)
        END DO
    END DO
!
    RETURN
!
    END SUBROUTINE GSCOND
