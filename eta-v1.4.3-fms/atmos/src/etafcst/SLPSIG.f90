    SUBROUTINE SLPSIG (PD, FIS, SM, TSIG, QSIG, CWMSIG, HBM2, U00, SPL, LSL, UL, DETA, PT, PSLP)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SLPSIG
!> 
!> SUBROUTINE: SLPSIG - SLP REDUCTION
!> PROGRAMMER: BLACK
!> ORG: W/NP22
!> DATE: 99-04-22
!>
!> ABSTRACT:
!> THIS ROUTINE COMPUTES THE SEA LEVEL PRESSURE REDUCTION USING EITHER THE MESINGER
!> RELAXATION METHOD OR THE STANDARD NCEP REDUCTION FOR SIGMA COORDINATES.
!> A BY-PRODUCT IS THE SET OF VALUES FOR THE UNDERGROUND TEMPERATURES ON THE SPECIFIED 
!> PRESSURE LEVELS.
!>
!> PROGRAM HISTORY LOG:
!> 99-09-23  T BLACK  - REWRITTEN FROM ROUTINE SLP (ETA COORDINATES)
!> 00-08-17  H CHUANG - MODIFIED THE ROUTINE TO BE INCLUDED IN THE QUILT INSTEAD OF POST
!> 18-01-15  LUCCI    - MODERNIZATION OF THE CODE, INCLUDING:
!>                      * F77 TO F90/F95
!>                      * INDENTATION & UNIFORMIZATION CODE
!>                      * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                      * DOCUMENTATION WITH DOXYGEN
!>                      * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> PD     - SFC PRESSURE MINUS PTOP
!> FIS    - SURFACE GEOPOTENTIAL
!> SM     -
!> TSIG   - 
!> QSIG   -
!> CWMSIG -
!> HBM2   -
!> U00    -
!> SPL    -
!> LSL    -
!> UL     - 
!> DETA   -
!> PT     - TOP PRESSURE OF DOMAIN
!>
!> OUTPUT ARGUMENT LIST:
!> PSLP - THE FINAL REDUCED SEA LEVEL PRESSURE ARRAY
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MPPCOM
!>              PARA
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : QUILT
!>
!> CALLS      : MPI_ALLREDUCE
!>              UPDATE
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARA
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    INCLUDE "mpif.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: LMP1   = LM  + 1
    INTEGER(KIND=I4KIND), PARAMETER :: LP1    = LSM + 1
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM   = IM  * JM - JM / 2
    INTEGER(KIND=I4KIND), PARAMETER :: IM_JM  = IM  * JM
    INTEGER(KIND=I4KIND), PARAMETER :: JM2    = JM  - 2
    INTEGER(KIND=I4KIND), PARAMETER :: NFILL  = 8
    INTEGER(KIND=I4KIND), PARAMETER :: NRLX1  = 500 
    INTEGER(KIND=I4KIND), PARAMETER :: NRLX2  = 100
    INTEGER(KIND=I4KIND), PARAMETER :: KSLPD  = 1 
!
    REAL   (KIND=R4KIND), PARAMETER :: OVERRC = 1.5      
    REAL   (KIND=R4KIND), PARAMETER :: AD05   = OVERRC * 0.05
    REAL   (KIND=R4KIND), PARAMETER :: CFT0   = OVERRC - 1.
    REAL   (KIND=R4KIND), PARAMETER :: RD     = 287.04
    REAL   (KIND=R4KIND), PARAMETER :: ROG    = RD / 9.8  
    REAL   (KIND=R4KIND), PARAMETER :: PQ0    = 379.90516  
    REAL   (KIND=R4KIND), PARAMETER :: A2     =  17.2693882
    REAL   (KIND=R4KIND), PARAMETER :: A3     = 273.16
    REAL   (KIND=R4KIND), PARAMETER :: A4     =  35.86
    REAL   (KIND=R4KIND), PARAMETER :: GAMMA  =   6.5E-3
    REAL   (KIND=R4KIND), PARAMETER :: RGAMOG = GAMMA * ROG
    REAL   (KIND=R4KIND), PARAMETER :: H1M12  = 1.E-12
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LSL
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & PT
!
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED)                    , INTENT(IN)          ::&
    & PD      , FIS     , SM      , HBM2    , U00
!
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED)                                          ::& 
    & HTM2D   , TG       
!
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED, LSM)                                     ::& 
    & HTM
!
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED, LSL)                                     ::&
    & T       , Q       , FI
!
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED, LM)                , INTENT(IN)          ::&
    & TSIG    , QSIG    , CWMSIG 
!
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED, LM)                                      ::&
    & IW
!
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED, LMP1)                                    ::& 
    & ALPINT  , ZINT    , PINT
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & IWU     , IWL
!
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED)                    , INTENT(INOUT)       ::&
    & PSLP
!
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED)                                          ::& 
    & TTV     , SLPX    , P1
!
    REAL   (KIND=R4KIND), DIMENSION(LSM)                                  , INTENT(IN)          ::&
    & SPL
!
    REAL   (KIND=R4KIND), DIMENSION(LSM)                                                        ::&
    & SPLI
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                   , INTENT(IN)          ::&
    & DETA
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                                         ::&
    & RDETA   , AETA    , F4Q2  
!
    REAL   (KIND=R4KIND), DIMENSION(LMP1)                                                       ::&
    & ETA     , DFL     
!
    REAL   (KIND=R4KIND), DIMENSION(2 * LM)                               , INTENT(IN)          ::&
    & UL
!  
    REAL   (KIND=R4KIND), DIMENSION(7, 7)                                                       ::&
    & TEMPT  
!
    INTEGER(KIND=I4KIND), DIMENSION(LSM)                                                        ::&
    & KMNTM
!
    INTEGER(KIND=I4KIND), DIMENSION(IMJM, LSM)                                                  ::&
    & IMNT    , JMNT  
!
    INTEGER(KIND=I4KIND), DIMENSION(IM, MY_JSD:MY_JED)                                          ::&
    & LMH     , NL1X
!
    INTEGER(KIND=I4KIND), DIMENSION(IM_JM)                                                      ::&
    & IHOLD   , JHOLD
!
    INTEGER(KIND=I4KIND), DIMENSION(JM)                                                         ::&
    & IHE     , IHW     , IVE     , IVW
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & SIGMA   , STDRD   
!
    LOGICAL(KIND=L4KIND), DIMENSION(IM, MY_JSD:MY_JED)                                          ::&
    & DONE
!
    STDRD = .FALSE. 
!------------------------------------------- 
! CALCULATE THE I-INDEX EAST-WEST INCREMENTS
!------------------------------------------- 
    PRINT*,'IN SLPSIG'
!
    DO J=1,JM
        IHE(J) = MOD(J+1, 2)
        IHW(J) = IHE(J) - 1
        IVE(J) = MOD(J  , 2)
        IVW(J) = IVE(J) - 1
    END DO
!--------------------------------------------------------- 
! INITIALIZE ARRAYS. LOAD SLP ARRAY WITH SURFACE PRESSURE.
!--------------------------------------------------------- 
!
!$omp parallel do
!
    DO J=JSTA_I,JEND_I
        DO I=1,IM
            PSLP(I,J) = PD(I,J) + PT
             TTV(I,J) = 0.
             LMH(I,J) = 0
        END DO
    END DO
!--------------------- 
! COMPUTE ETA AND AETA
!--------------------- 
    ETA(1) = 0.0
!
    DO L=2,LM+1
        IF (L == 2) WRITE(6,*) 'L, DETA(L-1), ETA(L-1): ',L, DETA(L-1), ETA(L-1)
        ETA(L) = ETA(L-1) + DETA(L-1)
    END DO
!
    DO L=1,LM
        AETA(L) = 0.5 * (ETA(L) + ETA(L+1))
    END DO
!---------------------------------------------------------------------------------------- 
! CALCULATE SEA LEVEL PRESSURE FOR PROFILES (AND POSSIBLY FOR POSTING BY POST PROCESSOR).
! "STDRD" REFERS TO THE "STANDARD" SLP REDUCTION SCHEME.
!----------------------------------------------------------------------------------------
    IF (STDRD) GOTO 400
!--------------------------------------------------------------------------------------------------
! CREATE A 3-D "HEIGHT MASK" FOR THE SPECIFIED PRESSURE LEVELS (1 => ABOVE GROUND) AND A 2-D
! INDICATOR ARRAY THAT SAYS WHICH PRESSURE LEVEL IS THE LOWEST ONE ABOVE THE GROUND
!--------------------------------------------------------------------------------------------------
    DO 100 L=1,LSM
        SPLL = SPL(L)
!    
        DO J=JSTA_I,JEND_I
            DO I=1,IM
                PSFC = PD(I,J) + PT
                PCHK = PSFC
                IF (NFILL > 0) THEN
                    DO LL=1,NFILL
                        PCHK = PCHK - DETA(LM+1-LL) * PD(I,J)
                    END DO
                END IF
!
                IF (SM(I,J) > 0.5 .AND. FIS(I,J) < 10.) PCHK = PSFC
!
                IF (SPLL < PCHK) THEN
                    HTM(I,J,L) = 1.
                ELSE
                    HTM(I,J,L) = 0.
                    IF (L > 1 .AND. HTM(I,J,L-1) > 0.5) LMH(I,J) = L-1
                END IF
!            
                IF (L == LSM .AND. HTM(I,J,L) > 0.5)    LMH(I,J) = LSM
            END DO
        END DO
!    
100 END DO
!--------------------------------------------------- 
!  INTERPOLATE T AND Q FROM SIGMA TO PRESSURE LEVELS
!--------------------------------------------------- 
!
!--------------------------------------------------------------------------------------------------
! VALUES OF T ON THE OUTPUT PRESSURE LEVELS ABOVE GROUND MUST BE KNOWN BEFORE THEY CAN BE FILLED 
! IN BELOW GROUND IN THE ROUTINE WHICH COMPUTES THE MESINGER SEA LEVEL PRESSURE. THEREFORE DO
! ALL INTERPOLATION ABOVE GROUND NOW.
!--------------------------------------------------------------------------------------------------
!
!----------------------- 
! DEFINE ALPINT AND PINT
!-----------------------
    DO L=1,LMP1
        DO J=JSTA_I,JEND_I
            DO I=1,IM
                  PINT(I,J,L) =     PD(I,J) * ETA(L) + PT
                ALPINT(I,J,L) = LOG(PD(I,J) * ETA(L) + PT)
            END DO
        END DO
    END DO
!------------- 
! COMPUTE ZINT
!-------------
!
!$omp parallel do 
!
    DO J=JSTA_I,JEND_I
        DO I=1,IM
            ZINT(I,J,LMP1) = FIS(I,J) / 9.8
        END DO
    END DO
!------------------------------------ 
! COMPUTE VALUES FROM THE SURFACE UP.
!------------------------------------ 
    DO L=LM,1,-1
!
!$omp parallel do 
!
        DO J=JSTA_I,JEND_I
            DO I=1,IM
                ZINT(I,J,L) =    ZINT(I,J,L+1) +   TSIG(I,J,L) * (QSIG(I,J,L) * .608 + 1.) * RD &
    &                       * (ALPINT(I,J,L+1) - ALPINT(I,J,L)) / 9.8
            END DO
        END DO
    END DO
!
    IF (I == 10 .AND. J == 86) THEN
        PRINT*,'WRITING SAMPLE INPUT HEIGHT, T'
        DO L=1,LMP1
            PRINT*,'I,J,L,P,T,ZINT = ', I, J, L, PINT(I,J,L), TSIG(I,J,L), ZINT(I,J,L)
        END DO
    END IF
!--------------------------------
!  SET UP UTIM FOR THIS TIME STEP
!--------------------------------
    UTIM   = 1.
    CLIMIT = 1.E-20
!
    DO 75 L=1,LM
        IF (L == 1) THEN
!
!$omp parallel do
!
            DO J=JSTA_I,JEND_I
                DO I=1,IM
                    IW(I,J,L) = 0.
                END DO
            END DO
!
            GOTO 75
!
        END IF
!
!$omp parallel do                                                                                 
!$omp private (CWMKL    , FIQ     , HH      , LML     , PP      , QI      , QKL     , QW      ,   &
!$omp          TKL      , TMT0    , TMT15   , U00KL   ) 
!                                                                                                             
        DO 70 J=JSTA_I,JEND_I
            DO 70 I=1,IM
                LML   = LM - LM
                HH    = 1 * HBM2(I,J)
!
                TKL   =   TSIG(I,J,L)
                QKL   =   QSIG(I,J,L)
                CWMKL = CWMSIG(I,J,L)
!
                TMT0  = (TKL - 273.16)   * HH
                TMT15 = AMIN1(TMT0,-15.) * HH
!
                PP    = PD(I,J) * AETA(L) + PT
                QW    = HH * PQ0 / PP * EXP(HH * A2 * (TKL-A3) / (TKL-A4))
                QI    = QW * (1. + 0.01 * AMIN1(TMT0,0.))
                U00KL = U00(I,J) + UL(L + LML) * (0.95 - U00(I,J)) * UTIM
!            
                IF (TMT0 < -15.0) THEN
                    FIQ = QKL - U00KL * QI
                    IF (FIQ > 0. .OR. CWMKL > CLIMIT) THEN
                        IW(I,J,L) = 1.
                    ELSE
                        IW(I,J,L) = 0.
                    END IF
                END IF
!            
                IF (TMT0 >= 0.0) IW(I,J,L) = 0.
!
                IF (TMT0 < 0.0 .AND. TMT0 >= -15.0) THEN
                    IW(I,J,L) = 0.
                    IF (IW(I,J,L-1) == 1. .AND. CWMKL > CLIMIT) IW(I,J,L) = 1.
                END IF
!            
     70 END DO
 75 END DO
!
    DO 228 LP=1,LSL
        ALSL  = LOG(SPL(LP))
        NHOLD = 0
!    
        DO 125 J=JSTA_I,JEND_I
            DO 125 I=1,IM
!            
                 T(I,J,LP) = -1.E6
                 Q(I,J,LP) = -1.E6
                FI(I,J,LP) = -1.E6
!------------------------------------------------------------------------------------------------------            
! LOCATE VERTICAL INDEX OF MODEL INTERFACE JUST BELOW THE PRESSURE LEVEL TO WHICH WE ARE INTERPOLATING.
!------------------------------------------------------------------------------------------------------              
                DO 115 L=2,LM+1
                    IF (ALPINT(I,J,L) >= ALSL) THEN
                        NL1X(I,J) = L
                        NHOLD = NHOLD + 1
                        IHOLD(NHOLD) = I
                        JHOLD(NHOLD) = J
                        GOTO 125
                    ELSE IF(ALPINT(I,J,LMP1) < ALSL) THEN
                        NL1X(I,J) = LMP1
                        NHOLD = NHOLD + 1
                        IHOLD(NHOLD) = I
                        JHOLD(NHOLD) = J
                        GOTO 125
                    END IF
            115 END DO
!            
    125 END DO
!
        IF (NHOLD == 0) GOTO 228
!    
        TRF = 2. * ALSL
!--------------------------------------------------------------------------------------------------  
! NL1X WAS SET BASED ON INTERFACE PRESSURE VALUES. HOWEVER, WE ARE READY TO INTERPOLATE VALUES FROM
! MIDLAYER LOCATIONS SO ADJUST NL1X UP ONE LAYER IF NEEDED SO THAT THE INDEX OF THE MIDLAYER
! LOCATION IS IMMEDIATELY BELOW THE PRESSURE LEVEL TO WHICH INTERPOLATION IS BEING DONE.
!--------------------------------------------------------------------------------------------------
!
!$omp parallel do                                                                                 &  
!$omp private (AHF      , AHFO    , AHFQ    , AHFQ2   , AHFQC   , AHFQI   , AI      , B       ,   &
!$omp          BI       , BOM     , BQ      , BQ2     , BQC     , BQC_2   , BQI_2   , BQI     ,   &
!$omp          FAC      , GMIW    , GMIW_2  , IWL     , IWU     , LMP1    , PL      , PNL1    ,   &
!$omp          PU       , Q2A     , Q2B     , QABV    , QI      , QINT    , QL      , QSAT    ,   &
!$omp          QU       , QW      , RHU     , TABV    , TBLO    , TL      , TMT0    , TMT15   ,   &
!$omp          TU       , TVRABV  , TVRBLO  , TVRL    , TVRU    , ZL      , ZU      )
!
        DO 220 NN=1,NHOLD
            I = IHOLD(NN)
            J = JHOLD(NN)
            DIFU = ALSL - ALPINT(I,J,NL1X(I,J)-1)
            DIFL =        ALPINT(I,J,NL1X(I,J)  ) - ALSL
!
            IF (DIFU < DIFL) NL1X(I,J) = NL1X(I,J) - 1
!
            PNL1 = PINT(I,J,NL1X(I,J))
!--------------------------------------------------------------------------------------------------
! VERTICAL INTERPOLATION OF GEOPOTENTIAL, TEMPERATURE, SPECIFIC HUMIDITY, CLOUD WATER/ICE, AND TKE.
!--------------------------------------------------------------------------------------------------
!
!---------------------------------------------------- 
! EXTRAPOLATE ABOVE THE TOPMOST MIDLAYER OF THE MODEL
!----------------------------------------------------        
            IF (NL1X(I,J) == 1) THEN
                PU = PINT(I,J,2)
                ZU = ZINT(I,J,2)
                TU = 0.5 * (TSIG(I,J,1) + TSIG(I,J,2))
                QU = 0.5 * (QSIG(I,J,1) + QSIG(I,J,2))
!            
                  IWU = 0.5 * (IW(I,J,1) + IW(I,J,2))
                 TMT0 = TU - 273.16
                TMT15 = AMIN1(TMT0,-15.)
                AI = 0.008855
                BI = 1.
!
                IF (TMT0 < -20.) THEN
                    AI = 0.007225
                    BI = 0.9674
                END IF
!
                  QW = PQ0 / PU * EXP(A2 * (TU-A3) / (TU-A4))
                  QI = QW * (BI + AI * AMIN1(TMT0,0.))
                QINT = QW * (1. - 0.00032 * TMT15 * (TMT15+15.))
!
                IF (TMT0 < -15.) THEN
                    QSAT = QI
                ELSE IF(TMT0 >= 0.) THEN
                    QSAT = QINT
                ELSE
                    IF (IWU > 0.) THEN
                        QSAT = QI
                    ELSE
                        QSAT = QINT
                    END IF
                END IF
!-------------------------------------------------------            
! USE RH FOR LIQUID WATER NO MATTER WHAT.
! DELETE THE FOLLOWING LINE TO SWITCH BACK TO RH FOR ICE
!------------------------------------------------------- 
                QSAT = QW
!            
                RHU = QU / QSAT
!            
                IF (RHU > 1.) THEN
                    RHU = 1.
                     QU = RHU * QSAT
                END IF
!            
                IF (RHU < 0.01) THEN
                    RHU = 0.01
                     QU = RHU * QSAT
                END IF
!            
                  TVRU = TU     * (1. + 0.608 * QU)
                TVRABV = TVRU   * (SPL(LP) / PU) ** RGAMOG
                  TABV = TVRABV / (1. + 0.608 * QU)
!            
                 TMT0 = TABV - 273.16
                TMT15 = AMIN1(TMT0,-15.)
                   AI = 0.008855
                   BI = 1.
!
                IF (TMT0 < -20.) THEN
                    AI = 0.007225
                    BI = 0.9674
                END IF
!
                  QW = PQ0 / SPL(LP) * EXP(A2 * (TABV-A3) / (TABV-A4))
                  QI = QW * (BI + AI * AMIN1(TMT0,0.))
                QINT = QW *(1. - 0.00032 * TMT15 * (TMT15+15.))
!
                IF (TMT0 < -15.) THEN
                    QSAT = QI
                ELSE IF (TMT0 >= 0.) THEN
                    QSAT = QINT
                ELSE
                    IF (IWU > 0.) THEN
                        QSAT = QI
                    ELSE
                        QSAT = QINT
                    END IF
                END IF
!-------------------------------------------------------             
! USE RH FOR LIQUID WATER NO MATTER WHAT.
! DELETE THE FOLLOWING LINE TO SWITCH BACK TO RH FOR ICE
!------------------------------------------------------- 
                QSAT = QW         
                QABV = RHU * QSAT
                QABV = AMAX1(1.E-12,QABV)
                   B = TABV
                  BQ = QABV
                GMIW = IW(I,J,1)
                 AHF = 0.
                AHFQ = 0.
                 FAC = 0.
!            
            ELSE IF (NL1X(I,J) == LMP1) THEN
!----------------------------------------------------------------- 
! EXTRAPOLATE BELOW LOWEST MODEL MIDLAYER (BUT STILL ABOVE GROUND)
!----------------------------------------------------------------- 
                PL  = PINT(I,J,LM-1)
                ZL  = ZINT(I,J,LM-1)
!
                TL  = 0.5  * (TSIG(I,J,LM-2) + TSIG(I,J,LM-1))
                QL  = 0.5  * (QSIG(I,J,LM-2) + QSIG(I,J,LM-1))
!
                IWL = 0.5 *    (IW(I,J,LM-2) +   IW(I,J,LM-1))
!
                 TMT0 = TL - 273.16
                TMT15 = AMIN1(TMT0,-15.)
                   AI = 0.008855
                   BI = 1.
!
                IF (TMT0 < -20.) THEN
                    AI = 0.007225
                    BI = 0.9674
                END IF
!
                  QW = PQ0 / PL * EXP(A2 * (TL-A3) / (TL-A4))
                  QI = QW * (BI + AI * AMIN1(TMT0,0.))
                QINT = QW * (1. - 0.00032 * TMT15 * (TMT15+15.))
!
                IF (TMT0 < -15.) THEN
                    QSAT = QI
                ELSE IF (TMT0 >= 0.) THEN
                    QSAT = QINT
                ELSE
                    IF (IWL > 0.) THEN
                        QSAT = QI
                    ELSE
                        QSAT = QINT
                    END IF
                END IF
!-------------------------------------------------------             
! USE RH FOR LIQUID WATER NO MATTER WHAT.
! DELETE THE FOLLOWING LINE TO SWITCH BACK TO RH FOR ICE
!-------------------------------------------------------
                QSAT = QW        
                RHL  = QL / QSAT
!            
                IF (RHL > 1.) THEN
                    RHL = 1.
                     QL = RHL * QSAT
                END IF
!            
                IF (RHL < 0.01) THEN
                    RHL = 0.01
                     QL = RHL * QSAT
                END IF
!            
                  TVRL = TL     * (1. + 0.608 * QL)
                TVRBLO = TVRL   * (SPL(LP) / PL) ** RGAMOG
                  TBLO = TVRBLO / (1. + 0.608 * QL)
!            
                 TMT0 = TBLO - 273.16
                TMT15 = AMIN1(TMT0,-15.)
                   AI = 0.008855
                   BI = 1.
!
                IF (TMT0 < -20.) THEN
                    AI = 0.007225
                    BI = 0.9674
                END IF
!
                  QW = PQ0 / SPL(LP) * EXP(A2 * (TBLO-A3) / (TBLO-A4))
                  QI = QW * (BI + AI * AMIN1(TMT0,0.))
                QINT = QW * (1. - 0.00032 * TMT15 * (TMT15+15.))
!
                IF (TMT0 < -15.) THEN
                    QSAT = QI
                ELSE IF (TMT0 >= 0.) THEN
                    QSAT = QINT
                ELSE
                    IF (IWL > 0.) THEN
                        QSAT = QI
                    ELSE
                        QSAT = QINT
                    END IF
                END IF
!------------------------------------------------------- 
! USE RH FOR LIQUID WATER NO MATTER WHAT.
! DELETE THE FOLLOWING LINE TO SWITCH BACK TO RH FOR ICE
!------------------------------------------------------- 
                QSAT = QW           
                QBLO = RHL*QSAT
                QBLO = AMAX1(1.E-12,QBLO)
                   B = TBLO
                  BQ = QBLO
                 AHF = 0.
                AHFQ = 0.
                 FAC = 0.
            ELSE
!---------------------------------------------------- 
! INTERPOLATION BETWEEN NORMAL LOWER AND UPPER BOUNDS
!----------------------------------------------------
                   B = TSIG(I,J,NL1X(I,J))
                  BQ = QSIG(I,J,NL1X(I,J))
!
                 FAC = 2. * ALOG(PT + PD(I,J) * AETA(NL1X(I,J)))
!
                 AHF = (B - TSIG(I,J,NL1X(I,J)-1)) / (ALPINT(I,J,NL1X(I,J)+1)                     &
    &                -    ALPINT(I,J,NL1X(I,J)-1))
!
                AHFQ = (BQ - QSIG(I,J,NL1X(I,J) - 1)) / (ALPINT(I,J,NL1X(I,J) + 1)                &
    &                -     ALPINT(I,J,NL1X(I,J) - 1))
            END IF
!        
             T(I,J,LP) = B + AHF *   (TRF - FAC)
             Q(I,J,LP) = BQ + AHFQ * (TRF - FAC)
             Q(I,J,LP) = AMAX1(Q(I,J,LP),H1M12)
            FI(I,J,LP) = (PNL1 - SPL(LP)) / (SPL(LP) + PNL1)                                      &
    &                  * ((ALSL + ALPINT(I,J,NL1X(I,J)) -  FAC) * AHF + B)                        &
    &                  * RD * 2. + ZINT(I,J,NL1X(I,J)) * 9.8
!        
            IF (I == 86 .AND. J == 10) THEN
                PRINT*,'LP,NL1X,T,Q,FI= ', LP, NL1X(I,J), T(I,J,LP), Q(I,J,LP), FI(I,J,LP)
            END IF
    220 END DO
!
228 END DO
!--------------------------------------------------------------------------------------------------
! WE REACH THIS LINE IF WE WANT THE MESINGER ETA SLP REDUCTION BASED ON RELAXATION TEMPERATURES. 
! THE FIRST STEP IS TO FIND THE HIGHEST LAYER CONTAINING MOUNTAINS.
!--------------------------------------------------------------------------------------------------
    DO 230 L=LSM,1,-1
!    
        DO J=JSTA_I,JEND_I
            DO I=1,IM
                IF (HTM(I,J,L) < 0.5) GOTO 230
            END DO
        END DO
!    
        LHMNT = L+1
!
        GOTO 235
!
230 END DO
!
235 CONTINUE
!
    CALL MPI_ALLREDUCE(LHMNT, LXXX, 1, MPI_INTEGER, MPI_MIN,MPI_COMM_COMP, IERR)
!
    LHMNT = LXXX
!
    IF (LHMNT == LP1) THEN
        GOTO 430
    END IF
!--------------------------------------------------------
! NOW GATHER THE ADDRESSES OF ALL THE UNDERGROUND POINTS.
!--------------------------------------------------------
!
!$omp parallel do                                                                                 & 
!$omp private (KMN, KOUNT)
!
    DO 250 L=LHMNT,LSM
        KMN = 0
        KMNTM(L) = 0
        KOUNT = 0
!
        DO 240 J=JSTA_IM2,JEND_IM2
            DO 240 I=2,IM-1
                KOUNT = KOUNT+1
                IMNT(KOUNT,L) = 0
                JMNT(KOUNT,L) = 0
!
                IF (HTM(I,J,L) > 0.5) GOTO 240
!
                KMN = KMN+1
                IMNT(KMN,L) = I
                JMNT(KMN,L) = J
    240 END DO
        KMNTM(L) = KMN
250 END DO
!--------------------------------------------------------------------
! AS THE FIRST GUESS, SET THE UNDERGROUND TEMPERATURES EQUAL TO 0.0C.
!--------------------------------------------------------------------
    KMM = KMNTM(LSM)
!
!$omp parallel do                                                                                 & 
!$omp private (I        , J       , LAMP1), SHARED (T)
!
    DO 260 KM=1,KMM
        I = IMNT(KM,LSM)
        J = JMNT(KM,LSM)
        LMAP1 = LMH(I,J) + 1
        DO 260 L=LMAP1,LSM
            T(I,J,L) = 273.15
260 END DO
!-----------------------------------------------------------------------------------------
! CREATE A TEMPORARY TV ARRAY, AND FOLLOW BY SEQUENTIAL OVERRELAXATION, DOING NRLX PASSES.
!-----------------------------------------------------------------------------------------
    NRLX = NRLX1
!
!$omp parallel do                                                                                 & 
!$omp private (I        , J       , TINIT   , TTV)
!
    DO 300 L=LHMNT,LSM
!    
        DO 270 J=JSTA_I,JEND_I
            DO 270 I=1,IM
                  TTV(I,J) =   T(I,J,L)
                HTM2D(I,J) = HTM(I,J,L)
270 END DO
!--------------------------------------------------------------------------------------------------    
! FOR GRID BOXES NEXT TO MOUNTAINS, COMPUTE TV TO USE AS BOUNDARY CONDITIONS FOR THE RELAXATION
! UNDERGROUND
!--------------------------------------------------------------------------------------------------     
        CALL UPDATE(HTM2D)
!
        DO J=JSTA_IM2,JEND_IM2
            DO I=2,IM-1
                IF (HTM2D(I,J) > 0.5 .AND. HTM2D(I+IHW(J),J-1) * HTM2D(I+IHE(J),J-1) *            &
    &                                      HTM2D(I+IHW(J),J+1) * HTM2D(I+IHE(J),J+1) *            &
    &                                      HTM2D(I-1     ,J  ) * HTM2D(I+1     ,J  ) *            &
    &                                      HTM2D(I       ,J-2) * HTM2D(I       ,J+2) < 0.5) THEN 
                    TTV(I,J) = T(I,J,L) * (1. + 0.608 * Q(I,J,L))
                END IF
            END DO
        END DO
!
        KMM = KMNTM(L)
!
        DO 285 N=1,NRLX
            CALL UPDATE(TTV)
!
            DO 280 KM=1,KMM
                I = IMNT(KM,L)
                J = JMNT(KM,L)
                TINIT = TTV(I,J)
                TTV(I,J) = AD05 * (4. * (TTV(I+IHW(J),J-1) + TTV(I+IHE(J),J-1)                    &
    &                    +               TTV(I+IHW(J),J+1) + TTV(I+IHE(J),J+1))                   &
    &                    +               TTV(I-1     ,J  ) + TTV(I+1     ,J  )                    &
    &                    +               TTV(I       ,J-2) + TTV(I       ,J+2))                   &
    &                    - CFT0 *        TTV(I       ,J  )
!
                IF (KMM == 1 .AND. N == 500) THEN
                    PRINT*,'L,N,TTV NEIGHBOR', L, N, TTV(I+IHW(J),J-1), TTV(I+IHE(J),J-1),        &
    &                                                TTV(I+IHW(J),J+1), TTV(I+IHE(J),J+1),        &
    &                                                TTV(I-1     ,J  ), TTV(I+1     ,J  ),        &
    &                                                TTV(I       ,J-2), TTV(I       ,J+2),        &
    &                                                TTV(I       ,J  )
                END IF
!
        280 END DO
!        
    285 END DO
!
        DO 290 KM=1,KMM
            I = IMNT(KM,L)
            J = JMNT(KM,L)
            T(I,J,L) = TTV(I,J)
    290 END DO
!
300 END DO
!--------------------------------------------------------------------------------------------------
! CALCULATE THE SEA LEVEL PRESSURE AS PER THE NEW SCHEME.
! INTEGRATE THE HYDROSTATIC EQUATION DOWNWARD FROM THE GROUND THROUGH EACH OUTPUT PRESSURE LEVEL 
! (WHERE TV IS NOW KNOWN) TO FIND GZ AT THE NEXT MIDPOINT BETWEEN PRESSURE LEVELS. WHEN GZ=0 
! IS REACHED, SOLVE FOR THE PRESSURE.
!--------------------------------------------------------------------------------------------------
!
!----------------------------------------------------------- 
! COUNT THE POINTS WHERE SLP IS DONE BELOW EACH OUTPUT LEVEL
!-----------------------------------------------------------
    KOUNT = 0
    DO J=JSTA_I,JEND_I
        DO I=1,IM
              P1(I,J) = SPL(LMH(I,J))
            DONE(I,J) = .FALSE. 
            IF (FIS(I,J) < 10.) THEN
                PSLP(I,J) = PD(I,J) + PT
                DONE(I,J) = .TRUE. 
                KOUNT = KOUNT+1
            END IF
        END DO
    END DO
!
    KMM = KMNTM(LSM)
!
!$omp parallel do                                                                                 & 
!$omp private (GZ1      , GZ2     , I       , J       , LMAP1   , P1      , P2), SHARED (PSLP)
!
    DO 320 KM=1,KMM
        I = IMNT(KM,LSM)
        J = JMNT(KM,LSM)
!
        LMHIJ = LMH(I,J)
!
        GZ1   =  FI(I,J,LMHIJ)
!
        P1(I,J) = SPL(LMHIJ)
!    
        LMAP1 = LMHIJ + 1
        DO L=LMAP1,LSM
            P2 = SPL(L)
            TLYR = 0.5 * (T(I,J,L) + T(I,J,L-1))
            GZ2 = GZ1 + RD * TLYR * ALOG(P1(I,J) / P2)
            FI(I,J,L) = GZ2
            IF (GZ2 <= 0.) THEN
                PSLP(I,J) = P1(I,J) / EXP(-GZ1 / (RD * T(I,J,L-1)))
                DONE(I,J) = .TRUE. 
                KOUNT = KOUNT + 1
                GO TO 320
            END IF
            P1(I,J) = P2
            GZ1 = GZ2
        END DO
320 END DO
!--------------------------------------------------------------------------------------------------
! WHEN SEA LEVEL IS BELOW THE LOWEST OUTPUT PRESSURE LEVEL, SOLVE THE HYDROSTATIC EQUATION BY
! CHOOSING A TEMPERATURE AT THE MIDPOINT OF THE LAYER BETWEEN THAT LOWEST PRESSURE LEVEL AND THE
! GROUND BY EXTRAPOLATING DOWNWARD FROM T ON THE LOWEST PRESSURE LEVEL USING THE DT/DFI BETWEEN 
! THE LOWEST PRESSURE LEVEL AND THE ONE ABOVE IT.
!--------------------------------------------------------------------------------------------------
    TOTAL = (IM-2) * (JM-4)
!
    DO 340 LP=LSM,1,-1
        IF (KOUNT == TOTAL) GO TO 350
        DO 330 J=JSTA_I,JEND_I
            DO 330 I=1,IM
                IF (FI(I,J,LP) < 0. .OR. DONE(I,J)) GO TO 330
                SLOPE = (T(I,J,LP) - T(I,J,LP-1)) / (FI(I,J,LP) - FI(I,J,LP-1))
!
                TLYR = T(I,J,LP) - 0.5 * FI(I,J,LP) * SLOPE
                PSLP(I,J) = P1(I,J) / EXP(-FI(I,J,LP) / (RD * TLYR))
                DONE(I,J)= .TRUE. 
                KOUNT = KOUNT + 1
    330 END DO
340 END DO
!
350 CONTINUE
!
    PRINT*, 'DONE WITH MESINGER SLP REDUCTION'
!-------------------------- 
! SKIP THE STANDARD SCHEME.
!-------------------------- 
    GOTO 430
!-------------------------------------------------------------------------
! IF YOU WANT THE "STANDARD" ETA/SIGMA REDUCTION THIS IS WHERE IT IS DONE.
!-------------------------------------------------------------------------
400 CONTINUE
!
    DO 410 J=JSTA_I,JEND_I
        DO 410 I=1,IM
            IF (FIS(I,J) >= 1.) THEN
                  LMA = LM
                ALPP1 = ALOG(PD(I,J) + PT)
                 SLOP = 0.0065 * ROG * TSIG(I,J,LM)
!
                IF (SLOP < 0.50) THEN
                    SLPP = ALPP1 + FIS(I,J) / (RD * TSIG(I,J,LMA))
                ELSE
                     TTT = -(ALOG(PD(I,J) + PT) + ALPP1) * SLOP * 0.50 + TSIG(I,J,LMA)
                    SLPP = (-TTT + SQRT(TTT * TTT + 2. * SLOP * (FIS(I,J) / RD                    &
    &                    + ( TTT + 0.50 * SLOP * ALPP1) * ALPP1))) / SLOP
                END IF
                PSLP(I,J) = EXP(SLPP)
            END IF
410 END DO
!--------------------------------------------------------------------------------------------------
! AT THIS POINT WE HAVE A SEA LEVEL PRESSURE FIELD BY EITHER METHOD. 5-POINT AVERAGE THE FIELD ON 
! THE E-GRID.
!--------------------------------------------------------------------------------------------------
430 CONTINUE
!--------------------------------------- 
! EXTRAPOLATE VALUES TO THE OUTER 2 ROWS
!---------------------------------------
    PRINT*, 'ME IN SLPSIGC= ', ME
!
    IF (ME == 0) THEN
        IF ((JEND_I-JSTA_I) < 5) THEN
            PRINT*, 'HERE!'
            DO J=1,2
                IEND = IM - 1 - MOD(J+1,2)
                DO I=2,IEND
                    PSLP(I,J) = PSLP(I,3)
                END DO
            END DO
        ELSE
            DO J=1,2
                IEND = IM - 1 - MOD(J+1,2)
                DO I=2,IEND
                    PSLP(I,J) = 1.5 * PSLP(I,J+2) - 0.5 * PSLP(I,J+4)
                END DO
            END DO
        END IF
    END IF
!
    IF (ME == (NUM_PROCS-1)) THEN
        IF ((JEND_I-JSTA_I) < 5) THEN
            DO J=JM-1,JM
                IEND = IM - 1 - MOD(J+1,2)
                DO I=2,IEND
                    PSLP(I,J) = PSLP(I,IM-2)
                END DO
            END DO
        ELSE
            DO J=JM-1,JM
                IEND = IM - 1 - MOD(J+1,2)
                DO I=2,IEND
                    PSLP(I,J) = 1.5 * PSLP(I,J-2) - 0.5 * PSLP(I,J-4)
                END DO
            END DO
        END IF
    END IF
!
    DO J=JSTA_I,JEND_I
        PSLP(1,J) = 1.5 * PSLP(2,J)   - 0.5 * PSLP(3,J)
    END DO
!
    DO J=JSTA_I,JEND_I
        I = IM - MOD(J+1,2)
        PSLP(I,J) = 1.5 * PSLP(I-1,J) - 0.5 * PSLP(I-2,J)
    END DO
!
    PRINT*, 'PSLP BEFORE SMOOTHING'
!
!$omp parallel do 
!
    DO 440 J=JSTA_I,JEND_I
        DO 440 I=1,IM
            SLPX(I,J) = PSLP(I,J)
440 END DO
!
    DO 480 KS=1,KSLPD
!    
        CALL UPDATE(PSLP)
!
!$omp parallel do                                                                                 &
!$omp private (IHH2)
!
        DO 460 J=JSTA_IM2,JEND_IM2
            IHH2 = IM - 1 - MOD(J+1,2)
            DO 460 I=2,IHH2
!--------------------------------------------------------
! EXTRA AVERAGING UNDER MOUNTAINS TAKEN OUT, FM, MARCH 96
!--------------------------------------------------------            
                SLPX(I,J) = 0.125 * (PSLP(I+IHW(J),J-1) + PSLP(I+IHE(J),J-1)                      &
    &                     +          PSLP(I+IHW(J),J+1) + PSLP(I+IHE(J),J+1)                      &
    &                     + 4. *     PSLP(I       ,J  ))
    460 END DO
!
!$omp paralle do   
!  
        DO J=JSTA_I,JEND_I
            DO I=1,IM
                PSLP(I,J) = SLPX(I,J)
            END DO
        END DO
!
        PRINT*,'SAMPLE PSLP', IM/2, JSTA_I, JEND_I, PSLP(IM/2,JSTA_I), PSLP(IM/2,JEND_I)
!    
480 END DO
!
    WRITE(6,*) 'LEAVING SLPSIG'
!
    RETURN
!
    END SUBROUTINE SLPSIG
