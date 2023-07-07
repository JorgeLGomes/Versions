    SUBROUTINE GSMDRIVE
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE GSMCONST
!> 
!> SUBPROGRAM: GSMDRIVE - GRID-SCALE MICROPHYSICAL PROCESSES - CONDENSATION AND PRECIPITATION
!> PROGRAMMER: FERRIER
!> ORG: W/NP22     
!> DATE: 01-02-??
!>
!> ABSTRACT:
!> MERGES ORIGINAL GSCOND & PRECPD SUBROUTINES.
!> ODE HAS BEEN SUBSTANTIALLY STREAMLINED AND RESTRUCTURED.
!> EXCHANGE BETWEEN WATER VAPOR & SMALL CLOUD CONDENSATE IS CALCULATED USING THE ORIGINAL 
!> ASAI (1965, J. JAPAN) ALGORITHM.
!> SEE ALSO REFERENCES TO YAU AND AUSTIN (1979, JAS), RUTLEDGE AND HOBBS (1983, JAS), AND 
!> TAO ET AL. (1989, MWR).
!> THIS ALGORITHM REPLACES THE SUNDQVIST ET AL. (1989, MWR) PARAMETERIZATION.
!> DRIVER OF THE NEW MICROPHYSICS 
!>
!> NOTE:  CODE IS CURRENTLY SET UP W/O THREADING
!>
!> PROGRAM HISTORY LOG:
!> 01-04-XX    FERRIER  - BETA-TESTED VERSION
!> 01-05-21    FERRIER  - ADDED GRADUAL LATENT HEATING TO REMOVE EXTERNAL WAVES
!> 01-05-30    FERRIER  - CHANGED DEFAULT TO UNIFORM MARITIME CONDITIONS FOR TESTING
!> 18-01-15    LUCCI    - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> PRIOR PROGRAM HISTORY LOG:
!> HERITAGE AS SUBROUTINE GSCOND:
!> 94-??-??  ZHAO       - ORIGINATOR
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 95-03-28  BLACK      - ADDED EXTERNAL EDGE
!> 98-11-02  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
!>
!> HERITAGE AS SUBROUTINE PRECPD:
!> 94-??-??  ZHAO       - ORIGINATOR
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 95-11-20  ABELES     - PARALLEL OPTIMIZATION
!> 96-03-29  BLACK      - REMOVED SCRCH COMMON
!> 96-07-18  ZHAO       - NEW WMIN CALCULATION
!> 96-09-25  BALDWIN    - NEW SR CALCULATION
!> 98-11-02  BLACK      - MODIFICATION FOR DISTRIBUTED MEMORY
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
!> USE MODULES:  ACMCLH
!>               CLDWTR
!>               CMICRO_START
!>               CMICRO_STATS
!>               CTLBLK
!>               C_FRACN
!>               DYNAM
!>               F77KINDS
!>               GLB_TABLE
!>               LOOPS
!>               MAPPINGS
!>               MASKS
!>               MPPCOM
!>               PARMETA
!>               PARM_TBL
!>               PHYS
!>               PPTASM
!>               PVRBLS
!>               TEMPCOM
!>               TOPO
!>               VRBLS 
!>
!> DRIVER      : EBU
!>
!> CALLS       : GSMCOLUMN
!>               GSMCONST
!>               MPI_REDUCE
!>--------------------------------------------------------------------------------------------------
    USE ACMCLH
    USE CLDWTR
    USE CMICRO_START
    USE CMICRO_STATS
    USE CTLBLK
    USE C_FRACN
    USE DYNAM
    USE F77KINDS
    USE GLB_TABLE
    USE LOOPS
    USE MAPPINGS
    USE MASKS
    USE CMICRO_CONS
    USE MPPCOM
    USE PARMETA
    USE PARM_TBL
    USE PHYS
    USE PPTASM
    USE PVRBLS
    USE TEMPCOM
    USE TOPO
    USE VRBLS

    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!
#include "sp.h"
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & NOZ
!----------------------------------------------------
! THIS VARIABLE IS FOR DEBUGGING PURPOSES (IF .TRUE.)
!---------------------------------------------------- 
    LOGICAL(KIND=L4KIND), PARAMETER :: PRINT_DIAG = .TRUE. 
!-------------------------------------------------------------------------
! THE FOLLOWING VARIABLES ARE FOR MICROPHYSICAL STATISTICS (NON-ESSENTIAL)
!-------------------------------------------------------------------------
    INTEGER(KIND=I4KIND), PARAMETER :: ITHILO    = ITHI   - ITLO + 1
    INTEGER(KIND=I4KIND), PARAMETER :: ITHILO_N  = ITHILO *  4
    INTEGER(KIND=I4KIND), PARAMETER :: ITHILO_QM = ITHILO *  5
    INTEGER(KIND=I4KIND), PARAMETER :: ITHILO_QT = ITHILO * 22
!
    INTEGER(KIND=I4KIND), DIMENSION(ITLO:ITHI,  4)                                              ::&
    & NSTATS_0
!
    REAL   (KIND=R4KIND), DIMENSION(ITLO:ITHI,  5)                                              ::&
    & QMAX_0
!
    REAL   (KIND=R4KIND), DIMENSION(ITLO:ITHI, 22)                                              ::&
    & QTOT_0
!
    REAL   (KIND=R4KIND), SAVE                                                                  ::&
    & THOUR_PRINT
!
    REAL   (KIND=R4KIND), DIMENSION(2), SAVE                                                    ::&
    & PRECMAX , PRECTOT , PRECMAX_0, PRECTOT_0
!
    REAL   (KIND=R4KIND), PARAMETER :: DTHOUR_PRINT = 3. ! PRINT STATISTICS EVERY 3 H
!--------------------------------------- 
! BEGIN SECTION ON HYDROMETEOR FRACTIONS
!--------------------------------------- 
!
!------------------------------------------------------------ 
! SAVED VALUES USE REAL (REAL*4) ARRAYS RATHER THAN INTEGER*2
!------------------------------------------------------------ 
    REAL   (KIND=R4KIND)                                                                        ::&
    & FICE    , FRAIN   , DUM
!----------------------------------------------------- 
! SAVED VALUES USE INTEGER*2 ARRAYS RATHER THAN REAL*4
!----------------------------------------------------- 
!
!--------------------------------------------------------------------------------------------------
! PARAMETERS USED IN CONVERTING FROM REAL TO INTEGER*2 FIELDS FOR FRACTION OF ICE (F_ICE), FRACTION 
! OF RAIN (F_RAIN), AND MASS RATIO OF RIMED ICE ("RIME FACTOR", F_RIMEF).
!--------------------------------------------------------------------------------------------------
!
!------------------------------------- 
! END SECTION ON HYDROMETEOR FRACTIONS
!------------------------------------- 
!
!---------------------------------------- 
! LOCAL ARRAYS AND PARAMETERS IN GSMDRIVE 
!---------------------------------------- 
!
!--------------------------------------------------------------------------------------------------
! COMMENTS ON 23 AUGUST 2001
! EPSQ=1.E-20 IS THE LOWER LIMIT FOR SPECIFIC HUMIDITY AND CLOUD CONDENSATE.
! THE VALUE OF EPSQ WILL NEED TO BE CHANGED IN THE OTHER SUBROUTINES IN ORDER TO MAKE IT CONSISTENT
! THROUGHOUT THE ETA CODE.
!--------------------------------------------------------------------------------------------------
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ   =    1.E-20
    REAL   (KIND=R4KIND), PARAMETER :: GRAV   =    9.81
    REAL   (KIND=R4KIND), PARAMETER :: RHOL   = 1000.
    REAL   (KIND=R4KIND), PARAMETER :: T0C    =  273.15
    REAL   (KIND=R4KIND), PARAMETER :: T_ICE  =  -10.
    REAL   (KIND=R4KIND), PARAMETER :: T_ICEK = T0C + T_ICE
    REAL   (KIND=R4KIND), PARAMETER :: RRHOL  = 1.  / RHOL
    REAL   (KIND=R4KIND), PARAMETER :: CLIMIT =    1.E-12
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ARAIN   , ASNOW
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                                         ::&
    & P_COL   , QI_COL  , QR_COL  , QV_COL    , QW_COL  , RIMEF_COL, T_COL   , THICK_COL ,        &
    & WC_COL
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    REAL   (KIND=R4KIND)                                                                        ::&
    & DTPH    , PDSL    , CONSTA  , TIME_MODEL, HDTPH   , DEL_HYD  , DEL_QT  , RIMEF_BULK,        &
    & TC      , WC      , QI      , QR        , QW
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K       , LSFC      , I_INDEX , J_INDEX  , IRTN    , II        ,        &
    & JJ      , KK
!---------------- 
! BEGIN EXECUTION 
!---------------- 
!
!------------------------  
! MICROPHYSICAL CONSTANTS 
!------------------------ 
    DTPH = DTQ2 ! PHYSICS TIME STEP (S)
!------------------------------------ 
! INITIALIZE CONSTANTS FOR STATISTICS 
!------------------------------------ 
    IF (MICRO_START) THEN
        MICRO_START = .FALSE. ! NO NEED TO CALCULATE THESE PARAMETERS AGAIN
!------------------------------------------------------------------------     
! BEGIN: INITIALIZE ARRAYS IF NOT INITIALIZED IN RESTART FILES (01/16/02)
!------------------------------------------------------------------------     
        DUM = 0.
        DO J=MYJS2,MYJE2
            DO I=MYIS,MYIE
                IF (HBM2(I,J) > 0.5) THEN
                    LSFC = LMH(I,J) ! "K" OF SURFACE
                    DO K=1,LSFC
                        DUM = MAX(DUM,F_RAIN(K,I,J), F_ICE(K,I,J), F_RIMEF(K,I,J))
                        IF (DUM > 0.) GO TO 110
                    END DO
                END IF
            END DO
        END DO
!
    110 IF (DUM <= 0.) THEN
            IF (MYPE == 0) WRITE(0, "(A,2I3)" ) 'INITIALIZE INTERNAL MICROPHYSICS ARRAYS'
            DO J=MYJS2,MYJE2
                DO I=MYIS,MYIE
                    IF (HBM2(I,J) > 0.5) THEN
                        LSFC = LMH(I,J) ! "K" OF SURFACE
                        DO K=1,LSFC
                               F_RAIN(K,I,J) = 0.
!
                            IF (T(I,J,K) <= T_ICEK) THEN
                                F_ICE(K,I,J) = 1.
                            ELSE
                                F_ICE(K,I,J) = 0.
                            END IF
!
                              F_RIMEF(K,I,J) = 1.
                        END DO
                    END IF
                END DO
            END DO
        END IF
!    
        CALL GSMCONST(DTPH) ! INITIALIZE LOOKUP TABLES & CONSTANTS
!
        THOUR_PRINT = -DTPH / 3600.
!
        IF (PRINT_DIAG) THEN
!-----------------------------         
! TOTAL AND MAXIMUM QUANTITIES
!-----------------------------          
            NSTATS  = 0     ! MICROPHYSICAL STATISTICS DEALING W/ GRID-POINT COUNTS
            QMAX    = 0.    ! MICROPHYSICAL STATISTICS DEALING W/ HYDROMETEOR MASS
            QTOT    = 0.    ! MICROPHYSICAL STATISTICS DEALING W/ HYDROMETEOR MASS
            PRECMAX = 0.    ! MAXIMUM PRECIP RATES (RAIN, SNOW) AT SURFACE (MM/H)
            PRECTOT = 0.    ! TOTAL PRECIPITATION (RAIN, SNOW) ACCUMULATION AT SURFACE
        END IF
    END IF
!---------------------------------  
! LOOP HORIZONTALLY THROUGH DOMAIN 
!--------------------------------- 
    DO 100 J=MYJS2,MYJE2
        DO 100 I=MYIS,MYIE
            IF (HBM2(I,J) < 0.5) GOTO 100    ! IGNORE COLUMNS NEAR LATERAL BOUNDARIES
            LSFC   = LMH(I,J)                ! "K" OF SURFACE
            PDSL   =  PD(I,J) * RES(I,J)     ! (PSFC-PTOP)/ETA_SFC
            CONSTA = PDSL / GRAV             ! (PSFC-PTOP)/(G*ETA_SFC), USED FOR THICK BELOW
!-----------------------------------         
! INITIALIZE COLUMN DATA (1D ARRAYS)
!-----------------------------------         
            IF (CWM(I,J,1) <= CLIMIT) CWM(I,J,1) = EPSQ
              F_ICE(1,I,J) = 1.
             F_RAIN(1,I,J) = 0.
            F_RIMEF(1,I,J) = 1.
!--------------------------------  
! INICIALIZACAO DA VARIAVEL T_COL
!--------------------------------  
            T_COL(1:LM) = 0.
!
            DO K=1,LSFC
!-----------------------------------------------             
! PRESSURE (PA) = (PSFC-PTOP)*(ETA/ETA_SFC)+PTOP
!-----------------------------------------------            
                P_COL(K) = PDSL * AETA(K) + PT
!-----------------------------------------------------------------            
! LAYER THICKNESS = RHO*DZ = -DP/G = (PSFC-PTOP)*D_ETA/(G*ETA_SFC)
!----------------------------------------------------------------- 
                THICK_COL(K) = CONSTA * DETA(K)
                    T_COL(K) = T(I,J,K)
                          TC = T_COL(K) - T0C
                   QV_COL(K) = MAX(EPSQ, Q(I,J,K))
!
                IF (CWM(I,J,K) <= CLIMIT) THEN
                    WC_COL(K) = 0.
!
                    IF (TC < T_ICE) THEN
                        F_ICE(K,I,J) = 1.
                    ELSE
                        F_ICE(K,I,J) = 0.
                    END IF
!
                       F_RAIN(K,I,J) = 0.
                      F_RIMEF(K,I,J) = 1.
                ELSE
                    WC_COL(K) = CWM(I,J,K)
                END IF
!---------------------------------------------------------------------------             
! DETERMINE COMPOSITION OF CONDENSATE IN TERMS OF CLOUD WATER, ICE, AND RAIN
 !---------------------------------------------------------------------------           
                WC = WC_COL(K)
                QI = 0.
                QR = 0.
                QW = 0.
!
                FICE =   F_ICE(K,I,J)
                FRAIN = F_RAIN(K,I,J)
!---------------------             
! REAL*4 ARRAY STORAGE
!---------------------             
                IF (FICE >= 1.) THEN
                    QI = WC
                ELSE IF (FICE <= 0.) THEN
                    QW = WC
                ELSE
                    QI = FICE * WC
                    QW = WC   - QI
                END IF
!
                IF (QW > 0. .AND. FRAIN > 0.) THEN
                    IF (FRAIN >= 1.) THEN
                        QR = QW
                        QW = 0.
                    ELSE
                        QR = FRAIN * QW
                        QW = QW    - QR
                    END IF
                END IF
!
                RIMEF_COL(K) = F_RIMEF(K,I,J) ! (REAL)
!
                QI_COL(K) = QI
                QR_COL(K) = QR
                QW_COL(K) = QW
            END DO
!------------------------------------------------------ 
! PERFORM THE MICROPHYSICAL CALCULATIONS IN THIS COLUMN
!------------------------------------------------------        
            I_INDEX = I
            J_INDEX = J
!-------------------------------------------------------- 
! DIEGO - NEW RHGRD CALCULATION USING SM(I,J) (23/10/2020
!--------------------------------------------------------
! when over ocean SM=1 and RHGRD=0.99
! when over land  SM=0 and RHGRD=0.95
!
            RHGRD = RHGRDS * SM(I,J) + RHGRDL * (1.-SM(I,J))
!	    	    
            CALL GSMCOLUMN(ARAIN , ASNOW , DTPH  , I_INDEX  , J_INDEX, LSFC     , P_COL , QI_COL, &
    &                      QR_COL, QV_COL, QW_COL, RIMEF_COL, T_COL  , THICK_COL, WC_COL)
!----------------------  
! UPDATE STORAGE ARRAYS
!----------------------        
            DO K=1,LSFC
                 TRAIN(I,J,K) = (T_COL(K) - T(I,J,K)) / DTPH
                     T(I,J,K) =  T_COL(K)
                     Q(I,J,K) = QV_COL(K)
                TLATGS(I,J,K) =  T_COL(K) - T(I,J,K)
                   CWM(I,J,K) = WC_COL(K)
!---------------------             
! REAL*4 ARRAY STORAGE
!---------------------            
                F_RIMEF(K,I,J) = MAX(1., RIMEF_COL(K))
!
                IF (QI_COL(K) <= CLIMIT) THEN
                    F_ICE(K,I,J) = 0.
                    IF (T_COL(K) < T_ICEK) F_ICE(K,I,J) = 1.
                ELSE
                    F_ICE(K,I,J) = MAX(0., MIN(1., QI_COL(K) / WC_COL(K)))
                END IF
!
                IF (QR_COL(K) <= CLIMIT) THEN
                    DUM = 0
                ELSE
                    DUM = QR_COL(K) / (QR_COL(K) + QW_COL(K))
                END IF
!
                   F_RAIN(K,I,J) = DUM
!
            END DO
!--------------------------------------------------------------------------------------------------        
! UPDATE ACCUMULATED PRECIPITATION STATISTICS
!       
! SURFACE PRECIPITATION STATISTICS; 
! SR IS FRACTION OF SURFACE PRECIPITATION (IF >0) ASSOCIATED WITH SNOW
!--------------------------------------------------------------------------------------------------       
             APREC(I,J) = (ARAIN + ASNOW) * RRHOL  ! ACCUMULATED SURFACE PRECIP (DEPTH IN M) - YING
!
              PREC(I,J) =   PREC(I,J) + APREC(I,J)
            ACPREC(I,J) = ACPREC(I,J) + APREC(I,J)
!
            IF (APREC(I,J) < 1.E-8) THEN
                SR(I,J) = 0.
            ELSE
                SR(I,J) = RRHOL * ASNOW / APREC(I,J)
            END IF
!-----------------         
! DEBUG STATISTICS
!-----------------         
            IF (PRINT_DIAG) THEN
                PRECTOT(1) = PRECTOT(1) + ARAIN
                PRECTOT(2) = PRECTOT(2) + ASNOW
!
                PRECMAX(1) = MAX(PRECMAX(1), ARAIN)
                PRECMAX(2) = MAX(PRECMAX(2), ASNOW)
            END IF
!       
            CONTINUE ! END "I" & "J" LOOPS
100 END DO
!------------------------------ 
! END OF MAIN MICROPHYSICS LOOP 
!------------------------------ 
    TIME_MODEL = FLOAT(NTSD-1) * DT / 3600.
!
    IF (PRINT_DIAG .AND. TIME_MODEL >= THOUR_PRINT) THEN
        CALL MPI_REDUCE(NSTATS , NSTATS_0 , ITHILO_N , MPI_INTEGER, MPI_SUM, 0, MPI_COMM_COMP,IRTN)
        CALL MPI_REDUCE(QMAX   , QMAX_0   , ITHILO_QM, MPI_REAL   , MPI_MAX, 0, MPI_COMM_COMP,IRTN)
        CALL MPI_REDUCE(PRECMAX, PRECMAX_0, 2        , MPI_REAL   , MPI_MAX, 0, MPI_COMM_COMP,IRTN)
        CALL MPI_REDUCE(QTOT   , QTOT_0   , ITHILO_QT, MPI_REAL   , MPI_SUM, 0, MPI_COMM_COMP,IRTN)
        CALL MPI_REDUCE(PRECTOT, PRECTOT_0, 2        , MPI_REAL   , MPI_SUM, 0, MPI_COMM_COMP,IRTN)
!
        IF (MYPE == 0) THEN
            HDTPH = 3600. / DTPH ! CONVERT PRECIP RATES TO MM/H
!
            DO K=ITLO,ITHI
                QMAX_0(K,1) = 1000. * QMAX_0(K,1)
                QMAX_0(K,2) = 1000. * QMAX_0(K,2)
                QMAX_0(K,3) = 1000. * QMAX_0(K,3)
                QMAX_0(K,4) = HDTPH * QMAX_0(K,4)
                QMAX_0(K,5) = HDTPH * QMAX_0(K,5)
            END DO
!
            PRECMAX_0(1) = HDTPH * PRECMAX_0(1)
            PRECMAX_0(2) = HDTPH * PRECMAX_0(2)
        
!GSM            WRITE(6,"(A,F5.2,4(A,G11.4))") '{ TIME(H)=',TIME_MODEL,'  TRAIN_SFC=',PRECTOT_0(1),   &
!GSM    &       '  TSNOW_SFC=',PRECTOT_0(2),'  RRMAX_SFC(MM/H)=',PRECMAX_0(1),                        &
!GSM    &       '  SRMAX_SFC(MM/H)=',PRECMAX_0(2)
!  
!GSM            WRITE(6,"(3A)") '{ (C) <--------- COUNTS ----------> ',                               &
!GSM    &       '<----------- G/KG ----------> <----- MM/H ------>'   ,                               &
!GSM            ' <---- KG/M**2 * # GRIDS ---->'
!
!GSM            WRITE(6,"(3A)") '{  T     NCICE  NCMIX  NCWAT NCRAIN  '       ,                       &
!GSM    &       'QIMAX     QWMAX     QRMAX     SRMAX     RRMAX     QITOT     ',                       &
!GSM    &       'QWTOT     QRTOT'
!
!GSM            DO K=ITLO,ITHI
!GSM                WRITE(6,"(A,I3,I9,3I7,8G10.4)") '{ ', K, (NSTATS_0(K,II), II=1,4),                &
!GSM    &                                                      (QMAX_0(K,JJ), JJ=1,5),                &
!GSM    &                                                      (QTOT_0(K,KK), KK=1,3)
!GSM            END DO
 !       
!GSM            WRITE(6,"(3A)") '{  T   TCOND     TICND     TIEVP     TIDEP     TREVP     ',          &
!GSM    &       'TRAUT     TRACW     TIMLT     TIACW     TIACWI    TIACWR    ','TIACR'
!
!GSM            DO K=ITLO,ITHI
!GSM                WRITE(6,"(A,I3,12G10.4)") '{ ',K,(QTOT_0(K,II), II=4,15)
!GSM            END DO
!        
!GSM            WRITE(6,"(2A)") '{  T   DEL_QT   TVDIF   DEL_HYD        TWDIF  TIDIF       ',         &
!GSM    &       'TRDIF    DARAIN   DASNOW    RIMEF'
!
            DO K=ITLO,ITHI
                DEL_HYD = 0.
!
                DO II=17,19
                    DEL_HYD = DEL_HYD + QTOT_0(K,II)
                END DO
!
                DEL_QT = 0.
!
                DO II=16,21
                    DEL_QT = DEL_QT + QTOT_0(K,II)
                END DO
!
                IF (QTOT_0(K,22) > 0.) THEN
                    RIMEF_BULK = QTOT_0(K,1) / QTOT_0(K,22)
                ELSE
                    RIMEF_BULK = 1.
                END IF
!
                WRITE(6,"(A,I3,9G10.4)") '{ ', K, DEL_QT, QTOT_0(K,16),            DEL_HYD   ,    &
    &                                                    (QTOT_0(K,II), II=17,21), RIMEF_BULK
            END DO
!       
        END IF
!
        NSTATS  = 0     ! MICROPHYSICAL STATISTICS DEALING W/ GRID-POINT COUNTS
        QMAX    = 0.    ! MICROPHYSICAL STATISTICS DEALING W/ HYDROMETEOR MASS
        QTOT    = 0.    ! MICROPHYSICAL STATISTICS DEALING W/ HYDROMETEOR MASS
        PRECMAX = 0.    ! MAXIMUM PRECIP RATES (RAIN, SNOW) AT SURFACE (MM/H)
        PRECTOT = 0.    ! TOTAL PRECIPITATION (RAIN, SNOW) ACCUMULATION AT SURFACE
!
        THOUR_PRINT = THOUR_PRINT + DTHOUR_PRINT
    END IF
!----------------------- 
! RETURN TO MAIN PROGRAM 
!----------------------- 
    RETURN
!
200 FORMAT(A2,I5,F6.2,4(1X,A10,G11.4))
210 FORMAT(A2,I5,F6.2,4(1X,A10,I7))
!
    END SUBROUTINE GSMDRIVE

