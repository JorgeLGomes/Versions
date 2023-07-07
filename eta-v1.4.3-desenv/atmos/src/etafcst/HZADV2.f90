    SUBROUTINE HZADV2
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE HZADV2
!>
!> SUBPROGRAM: HZADV2 - HORIZONTAL ADVECTION OF VAPOR AND CLOUD
!> PROGRAMMER: JANJIC
!> ORG: W/NP22
!> DATE: 96-07-19
!>
!> ABSTRACT:
!> HZADV2 CALCULATES THE CONTRIBUTION OF THE HORIZONTAL ADVECTION TO THE TENDENCIES OF SPECIFIC 
!> HUMIDITY AND CLOUD WATER AND THEN UPDATES THOSE VARIABLES.  AN ANTI-FILTERING TECHNIQUE IS USED.
!>
!> PROGRAM HISTORY LOG:
!> 96-07-19  JANJIC   - ORIGINATOR
!> 98-11-02  BLACK    - MODIFIED FOR DISTRIBUTED MEMORY
!> 99-03-17  TUCCILLO - INCORPORATED MPI_ALLREDUCE FOR GLOBAL SUM
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
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: CLDWTR
!>              CONTIN
!>              CTLBLK
!>              DYNAM
!>              F77KINDS
!>              GLB_TABLE
!>              INDX
!>              LOOPS
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              PARMETA
!>              PVRBLS
!>              VRBLS
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : EBU
!>              NEWFLT
!>
!> CALLS      : MPI_ALLREDUCE                
!>--------------------------------------------------------------------------------------------------
    USE CLDWTR
    USE CONTIN
    USE CTLBLK
    USE DYNAM
    USE F77KINDS
    USE GLB_TABLE
    USE INDX
    USE LOOPS
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE PARMETA
    USE PVRBLS
    USE VRBLS
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ   =  2.E-12
    REAL   (KIND=R4KIND), PARAMETER :: CLIMIT =  1.E-20
    REAL   (KIND=R4KIND), PARAMETER :: FF1    =  0.52500
    REAL   (KIND=R4KIND), PARAMETER :: FF2    = -0.64813
    REAL   (KIND=R4KIND), PARAMETER :: FF3    =  0.24520
    REAL   (KIND=R4KIND), PARAMETER :: FF4    = -0.12189
!
    INCLUDE "mpif.h"
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: IM1    = IM - 1
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM   = IM * JM - JM / 2
    INTEGER(KIND=I4KIND), PARAMETER :: JAMD   = (JAM * 2 - 10) * 3
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & IFPA    , IFQA    ,                                                                         &
    & IFPF    , IFQF    ,                                                                         &
    & JFPA    , JFQA    ,                                                                         &
    & JFPF    , JFQF
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & AFP     , AFQ     ,                                                                         &
    & Q1      , DQST    ,                                                                         &
    & W1      , DWST    ,                                                                         &
    & DVOL
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & DARE    , EMH 
!
    REAL   (KIND=R4KIND), DIMENSION(4,LM)                                                       ::&
    & GSUMS   , XSUMS
!
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K       , JFP     , JFQ     , IRECV
!
    REAL(KIND=R4KIND)                                                                           ::&
    & HTMIJL  , SUMPQ   , SUMNQ   , SUMPW   , SUMNW   , DVOLP   , TTA     , TTB     , PP      ,   &
    & QP      , DQSTIJ  , DWSTIJ  , RFACQ   , RFACW   , RFQIJ   , RFWIJ   , Q1IJ    , W1IJ    ,   &
    & ENH     , D2PQQ   , D2PQW   , QSTIJ   , WSTIJ   , Q00     , QP0     , Q0Q     , W00     ,   &
    & WP0     , W0Q
!
    ENH = FLOAT(IDTAD) * DT / (08. * DY)
!
    DO J=MYJS_P2,MYJE_P2
        DO I=MYIS_P1,MYIE_P1
            EMH (I,J) = FLOAT(IDTAD) * DT / (08. * DX(I,J))
            DARE(I,J) = HBM2(I,J) * DX(I,J) * DY
        END DO
    END DO
!
!$omp parallel do
!$omp private (DQSTIJ   , DVOLP   , DWSTIJ  , HTMIJL  , JFP     , JFQ     , PP      , QP      ,   &
!$omp          SUMNQ    , SUMNW   , SUMPQ   , SUMPW   , TTA     , TTB     )
!
    DO K=1,LM
!
        DO 200 J=MYJS_P2,MYJE_P2
            DO 200 I=MYIS_P1,MYIE_P1
                DVOL(I,J,K) = DARE(I,J) * PDSL(I,J) * DETA(K)
                HTMIJL = HTM(I,J,K)
!
                 Q  (I,J,K) = AMAX1(Q  (I,J,K), EPSQ  ) * HTMIJL
                 CWM(I,J,K) = AMAX1(CWM(I,J,K), CLIMIT) * HTMIJL
                Q1  (I,J,K) = Q  (I,J,K)
                W1  (I,J,K) = CWM(I,J,K)
    200 END DO
!
        SUMPQ = 0.
        SUMNQ = 0.
        SUMPW = 0.
        SUMNW = 0.
!    
        DO 300 J=MYJS2_P1,MYJE2_P1
            DO 300 I=MYIS1_P1,MYIE1_P1
!            
                DVOLP = DVOL(I,J,K) * HBM3(I,J)
                TTA = (U(I       ,J-1,K) + U(I+IHW(J),J  ,K)                                      &
    &               +  U(I+IHE(J),J  ,K) + U(I       ,J+1,K)) * HBM2(I,J) * EMH(I,J)
!
                TTB = (V(I       ,J-1,K) + V(I+IHW(J),J  ,K)                                      &
    &               +  V(I+IHE(J),J  ,K) + V(I       ,J+1,K)) * HBM2(I,J) * ENH
!            
                PP = - TTA - TTB
                QP =   TTA - TTB
!            
                JFP = INT(SIGN(1.,PP))
                JFQ = INT(SIGN(1.,QP))
!            
                IFPA(I,J,K) = IHE(J) + I + ( JFP-1) / 2
                IFQA(I,J,K) = IHE(J) + I + (-JFQ-1) / 2
!            
                JFPA(I,J,K) = J + JFP
                JFQA(I,J,K) = J + JFQ
!            
                IFPF(I,J,K) = IHE(J) + I + (-JFP-1) / 2
                IFQF(I,J,K) = IHE(J) + I + ( JFQ-1) / 2
!            
                JFPF(I,J,K) = J - JFP
                JFQF(I,J,K) = J - JFQ
!            
                PP = ABS(PP) * HTM(I,J,K) * HTM(IFPA(I,J,K),JFPA(I,J,K),K)
                QP = ABS(QP) * HTM(I,J,K) * HTM(IFQA(I,J,K),JFQA(I,J,K),K)
!            
                AFP (I,J,K) = (((FF4*PP+FF3) * PP + FF2) * PP + FF1) * PP
                AFQ (I,J,K) = (((FF4*QP+FF3) * QP + FF2) * QP + FF1) * QP
!            
                DQSTIJ = (Q  (IFPA(I,J,K),JFPA(I,J,K),K) - Q  (I,J,K)) * PP                       &
    &                  + (Q  (IFQA(I,J,K),JFQA(I,J,K),K) - Q  (I,J,K)) * QP
!
                DWSTIJ = (CWM(IFPA(I,J,K),JFPA(I,J,K),K) - CWM(I,J,K)) * PP                       &
    &                  + (CWM(IFQA(I,J,K),JFQA(I,J,K),K) - CWM(I,J,K)) * QP
!            
                DQST(I,J,K) = DQSTIJ
                DWST(I,J,K) = DWSTIJ
!            
    300 END DO
!---------------------------- 
! GLOBAL SUM FOR CONSERVATION
!----------------------------
        DO 310 J=MYJS2,MYJE2
            DO 310 I=MYIS1,MYIE1
!            
                DVOLP =  DVOL(I,J,K) * HBM3(I,J)
                DQSTIJ = DQST(I,J,K) * DVOLP
                DWSTIJ = DWST(I,J,K) * DVOLP
!            
                IF (DQSTIJ > 0.) THEN
                    SUMPQ = SUMPQ + DQSTIJ
                ELSE
                    SUMNQ = SUMNQ + DQSTIJ
                END IF
!            
                IF (DWSTIJ > 0.) THEN
                    SUMPW = SUMPW + DWSTIJ
                ELSE
                    SUMNW = SUMNW + DWSTIJ
                END IF
!            
    310 END DO
!
        XSUMS(1,K) = SUMPQ
        XSUMS(2,K) = SUMNQ
        XSUMS(3,K) = SUMPW
        XSUMS(4,K) = SUMNW
!--------------- 
! END OF LM LOOP
!--------------- 
    END DO               
!----------------- 
! GLOBAL REDUCTION
!-----------------
    CALL MPI_ALLREDUCE(XSUMS, GSUMS, 4*LM, MPI_REAL, MPI_SUM, MPI_COMM_COMP, IRECV)
!------------------------
! END OF GLOBAL REDUCTION
!------------------------
!
!$omp parallel do
!$omp private (D2PQQ    , D2PQW   , DQSTIJ  , DVOLP   , DWSTIJ  , Q00     , Q0Q     , Q1IJ    ,   &
!$omp          QP0      , QSTIJ   , RFACQ   , RFACW   , RFQIJ   , RFWIJ   , SUMNQ   , SUMNW   ,   &
!$omp          SUMPQ    , SUMPW   , W00     , W0Q     , W1IJ    , WP0     , WSTIJ   )
!
    DO K=1,LM
!
        SUMPQ = GSUMS(1,K)
        SUMNQ = GSUMS(2,K)
        SUMPW = GSUMS(3,K)
        SUMNW = GSUMS(4,K)
!-------------------------------    
! FIRST MOMENT CONSERVING FACTOR 
!------------------------------- 
        IF (SUMPQ > 1.) THEN
            RFACQ = -SUMNQ / SUMPQ
        ELSE
            RFACQ = 1.
        END IF
!    
        IF (SUMPW > 1.) THEN
            RFACW = -SUMNW / SUMPW
        ELSE
            RFACW = 1.
        END IF
!    
        IF (RFACQ < 0.9 .OR. RFACQ > 1.1) RFACQ = 1.
        IF (RFACW < 0.9 .OR. RFACW > 1.1 )RFACW = 1.
!--------------------------------- 
! IMPOSE CONSERVATION ON ADVECTION
!---------------------------------
        IF (RFACQ < 1.) THEN
            DO J=MYJS2_P1,MYJE2_P1
                DO I=MYIS1_P1,MYIE1_P1
                    DQSTIJ = DQST(I,J,K)
                    RFQIJ   = HBM3(I,J) * (RFACQ-1.) + 1.
                    IF (DQSTIJ < 0.) DQSTIJ = DQSTIJ / RFQIJ
                    Q1(I,J,K) = Q(I,J,K) + DQSTIJ
                END DO
            END DO
        ELSE
            DO J=MYJS2_P1,MYJE2_P1
                DO I=MYIS1_P1,MYIE1_P1
                    DQSTIJ = DQST(I,J,K)
                    RFQIJ = HBM3(I,J) * (RFACQ-1.) + 1.
                    IF (DQSTIJ >= 0.) DQSTIJ = DQSTIJ * RFQIJ
                    Q1(I,J,K) = Q(I,J,K) + DQSTIJ
                END DO
            END DO
        END IF
!
        IF (RFACW < 1.) THEN
            DO J=MYJS2_P1,MYJE2_P1
                DO I=MYIS1_P1,MYIE1_P1
                    DWSTIJ = DWST(I,J,K)
                    RFWIJ  = HBM3(I,J) * (RFACW-1.) + 1.
                    IF (DWSTIJ < 0.) DWSTIJ = DWSTIJ / RFWIJ
                    W1(I,J,K) = CWM(I,J,K) + DWSTIJ
                END DO
            END DO
        ELSE
            DO J=MYJS2_P1,MYJE2_P1
                DO I=MYIS1_P1,MYIE1_P1
                    DWSTIJ = DWST(I,J,K)
                    RFWIJ  = HBM3(I,J) * (RFACW-1.) + 1.
                    IF (DWSTIJ >= 0.) DWSTIJ = DWSTIJ * RFWIJ
                    W1(I,J,K) = CWM(I,J,K) + DWSTIJ
                END DO
            END DO
        END IF
!-------------------- 
! ANTI-FILTERING STEP
!-------------------- 
        SUMPQ = 0.
        SUMNQ = 0.
        SUMPW = 0.
        SUMNW = 0.
!------------------------ 
! ANTI-FILTERING LIMITERS 
!------------------------
        DO 330 J=MYJS2,MYJE2
            DO 330 I=MYIS1,MYIE1
!            
                DVOLP = DVOL(I,J,K)
                Q1IJ  =   Q1(I,J,K)
                W1IJ  =   W1(I,J,K)
!            
                D2PQQ = (        (Q1(IFPA(I,J,K),JFPA(I,J,K),K)-Q1IJ)                             &
    &                 -  (Q1IJ -  Q1(IFPF(I,J,K),JFPF(I,J,K),K))                                  &
    &                 *          HTM(IFPF(I,J,K),JFPF(I,J,K),K)) * AFP(I,J,K)                     &
    &                 +         ((Q1(IFQA(I,J,K),JFQA(I,J,K),K)-Q1IJ)                             &
    &                 -  (Q1IJ -  Q1(IFQF(I,J,K),JFQF(I,J,K),K))                                  &
    &                 *          HTM(IFQF(I,J,K),JFQF(I,J,K),K)) * AFQ(I,J,K)
!            
                D2PQW = (    (W1(IFPA(I,J,K),JFPA(I,J,K),K)-W1IJ)                                 &
    &                 - (W1IJ-W1(IFPF(I,J,K),JFPF(I,J,K),K))                                      &
    &                 *      HTM(IFPF(I,J,K),JFPF(I,J,K),K)) * AFP(I,J,K)                         &
    &                 +     ((W1(IFQA(I,J,K),JFQA(I,J,K),K)-W1IJ)                                 &
    &                 - (W1IJ-W1(IFQF(I,J,K),JFQF(I,J,K),K))                                      &
    &                 *      HTM(IFQF(I,J,K),JFQF(I,J,K),K)) * AFQ(I,J,K)
!            
                QSTIJ = Q1IJ - D2PQQ
                WSTIJ = W1IJ - D2PQW
!            
                Q00 = Q(I          ,J          ,K)
                QP0 = Q(IFPA(I,J,K),JFPA(I,J,K),K)
                Q0Q = Q(IFQA(I,J,K),JFQA(I,J,K),K)
!            
                W00 = CWM(I          ,J          ,K)
                WP0 = CWM(IFPA(I,J,K),JFPA(I,J,K),K)
                W0Q = CWM(IFQA(I,J,K),JFQA(I,J,K),K)
!            
                QSTIJ = AMAX1(QSTIJ, AMIN1(Q00,QP0,Q0Q))
                QSTIJ = AMIN1(QSTIJ, AMAX1(Q00,QP0,Q0Q))
                WSTIJ = AMAX1(WSTIJ, AMIN1(W00,WP0,W0Q))
                WSTIJ = AMIN1(WSTIJ, AMAX1(W00,WP0,W0Q))
!            
                DQSTIJ = QSTIJ - Q1IJ
                DWSTIJ = WSTIJ - W1IJ
!            
                DQST(I,J,K) = DQSTIJ
                DWST(I,J,K) = DWSTIJ
!            
                DQSTIJ = DQSTIJ * DVOLP
                DWSTIJ = DWSTIJ * DVOLP
!            
                IF (DQSTIJ > 0.) THEN
                    SUMPQ = SUMPQ + DQSTIJ
                ELSE
                    SUMNQ = SUMNQ + DQSTIJ
                END IF
!            
                IF (DWSTIJ > 0.) THEN
                    SUMPW = SUMPW + DWSTIJ
                ELSE
                    SUMNW = SUMNW + DWSTIJ
                END IF
!            
    330 END DO
!
        XSUMS(1,K) = SUMPQ
        XSUMS(2,K) = SUMNQ
        XSUMS(3,K) = SUMPW
        XSUMS(4,K) = SUMNW
!---------------  
! END OF LM LOOP
!---------------
    END DO               
!----------------- 
! GLOBAL REDUCTION
!----------------- 
    CALL MPI_ALLREDUCE(XSUMS, GSUMS, 4*LM, MPI_REAL, MPI_SUM, MPI_COMM_COMP, IRECV)
!------------------------ 
! END OF GLOBAL REDUCTION
!------------------------
!
!$omp parallel do
!$omp private (DQSTIJ   , DWSTIJ  , HTMIJL  , RFACQ   , RFACW   , RFQIJ   , RFWIJ   , SUMNW   ,   &
!$omp          SUMNQ    , SUMPQ   , SUMPW   )
!
    DO K=1,LM
!    
        SUMPQ = GSUMS(1,K)
        SUMNQ = GSUMS(2,K)
        SUMPW = GSUMS(3,K)
        SUMNW = GSUMS(4,K)
!-------------------------------    
! FIRST MOMENT CONSERVING FACTOR
!-------------------------------
        IF (SUMPQ > 1.) THEN
            RFACQ = -SUMNQ / SUMPQ
        ELSE
            RFACQ = 1.
        END IF
!    
        IF (SUMPW > 1.) THEN
            RFACW = -SUMNW / SUMPW
        ELSE
            RFACW = 1.
        END IF
!    
        IF (RFACQ < 0.9 .OR. RFACQ > 1.1) RFACQ = 1.
        IF (RFACW < 0.9 .OR. RFACW > 1.1) RFACW = 1.
!-------------------------------------- 
! IMPOSE CONSERVATION ON ANTI-FILTERING 
!--------------------------------------
        IF (RFACQ < 1.) THEN
            DO J=MYJS2,MYJE2
                DO I=MYIS1,MYIE1
                    DQSTIJ = DQST(I,J,K)
                    RFQIJ  = HBM2(I,J) * (RFACQ-1.) + 1.
                    IF (DQSTIJ >= 0.) DQSTIJ = DQSTIJ * RFQIJ
                    Q(I,J,K) = Q1(I,J,K) + DQSTIJ
                END DO
            END DO
        ELSE
            DO J=MYJS2,MYJE2
                DO I=MYIS1,MYIE1
                    DQSTIJ = DQST(I,J,K)
                    RFQIJ  = HBM2(I,J) * (RFACQ-1.) + 1.
                    IF (DQSTIJ < 0.) DQSTIJ = DQSTIJ / RFQIJ
                    Q(I,J,K) = Q1(I,J,K) + DQSTIJ
                END DO
            END DO
        END IF
!
        IF (RFACW < 1.) THEN
            DO J=MYJS2,MYJE2
                DO I=MYIS1,MYIE1
                    DWSTIJ = DWST(I,J,K)
                    RFWIJ = HBM2(I,J) * (RFACW-1.) + 1.
                    IF (DWSTIJ >= 0.) DWSTIJ = DWSTIJ * RFWIJ
                    CWM(I,J,K) = W1(I,J,K) + DWSTIJ
                END DO
            END DO
        ELSE
            DO J=MYJS2,MYJE2
                DO I=MYIS1,MYIE1
                    DWSTIJ = DWST(I,J,K)
                    RFWIJ = HBM2(I,J) * (RFACW-1.) + 1.
                    IF (DWSTIJ < 0.) DWSTIJ = DWSTIJ / RFWIJ
                    CWM(I,J,K) = W1(I,J,K) + DWSTIJ
                END DO
            END DO
        END IF
!
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                HTMIJL = HTM(I,J,K)
                  Q(I,J,K) = AMAX1(  Q(I,J,K),   EPSQ) * HTMIJL
                CWM(I,J,K) = AMAX1(CWM(I,J,K), CLIMIT) * HTMIJL
            END DO
        END DO
!---------------  
! END OF LM LOOP 
!--------------- 
    END DO       
!
    RETURN
!
    END SUBROUTINE HZADV2
