    SUBROUTINE EPS
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE EPS
!> 
!> SUBPROGRAM: EPS - ?????
!> PROGRAMMER: JANJIC
!> ORG: W/NP22
!> DATE: 9?-??-??
!>
!> ABSTRACT:
!> EPS COMPUTES THE VERTICAL AND HORIZONTAL ADVECTION OF DZ / DT.
!>
!> PROGRAM HISTORY LOG:
!> 9?-??-??  JANJIC      - ORIGINATOR
!> 00-01-05  BLACK       - DISTRIBUTED MEMORY AND THREADS
!> 26-11-12  RISTIC IVAN - LOT OF CHANGES AND BUG FIX
!> 18-01-15  LUCCI       - MODERNIZATION OF THE CODE, INCLUDING:
!>                         * F77 TO F90/F95
!>                         * INDENTATION & UNIFORMIZATION CODE
!>                         * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                         * DOCUMENTATION WITH DOXYGEN
!>                         * OPENMP FUNCTIONALITY
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
!> USE MODULES: CONTIN
!>              CTLBLK
!>              DYNAM
!>              F77KINDS
!>              GLB_TABLE
!>              INDX
!>              LOOPS
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              NHYDRO
!>              PARMETA
!>              PVRBLS
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>  
!> DRIVER     : EBU
!>
!> CALLS      : EXCH
!>              MPI_ALLREDUCE
!--------------------------------------------------------------------------------------------------
    USE CONTIN
    USE CTLBLK
    USE DYNAM
    USE EXCHM
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
    USE PVRBLS
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    INCLUDE "EXCHM.h"
    INCLUDE "mpif.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: IM1   = IM - 1 
    INTEGER(KIND=I4KIND), PARAMETER :: JAMD  = (JAM * 2 - 10) * 3
    INTEGER(KIND=I4KIND), PARAMETER :: NTSHY =  3
!
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ  =    1.E-12 
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ2 =    0.12
    REAL   (KIND=R4KIND), PARAMETER :: FF1   =    0.52500
    REAL   (KIND=R4KIND), PARAMETER :: FF2   =   -0.64813
    REAL   (KIND=R4KIND), PARAMETER :: FF3   =    0.24520
    REAL   (KIND=R4KIND), PARAMETER :: FF4   =   -0.12189
    REAL   (KIND=R4KIND), PARAMETER :: EPSFC =    1. / 1.05
    REAL   (KIND=R4KIND), PARAMETER :: EPSN  = -EPSFC   
    REAL   (KIND=R4KIND), PARAMETER :: EPSP  =  EPSFC
    REAL   (KIND=R4KIND), PARAMETER :: ZERO  =    1.E-06
    REAL   (KIND=R4KIND), PARAMETER :: G     =    9.8
    REAL   (KIND=R4KIND), PARAMETER :: CP    = 1004.6
    REAL   (KIND=R4KIND), PARAMETER :: CAPA  =  287.04 / CP
    REAL   (KIND=R4KIND), PARAMETER :: GMA   = -287.04 * (1. - CAPA) / 2.
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & TOP     , BOT
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K       , LMP     , LAP     , LLAP    , JFP     , JFQ     , IRECV   ,   &
    & JHL     , JHH     , IHL     , IHH     , JX      , IX      , KN      , KP   
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ADDT    , RDT     , RR      , ARR     , DWP     , DETAL   , RFACW   , W4P     , AFRP    ,   &
    & D2PQW   , WP      , W00     , WP0     , ENH     , HM      , DVOLP   , TTA     , TTB     ,   &
    & PP      , QP      , DWSTIJ  , RFWIJ   , W1IJ    , WSTIJ   , W0Q     , WA      , DWDTMX  ,   &
    & DWDTMN  , DWDTT   , GDT     , GDT2    , WGHT    , FFC     , PDP     , RDP     , DPPL    ,   &
    & DPSTR   , PP1     , TFC     , TTFC    , PSTRUP  , RDPDN   , RDPUP   , RPD     , PSTRDN  ,   &
    & TMP     , HBM2IJ  , DPTU    , FCT     , DPTL    , DELP       
!--------------------------------------
! LOCAL VARIABLES REQUIRING DECLARATION 
!-------------------------------------- 
    REAL   (KIND=R4KIND), DIMENSION(LM+1)                                                       ::&
    & PONE    , PSTR    , PNP1    , COFF    , DFRC    , CHI     , CHIN    , WRES    , CLIM
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                                         ::&
    & RDPP    , C0      , B1      , B2      , B3      , W3      , W4
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                                         ::&
    & ETADTL  , DWL     , AFR
!
    INTEGER(KIND=I4KIND), DIMENSION(LM)                                                         ::&
    & LA
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & AFP     , AFQ     , W1      , DWST    , DVOL 
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & DARE    , EMH     , ANE     , ASE
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & IFPA    , IFQA    ,                                                                         &
    & IFPF    , IFQF    ,                                                                         &
    & JFPA    , JFQA    ,                                                                         &
    & JFPF    , JFQF
!
    REAL   (KIND=R4KIND), DIMENSION(2, LM)                                                      ::&
    & GSUMS   , XSUMS 
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & SUMPW   , SUMNW
!
    IF (NTSD <= NTSHY .OR. HYDRO) THEN
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO J=MYJS_P2,MYJE_P2
            DO I=MYIS_P1,MYIE_P1
                PINT(I,J,1) = PT
            END DO
        END DO
    
        DO K=1,LM
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
            DO J=MYJS_P2,MYJE_P2
                DO I=MYIS_P1,MYIE_P1
                     DWDT(I,J,K)  = 1.
                    PDWDT(I,J,K)  = 1.
                    PINT(I,J,K+1) = PDSL(I,J) * DETA(K) + PINT(I,J,K)
                END DO
            ENDDO
        END DO
!
        RETURN
!
    END IF
!
    ADDT = DT
    RDT  = 1. / ADDT
!--------------     
! TIME TENDENCY 
!--------------
!------- 
! OPENMP
!-------
! 
!$omp parallel do
!
    DO K=1,LM
        DO J=MYJS_P1,MYJE_P1
            DO I=MYIS_P1,MYIE_P1
                DWDT(I,J,K) = (W(I,J,K) - DWDT(I,J,K)) * HTM(I,J,K) * HBM2(I,J) * RDT
            END DO
        END DO
    END DO
!------------------- 
! VERTICAL ADVECTION 
!------------------- 
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!$omp private (AFR      , AFRP    , ARR     , BOT     , D2PQW   , DETAL   , DWL     , DWP     ,   &
!$omp          ETADTL   , LA      , LAP     , LLAP    , LMP     , RFACW   , RR      , SUMNW   ,   &
!$omp          SUMPW    , TOP     , W00     , W3      , W4      , W4P     , WP      , WP0     )
!
    DO 200 J=MYJS_P1,MYJE_P1
        DO 200 I=MYIS_P1,MYIE_P1
!
            LMP=LMH(I,J)
!
            DO K=1,LMP
                W3(K) = W(I,J,K)
                W4(K) = W3(K)
            END DO
        
            ETADTL(1) = ETADT(I,J,1) * 0.5
        
            DO K=2,LMP-1
                ETADTL(K) = (ETADT(I,J,K-1) + ETADT(I,J,K)) * 0.5
            END DO
!        
            ETADTL(LMP) = ETADT(I,J,LMP-1) * 0.5
!
            SUMPW = 0.
            SUMNW = 0.
        
            DO K=1,LMP
                RR  = ETADTL(K) * (-ADDT)
                ARR = ABS(RR)
            
                IF (ARR > 0.) THEN
                    LAP = RR / ARR
                ELSE
                    LAP = 0
                END IF
!            
                LA(K) = LAP
                LLAP  = K + LAP
!            
                TOP = .FALSE. 
                BOT = .FALSE. 
!            
                IF (LLAP > 0 .AND. LLAP < LMH(I,J)+1 .AND. LAP /= 0) THEN
                    RR     = ARR / ABS(AETA(LLAP) - AETA(K))
                    AFR(K) = (((FF4 * RR + FF3) * RR + FF2) * RR + FF1) * RR
                    DWP    = (W3(LLAP) - W3(K)) * RR
                    DWL(K) = DWP
                ELSE
                    TOP = LLAP == 0
                    BOT = LLAP == LMH(I,J)+1
!                
                    RR     = 0.
                    AFR(K) = 0.
                    DWL(K) = 0.
                END IF
            END DO
!        
            IF (TOP) THEN
                IF (LA(2) < 0) THEN
                    DWL(1)   = -DWL(2)     * DETA(2)     / DETA(1)
                END IF
            END IF
!        
            IF (BOT) THEN
                IF (LA(LMP-1) > 0)THEN
                    DWL(LMP) = -DWL(LMP-1) * DETA(LMP-1) / DETA(LMP)
                END IF
            END IF
!        
            DO K=1,LMP
                DETAL = DETA(K)
                DWP   = DWL (K) * DETAL
                IF (DWP > 0.) THEN
                    SUMPW = SUMPW + DWP
                ELSE
                    SUMNW = SUMNW + DWP
                END IF
            END DO
!-------------------------------
! FIRST MOMENT CONSERVING FACTOR
!-------------------------------
            IF (SUMPW > 1.E-9) THEN
                RFACW = -SUMNW / SUMPW
            ELSE
                RFACW = 1.
            END IF
!        
            IF(RFACW < 0.9 .OR. RFACW > 1.1) RFACW = 1.
!---------------------------------
! IMPOSE CONSERVATION ON ADVECTION
!---------------------------------
            IF (RFACW < 1.) THEN
                DO K=1,LMP
                    DWP = DWL(K)
                    IF (DWP < 0.)  DWP = DWP / RFACW
                    W4(K) = W3(K) + DWP
                END DO
            ELSE
                DO K=1,LMP
                    DWP = DWL(K)
                    IF (DWP >= 0.) DWP = DWP * RFACW
                    W4(K) = W3(K) + DWP
                END DO
            END IF
!-------------------- 
! ANTI-FILTERING STEP
!-------------------- 
            SUMPW = 0.
            SUMNW = 0.
!------------------------ 
! ANTI-FILTERING LIMITERS
!------------------------ 
            DO 50 K=2,LMP-1
                DETAL = DETA(K)
!            
                W4P = W4(K)
!            
                LAP = LA(K)
!            
                IF (LAP /= 0) THEN
                    AFRP  = 2. * AFR(K) * (AETA(K+LAP) - AETA(K)) ** 2 / (AETA(K+LAP) - AETA(K-LAP))
                    D2PQW = ((W4(K+LAP) - W4P)       / (AETA(K+LAP) - AETA(K))                          &
    &                     -  (W4P       - W4(K-LAP)) / (AETA(K)     - AETA(K-LAP)))                     &
    &                     *   AFRP
                ELSE
                    D2PQW = 0.
                END IF
!            
                WP = W4P - D2PQW
!            
                W00 = W3(K)
                WP0 = W3(K+LAP)
!            
                WP = MAX(WP,MIN(W00,WP0))
                WP = MIN(WP,MAX(W00,WP0))
!            
                DWP = WP - W4P
!            
                DWL(K) = DWP
!            
                DWP = DWP * DETAL
!            
                IF (DWP > 0.) THEN
                    SUMPW = SUMPW + DWP
                ELSE
                    SUMNW = SUMNW + DWP
                END IF
!            
         50 END DO
!
            DWL(1)   = 0.
        
            DWL(LMP) = 0.
!------------------------------- 
! FIRST MOMENT CONSERVING FACTOR
!-------------------------------
            IF (SUMPW > 1.E-9) THEN
                RFACW = -SUMNW / SUMPW
            ELSE
                RFACW = 1.
            END IF
!        
            IF (RFACW < 0.9 .OR. RFACW > 1.1) RFACW = 1.
!--------------------------------------
! IMPOSE CONSERVATION ON ANTI-FILTERING 
!--------------------------------------
            IF (RFACW < 1.) THEN
                DO K=1,LMP
                    DWP = DWL(K)
                    IF (DWP >= 0.) DWP = DWP * RFACW
                    DWDT(I,J,K) = DWDT(I,J,K) - (W4(K) + DWP - W3(K)) * RDT
                END DO
            ELSE
                DO K=1,LMP
                    DWP = DWL(K)
                    IF (DWP < 0.)  DWP = DWP / RFACW
                    DWDT(I,J,K) = DWDT(I,J,K) - (W4(K) + DWP - W3(K)) * RDT
                END DO
            END IF
!
200 END DO
!-------------------------- 
! END OF VERTICAL ADVECTION
!-------------------------- 
!
    ENH = ADDT / (08. * DY)
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
    DO J=MYJS_P2,MYJE_P2
        DO I=MYIS_P1,MYIE_P1
            EMH(I,J)  = ADDT / (08. * DX(I,J))
            DARE(I,J) = HBM2(I,J)   * DX(I,J) * DY
        END DO
    END DO
!--------------------- 
! HORIZONTAL ADVECTION 
!--------------------- 
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!$omp private (DVOLP    , DWSTIJ  , HM      , JFP     , JFQ     , PP      , QP      , SUMNW   ,   &
!$omp          SUMPW    , TTA     , TTB     )
!
    DO 300 K=1,LM
!
        DO J=MYJS_P2,MYJE_P2
            DO I=MYIS_P1,MYIE_P1
                DVOL(I,J,K) = DARE(I,J)   * PDSL(I,J) * DETA(K)
                HM          =  HTM(I,J,K)
                  W1(I,J,K) =    W(I,J,K) * HM
            END DO
        END DO
!
        SUMPW = 0.
        SUMNW = 0.
!    
        DO 225 J=MYJS2_P1,MYJE2_P1
            DO 225 I=MYIS1_P1,MYIE1_P1
!            
                DVOLP = DVOL(I,J,K) * HBM3(I,J)
                  TTA = (U(I,J-1,K) + U(I+IHW(J),J,K) + U(I+IHE(J),J,K) + U(I,J+1,K))             &
    &                 *  HBM2(I,J)  * EMH(I,J)
!
                  TTB = (V(I,J-1,K) + V(I+IHW(J),J,K) + V(I+IHE(J),J,K) + V(I,J+1,K))             &
    &                 *  HBM2(I,J)  * ENH
!            
                PP = -TTA - TTB
                QP =  TTA - TTB
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
                PP = ABS(PP) * HTM(I,J,K) * HTM(IFPA(I,J,K), JFPA(I,J,K),K)
                QP = ABS(QP) * HTM(I,J,K) * HTM(IFQA(I,J,K), JFQA(I,J,K),K)
!            
                AFP(I,J,K) = (((FF4 * PP + FF3) * PP + FF2) * PP + FF1) * PP
                AFQ(I,J,K) = (((FF4 * QP + FF3) * QP + FF2) * QP + FF1) * QP
!            
                DWSTIJ = (W(IFPA(I,J,K), JFPA(I,J,K),K) - W(I,J,K)) * PP                          &
    &                  + (W(IFQA(I,J,K), JFQA(I,J,K),K) - W(I,J,K)) * QP
!            
                DWST(I,J,K) = DWSTIJ
!            
    225 END DO
!---------------------------- 
! GLOBAL SUM FOR CONSERVATION
!---------------------------- 
        DO 230 J=MYJS2,MYJE2
            DO 230 I=MYIS1,MYIE1
!           
                DVOLP  = DVOL(I,J,K) * HBM3(I,J)
                DWSTIJ = DWST(I,J,K) * DVOLP
!            
                IF (DWSTIJ > 0.) THEN
                    SUMPW = SUMPW + DWSTIJ
                ELSE
                    SUMNW = SUMNW + DWSTIJ
                END IF
!            
    230 END DO
!
        XSUMS(1,K) = SUMPW
        XSUMS(2,K) = SUMNW
    
        CONTINUE !END OF VERTICAL LOOP
!
300 END DO
!----------------- 
! GLOBAL REDUCTION
!----------------- 
    CALL MPI_ALLREDUCE(XSUMS, GSUMS, 2*LM, MPI_REAL, MPI_SUM, MPI_COMM_COMP, IRECV)
!
    IF (IRECV /= 0) THEN
        PRINT*, ' RETURN CODE FROM 1ST ALLREDUCE IN EPS = ', IRECV
        STOP
    END IF
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!$omp private (D2PQW    , DVOLP   , DWSTIJ  , RFACW   , RFWIJ   , SUMNW   , SUMPW   , W00     ,   &
!$omp          W0Q      , W1IJ    , WP0     , WSTIJ   )
!
    DO 400 K=1,LM
!   
        SUMPW = GSUMS(1,K)
        SUMNW = GSUMS(2,K)
!-------------------------------    
! FIRST MOMENT CONSERVING FACTOR 
!-------------------------------
        IF (SUMPW > 1.) THEN
            RFACW = -SUMNW / SUMPW
        ELSE
            RFACW = 1.
        END IF
!
        IF (RFACW < 0.9 .OR. RFACW > 1.1) RFACW = 1.
!---------------------------------
! IMPOSE CONSERVATION ON ADVECTION 
!---------------------------------
        IF (RFACW < 1.) THEN
            DO J=MYJS2_P1,MYJE2_P1
                DO I=MYIS1_P1,MYIE1_P1
                    DWSTIJ = DWST(I,J,K)
                    RFWIJ = HBM3(I,J) * (RFACW-1.) + 1.
                    IF (DWSTIJ < 0.) DWSTIJ = DWSTIJ / RFWIJ
                    W1(I,J,K) = W(I,J,K) + DWSTIJ
                END DO
            END DO
        ELSE
            DO J=MYJS2_P1,MYJE2_P1
                DO I=MYIS1_P1,MYIE1_P1
                    DWSTIJ = DWST(I,J,K)
                    RFWIJ  = HBM3(I,J) * (RFACW-1.) + 1.
                    IF (DWSTIJ >= 0.) DWSTIJ = DWSTIJ * RFWIJ
                    W1(I,J,K) = W(I,J,K) + DWSTIJ
                END DO
            END DO
        END IF
!-------------------- 
! ANTI-FILTERING STEP
!-------------------- 
        SUMPW = 0.
        SUMNW = 0.
!------------------------ 
! ANTI-FILTERING LIMITERS
!------------------------   
        DO 350 J=MYJS2,MYJE2
            DO 350 I=MYIS1,MYIE1
!            
                DVOLP = DVOL(I,J,K)
                W1IJ  =   W1(I,J,K)
!            
                D2PQW = ((             W1(IFPA(I,J,K), JFPA(I,J,K), K)  - W1IJ)                   &
    &                 -  (W1IJ     -   W1(IFPF(I,J,K), JFPF(I,J,K), K))                           &
    &                 *               HTM(IFPF(I,J,K), JFPF(I,J,K), K))                           &
    &                 * AFP(I,J,K) + ((W1(IFQA(I,J,K), JFQA(I,J,K), K)  - W1IJ)                   &
    &                 -  (W1IJ     -   W1(IFQF(I,J,K), JFQF(I,J,K), k))                           &
    &                 *               HTM(IFQF(I,J,K), JFQF(I,J,K), K)) * AFQ(I,J,K)
!            
                WSTIJ = W1IJ - D2PQW
!            
                W00 = W(I,J,K)
                WP0 = W(IFPA(I,J,K), JFPA(I,J,K), K)
                W0Q = W(IFQA(I,J,K), JFQA(I,J,K), K)
!            
                WSTIJ = AMAX1(WSTIJ, AMIN1(W00,WP0,W0Q))
                WSTIJ = AMIN1(WSTIJ, AMAX1(W00,WP0,W0Q))
!            
                DWSTIJ = WSTIJ  - W1IJ
!            
                DWST(I,J,K) = DWSTIJ
!            
                DWSTIJ = DWSTIJ * DVOLP
!            
                IF (DWSTIJ > 0.) THEN
                    SUMPW = SUMPW + DWSTIJ
                ELSE
                    SUMNW = SUMNW + DWSTIJ
                END IF
!            
    350 END DO
!
        XSUMS(1,K) = SUMPW
        XSUMS(2,K) = SUMNW
!    
        CONTINUE !END OF VERTICAL LOOP
400 END DO
!----------------- 
! GLOBAL REDUCTION
!----------------- 
    CALL MPI_ALLREDUCE(XSUMS, GSUMS, 2*LM, MPI_REAL, MPI_SUM, MPI_COMM_COMP, IRECV)
!
    IF (IRECV /= 0) THEN
        PRINT*, ' RETURN CODE FROM 2ND ALLREDUCE IN EPS = ', IRECV
        STOP
    END IF
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!$omp private (DWSTIJ   , RFACW   , RFWIJ   , SUMNW   , SUMPW)
!
    DO 425 K=1,LM
!    
        SUMPW = GSUMS(1,K)
        SUMNW = GSUMS(2,K)
!-------------------------------    
! FIRST MOMENT CONSERVING FACTOR
!-------------------------------
        IF (SUMPW > 1.) THEN
            RFACW = -SUMNW / SUMPW
        ELSE
            RFACW = 1.
        END IF
!
        IF (RFACW < 0.9 .OR. RFACW > 1.1) RFACW = 1.
!--------------------------------------
! IMPOSE CONSERVATION ON ANTI-FILTERING
!--------------------------------------
        IF (RFACW < 1.) THEN
            DO J=MYJS2,MYJE2
                DO I=MYIS1,MYIE1
                    DWSTIJ = DWST(I,J,K)
                    RFWIJ = HBM2(I,J) * (RFACW-1.) + 1.
                    IF (DWSTIJ >= 0.) DWSTIJ = DWSTIJ * RFWIJ
                    DWDT(I,J,K) = (DWDT(I,J,K) - (W1(I,J,K) + DWSTIJ - W(I,J,K)) * RDT) * HBM3(I,J)
                END DO
            END DO
        ELSE
            DO J=MYJS2,MYJE2
                DO I=MYIS1,MYIE1
                    DWSTIJ = DWST(I,J,K)
                    RFWIJ = HBM2(I,J) * (RFACW-1.) + 1.
                    IF (DWSTIJ < 0.) DWSTIJ = DWSTIJ / RFWIJ
                    DWDT(I,J,K) = (DWDT(I,J,K) - (W1(I,J,K) + DWSTIJ - W(I,J,K)) * RDT) * HBM3(I,J)
                END DO
            END DO
        END IF
!
425 END DO
!--------------------------------------------------
! RESTRICTING THE ACCELERATION ALONG THE BOUNDARIES
!--------------------------------------------------
    JHL = 0
!
    IF (JHL > 1) THEN
        JHH = JM  - JHL + 1
!    
        IHL = JHL / 2   + 1
        IHH = IM  - IHL + MOD(J,2)
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (IX     , JX)
!
        DO 450 K=1,LM
!        
            DO J=1,JHL
                IF (J >= MY_JS_GLB-JBPAD2 .AND. J <= MY_JE_GLB+JTPAD2) THEN
                    JX = J - MY_JS_GLB + 1
                    DO I=1,IM
                        IF (I >= MY_IS_GLB-ILPAD2 .AND. I <= MY_IE_GLB+IRPAD2) THEN
                            IX = I - MY_IS_GLB + 1
                            DWDT(IX,JX,K) = MAX(DWDT(IX,JX,k) * HBM3(IX,JX), -0.03)
                            DWDT(IX,JX,K) = MIN(DWDT(IX,JX,K) * HBM3(IX,JX),  0.03)
                        END IF
                    END DO
                END IF
            END DO
!        
            DO J=JHH,JM
                IF (J >= MY_JS_GLB-JBPAD2 .AND. J <= MY_JE_GLB+JTPAD2) THEN
                    JX = J - MY_JS_GLB + 1
                    DO I=1,IM
                        IF (I >= MY_IS_GLB-ILPAD2 .AND. I <= MY_IE_GLB+IRPAD2) THEN
                            IX = I - MY_IS_GLB + 1
                            DWDT(IX,JX,K) = MAX(DWDT(IX,JX,K) * HBM3(IX,JX), -0.03)
                            DWDT(IX,JX,K) = MIN(DWDT(IX,JX,K) * HBM3(IX,JX),  0.03)
                        END IF
                    END DO
                END IF
            END DO
!        
            DO J=1,JM
                IF (J >= MY_JS_GLB-JBPAD2 .AND. J <= MY_JE_GLB+JTPAD2) THEN
                    JX = J - MY_JS_GLB + 1
                    DO I=1,IHL
                        IF (I >= MY_IS_GLB-ILPAD2 .AND. I <= MY_IE_GLB+IRPAD2) THEN
                            IX = I - MY_IS_GLB + 1
                            DWDT(IX,JX,L) = MAX(DWDT(IX,JX,L) * HBM3(IX,JX), -0.03)
                            DWDT(IX,JX,L) = MIN(DWDT(IX,JX,L) * HBM3(IX,JX),  0.03)
                        END IF
                    END DO
                END IF
            END DO
!        
            DO J=1,JM
                IF (J >= MY_JS_GLB-JBPAD2 .AND. J <= MY_JE_GLB+JTPAD2) THEN
                    JX = J - MY_JS_GLB + 1
                    DO I=IHH,IM
                        IF (I >= MY_IS_GLB-ILPAD2 .AND. I <= MY_IE_GLB+IRPAD2) THEN
                            IX = I - MY_IS_GLB + 1
                            DWDT(IX,JX,K) = MAX(DWDT(IX,JX,K) * HBM3(IX,JX), -0.03)
                            DWDT(IX,JX,K) = MIN(DWDT(IX,JX,K) * HBM3(IX,JX),  0.03)
                        END IF
                    END DO
                END IF
            END DO
!        
    450 END DO
!
    END IF
!
    CALL EXCH(DWDT,LM,1,1)
!
    WA = 0.15
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (ANE    , ASE)
!
    DO 500 K=1,LM
!    
        DO J=MYJS1_P1,MYJE2_P1
            DO I=MYIS_P1,MYIE1_P1
                ANE(I,J) = (DWDT(I+IHE(J),J+1,K) - DWDT(I,J,K)) * HTM(I,J,K)                      &
    &                    *   HTM(I+IHE(J),J+1,K)
            END DO
        END DO
!    
        DO J=MYJS2_P1,MYJE1_P1
            DO I=MYIS_P1,MYIE1_P1
                ASE(I,J) = (DWDT(I+IHE(J),J-1,K) - DWDT(I,J,K)) * HTM(I+IHE(J),J-1,K)             &
    &                    *   HTM(I,J,K)
            END DO
        END DO
!    
        DO J=MYJS2,MYJE2
            DO I=MYIS1,MYIE1
                DWDT(I,J,K) = (ANE(I,J) - ANE(I+IHW(J),J-1) + ASE(I,J) - ASE(I+IHW(J),J+1))       &
    &                       *  WA * HBM2(I,J) + DWDT(I,J,K)
            END DO
        END DO
!
500 END DO
!
    WP     = 0.075
    KN     = 0
    KP     = 0
    DWDTMX = 0.
    DWDTMN = 0.
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!$omp private (DWDTMN   , DWDTMX  , DWDTT   , HM      , KN      , KP)
!
    DO 525 K=1,LM
!
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                HM     =  HTM(I,J,K) * HBM2(I,J)
                DWDTT  = DWDT(I,J,K) * HM
!            
                DWDTMX = AMAX1(DWDTT,DWDTMX)
                DWDTMN = AMIN1(DWDTT,DWDTMN)
!            
                IF (DWDTT < EPSN) THEN
                    DWDTT = EPSN
                    KN    = KN + 1
                END IF
!            
                IF (DWDTT > EPSP) THEN
                    DWDTT = EPSP
                    KP    = KP + 1
                END IF
!            
                DWDT(I,J,K) = (DWDTT/G+1.) * (1.-WP) + PDWDT(I,J,K) * WP
            END DO
        END DO
!
525 END DO
!
    GDT  = G * DT
    GDT2 = GDT * GDT
    WGHT = 0.35
    FFC  = -4. * R / GDT2
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!$omp private (CHI      , CHIN    , CHMOD   , CLIM    , COFF    , DELP    , DFRC    , DP      ,   &
!$omp          DPSTR    , DPTL    , DPTU    , DWDTT   , HBM2IJ  , IMAX    , INOT    , ITER    ,   &
!$omp          ITMX     , JMAX    , K       , LMP     , PDP     , PNP1    , PONE    , PP1     ,   &
!$omp          PSTR     , RDP     , RESDL   , WIL     , WRES    )
!
    DO 600 J=MYJS2,MYJE2
        DO 600 I=MYIS1,MYIE1
!        
            LMP =  LMH(I,J)
            PDP = PDSL(I,J)
            RPD = 1. / PDP
!
            PONE(1) = PT
            PSTR(1) = PT
            PNP1(1) = PT
            CHI(1)  = 0.
!        
            DO K=2,LMP+1
                CHI(K)    = 0.
                DPPL      = DETA(K-1) * PDP
                RDPP(K-1) = 1. / DPPL
                PONE(K)   = PINT(I,J,K)
                DPSTR     = DWDT(I,J,K-1) * DPPL
                PSTR(K)   = PSTR(K-1) + DPSTR
                PP1       = PNP1(K-1) + DPSTR
                PNP1(K)   = (PP1-PONE(K)) * WGHT + PONE(K)
                TFC       = Q(I,J,K-1) * 0.608 + 1.
                TTFC      = -CAPA * TFC + 1.
                COFF(K-1) = T(I,J,K-1) * TTFC * TFC * DPPL * FFC / ((PNP1(K-1)+PNP1(K)) *         &
    &                                                               (PNP1(K-1)+PNP1(K)))
            END DO
!        
            PSTRUP = -(PSTR(1) + PSTR(2) - PONE(1) - PONE(2)) * COFF(1)
!        
            DO K=2,LMP
                RDPDN = RDPP(K)
                RDPUP = RDPP(K-1)
!            
                PSTRDN = -(PSTR(K) + PSTR(K+1) - PONE(K) - PONE(K+1)) * COFF(K)
!            
                B1(K)  =  COFF(K-1) + RDPUP
                B2(K)  = (COFF(K-1) + COFF(K)) - (RDPUP+RDPDN)
                B3(K)  =  COFF(K)   + RDPDN
                C0(K)  = PSTRUP     + PSTRDN
!            
                PSTRUP = PSTRDN
            END DO
!        
            B1(2)   = 0.
            B2(LMP) = B2(LMP) + B3(LMP)
!------------ 
! ELIMINATION 
!------------ 
            DO K=3,LMP
                TMP   = -B1(K)   / B2(K-1)
                B2(K) =  B3(K-1) * TMP + B2(K)
                C0(K) =  C0(K-1) * TMP + C0(K)
            END DO
!        
            CHI(1) = 0.
!------------------
! BACK SUBSTITUTION 
!------------------
            CHI(LMP)   =  C0(LMP) / B2(LMP)
            CHI(LMP+1) = CHI(LMP)
!        
            DO K=LMP-1,2,-1
                CHI(K) = (-B3(k) * CHI(K+1) + C0(K)) / B2(K)
            END DO
!
            HBM2IJ = HBM2(I,J)
            DPTU   = 0.
            FCT    = 0.5 / CP * HBM2IJ
!        
            DO K=1,LMP
                DPTL          = (CHI(K+1)      + PSTR(K+1)    - PINT(I,J,K+1)) * HBM2IJ
                PINT(I,J,K+1) = PINT(I,J,K+1)  + DPTL
                T   (I,J,K)   = (DPTU+DPTL)    * RTOP(I,J,K)  * FCT            + T(I,J,K)
                DELP          = (PINT(I,J,K+1) - PINT(I,J,K)) * RDPP(K)
                W(I,J,K)      = ((DELP         - DWDT(I,J,K)) * GDT            + W(I,J,K)) * HBM2IJ
                DWDT(I,J,K)   = (DELP-1.)      * HBM2IJ       + 1.
!
                DPTU = DPTL
            END DO
!        
            DO K=LMP+1,LM
                PINT(I,J,K+1) = PDSL(I,J) * DETA(K) + PINT(I,J,K)
            END DO
!
600 END DO
!
    RETURN
!
    END SUBROUTINE EPS
