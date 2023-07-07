    SUBROUTINE VTADV
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE VTADV
!>
!> SUBPROGRAM: VTADV - VERTICAL ADVECTION
!> PROGRAMMER: JANJIC
!> ORG: W/NP22
!> DATE: 93-11-17
!>
!> ABSTRACT:
!> VTADV CALCULATES THE CONTRIBUTION OF THE VERTICAL ADVECTION TO THE TENDENCIES OF TEMPERATURE,
!> SPECIFIC HUMIDITY, WIND COMPONENTS, AND TURBULENT KINETIC ENERGY AND THEN UPDATES THOSE
!> VARIABLES. FOR ALL VARIABLES EXCEPT SPECIFIC HUMIDITY A SIMPLE CENTERED DIFFERENCE SCHEME IN 
!> SPACE IS USED IN CONJUNCTION WITH THE PURE EULER-BACKWARD TIME SCHEME.
!> A PIECEWISE LINEAR SCHEME IS USED TO CALCULATE THE VERTICAL ADVECTION OF SPECIFIC HUMIDITY SO 
!> THAT NO FALSE MAXIMA OR MINIMA ARE PRODUCED.
!>
!> PIECEWISE LINEAR SCHEME (MESINGER AND JOVIC, NCEP OFFICE NOTE #439) HERE USED FOR ALL VARIABLES,
!> AVOIDING FALSE ADVECTION FROM BELOW GROUND, AND AT THE SAME TIME MAKING THE ETA VERY NEARLY A 
!> FINITE VOLUME MODEL.
!>
!> PROGRAM HISTORY LOG:
!> 87-06-??  JANJIC     - ORIGINATOR
!> 90-??-??  MESINGER   - INSERTED PIECEWISE LINEAR SCHEME FOR SPECIFIC HUMIDITY
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 95-11-20  ABELES     - PARALLEL OPTIMIZATION
!> 96-03-29  BLACK      - ADDED EXTERNAL EDGE; REMOVED SCRCH COMMON
!> 98-10-30  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
!> 00-02-04  BLACK      - ADDED CLOUD WATER/ICE
!> 01-12-11  BLACK      - SMOOTHING FOR CFL VIOLATION
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
!>              NHYDRO
!>              PARMETA
!>              PVRBLS
!>              TEMPCOM
!>              TOPO
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!v
!> DRIVER     : EBU
!>              NEWFLT
!>
!> CALLS      : -----                 
!>--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------
! THIS IS THE CODE AS E-MAILED APRIL 10, 07, BUT WITH BUGS CORRECTED:
! 1)  FFD1, TWICE, 1 IS REMOVED,
! 2)  F4D, IN LOOP 1320, THE SECOND PAIR CHANGED INTO F4Q;
! 3)  ..TEND1.. 2ND TIME IN LOOP 1330 REPL. BY ..TEND2..
! AND ALSO:
! A) DEAD PIECES OF THE CODE - ABANDONED FINITE DIFFERENCE VERTICAL ADVECTION OF T AND U, V 
!    (TWICE "GO TO ..." SEGMENTS) - DELETED;
! B) SWITCHING TO CALLS EVERY TIME STEP COMMENTED OUT, AS THIS TAKES NON-NEGLIGIBLE TIME, AND
!    IN THE ETA HORIZONTAL ADVECTION IS CALLED EVERY SECOND TIME STEP ALSO (CONSISTENCY)
!--------------------------------------------------------------------------------------------
    USE CLDWTR
    USE CONTIN
    USE CTLBLK
    USE DYNAM
    USE F77KINDS
    USE GLB_TABLE
    USE INDX
    USE LOOPS    , ONLY : JAM
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
    REAL   (KIND=R4KIND), PARAMETER :: EDQMX   =  2.E-5
    REAL   (KIND=R4KIND), PARAMETER :: EDQMN   = -2.E-5
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ    =  1.E-12
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ2   =  0.12
    REAL   (KIND=R4KIND), PARAMETER :: CFL_MAX =  0.97
!
    INTEGER(KIND=I4KIND), PARAMETER :: KSMUD   =  0
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM    = IM * JM - JM / 2
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & NOSLA   , NOSLAW 
!
    REAL   (KIND=R4KIND), DIMENSION(LM1)                                                        ::&
    & WFA     , WFB
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & ETADTL  ,                                                                                   &
    & TQ2B    ,                                                                                   &
    & DQTI    , DQBI    ,                                                                         &
    & RPDX    , RPDY    ,                                                                         &
    & QDEDB   , QDEUB   ,                                                                         &
    & EDBD    ,                                                                                   &
    & DQDEB   ,                                                                                   &
    & DUTI    , DUBI    ,                                                                         &
    & DVTI    , DVBI    ,                                                                         &
    & UDEDB   , UDEUB   ,                                                                         &
    & VDEDB   , VDEUB   ,                                                                         &
    & EDBDU   , EDBDV   ,                                                                         &
    & DUDEB   , DVDEB   ,                                                                         & 
    & FNE     , FSE 
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & SAM     ,                                                                                   &
    & QBI     ,                                                                                   &
    & Q2ST    ,                                                                                   &
    & SAMU    ,                                                                                   &
    & UBI     ,                                                                                   &
    & SAMV    ,                                                                                   &
    & VBI     
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM1)                              ::&
    & ARRAY1  ,                                                                                   &
    & ARRAY2  ,                                                                                   &
    & ARRAYU1 ,                                                                                   &
    & ARRAYU2 ,                                                                                   &
    & ARRAYV1 ,                                                                                   &
    & ARRAYV2 
!
    REAL   (KIND=R4KIND), DIMENSION(:,:,:), ALLOCATABLE                                         ::&
    & S
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & VAD_TEND1         , VAD_TEND2
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                                         ::&
    & VAD_TNDX1         , VAD_TNDX2
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & LBOT_CFL_T        , LTOP_CFL_T        ,                                                     &
    & LBOT_CFL_U        , LTOP_CFL_U        ,                                                     &
    & LBOT_CFL_V        , LTOP_CFL_V
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & K       , NMSAP   , NMSAPW  , I       , J       , NSMUD   , KS      , NS      , MSA     ,   &
    & IER     , LSTART  , LSTOP
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & DTAD    , CFL     , EXTREM  , DQTIK   , ASTIK   , ASBIK   , QDEDTK  , QDEUTK  , SEDBK   ,   &
    & DQDEK   , EDBFK   , EDTDK   , TQ2AK   , EXTREMU , EXTREMV , DUTIK   , DVTIK   , ASTIKU  ,   &
    & ASTIKV  , ASBIKU  , ASBIKV  , VMIJ    , UDEDTK  , VDEDTK  , VDEUTK  , UDEUTK  , SEDBKU  ,   &
    & SEDBKV  , DUDEK   , DVDEK   , EDBFKU  , EDBFKV  , EDTDKU  , EDTDKV
!
    DTAD = IDTAD * DT
!------------------------------------------ 
! DEFINE ADDED UPSTREAM ADVECTION CONSTANTS 
!------------------------------------------
!    DO 25 K=1,LM1
!        WFA(K) = DETA(K  ) / (DETA(K) + DETA(K+1))
!        WFB(K) = DETA(K+1) / (DETA(K) + DETA(K+1))
! 25 END DO
       WFA(1:LM1) = DETA(1:LM1  ) / (DETA(1:LM1) + DETA(2:LM))
       WFB(1:LM1) = DETA(2:LM) / (DETA(1:LM1) + DETA(2:LM))  
!------------------------------------------- 
! NO MOISTURE SLOPE ADJUSTMENT IF NOT WANTED
!-------------------------------------------
    NOSLA  = .FALSE. 
    NOSLAW = .FALSE. 
!----------------------------------------------------- 
! IF FALSE, NUMBER OF MOISTURE SLOPE ADJUSTMENT PASSES
!----------------------------------------------------- 
    NMSAP  = 3
    NMSAPW = 3
!----------------------------------------  
! SMOOTHING VERTICAL VELOCITY AT H POINTS 
!----------------------------------------
    IF (KSMUD > 0) THEN
!
!$omp parallel do 
!
        DO 90 K=1,LM1
            DO 50 J=MYJS_P4,MYJE_P4
                DO 50 I=MYIS_P4,MYIE_P4
                    ETADT(I,J,K) = ETADT(I,J,K) * HBM2(I,J)
         50 END DO
!
            NSMUD = KSMUD
!-------------------------------------------------------------------------
! HE FNE, FSE, ETADTL, AND ETADT ARRAYS ARE ON OR ASSOCIATED WITH H POINTS
!-------------------------------------------------------------------------
            DO 90 KS=1,NSMUD
                DO 80 J=MYJS_P3,MYJE1_P3
                    DO 80 I=MYIS_P3,MYIE_P3
                        FNE(I,J) = (ETADT(I+IHE(J),J+1,K  ) - ETADT(I       ,J  ,K  ))            & 
    &                            *    HTM(I       ,J  ,K+1) *   HTM(I+IHE(J),J+1,K+1)
             80 END DO
!
                DO 82 J=MYJS1_P3,MYJE_P3
                    DO 82 I=MYIS_P3,MYIE_P3
                        FSE(I,J) = (ETADT(I+IHE(J),J-1,K  ) - ETADT(I,J,K  ))                     &
    &                            *    HTM(I+IHE(J),J-1,K+1) *   HTM(I,J,K+1)
             82 END DO
!
                DO 84 J=MYJS2_P1,MYJE2_P1
                    DO 84 I=MYIS_P1,MYIE_P1
                        ETADTL(I,J) = (FNE(I,J) - FNE(I+IHW(J),J-1)                               &
    &                               +  FSE(I,J) - FSE(I+IHW(J),J+1))                              &
    &                               * HBM2(I,J)
             84 END DO
!
                DO 86 J=MYJS2_P1,MYJE2_P1
                    DO 86 I=MYIS_P1,MYIE_P1
                        ETADT(I,J,K) = ETADTL(I,J) * 0.125 + ETADT(I,J,K)
             86 END DO
!
     90 END DO
!
    END IF
!------------------------------------------------------------------------------------------
!  IF THE CFL CRITERION IS VIOLATED THEN LOCATE VERTICAL LIMITS BETWEEN WHICH TO SMOOTH THE
!  TENDENCIES
!------------------------------------------------------------------------------------------
!
!$omp parallel do
! 
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            LTOP_CFL_T(I,J) = 0
            LBOT_CFL_T(I,J) = 0
            LTOP_CFL_U(I,J) = 0
            LBOT_CFL_U(I,J) = 0
            LTOP_CFL_V(I,J) = 0
            LBOT_CFL_V(I,J) = 0
        END DO
    END DO
!
    DO K=1,LM1
!
!$omp parallel do private (CFL)   
! 
        DO J=MYJS2,MYJE2
            DO I=MYIS,MYIE
!------------             
! MASS POINTS
!------------           
                CFL = ETADT(I,J,K) * DTAD * HBM2(I,J)/(0.5*(DETA(K)+DETA(K+1)))
                IF (ABS(CFL) > CFL_MAX) THEN
                    IF (LTOP_CFL_T(I,J) == 0) LTOP_CFL_T(I,J) = MAX(K,2)
                    IF (LBOT_CFL_T(I,J) < K ) LBOT_CFL_T(I,J) = MAX(K,2)
                END IF
!------------            
! U COMPONENT
!------------            
                CFL = (ETADT(I+IVW(J),J,K) + ETADT(I+IVE(J),J,K))                                 &
    &                * DTAD * VBM2(I,J) / (DETA(K)+DETA(K+1))
!            
                IF (ABS(CFL) > CFL_MAX) THEN
                    IF (LTOP_CFL_U(I,J) == 0) LTOP_CFL_U(I,J) = MAX(K,2)
                    IF (LBOT_CFL_U(I,J) < K ) LBOT_CFL_U(I,J) = MAX(K,2)
                END IF
!------------           
! V COMPONENT
!------------            
                CFL = (ETADT(I,J-1,K) + ETADT(I,J+1,K)) * DTAD * VBM2(I,J) / (DETA(K) + DETA(K+1))
!            
                IF (ABS(CFL) > CFL_MAX) THEN
                    IF (LTOP_CFL_V(I,J) == 0) LTOP_CFL_V(I,J) = MAX(K,2)
                    IF (LBOT_CFL_V(I,J) < K ) LBOT_CFL_V(I,J) = MAX(K,2)
                END IF
!            
            END DO
        END DO
!    
    END DO
!------------------------------------------------------------
! PIECEWISE LINEAR UPSTREAM VERTICAL ADVECTION OF Q AND CLOUD 
!------------------------------------------------------------
    ALLOCATE(S(IDIM1:IDIM2, JDIM1:JDIM2,LM), STAT=I)
!-------------------------------------------------------------------------------------------------
! INTIALIZE Q AT THE BOTTOM INTERFACE AND THE SLOPE ADJUSTMENT MASK (SAM=1 FOR SA PERMITTED, 0 FOR
! NOT PERMITTED)
!-------------------------------------------------------------------------------------------------
!
!--------------------------------
! LOOP OVER S VARIABLES (Q,CWM,T)
!--------------------------------
    DO 400 NS=1,3
!
        IF (NS == 1) THEN
!
!$omp parallel do
!
            DO K=1,LM
                DO J=JDIM1,JDIM2
                    DO I=IDIM1,IDIM2
                        S(I,J,K) = Q(I,J,K)
                    END DO
                END DO
            END DO
        ELSE IF (NS == 2) THEN
!
!$omp parallel do
! 
            DO K=1,LM
                DO J=JDIM1,JDIM2
                    DO I=IDIM1,IDIM2
                        S(I,J,K) = CWM(I,J,K)
                    END DO
                END DO
            END DO
        ELSE IF (NS == 3) THEN
!
!$omp parallel do 
!
            DO K=1,LM
                DO J=JDIM1,JDIM2
                    DO I=IDIM1,IDIM2
                        S(I,J,K) = T(I,J,K)
                    END DO
                END DO
            END DO
        ELSE
            WRITE(0,*)'ERROR IN VTADV. WILL STOP NOW.'
            STOP
        END IF
!
!$omp parallel do 
!
        DO 175 K=1,LM
            DO 175 J=MYJS2,MYJE2
                DO 175 I=MYIS,MYIE
                    QBI(I,J,K) = S(I,J,K)
                    SAM(I,J,K) = 1.
    175 END DO
!
        IF (NOSLA) GOTO 290
!------------------------------------------------------
! THE SLOPE ADJUSTMENT CODE
! NO SLOPE PERMITTED AT THE TOP AND AT THE BOTTOM LAYER
!------------------------------------------------------ 
!
!$omp parallel do 
!
        DO 190 J=MYJS2,MYJE2
            DO 190 I=MYIS,MYIE
                SAM(I,J, 1) = 0.
                SAM(I,J,LM) = 0.
    190 END DO
!
!$omp parallel do 
!
        DO 200 K=1,LM1
            DO 200 J=MYJS2,MYJE2
                DO 200 I=MYIS,MYIE
                    SAM(I,J,K) = SAM(I,J,K) * HTM(I,J,K+1)
    200 END DO
!----------------------------------------------------------------------------------------
! NOW, SEARCH FOR THE MAXIMA AND MINIMA OF Q (AT THE FIRST PASS) AND FOR LAYERS WHICH HAD
! OVERADJUSTED (AT SUBSEQUENT PASSES) DUE TO ROUND-OFF ERRORS
!----------------------------------------------------------------------------------------
!
!$omp parallel do private (DQBI, DQTI, EXTREM)
!
        DO 220 K=2,LM1
            DO 220 J=MYJS2,MYJE2
                DO 220 I=MYIS,MYIE
                    DQTI(I,J) = S(I,J,K  ) - S(I,J,K-1)
                    DQBI(I,J) = S(I,J,K+1) - S(I,J,K  )
!
                    EXTREM    = DQTI(I,J) * DQBI(I,J)
!
                    IF (EXTREM <= 0.) SAM(I,J,K) = 0.
    220 END DO
!
!$omp parallel do 
!
        DO 230 K=2,LM1
            DO 230 J=MYJS2,MYJE2
                DO 230 I=MYIS,MYIE
                    ARRAY1(I,J,K) = WFA(K-1) * (1.-SAM(I,J,K-1)) + WFB(K-1)
                    ARRAY2(I,J,K) = WFA(K  ) + WFB(K) * (1.-SAM(I,J,K+1))
    230 END DO
!
        DO 260 MSA=1,NMSAP
!-------------------------------------------------------------------------------------------------
! CALCULATE DQ AT INTERFACES AND ADJUST THE SLOPES WHERE AND TO THE EXTENT PERMITTED OBSERVING THE
! MONOTONICITY CONDITION (E.G. VAN LEER, J. COMP. PHYS. 1977, 276-299)
!-------------------------------------------------------------------------------------------------
!
!$omp parallel do 
!
            DO 240 J=MYJS2,MYJE2
                DO 240 I=MYIS,MYIE
                    DQBI(I,J) = 2. * S(I,J,2) - QBI(I,J,2) - QBI(I,J,1)
        240 END DO
!        
            DO 250 K=2,LM1
!
!$omp parallel do private (ASBIK, ASTIK, DQTIK)
!
                DO 250 J=MYJS2,MYJE2
                    DO 250 I=MYIS,MYIE
                        DQTIK =   DQBI(I,J)
                        ASTIK = ARRAY1(I,J,K) * DQTIK
                        DQBI(I,J) = 2. * S(I,J,K+1) - QBI(I,J,K+1) - QBI(I,J,K)
                        ASBIK = ARRAY2(I,J,K) * DQBI(I,J)
                         QBI(I,J,K) = QBI(I,J,K)                                                  &
    &                               + (ASTIK-SIGN(1.,ASTIK) * DIM(ABS(ASTIK), ABS(ASBIK)))        &
    &                               *  SAM(I,J,K)
            250 END DO
!
    260 END DO
!------------------------------------------------------------------------------------------------
! SLOPE ADJUSTMENT OF THE LAYERS ABOVE THAT NEXT TO THE SURFACE IS DONE; NOW ADJUST THE LOWERMOST
! LAYER
!------------------------------------------------------------------------------------------------
        DO 270 K=9,LM1
!
!$omp parallel do 
!
            DO 270 J=MYJS2,MYJE2
                DO 270 I=MYIS,MYIE
                    IF (HTM(I,J,K+1) == 0.) QBI(I,J,K) = 2. * S(I,J,K) - QBI(I,J,K-1)
        270 END DO
!
!$omp parallel do 
!
        DO 280 J=MYJS2,MYJE2
            DO 280 I=MYIS,MYIE
                QBI(I,J,LM) = 2. * S(I,J,LM) - QBI(I,J,LM1)
    280 END DO
!--------------------------------- 
! END OF THE SLOPE ADJUSTMENT CODE 
!--------------------------------- 
    290 CONTINUE
!
!$omp parallel do 
!
        DO 300 J=MYJS2,MYJE2
            DO 300 I=MYIS,MYIE
                QDEDB(I,J) = 0.
                QDEUB(I,J) = 0.
                DQDEB(I,J) = 2. * (QBI(I,J,1) - S(I,J,1)) * RDETA(1)
                EDBD (I,J) = 0.
    300 END DO
!    
        DO 320 K=1,LM1
!
!$omp parallel do private (DQDEK, EDBFK, EDTDK, QDEDTK, QDEUTK, SEDBK)
!
            DO 320 J=MYJS2,MYJE2
                DO 320 I=MYIS,MYIE
                    QDEDTK     = QDEDB(I,J)
                    QDEUTK     = QDEUB(I,J)
                    SEDBK      = SIGN(1.,ETADT(I,J,K))
                    DQDEK      = DQDEB(I,J)
                    DQDEB(I,J) = 2. * (QBI(I,J,K+1) - S(I,J,K+1)) * RDETA(K+1)
                    EDBFK      = ETADT(I,J,K) * F4D
                    QDEDB(I,J) = (1.+SEDBK) * (QBI(I,J,K) + DQDEK * EDBFK) * (-EDBFK)
                    QDEUB(I,J) = (1.-SEDBK) * (2.*S(I,J,K+1) - QBI(I,J,K+1) + DQDEB(I,J)*EDBFK)   &
    &                          *  EDBFK
!
                    EDTDK      = EDBD(I,J)
                    EDBD (I,J) = ETADT(I,J,K) * (-F4Q)
                    S(I,J,K)   = S(I,J,K) + (QDEDTK - QDEUTK - QDEDB(I,J) + QDEUB(I,J)            &
    &                          + S(I,J,K) * (EDBD(I,J)-EDTDK)) * RDETA(K)
    320 END DO
!
!$omp parallel do 
!
        DO 330 J=MYJS2,MYJE2
            DO 330 I=MYIS,MYIE
                S(I,J,LM) = S(I,J,LM) + (QDEDB(I,J)  -  QDEUB(I,J)                                &
    &                     + S(I,J,LM) * (-EDBD(I,J))) * RDETA(LM)
    330 END DO
!-------------------------------------------------------- 
! NEGATIVE MOISTURE MAY OCCUR DUE TO VIOLATION OF THE CFL
!--------------------------------------------------------
        DO 350 K=1,LM1
!
!$omp parallel do 
!
            DO 350 J=MYJS2,MYJE2
                DO 350 I=MYIS,MYIE
                    IF (S(I,J,K) < EPSQ) THEN
                        DQBI(I,J)  = S(I,J,K  )
                        S(I,J,K  ) = EPSQ
                        S(I,J,K+1) = S(I,J,K+1) + DETA(K) * RDETA(K+1) * DQBI(I,J)
                    END IF
        350 END DO
!
!$omp parallel do 
!
        DO 360 J=MYJS2,MYJE2
            DO 360 I=MYIS,MYIE
                IF (S(I,J,LM) < EPSQ) S(I,J,LM) = EPSQ
    360 END DO
!
        IF (NS == 1) THEN
!
!$omp parallel do 
!
            DO K=1,LM
                DO J=JDIM1,JDIM2
                    DO I=IDIM1,IDIM2
                        Q(I,J,K) = S(I,J,K)
                    END DO
                END DO
            END DO
        ELSE IF (NS == 2) THEN
!
!$omp parallel do 
!
            DO K=1,LM
                DO J=JDIM1,JDIM2
                    DO I=IDIM1,IDIM2
                        CWM(I,J,K) = S(I,J,K)
                    END DO
                END DO
            END DO
        ELSE IF (NS == 3) THEN
!
!$omp parallel do 
!
            DO K=1,LM
                DO J=JDIM1,JDIM2
                    DO I=IDIM1,IDIM2
                        T(I,J,K) = S(I,J,K)
                    END DO
                END DO
            END DO
        ELSE
            WRITE(0,*)'ERROR IN VTADV. WILL STOP NOW.'
            STOP
        END IF
!
400 END DO
!
    DEALLOCATE(S,STAT=IER)
!----------------------------------- 
! VERTICAL (MATSUNO) ADVECTION OF Q2 
!----------------------------------- 
!
!$omp parallel do 
!
    DO 420 J=MYJS2,MYJE2
        DO 420 I=MYIS,MYIE
            TQ2B(I,J) = Q2(I,J,1) * ETADT(I,J,1) * F4Q2(1)
420 END DO
!
    DO 425 K=1,LM2
!
!$omp parallel do private (TQ2AK)
!
        DO 425 J=MYJS2,MYJE2
            DO 425 I=MYIS,MYIE
                TQ2AK = (Q2(I,J,K+1) - Q2(I,J,K)) * (ETADT(I,J,K) + ETADT(I,J,K+1)) * F4Q2(K+1)
                Q2ST(I,J,K) = TQ2AK + TQ2B(I,J) + Q2(I,J,K)
                TQ2B(I,J)   = TQ2AK
    425 END DO
!
!$omp parallel do private (TQ2AK)
!
    DO 440 J=MYJS2,MYJE2
        DO 440 I=MYIS,MYIE
            TQ2AK = (Q2(I,J,LM) - Q2(I,J,LM1)) * ETADT(I,J,LM1) * F4Q2(LM)
            Q2ST(I,J,LM1) = TQ2AK + TQ2B(I,J) + Q2(I,J,LM1)
            Q2ST(I,J,LM ) =                     Q2(I,J,LM)
440 END DO
!------------------------------- 
! SECOND (BACKWARD) MATSUNO STEP 
!------------------------------- 
!
!$omp parallel do 
!
    DO 450 J=MYJS2,MYJE2
        DO 450 I=MYIS,MYIE
            TQ2B(I,J) = Q2ST(I,J,1) * ETADT(I,J,1) * F4Q2(1)
450 END DO
!
    DO K=1,LM2
!
!$omp parallel do private (TQ2AK)
!
        DO J=MYJS2,MYJE2
            DO I=MYIS,MYIE
                TQ2AK = (Q2ST(I,J,K+1) - Q2ST(I,J,K)) * (ETADT(I,J,K) + ETADT(I,J,K+1)) * F4Q2(K+1)
                VAD_TEND1(I,J,K) = TQ2AK + TQ2B(I,J)
                TQ2B(I,J) = TQ2AK
            END DO
        END DO
!
    END DO
!
    DO J=MYJS2,MYJE2
        DO I=MYIS,MYIE
            TQ2AK = (Q2ST(I,J,LM) - Q2ST(I,J,LM1)) * ETADT(I,J,LM1) * F4Q2(LM)
            VAD_TEND1(I,J,LM1) = TQ2AK + TQ2B(I,J)
        END DO
    END DO
!--------------------------------------------------------------------- 
! IF THE CFL CRITERION IS VIOLATED THEN VERTICALLY SMOOTH THE TENDENCY
!--------------------------------------------------------------------- 
!
!$omp parallel do private (CFL, LBOT_CFL, LSTART, LSTOP, LTOP_CFL, VAD_TEND1, VAD_TNDX1)
!
    DO J=MYJS2,MYJE2
        DO I=MYIS,MYIE
!
            IF (LTOP_CFL_T(I,J) > 0) THEN
                LSTART =     LTOP_CFL_T(I,J)
                LSTOP  = MIN(LBOT_CFL_T(I,J),LM-2)
!            
                DO K=LSTART,LSTOP
                    VAD_TNDX1(K) = (VAD_TEND1(I,J,K-1)  + VAD_TEND1(I,J,K+1) + 2.                 &
    &                            *  VAD_TEND1(I,J,K  )) * 0.25
                END DO
                DO K=LSTART,LSTOP
                    VAD_TEND1(I,J,K) = VAD_TNDX1(K)
                END DO
            END IF
!        
        END DO
    END DO
!
    DO 470 K=1,LM2
!
!$omp parallel do 
!
        DO 470 J=MYJS2,MYJE2
            DO 470 I=MYIS,MYIE
                Q2(I,J,K) = VAD_TEND1(I,J,K) + Q2(I,J,K)
                Q2(I,J,K) =  AMAX1(Q2(I,J,K), EPSQ2)
    470 END DO
!
!$omp parallel do 
!
    DO 480 J=MYJS2,MYJE2
        DO 480 I=MYIS,MYIE
            Q2(I,J,LM1) = VAD_TEND1(I,J,LM1) + Q2(I,J,LM1)
            Q2(I,J,LM1) =  AMAX1(Q2(I,J,LM1), EPSQ2)
480 END DO
!------------------------------------------- 
! DEFINITION OF VARIABLES NEEDED AT V POINTS 
!------------------------------------------- 
!
!$omp parallel do 
!
    DO 500 K=1,LM1
        DO 500 J=MYJS_P1,MYJE_P1
            DO 500 I=MYIS_P1,MYIE_P1
                ETADT(I,J,K) = ETADT(I,J,K) * PDSL(I,J) * HBM2(I,J)
    500 END DO
!
!$omp parallel do 
!
    DO 510 J=MYJS2,MYJE2
        DO 510 I=MYIS,MYIE
            RPDX(I,J) = 1. / (PDSL(I+IVW(J),J  ) + PDSL(I+IVE(J),J  ))
            RPDY(I,J) = 1. / (PDSL(I       ,J-1) + PDSL(I       ,J+1))
510 END DO
!-------------------------------------------------------------------------------------------------- 
! PIECEWISE LINEAR UPSTREAM VERTICAL ADVECTION OF U AND V  
! INTIALIZE U AND V AT THE BOTTOM INTERFACE AND THE SLOPE ADJUSTMENT MASK (SAMU=1 AND SAMV=1 FOR SA
! PERMITTED, 0 FOR NOT PERMITTED)
!-------------------------------------------------------------------------------------------------- 
!
!$omp parallel do 
!
    DO 999 K=1,LM
        DO 999 J=MYJS2,MYJE2
            DO 999 I=MYIS,MYIE
                 UBI(I,J,K) = U(I,J,K)
                 VBI(I,J,K) = V(I,J,K)
                SAMU(I,J,K) = 1.
                SAMV(I,J,K) = 1.
999 END DO
!
    IF (NOSLAW) GOTO 998
!------------------------------------------------------
! THE SLOPE ADJUSTMENT CODE 
! NO SLOPE PERMITTED AT THE TOP AND AT THE BOTTOM LAYER
!------------------------------------------------------ 
!
!$omp parallel do 
!
    DO 997 J=MYJS2,MYJE2
        DO 997 I=MYIS,MYIE
            SAMU(I,J, 1) = 0.
            SAMU(I,J,LM) = 0.
            SAMV(I,J, 1) = 0.
            SAMV(I,J,LM) = 0.
997 END DO
!
!$omp parallel do
! 
    DO 996 K=1,LM1
        DO 996 J=MYJS2,MYJE2
            DO 996 I=MYIS,MYIE
                SAMU(I,J,K) = SAMU(I,J,K) * VTM(I,J,K+1)
                SAMV(I,J,K) = SAMV(I,J,K) * VTM(I,J,K+1)
996 END DO
!--------------------------------------------------------------------------------------------
! NOW, SEARCH FOR THE MAXIMA AND MINIMA OF U & V (AT THE FIRST PASS) AND FOR LAYERS WHICH HAD
! OVERADJUSTED (AT SUBSEQUENT PASSES) DUE TO ROUND-OFF ERRORS
!--------------------------------------------------------------------------------------------
!
!$omp parallel do private (DUBI, DVBI, DUTI, DVTI, EXTREMU, EXTREMV)
!
    DO 995 K=2,LM1
        DO 995 J=MYJS2,MYJE2
            DO 995 I=MYIS,MYIE
                DUTI(I,J) = U(I,J,K  ) - U(I,J,K-1)
                DVTI(I,J) = V(I,J,K  ) - V(I,J,K-1)
                DUBI(I,J) = U(I,J,K+1) - U(I,J,K  )
                DVBI(I,J) = V(I,J,K+1) - V(I,J,K  )
!
                EXTREMU = DUTI(I,J) * DUBI(I,J)
                EXTREMV = DVTI(I,J) * DVBI(I,J)
!
                IF (EXTREMU < 0.) SAMU(I,J,K) = 0.
                IF (EXTREMV < 0.) SAMV(I,J,K) = 0.
995 END DO
!
!$omp parallel do 
!
    DO 994 K=2,LM1
        DO 994 J=MYJS2,MYJE2
            DO 994 I=MYIS,MYIE
                ARRAYU1(I,J,K) = WFA(K-1) * (1.-SAMU(I,J,K-1)) + WFB(K-1)
                ARRAYV1(I,J,K) = WFA(K-1) * (1.-SAMV(I,J,K-1)) + WFB(K-1)
                ARRAYU2(I,J,K) = WFA(K  ) + WFB(K) * (1.-SAMU(I,J,K+1))
                ARRAYV2(I,J,K) = WFA(K  ) + WFB(K) * (1.-SAMV(I,J,K+1))
994 END DO
!
    DO 993 MSA=1,NMSAPW
!--------------------------------------------------------------------------------------------------
! CALCULATE DU & DV AT INTERFACES AND ADJUST THE SLOPES WHERE AND TO THE EXTENT PERMITTED OBSERVING
! THE MONOTONICITY CONDITION (E.G. VAN LEER, J. COMP. PHYS. 1977, 276-299, SCHEME USED HERE OF 
! MESINGER AND JOVIC, NCEP OFFICE NOTE #439)
!--------------------------------------------------------------------------------------------------
!
!$omp parallel do 
!
        DO 992 J=MYJS2,MYJE2
            DO 992 I=MYIS,MYIE
                DUBI(I,J) = 2. * U(I,J,2) - UBI(I,J,2) - UBI(I,J,1)
                DVBI(I,J) = 2. * V(I,J,2) - VBI(I,J,2) - VBI(I,J,1)
    992 END DO
!
!$omp parallel do private (ASBIKU, ASBIKV, ASTIKU, ASTIKV, DUTIK, DVTIK)
!
        DO 991 K=2,LM1
            DO 991 J=MYJS2,MYJE2
                DO 991 I=MYIS,MYIE
                    DUTIK = DUBI(I,J)
                    DVTIK = DVBI(I,J)
                    ASTIKU = ARRAYU1(I,J,K) * DUTIK
                    ASTIKV = ARRAYV1(I,J,K) * DVTIK
!
                    DUBI(I,J) = 2.0 * U(I,J,K+1) - UBI(I,J,K+1) - UBI(I,J,K)
                    DVBI(I,J) = 2.0 * V(I,J,K+1) - VBI(I,J,K+1) - VBI(I,J,K)
!
                    ASBIKU = ARRAYU2(I,J,K) * DUBI(I,J)
                    ASBIKV = ARRAYV2(I,J,K) * DVBI(I,J)
!
                    UBI(I,J,K) = UBI(I,J,K) + (ASTIKU-SIGN(1.,ASTIKU)                             &
    &                          * DIM(ABS(ASTIKU),ABS(ASBIKU))) * SAMU(I,J,K)
!
                    VBI(I,J,K) = VBI(I,J,K) + (ASTIKV-SIGN(1.,ASTIKV)                             &
    &                          * DIM(ABS(ASTIKV),ABS(ASBIKV))) * SAMV(I,J,K)
    991 END DO
993 END DO
!------------------------------------------------------------------------------------------------
! SLOPE ADJUSTMENT OF THE LAYERS ABOVE THAT NEXT TO THE SURFACE IS DONE; NOW ADJUST THE LOWERMOST
! LAYER
!------------------------------------------------------------------------------------------------
    DO 990 K=9,LM1
!
!$omp parallel do 
!
        DO 990 J=MYJS2,MYJE2
            DO 990 I=MYIS,MYIE
                IF (VTM(I,J,K+1) == 0.) UBI(I,J,K) = 2. * U(I,J,K) - UBI(I,J,K-1)
                IF (VTM(I,J,K+1) == 0.) VBI(I,J,K) = 2. * V(I,J,K) - VBI(I,J,K-1)
990 END DO
!
!$omp parallel do 
!
    DO 989 J=MYJS2,MYJE2
        DO 989 I=MYIS,MYIE
            UBI(I,J,LM) = 2. * U(I,J,LM) - UBI(I,J,LM1)
            VBI(I,J,LM) = 2. * V(I,J,LM) - VBI(I,J,LM1)
989 END DO
!--------------------------------- 
! END OF THE SLOPE ADJUSTMENT CODE 
!--------------------------------- 
998 CONTINUE
!
!$omp parallel do 
!
    DO 988 J=MYJS2,MYJE2
        DO 988 I=MYIS,MYIE
            UDEDB(I,J) = 0.
            VDEDB(I,J) = 0.
            UDEUB(I,J) = 0.
            VDEUB(I,J) = 0.
            DUDEB(I,J) = 2. * (UBI(I,J,1) - U(I,J,1)) * RDETA(1)
            DVDEB(I,J) = 2. * (VBI(I,J,1) - V(I,J,1)) * RDETA(1)
            EDBDU(I,J) = 0.
            EDBDV(I,J) = 0.
988 END DO
!
    DO 987 K=1,LM1
!
!$omp parallel do private (DQDEK, EDBFK, EDTDK, QDEDTK, QDEUTK, SEDBK)
!
        DO 987 J=MYJS2,MYJE2
            DO 987 I=MYIS,MYIE
                VMIJ   =   VTM(I,J,K) * VBM2(I,J)
!
                UDEDTK = UDEDB(I,J)
                VDEDTK = VDEDB(I,J)
                UDEUTK = UDEUB(I,J)
                VDEUTK = VDEUB(I,J)
!
                SEDBKU = SIGN(1.,(ETADT(I+IVW(J),J  ,K) + ETADT(I+IVE(J),J  ,K)))
                SEDBKV = SIGN(1.,(ETADT(I       ,J-1,K) + ETADT(I       ,J+1,K)))
!
                DUDEK = DUDEB(I,J)
                DVDEK = DVDEB(I,J)
!
                DUDEB(I,J) = 2. * (UBI(I,J,K+1) - U(I,J,K+1)) * RDETA(K+1)
                DVDEB(I,J) = 2. * (VBI(I,J,K+1) - V(I,J,K+1)) * RDETA(K+1)
!            
                EDBFKU = (ETADT(I+IVW(J),J  ,K) + ETADT(I+IVE(J),J  ,K)) * RPDX(I,J) * F4D
                EDBFKV = (ETADT(I       ,J-1,K) + ETADT(I       ,J+1,K)) * RPDY(I,J) * F4D
!            
                UDEDB(I,J) = (1.+SEDBKU) * (UBI(I,J,K) + DUDEK * EDBFKU) * (-EDBFKU)
                VDEDB(I,J) = (1.+SEDBKV) * (VBI(I,J,K) + DVDEK * EDBFKV) * (-EDBFKV)
!
                UDEUB(I,J) = (1.-SEDBKU) * (U(I,J,K+1) +     U(I,J,K+1)                           &
    &                      -              UBI(I,J,K+1) + DUDEB(I,J)                               &
    &                      *              EDBFKU)      * EDBFKU
!
                VDEUB(I,J) = (1.-SEDBKV) * (V(I,J,K+1) +     V(I,J,K+1)                           &
    &                      -              VBI(I,J,K+1) + DVDEB(I,J)                               &
    &                      *              EDBFKV)      * EDBFKV
                EDTDKU = EDBDU(I,J)
                EDTDKV = EDBDV(I,J)
!
                EDBDU(I,J) = (ETADT(I+IVW(J),J  ,K) + ETADT(I+IVE(J),J  ,K)) * RPDX(I,J) * (-F4Q)
                EDBDV(I,J) = (ETADT(I       ,J-1,K) + ETADT(I       ,J+1,K)) * RPDY(I,J) * (-F4Q)
!
                VAD_TEND1(I,J,K) = (UDEDTK - UDEUTK - UDEDB(I,J) + UDEUB(I,J) + U(I,J,K)          &
    &                            * (EDBDU(I,J) - EDTDKU)) * RDETA(K) * VMIJ
!
                VAD_TEND2(I,J,K) = (VDEDTK - VDEUTK - VDEDB(I,J) + VDEUB(I,J) + V(I,J,K)          &
    &                            * (EDBDV(I,J) - EDTDKV)) * RDETA(K) * VMIJ
987 END DO
!
!$omp parallel do 
!
    DO 986 J=MYJS2,MYJE2
        DO 986 I=MYIS,MYIE
            VMIJ = VTM(I,J,LM) * VBM2(I,J)
            VAD_TEND1(I,J,LM) = (UDEDB(I,J) - UDEUB(I,J) + U(I,J,LM) * (-EDBDU(I,J)))             &
    &                         * RDETA(LM) * VMIJ
            VAD_TEND2(I,J,LM) = (VDEDB(I,J) - VDEUB(I,J) + V(I,J,LM) * (-EDBDV(I,J)))             &
    &                         * RDETA(LM) * VMIJ
986 END DO
!-----------------------------------------------------------------------
! IF THE CFL CRITERION IS VIOLATED THEN VERTICALLY SMOOTH THE TENDENCIES
!-----------------------------------------------------------------------
!
!$omp parallel do private (LSTART, LSTOP, VAD_TNDX1, VAD_TNDX2)
!
    DO J=MYJS2,MYJE2
        DO I=MYIS,MYIE
!----------         
! COMPONENT
!----------        
            IF (LTOP_CFL_U(I,J) > 0) THEN
                LSTART =     LTOP_CFL_U(I,J)
                LSTOP  = MIN(LBOT_CFL_U(I,J), LM-1)
!
                DO K=LSTART,LSTOP
                    VAD_TNDX1(K) = (VAD_TEND1(I,J,K-1)  + VAD_TEND1(I,J,K+1) + 2.                 &
    &                            *  VAD_TEND1(I,J,K  )) * 0.25
                END DO
!
                DO K=LSTART,LSTOP
                    VAD_TEND1(I,J,K) = VAD_TNDX1(K)
                END DO
!
            END IF
!------------        
! V COMPONENT
!------------         
            IF (LTOP_CFL_V(I,J) > 0) THEN
                LSTART =     LTOP_CFL_V(I,J)
                LSTOP  = MIN(LBOT_CFL_V(I,J),LM-1)
!
                DO K=LSTART,LSTOP
                    VAD_TNDX2(K) = (VAD_TEND2(I,J,K-1)  + VAD_TEND2(I,J,K+1) + 2.                 &
    &                            *  VAD_TEND2(I,J,K  )) * 0.25
                END DO
                DO K=LSTART,LSTOP
                    VAD_TEND2(I,J,K) = VAD_TNDX2(K)
                END DO
!
            END IF
!        
        END DO
    END DO
!
    DO 580 K=1,LM
!
!$omp parallel do
! 
        DO 580 J=MYJS2,MYJE2
            DO 580 I=MYIS,MYIE
                U(I,J,K) = U(I,J,K) + VAD_TEND1(I,J,K)
                V(I,J,K) = V(I,J,K) + VAD_TEND2(I,J,K)
    580 END DO
!
    IF (.NOT. HYDRO) THEN
!
!$omp parallel do 
!
        DO K=1,LM1
            DO J=MYJS_P1,MYJE_P1
                DO I=MYIS_P1,MYIE_P1
                    ETADT(I,J,K) = ETADT(I,J,K) / PDSL(I,J)
                END DO
            END DO
        END DO
!
    END IF
!
    RETURN
!
    END SUBROUTINE VTADV
