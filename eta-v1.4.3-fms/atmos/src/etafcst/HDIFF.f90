    SUBROUTINE HDIFF
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE HDIFF
!>
!> SUBPROGRAM: HDIFF - HORIZONTAL DIFFUSION
!> PROGRAMMER: JANJIC
!> ORG: W/NP22
!> DATE: 93-11-17
!>
!> ABSTRACT:
!> HDIFF CALCULATES THE CONTRIBUTION OF THE HORIZONTAL DIFFUSION TO THE TENDENCIES OF TEMPERATURE,
!> SPECIFIC HUMIDITY, WIND COMPONENTS AND TURBULENT KINETIC ENERGY AND THEN UPDATES THOSE VARIABLES.
!> A SECOND-ORDER NONLINEAR SCHEME SIMILAR TO SMAGORINSKYS IS USED WHERE THE DIFFUSION COEFFICIENT 
!> IS A FUNCTION OF THE DEFORMATION FIELD AND OF THE TURBULENT KINETIC ENERGY.
!>
!> PROGRAM HISTORY LOG:
!> 87-06-??  JANJIC     - ORIGINATOR
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 96-03-28  BLACK      - ADDED EXTERNAL EDGE
!> 98-10-30  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
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
!> USE MODULES: CLDWTR
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
!>              PARM_TBL
!>              PHYS
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>
!> DRIVER     : EBU
!>              NEWFLT
!>
!> CALLS      : ZERO2
!>--------------------------------------------------------------------------------------------------
    USE CLDWTR
    USE CTLBLK
    USE DYNAM
    USE F77KINDS
    USE GLB_TABLE
    USE INDX
    USE LOOPS    , ONLY : JAM
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE PARMETA
    USE PARM_TBL
    USE PHYS
    USE PVRBLS
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND), PARAMETER :: DEFC  =  8.0
    REAL   (KIND=R4KIND), PARAMETER :: DEFM  = 32.0
    REAL   (KIND=R4KIND), PARAMETER :: SCQ2  = 50.0
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ2 =  0.12
    REAL   (KIND=R4KIND), PARAMETER :: FCDIF =  1.0
    REAL   (KIND=R4KIND), PARAMETER :: RFCP  =   .25 / 1004.6
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM  = IM * JM - JM / 2
    INTEGER(KIND=I4KIND), PARAMETER :: KSMUD = 1
    INTEGER(KIND=I4KIND), PARAMETER :: JAMD  = (JAM * 2 - 10) * 3
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & SECOND  , HEAT    , STTDF
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & Q2L     , UT      ,                                                                         &
    & HKNE    , HKSE    ,                                                                         &
    & VKNE    , VKSE    ,                                                                         &
    & HMASK   , HMSKL
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & TNE     , TSE     ,                                                                         &
    & QNE     , QSE     ,                                                                         &
!CHOU
    & CWMNE   , CWMSE   ,                                                                         & 
    & Q2NE    , Q2SE    ,                                                                         &
    & UNE     , USE     ,                                                                         &
    & VNE     , VSE     ,                                                                         &
!Chou    & TDIF    , QDIF    ,                                                                         &
    & TDIF    , QDIF    ,  CWMDIF    ,                                                            &
    & UDIF    , VDIF    ,                                                                         &
    & Q2DIF   ,                                                                                   &
    & DEF     , CKE  
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & LUL     , I       , J       , K       , KS
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & DH1     , DH2     , DH3     , DH4     , DHM     , DEFTK   , DEFSK   , UTK     , VTK      
!-------------------------------------------------------------------
! DIFFUSING Q2 AT GROUND LEVEL DOESNT MATTER, USTAR2 IS RECALCULATED
!-------------------------------------------------------------------
    SECOND = .TRUE. 
    HEAT   = .FALSE. 
    LUL    = UL(1)
!
    DO J=MYJS_P1,MYJE_P2
        DO I=MYIS_P2,MYIE_P2
            HMASK(I,J) = 1.
            HMSKL(I,J) = 1.
        END DO
    END DO
!
    IF (SIGMA) THEN
        DO 100 J=MYJS2_P1,MYJE2_P1
!        
            DO I=MYIS1_P1,MYIE1_P1
                DH1 = ABS(FIS(I+IHW(J), J-1) - FIS(I,J))
                DH2 = ABS(FIS(I+IHE(J), J-1) - FIS(I,J))
                DH3 = ABS(FIS(I+IHW(J), J+1) - FIS(I,J))
                DH4 = ABS(FIS(I+IHE(J), J+1) - FIS(I,J))
!            
                DHM = AMAX1(DH1, DH2, DH3, DH4) / DY
!            
                IF (DHM > 0.100) THEN
                    HMASK(I,J) = 0.
                    HMSKL(I,J) = 0.
                END IF 
            END DO
!        
    100 END DO
    END IF
!
    DO 600 KS=1,KSMUD
!------------------------------- 
! MAIN VERTICAL INTEGRATION LOOP 
!------------------------------- 
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!$omp private (CKE      , DEF     , DEFSK   , DEFTK   , HKNE    , HKSE    , HMSKL   , Q2DIF   ,   &
!$omp          Q2L      , Q2NE    , Q2SE    , QDIF    , QNE     , QSE     , TDIF    , TNE     ,   &
!$omp          TSE      , UDIF    , UNE     , USE     , UTK     , VDIF    , VKNE    , VKSE    ,   &
!$omp          VNE      , VSE     , VTK     )
!
        DO 500 K=1,LM
!
            CALL ZERO2(DEF)
            CALL ZERO2(Q2NE)
            CALL ZERO2(Q2SE)
            CALL ZERO2(QNE)
            CALL ZERO2(CWMNE)  !CHOU
            CALL ZERO2(CWMSE)  !CHOU 
            CALL ZERO2(QSE)
            CALL ZERO2(TNE)
            CALL ZERO2(TSE)
            CALL ZERO2(UNE)
            CALL ZERO2(USE)
            CALL ZERO2(VSE)
            CALL ZERO2(VNE)
            CALL ZERO2(VSE)
            CALL ZERO2(TDIF)
            CALL ZERO2(QDIF)
            CALL ZERO2(CWMDIF)   !chou
            CALL ZERO2(UDIF)
            CALL ZERO2(VDIF)
            CALL ZERO2(Q2DIF)
!
            DO 210 J=MYJS_P1,MYJE_P1
                DO 210 I=MYIS_P1,MYIE_P1
                    Q2L(I,J) = AMAX1(Q2(I,J,K),EPSQ2)
        210 END DO
!------------- 
! DEFORMATIONS 
!------------- 
            DO 220 J=MYJS1_P1,MYJE1_P1
                DO 220 I=MYIS_P1,MYIE1_P1
!                
                    IF (K < LUL) THEN
                        HMSKL(I,J) = HMASK(I,J)
                    ELSE
                        HMSKL(I,J) = HMASK(I,J)
                    END IF
!                
                    DEFTK = U(I+IHE(J),J  ,K) - U(I+IHW(J),J  ,K)                                 &
    &                     - V(I       ,J+1,K) + V(I       ,J-1,K)
!
                    DEFSK = U(I       ,J+1,K) - U(I       ,J-1,K)                                 &
    &                     + V(I+IHE(J),J  ,K) - V(I+IHW(J),J  ,K)
!
                    DEF(I,J) = DEFTK * DEFTK + DEFSK * DEFSK
                    DEF(I,J) = SQRT(DEF(I,J) + DEF(I,J)) * HBM2(I,J)
                    DEF(I,J) = AMAX1(DEF(I,J),DEFC)
                    DEF(I,J) = DEF(I,J) * HMSKL(I,J)
        220 END DO
!-------------------------------- 
! T, Q, Q2 DIAGONAL CONTRIBUTIONS 
!-------------------------------- 
            DO 250 J=MYJS_P1,MYJE1_P1
                DO 250 I=MYIS_P1,MYIE1_P1
                    HKNE(I,J) =  (DEF(I,J)   +   DEF(I+IHE(J),J+1))                               &
    &                         *   HTM(I,J,K) *   HTM(I+IHE(J),J+1,K)                              &
    &                         * HMSKL(I,J)   * HMSKL(I+IHE(J),J+1)
!
                    TNE (I,J) = (T (I+IHE(J),J+1,K) - T (I,J,K)) * HKNE(I,J)
                    QNE (I,J) = (Q (I+IHE(J),J+1,K) - Q (I,J,K)) * HKNE(I,J)
                    Q2NE(I,J) = (Q2(I+IHE(J),J+1,K) - Q2(I,J,K)) * HKNE(I,J)
                    CWMNE(I,J)= (CWM(I+IHE(J),J+1,K)- CWM(I,J,K))* HKNE(I,J)   !CHOU 20210206
        250 END DO
!        
            DO 260 J=MYJS1_P1,MYJE_P1
                DO 260 I=MYIS_P1,MYIE1_P1
                    HKSE(I,J) =  (DEF(I+IHE(J),J-1)   +   DEF(I,J))                               &
    &                         *   HTM(I+IHE(J),J-1,K) *   HTM(I,J,K)                              &
    &                         * HMSKL(I+IHE(J),J-1)   * HMSKL(I,J)
!
                    TSE (I,J) = (T (I+IHE(J),J-1,K) - T (I,J,K)) * HKSE(I,J)
                    QSE (I,J) = (Q (I+IHE(J),J-1,K) - Q (I,J,K)) * HKSE(I,J)
                    Q2SE(I,J) = (Q2(I+IHE(J),J-1,K) - Q2(I,J,K)) * HKSE(I,J)
                    CWMSE(I,J)= (CWM(I+IHE(J),J-1,K)- CWM(I,J,K))* HKSE(I,J)   !CHOU 20210206
        260 END DO
!
            DO 270 J=MYJS1,MYJE1
                DO 270 I=MYIS1,MYIE
                    TDIF (I,J) = (TNE (I,J) - TNE (I+IHW(J),J-1)                                  &
    &                          +  TSE (I,J) - TSE (I+IHW(J),J+1))                                 &
    &                          *  HDAC(I,J)
!
                    QDIF (I,J) = (QNE (I,J) - QNE (I+IHW(J),J-1)                                  &
    &                          +  QSE (I,J) - QSE (I+IHW(J),J+1))                                 &
    &                          *  HDAC(I,J) * FCDIF
!
                    Q2DIF(I,J) = (Q2NE(I,J) - Q2NE(I+IHW(J),J-1)                                  &
    &                          +  Q2SE(I,J) - Q2SE(I+IHW(J),J+1))                                 &
    &                          *  HDAC(I,J)
!
                    CWMDIF(I,J)= (CWMNE(I,J)- CWMNE(I+IHW(J),J-1)                                 &
    &                          +  CWMSE(I,J)- CWMSE(I+IHW(J),J+1))                                &
    &                          *  HDAC(I,J) * FCDIF
        270 END DO
!--------------------- 
! 2-ND ORDER DIFFUSION
 !-------------------- 
            IF (SECOND) THEN
                DO 280 J=MYJS2,MYJE2
                    DO 280 I=MYIS1,MYIE1
                        T(I,J,K) = T(I,J,K) + TDIF(I,J)
                        Q(I,J,K) = Q(I,J,K) + QDIF(I,J)
                        CWM(I,J,K)= CWM(I,J,K) + CWMDIF(I,J)    !CHOU ADD cloud (20210206)
            280 END DO
!
                IF (K /= LM) THEN
                    DO 290 J=MYJS2,MYJE2
                        DO 290 I=MYIS1,MYIE1
!DCR                            Q2(I,J,K) = Q2(I,J,K) + Q2DIF(I,J) * HTM(I,J,K+1)
                            Q2(I,J,K) = Q2(I,J,K) + Q2DIF(I,J) * HTM(I,J,K)
                            Q2(I,J,K) = AMAX1(Q2(I,J,K) + Q2DIF(I,J) * HTM(I,J,K),EPSQ2)
                290 END DO
                END IF
!            
                GOTO 360
            END IF
!---------------------------------- 
! 4-TH ORDER DIAGONAL CONTRIBUTIONS
!---------------------------------- 
            DO 310 J=MYJS,MYJE1
                DO 310 I=MYIS,MYIE1
                    TNE (I,J) = (TDIF (I+IHE(J),J+1) - TDIF (I,J)) * HKNE(I,J)
                    QNE (I,J) = (QDIF (I+IHE(J),J+1) - QDIF (I,J)) * HKNE(I,J)
                    Q2NE(I,J) = (Q2DIF(I+IHE(J),J+1) - Q2DIF(I,J)) * HKNE(I,J)
                    CWMNE(I,J)= (CWMDIF(I+IHE(J),J+1)- CWMDIF(I,J))* HKNE(I,J)   !CHOU 20210206
        310 END DO
!        
            DO 320 J=MYJS1,MYJE
                DO 320 I=MYIS,MYIE1
                    TSE (I,J) = (TDIF (I+IHE(J),J-1) - TDIF (I,J)) * HKSE(I,J)
                    QSE (I,J) = (QDIF (I+IHE(J),J-1) - QDIF (I,J)) * HKSE(I,J)
                    Q2SE(I,J) = (Q2DIF(I+IHE(J),J-1) - Q2DIF(I,J)) * HKSE(I,J)
                    CWMSE(I,J)= (CWMDIF(I+IHE(J),J-1)- CWMDIF(I,J))* HKSE(I,J)    !CHOU 20210206
        320 END DO
!
            DO 330 J=MYJS2,MYJE2
                DO 330 I=MYIS1,MYIE1
                    T(I,J,K) = T(I,J,K) - (TNE (I,J) - TNE(I+IHW(J),J-1)                          &
    &                        +             TSE (I,J) - TSE(I+IHW(J),J+1))                         &
    &                        *             HDAC(I,J)
!
                    Q(I,J,K) = Q(I,J,K) - (QNE (I,J) - QNE(I+IHW(J),J-1)                          &
    &                        +             QSE (I,J) - QSE (I+IHW(J),J+1))                        &
    &                        *             HDAC(I,J) * FCDIF
!CHOU 20210206
                    CWM(I,J,K) = CWM(I,J,K) - (CWMNE(I,J)-CWMNE(I+IHW(J),J-1)                     &
    &                         +               CWMSE(I,J)-CWMSE(I+IHW(J),J+1))                     &
    &                         *               HDAC(I,J) * FCDIF
       330 END DO
!
            IF (K /= LM) THEN
                DO 340 J=MYJS2,MYJE2
                    DO 340 I=MYIS1,MYIE1
                        Q2(I,J,K) = Q2(I,J,K) - (Q2NE(I,J) - Q2NE(I+IHW(J),J-1)                   &
    &                             +              Q2SE(I,J) - Q2SE(I+IHW(J),J+1))                  &
    &                             * HDAC(I,J) * HTM(I,J,K)
!dcr    &                             * HDAC(I,J) * HTM(I,J,K+1)
            340 END DO
            END IF
!----------------------------- 
! U, V, DIAGONAL CONTRIBUTIONS
!----------------------------- 
        360 DO 410 J=MYJS_P1,MYJE1_P1
            DO 410 I=MYIS_P1,MYIE1_P1
                VKNE(I,J) =  (DEF(I+IVE(J),J) +   DEF(I,J+1))                                     &
    &                     *   VTM(I,J,K)      *   VTM(I+IVE(J),J+1,K)                             &
    &                     * HMASK(I+IVE(J),J) * HMASK(I,J+1)
!
                UNE(I,J)  = (U(I+IVE(J),J+1,K) - U(I,J,K)) * VKNE(I,J)
                VNE(I,J)  = (V(I+IVE(J),J+1,K) - V(I,J,K)) * VKNE(I,J)
        410 END DO
!        
            DO 420 J=MYJS1_P1,MYJE_P1
                DO 420 I=MYIS_P1,MYIE1_P1
                    VKSE(I,J) =  (DEF(I,J-1)          +   DEF(I+IVE(J),J))                        &
    &                         *   VTM(I+IVE(J),J-1,K) *   VTM(I,J,K)                              &
    &                         *  HMASK(I,J-1)         * HMASK(I+IVE(J),J)
!
                    USE(I,J) = (U(I+IVE(J),J-1,K) - U(I,J,K)) * VKSE(I,J)
                    VSE(I,J) = (V(I+IVE(J),J-1,K) - V(I,J,K)) * VKSE(I,J)
        420 END DO
!
            DO 430 J=MYJS1,MYJE1
                DO 430 I=MYIS,MYIE1
                    UDIF(I,J) =  (UNE(I,J) - UNE(I+IVW(J),J-1)                                    &
    &                         +   USE(I,J) - USE(I+IVW(J),J+1))                                   &
    &                         * HDACV(I,J)
!
                    VDIF(I,J) =  (VNE(I,J) - VNE(I+IVW(J),J-1)                                    &
    &                         +   VSE(I,J) - VSE(I+IVW(J),J+1))                                   &
    &                         * HDACV(I,J)
        430 END DO
!--------------------- 
! 2-ND ORDER DIFFUSION
!--------------------- 
            IF (SECOND) THEN
                DO 440 J=MYJS2,MYJE2
                    DO 440 I=MYIS1,MYIE1
                        U(I,J,K) = U(I,J,K) + UDIF(I,J)
                        V(I,J,K) = V(I,J,K) + VDIF(I,J)
            440 END DO
            ELSE
!---------------------------------- 
! 4-TH ORDER DIAGONAL CONTRIBUTIONS 
!---------------------------------- 
                DO 450 J=MYJS,MYJE1
                    DO 450 I=MYIS,MYIE1
                        UNE(I,J) = (UDIF(I+IVE(J),J+1) - UDIF(I,J)) * VKNE(I,J)
                        VNE(I,J) = (VDIF(I+IVE(J),J+1) - VDIF(I,J)) * VKNE(I,J)
            450 END DO
!            
                DO 460 J=MYJS1,MYJE
                    DO 460 I=MYIS,MYIE1
                        USE(I,J) = (UDIF(I+IVE(J),J-1) - UDIF(I,J)) * VKSE(I,J)
                        VSE(I,J) = (VDIF(I+IVE(J),J-1) - VDIF(I,J)) * VKSE(I,J)
            460 END DO
!
                DO 470 J=MYJS2,MYJE2
                    DO 470 I=MYIS1,MYIE1
                        UTK = U(I,J,K)
                        VTK = V(I,J,K)
                        U(I,J,K) = U(I,J,K) - (UNE(I,J) - UNE(I+IVW(J),J-1)                       &
    &                            +             USE(I,J) - USE(I+IVW(J),J+1))                      &
    &                            *           HDACV(I,J)
!
                        V(I,J,K) = V(I,J,K) - (VNE(I,J) - VNE(I+IVW(J),J-1)                       &
    &                            +             VSE(I,J) - VSE(I+IVW(J),J+1))                      &
    &                            *           HDACV(I,J)
!
                        CKE(I,J) = 0.5 * (U(I,J,K) * U(I,J,K) - UTK * UTK                         &
    &                            +        V(I,J,K) * V(I,J,K) - VTK * VTK)
            470 END DO
!
                IF (HEAT) THEN
                    DO 480 J=MYJS2,MYJE2
                        DO 480 I=MYIS1,MYIE1
                            T(I,J,K) = -RFCP * (CKE(I+IHE(J),J) + CKE(I,J+1)                      &
    &                                +          CKE(I+IHW(J),J) + CKE(I,J-1))                     &
    &                                * HBM2(I,J) + T(I,J,K)
                480 END DO
                END IF
!            
            END IF
!
    500 END DO
!
600 END DO
!
    RETURN
!
    END SUBROUTINE HDIFF
