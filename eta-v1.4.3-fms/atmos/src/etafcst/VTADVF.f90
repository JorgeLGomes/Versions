    SUBROUTINE VTADVF
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE VTADVF
!> 
!> SUBROUTINE: VTADVF - VERTICAL ADVECTION
!> PROGRAMMER: JANJIC
!> ORG: W/NP22
!> DATE: 93-11-17
!>
!> ABSTRACT:
!> VTADVF CALCULATES THE CONTRIBUTION OF THE VERTICAL ADVECTION TO THE TENDENCIES OF TEMPERATURE, 
!> WIND COMPONENTS, AND TURBULENT KINETIC ENERGY AND THEN UPDATES THOSE VARIABLES. FOR ALL THESE 
!> VARIABLES A SIMPLE CENTERED DIFFERENCE SCHEME IN SPACE IS USED IN CONJUNCTION WITH THE PURE 
!> EULER-BACKWARD TIME SCHEME.
!>
!> PROGRAM HISTORY LOG:
!> 87-06-??  JANJIC     - ORIGINATOR
!> 90-??-??  MESINGER   - INSERTED PIECEWISE LINEAR SCHEME FOR SPECIFIC HUMIDITY
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 95-11-20  ABELES     - PARALLEL OPTIMIZATION
!> 96-03-29  BLACK      - ADDED EXTERNAL EDGE; REMOVED SCRCH COMMON
!> 98-11-24  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
!> 18-01-15  LUCCI      - INDENTATION, UNIFORMIZATION CODE AND FREE FORM WITH OPENMP
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
!>              PARMETA
!>              PVRBLS
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>
!> DRIVER     : DIGFLT
!>              NEWFLT
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
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
!
    INTEGER(KIND=I4KIND), PARAMETER :: KSMUD   =  0
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM    = IM * JM - JM / 2
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & NOSLA  
!
    REAL   (KIND=R4KIND), DIMENSION(LM1)                                                        ::&
    & WFA     , WFB
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & ETADTL  ,                                                                                   &
    & TTA     , TQ2A    ,                                                                         &
    & TUA     , TVA     ,                                                                         &
    & TTB     , TQ2B    ,                                                                         &
    & TUB     , TVB     ,                                                                         &
    & VM      ,                                                                                   &
    & RPDX    , RPDY    ,                                                                         &
    & FNE     , FSE   
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & TSTL    ,                                                                                   &
    & USTL    ,                                                                                   &
    & VSTL    ,                                                                                   &
    & Q2ST 
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K       , KS      , NMSAP   , NSMUD   
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & TTAK    , TQ2AK   , VMK     , TUAK    , TVAK 
!------------------------------------------ 
! DEFINE ADDED UPSTREAM ADVECTION CONSTANTS 
!------------------------------------------ 
    DO 25 K=1,LM1
        WFA(K) = DETA(K  ) / (DETA(K) + DETA(K+1))
        WFB(K) = DETA(K+1) / (DETA(K) + DETA(K+1))
 25 END DO
!------------------------------------------- 
! NO MOISTURE SLOPE ADJUSTMENT IF NOT WANTED
!------------------------------------------- 
    NOSLA = .FALSE.
!----------------------------------------------------- 
! IF FALSE, NUMBER OF MOISTURE SLOPE ADJUSTMENT PASSES
!----------------------------------------------------- 
    NMSAP = 3
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
!-------------------------------------------------------------------------- 
! THE FNE, FSE, ETADTL, AND ETADT ARRAYS ARE ON OR ASSOCIATED WITH H POINTS
!-------------------------------------------------------------------------- 
            DO 90 KS=1,NSMUD
                DO 80 J=MYJS_P3,MYJE1_P3
                    DO 80 I=MYIS_P3,MYIE_P3
                        FNE(I,J) = (ETADT(I+IHE(J),J+1,K  ) - ETADT(I     ,J  ,K  ))              &
    &                            *    HTM(I       ,J  ,K+1) * HTM(I+IHE(J),J+1,K+1)
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
                        ETADTL(I,J) = (FNE(I,J) - FNE(I + IHW(J), J-1)                            &
    &                               +  FSE(I,J) - FSE(I + IHW(J), J+1))                           &
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
!----------------------------------  
! VERTICAL (MATSUNO) ADVECTION OF T
!---------------------------------- 
!
!$omp parallel do
!
    DO 100 J=MYJS,MYJE
        DO 100 I=MYIS,MYIE
            TTB(I,J) = 0.
100 END DO
!
    DO 110 K=1,LM1
!
!$omp parallel do private (TTAK)
!
        DO 110 J=MYJS2,MYJE2
            DO 110 I=MYIS,MYIE
                TTAK        = (T(I,J,K+1) -   T(I,J,K)) * ETADT(I,J,K) * F4D
                TSTL(I,J,K) = (TTAK       + TTB(I,J))   * RDETA(K)     + T(I,J,K)
                TTB(I,J)    =  TTAK
110 END DO
!
!$omp parallel do
!
    DO 120 J=MYJS2,MYJE2
        DO 120 I=MYIS,MYIE
            TSTL(I,J,LM) = T(I,J,LM) + TTB(I,J) * RDETA(LM)
120 END DO
!------------------------------- 
! SECOND (BACKWARD) MATSUNO STEP
!-------------------------------  
!
!$omp parallel do
!
    DO 125 J=MYJS,MYJE
        DO 125 I=MYIS,MYIE
            TTB(I,J) = 0.
125 END DO
!
    DO 140 K=1,LM1
!
!$omp parallel do private (TTAK)
!
        DO 140 J=MYJS2,MYJE2
            DO 140 I=MYIS,MYIE
                TTAK     = (TSTL(I,J,K+1)    - TSTL(I,J,K)) * ETADT(I,J,K) * F4D
                T(I,J,K) = (TTAK             +  TTB(I,J))   * RDETA(K)     + T(I,J,K)
                TTB(I,J) =  TTAK
    140 END DO
!
!$omp parallel do 
!
    DO 150 J=MYJS2,MYJE2
        DO 150 I=MYIS,MYIE
            T(I,J,LM) = T(I,J,LM) + TTB(I,J) * RDETA(LM)
150 END DO
!-----------------------------------
! VERTICAL (MATSUNO) ADVECTION OF Q2 
!-----------------------------------
!
!$omp parallel do 
!
    DO 400 J=MYJS2,MYJE2
        DO 400 I=MYIS,MYIE
            TQ2B(I,J) = Q2(I,J,1) * ETADT(I,J,1) * F4Q2(1)
400 END DO
!
    DO 425 K=1,LM2
!
!$omp parallel do private(TQ2AK)
!
        DO 425 J=MYJS2,MYJE2
            DO 425 I=MYIS,MYIE
                TQ2AK       = (   Q2(I,J,K+1) -    Q2(I,J,K  ))                                   &
    &                       * (ETADT(I,J,K  ) + ETADT(I,J,K+1))                                   &
    &                       *       F4Q2(K+1)
!
                Q2ST(I,J,K) = TQ2AK + TQ2B(I,J) + Q2(I,J,K)
                TQ2B(I,J)   = TQ2AK
425 END DO
!
!$omp parallel do private(TQ2AK)
!
    DO 440 J=MYJS2,MYJE2
        DO 440 I=MYIS,MYIE
            TQ2AK         = (Q2(I,J,LM) -   Q2(I,J,LM1)) * ETADT(I,J,LM1) * F4Q2(LM)
            Q2ST(I,J,LM1) = TQ2AK       + TQ2B(I,J)      +    Q2(I,J,LM1)
            Q2ST(I,J,LM ) =  Q2(I,J,LM)
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
    DO 470 K=1,LM2
!
!$omp parallel do private(TQ2AK)
!
        DO 470 J=MYJS2,MYJE2
            DO 470 I=MYIS,MYIE
                TQ2AK       = ( Q2ST(I,J,K+1) -  Q2ST(I,J,K  ))                                   &
    &                       * (ETADT(I,J,K  ) + ETADT(I,J,K+1))                                   &
    &                       *       F4Q2(K+1)
!
                  Q2(I,J,K) = TQ2AK + TQ2B(I,J) + Q2(I,J,K)
                TQ2B(I,J)   = TQ2AK
470 END DO
!
!$omp parallel do private(TQ2AK)
!
    DO 480 J=MYJS2,MYJE2
        DO 480 I=MYIS,MYIE
            TQ2AK       =      (Q2ST(I,J,LM) - Q2ST(I,J,LM1)) * ETADT(I,J,LM1) * F4Q2(LM)
            Q2(I,J,LM1) = TQ2AK+TQ2B(I,J)    +   Q2(I,J,LM1)
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
            RPDX(I,J) =1. / (PDSL(I+IVW(J),J  ) + PDSL(I+IVE(J),J  ))
            RPDY(I,J) =1. / (PDSL(I       ,J-1)      + PDSL(I  ,J+1))
510 END DO
!---------------------------------------- 
! VERTICAL (MATSUNO) ADVECTION OF U AND V
!----------------------------------------
!
!$omp parallel do
!
    DO 520 J=MYJS,MYJE
        DO 520 I=MYIS,MYIE
            TUB(I,J) = 0.
            TVB(I,J) = 0.
520 END DO
!
    DO 540 K=1,LM1
!
!$omp parallel do private (TUAK, TVAK, VMK)
!
        DO 540 J=MYJS2,MYJE2
            DO 540 I=MYIS,MYIE
                VMK         =    VTM(I       ,J,K+1) *  VBM2(I       ,J)
!
                TUAK        = (ETADT(I+IVW(J),J,K  ) + ETADT(I+IVE(J),J,K))                       &
    &                       *     (U(I       ,J,K+1) -     U(I       ,J,K))                       &
    &                       *   RPDX(I       ,J    ) * (VMK * F4D)
!
                USTL(I,J,K) = (TUAK + TUB(I,J))      * RDETA(K)             +  U(I,J,K)
                 TUB(I,J)   =  TUAK
!
                TVAK        = (ETADT(I,J-1,K)        + ETADT(I,J+1,K))                            &
    &                       * (    V(I,J,K+1)        -      V(I,J ,K))                            &
    &                       *   RPDY(I,J)            * (VMK * F4D)
!
                VSTL(I,J,K) = (TVAK + TVB(I,J))      * RDETA(K)             +  V(I,J,K)
!
                 TVB(I,J)   = TVAK
540 END DO
!
!$omp parallel do
!
    DO 550 J=MYJS2,MYJE2
        DO 550 I=MYIS,MYIE
            USTL(I,J,LM) = U(I,J,LM) + TUB(I,J) * RDETA(LM)
            VSTL(I,J,LM) = V(I,J,LM) + TVB(I,J) * RDETA(LM)
550 END DO
!------------------------------- 
! SECOND (BACKWARD) MATSUNO STEP
!------------------------------- 
!
!$omp parallel do
!
    DO 560 J=MYJS,MYJE
        DO 560 I=MYIS,MYIE
            TUB(I,J) = 0.
            TVB(I,J) = 0.
560 END DO
!
    DO 580 K=1,LM1
!
!$omp parallel do private (TUAK, TVAK, VMK)
!
        DO 580 J=MYJS2,MYJE2
            DO 580 I=MYIS,MYIE
                VMK      =    VTM(I       ,J,K+1) *  VBM2(I       ,J)
                TUAK     = (ETADT(I+IVW(J),J,K  ) + ETADT(I+IVE(J),J,K))                          &
    &                    * (USTL(I        ,J,K+1) -  USTL(I       ,J,K))                          &
    &                    *  RPDX(I,J) * (VMK * F4D)
!
                U(I,J,K) = (TUAK + TUB(I,J)) * RDETA(K) +  U(I,J,K)
              TUB(I,J)   =  TUAK
              TVAK       = (ETADT(I,J-1,K  ) + ETADT(I,J+1,K))                                  &
    &                    * ( VSTL(I,J  ,K+1) -  VSTL(I,J  ,K))                                  &
    &                    * RPDY(I,J) * (VMK * F4D)
!
                V(I,J,K) = (TVAK + TVB(I,J)) * RDETA(K) + V(I,J,K)
                TVB(I,J) =  TVAK											  
580 END DO
!
!$omp parallel do
!
    DO 590 J=MYJS2,MYJE2
        DO 590 I=MYIS,MYIE
            U(I,J,LM) = U(I,J,LM) + TUB(I,J) * RDETA(LM)
            V(I,J,LM) = V(I,J,LM) + TVB(I,J) * RDETA(LM)
590 END DO
!
    RETURN
!
    END SUBROUTINE VTADVF
