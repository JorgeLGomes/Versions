    SUBROUTINE HZADV_LM1
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE HZADV_LM1
!>
!> SUBPROGRAM: HZADV_LM1 - HORIZONTAL ADVECTION
!> PROGRAMMER: JANJIC
!> ORG: W/NP22
!> DATE: 93-10-28
!>
!> ABSTRACT:
!> HZADV CALCULATES THE CONTRIBUTION OF THE HORIZONTAL ADVECTION TO THE TENDENCIES OF TEMPERATURE,
!> WIND COMPONENTS, AND TURBULENT KINETIC ENERGY AND THEN UPDATES THOSE VARIABLES.
!> THE JANJIC ADVECTION SCHEME FOR THE ARAKAWA E GRID IS USED FOR ALL VARIABLES INSIDE THE FIFTH ROW
!> AN UPSTREAM SCHEME IS USED ON ALL VARIABLES IN THE THIRD, FOURTH, AND FIFTH OUTERMOST ROWS. 
!> A MODIFIED EULER-BACKWARD TIME SCHEME (HEUN) IS USED. UNDERGROUND WINDS MUST BE EQUAL TO ZERO 
!> SINCE THEY ARE USED EXPLICITLY WITHOUT THE VELOCITY MASK IN THE FLUX CALCULATIONS.
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
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>
!> DRIVER     : EBU
!>
!> CALLS      : ZERO2
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
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND), PARAMETER :: TLC  = 2. * 0.703972477
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: IM1  = IM - 1
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM = IM * JM - JM / 2
    INTEGER(KIND=I4KIND), PARAMETER :: JAMD = (JAM * 2 - 10) * 3
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & ITER2
! 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & HM      , VM      ,                                                                         &
    & RDPD    ,                                                                                   &
    & ADPDX   , ADPDY   ,                                                                         &
    & RDPDX   , RDPDY   ,                                                                         &
    & ADT     ,                                                                                   &
    & ADU     , ADV     ,                                                                         &
    & ADQ2M   , ADQ2L   ,                                                                         &
    & Q2MNS   , Q2LNS   ,                                                                         &
    & UDY     , VDX   
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & DPDE    ,                                                                                   &
    & TEMPA   , TEMPB   ,                                                                         &
    & TST     ,                                                                                   &
    & UST     , VST     ,                                                                         &
    & Q2M     , Q2L     ,                                                                         &
    & TEW     , TNS     ,                                                                         &
    & Q2MEW   , Q2LEW 
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & TNE     , TSE     ,                                                                         &
    & Q2MNE   , Q2MSE   ,                                                                         &
    & Q2LNE   , Q2LSE   ,                                                                         &
    & UEW     , UNS     ,                                                                         &
    & VEW     , VNS     ,                                                                         &
    & UNE     , USE     ,                                                                         &
    & VNE     , VSE     ,                                                                         &
    & FEW     , FNS     ,                                                                         &
    & FNE     , FSE   
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2 , LM)                              ::&
    & ADQ2HL
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2 , LM+1)                            ::&
    & Q2ML
!
    REAL   (KIND=R4KIND), DIMENSION(JAMD)                                                       ::&
    & ARRAY0  , ARRAY1  , ARRAY2  , ARRAY3 
!
    INTEGER(KIND=R4KIND), DIMENSION(JAMD)                                                       ::&
    & KHHAS   , IHLAS   , JHLAS   , KVHAS   , IVLAS   , JVLAS   , ISPA    , ISQA
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & UPSTRM
!
    LOGICAL(KIND=L4KIND), DIMENSION(JAM)                                                        ::&
    & LJRA
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IEND    , I       , J       , K       , JAKONE  , JA      , IHL     , IHH     , JAKTWO  ,   &
    & IVL     , IVH     , JAK     , IX      , JX      , ISP     , ISQ     , IFP     , IFQ     ,   &
    & IPQ        
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & TTA     , TTB     , PP      , QP      , F0      , F1     , F2       , F3
!
    INCLUDE "CHECKIN"
!-------------------------------------------- 
! FIGURE OUT IF WE ARE IN THE UPSTREAM REGION
!-------------------------------------------- 
    UPSTRM = .FALSE. 
!
    IF (MYPE <= INPES-1)        UPSTRM = .TRUE. 
    IF (MYPE >= NPES-INPES)     UPSTRM = .TRUE. 
    IF (MOD(MYPE,INPES) == 0)   UPSTRM = .TRUE. 
    IF (MOD(MYPE+1,INPES) == 0) UPSTRM = .TRUE. 
!
    JAKONE = 0
!
    DO 25 JA=1,JAM
        IHL = IHLA(JA)
        IHH = IHHA(JA)
        J   =  JRA(JA)
!
        LJRA(JA) = .FALSE. 
!    
        IF (J >= MY_JS_GLB-JBPAD2 .AND. J <= MY_JE_GLB+JTPAD2) THEN
            LJRA(JA) = .TRUE. 
            DO I=IHL,IHH
                IF (I >= MY_IS_GLB-ILPAD2 .AND. I <= MY_IE_GLB+IRPAD2) THEN
                    JAKONE = JAKONE + 1
!
                    KHHAS(JAKONE) = JA
                    IHLAS(JAKONE) = I
                    JHLAS(JAKONE) = J
                END IF
            END DO
        END IF
!    
 25 END DO
!
    JAKTWO = 0
!
    DO 50 JA=1,JAM
        IVL = IVLA(JA)
        IVH = IVHA(JA)
        J   =  JRA(JA)
!    
        DO 50 I=IVL,IVH
            IF (I >= MY_IS_GLB-ILPAD2 .AND. I <= MY_IE_GLB+IRPAD2                                 &
    &                                 .AND.                                                       &
    &           J >= MY_JS_GLB-JBPAD2 .AND. J <= MY_JE_GLB+JTPAD2) THEN
                JAKTWO = JAKTWO + 1
!
                KVHAS(JAKTWO) = JA
                IVLAS(JAKTWO) = I
                JVLAS(JAKTWO) = J
            END IF
 50 END DO
!
    DO 70 J=MYJS_P5,MYJE_P5
        DO 70 I=MYIS_P4,MYIE_P4
            Q2ML(I,J,1) = 0.
 70 END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
    DO 80 K=2,LM+1
        DO 80 J=MYJS_P5,MYJE_P5
            DO 80 I=MYIS_P4,MYIE_P4
                Q2ML(I,J,K) = Q2(I,J,K-1)
 80 END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!$omp private (ADPDX    , ADPDY   , ADQ     , ADQ2L   , ADQ2M   , ADT     , ADU     , ADV     ,   &
!$omp          ARRAY0   , ARRAY1  , ARRAY2  , ARRAY3  , DPDE    , F0      , F1      , F2      ,   &
!$omp          F3       , FEW     , FNE     , FNS     , FSE     , HM      , I       , IFP     ,   &
!$omp          IFQ      , IHH     , IHL     , IPQ     , ISP     , ISPA    , ISQ     , ISQA    ,   &
!$omp          ITER2    , IX      , IY      , J       , JA      , JAK     , K       , PP      ,   &
!$omp          Q2L      , Q2LEW   , Q2LNE   , Q2LNS   , Q2LSE   , Q2M     , Q2MEW   , Q2MNE   ,   &
!$omp          Q2MNS    , Q2MSE   , QEW     , QNE     , QNS     , QP      , QSE     , QST     ,   &
!$omp          RDPD     , RDPDX   , RDPDY   , TEMPA   , TEMPB   , TEW     , TNE     , TNS     ,   &
!$omp          TSE      , TST     , TTA     , TTB     , UDY     , UEW     , UNE     , UNS     ,   &
!$omp          USE      , UST     , VDX     , VEW     , VM      , VNE     , VNS     , VSE     ,   &
!$omp          VST      )
!
    DO 500 K=LM,1,-1
!
        CALL ZERO2(ADT)
        CALL ZERO2(ADU)
        CALL ZERO2(ADV)
        CALL ZERO2(ADQ2M)
        CALL ZERO2(ADQ2L)
        CALL ZERO2(DPDE)
        CALL ZERO2(FEW)
        CALL ZERO2(FNE)
        CALL ZERO2(FNS)
        CALL ZERO2(FSE)
        CALL ZERO2(Q2L)
        CALL ZERO2(Q2LEW)
        CALL ZERO2(Q2LNE)
        CALL ZERO2(Q2LSE)
        CALL ZERO2(Q2M)
        CALL ZERO2(Q2MEW)
        CALL ZERO2(Q2MNE)
        CALL ZERO2(Q2MSE)
        CALL ZERO2(RDPD)
        CALL ZERO2(TEMPA)
        CALL ZERO2(TEMPB)
        CALL ZERO2(TEW)
        CALL ZERO2(TNE)
        CALL ZERO2(TNS)
        CALL ZERO2(TSE)
        CALL ZERO2(TST)
        CALL ZERO2(UDY)
        CALL ZERO2(UEW)
        CALL ZERO2(UNE)
        CALL ZERO2(UNS)
        CALL ZERO2(USE)
        CALL ZERO2(UST)
        CALL ZERO2(VEW)
        CALL ZERO2(VNE)
        CALL ZERO2(VNS)
        CALL ZERO2(VSE)
        CALL ZERO2(VST)
        CALL ZERO2(VM)
!
        ITER2 = .FALSE. 
!
        DO J=MYJS_P4,MYJE_P4
            DO I=MYIS_P4,MYIE_P4
                Q2M(I,J) = Q2ML(I,J,K)
            END DO
        END DO
!    
        DO 110 J=MYJS_P5,MYJE_P5
            DO 110 I=MYIS_P4,MYIE_P4
                  HM(I,J) =  HTM(I,J,K) * HBM2(I,J)
                DPDE(I,J) = PDSL(I,J)   * DETA(K)
                RDPD(I,J) = 1. / DPDE(I,J)
                 UST(I,J) =    U(I,J,K  )
                 VST(I,J) =    V(I,J,K  )
                 TST(I,J) =    T(I,J,K  )
                 Q2L(I,J) = Q2ML(I,J,K+1)
    110 END DO
!
        DO 120 J=MYJS1_P4,MYJE1_P4
            DO 120 I=MYIS_P4,MYIE_P4
                   VM(I,J) = VTM(I,J,K) * VBM2(I,J)
                ADPDX(I,J) = DPDE(I+IVW(J),J  ) + DPDE(I+IVE(J),J  )
                ADPDY(I,J) = DPDE(I       ,J-1) + DPDE(I       ,J+1)
                RDPDX(I,J) = 1. / ADPDX(I,J)
                RDPDY(I,J) = 1. / ADPDY(I,J)
    120 END DO
!----------------------------------------------------------- 
! MASS FLUXES AND MASS POINTS ADVECTION COMPONENTS 
! 
! THE NS AND EW FLUXES IN THE FOLLOWING LOOP ARE ON V POINTS
!-----------------------------------------------------------
    125 DO 130 J=MYJS1_P4,MYJE1_P4
            DO 130 I=MYIS_P4,MYIE_P4
                  UDY(I,J) = UST(I,J) * DY
                  FEW(I,J) = UDY(I,J) * ADPDX(I,J)
                  TEW(I,J) = FEW(I,J) *  (TST(I+IVE(J),J  ) - TST(I+IVW(J),J  ))
                Q2MEW(I,J) = FEW(I,J) *  (Q2M(I+IVE(J),J  ) - Q2M(I+IVW(J),J  ))
                Q2LEW(I,J) = FEW(I,J) *  (Q2L(I+IVE(J),J  ) - Q2L(I+IVW(J),J  ))
                  VDX(I,J) = VST(I,J) *    DX(I       ,J  )
                  FNS(I,J) = VDX(I,J) * ADPDY(I       ,J  )
                  TNS(I,J) = FNS(I,J) *  (TST(I       ,J+1) - TST(I       ,J-1))
                Q2MNS(I,J) = FNS(I,J) *  (Q2M(I       ,J+1) - Q2M(I       ,J-1))
                Q2LNS(I,J) = FNS(I,J) *  (Q2L(I       ,J+1) - Q2L(I       ,J-1))
    130 END DO
!-------------------------------------------------------------------------------------- 
! DIAGONAL FLUXES AND DIAGONALLY AVERAGED WIND 
! 
! THE NE AND SE FLUXES ARE ON H POINTS (ACTUALLY JUST TO THE NE AND SE OF EACH H POINT)
!--------------------------------------------------------------------------------------
        DO 145 J=MYJS2_P4,MYJE2_P4
            DO 145 I=MYIS_P4,MYIE_P4
                TEMPA(I,J) = UDY(I,J) + VDX(I,J)
                TEMPB(I,J) = UDY(I,J) - VDX(I,J)
    145 END DO
!    
        DO 150 J=MYJS2_P4,MYJE2_P4
            DO 150 I=MYIS_P4,MYIE_P4
                  FNE(I,J) = (TEMPA(I+IHE(J),J) + TEMPA(I,J+1)) * (DPDE(I,J) + DPDE(I+IHE(J),J+1))
                  TNE(I,J) = FNE(I,J) * (TST(I+IHE(J),J+1) - TST(I,J))
                Q2MNE(I,J) = FNE(I,J) * (Q2M(I+IHE(J),J+1) - Q2M(I,J))
                Q2LNE(I,J) = FNE(I,J) * (Q2L(I+IHE(J),J+1) - Q2L(I,J))
                  FSE(I,J) = (TEMPB(I+IHE(J),J) + TEMPB(I,J-1)) * (DPDE(I,J) + DPDE(I+IHE(J),J-1))
                  TSE(I,J) = FSE(I,J) * (TST(I+IHE(J),J-1) - TST(I,J))
                Q2MSE(I,J) = FSE(I,J) * (Q2M(I+IHE(J),J-1) - Q2M(I,J))
                Q2LSE(I,J) = FSE(I,J) * (Q2L(I+IHE(J),J-1) - Q2L(I,J))
    150 END DO
!---------------------------------------------- 
! THERMODYNAMIC EQUATION AND MOISTURE 
!
! THE AD ARRAYS IN THE 170 LOOP ARE ON H POINTS
!----------------------------------------------
        DO 170 J=MYJS5_P2,MYJE5_P2
            DO 170 I=MYIS_P2,MYIE_P2
                  ADT(I,J) = (TEW(I+IHW(J),J  ) + TEW(I+IHE(J),J  )                               &
    &                      +  TNS(I       ,J-1) + TNS(I       ,J+1)                               &
    &                      +  TNE(I+IHW(J),J-1) + TNE(I       ,J  )                               &
    &                      +  TSE(I       ,J  ) + TSE(I+IHW(J),J+1))                              &
    &                      * RDPD(I       ,J  ) * FAD(I       ,J   )
!
                ADQ2M(I,J) = (Q2MEW(I+IHW(J),J  ) + Q2MEW(I+IHE(J),J  )                           &
    &                      +  Q2MNS(I       ,J-1) + Q2MNS(I       ,J+1)                           &
    &                      +  Q2MNE(I+IHW(J),J-1) + Q2MNE(I       ,J  )                           &
    &                      +  Q2MSE(I       ,J  ) + Q2MSE(I+IHW(J),J+1))                          &
    &                      *   RDPD(I       ,J  ) *   FAD(I       ,J  )
!
                ADQ2L(I,J) = (Q2LEW(I+IHW(J),J  ) + Q2LEW(I+IHE(J),J  )                           &
    &                      +  Q2LNS(I       ,J-1) + Q2LNS(I       ,J+1)                           &
    &                      +  Q2LNE(I+IHW(J),J-1) + Q2LNE(I       ,J  )                           &
    &                      +  Q2LSE(I       ,J  ) + Q2LSE(I+IHW(J),J+1))                          &
    &                      *   RDPD(I       ,J  ) *   FAD(I       ,J  )
    170 END DO
!---------------------------------- 
! UPSTREAM ADVECTION OF T, Q AND Q2 
!---------------------------------- 
        IF (UPSTRM) THEN
            DO 171 JAK=1,JAKONE
                JA = KHHAS(JAK)
                I =  IHLAS(JAK)
                J =  JHLAS(JAK)
!
                IX = I - MY_IS_GLB + 1
                JX = J - MY_JS_GLB + 1
!
                TTA = EMT(JA) * (UST(IX        ,JX-1) + UST(IX+IHW(JX),JX  )                      &
    &               +            UST(IX+IHE(JX),JX  ) + UST(IX        ,JX+1))
!
                TTB = ENT     * (VST(IX        ,JX-1) + VST(IX+IHW(JX),JX  )                      &
    &               +            VST(IX+IHE(JX),JX  ) + VST(IX        ,JX+1))
!
                PP = -TTA - TTB
                QP =  TTA - TTB
!            
                IF (PP < 0.) THEN
                    ISPA(JAK) = -1
                ELSE
                    ISPA(JAK) =  1
                END IF
!            
                IF (QP < 0.) THEN
                    ISQA(JAK) = -1
                ELSE
                    ISQA(JAK) =  1
                END IF
!            
                PP = ABS(PP)
                QP = ABS(QP)
!
                ARRAY3(JAK) = PP * QP
                ARRAY0(JAK) = ARRAY3(JAK) - PP - QP
                ARRAY1(JAK) = PP - ARRAY3(JAK)
                ARRAY2(JAK) = QP - ARRAY3(JAK)
        171 END DO
!        
            JAK = 0
!
            DO 173 JA=1,JAM
                IHL = IHLA(JA)
                IHH = IHHA(JA)
                J   =  JRA(JA)
!
                IF ( .NOT. LJRA(JA)) GOTO 173
!            
                DO I=IHL,IHH
                    IF (I >= MY_IS_GLB-ILPAD2 .AND. I <= MY_IE_GLB+IRPAD2) THEN
                        JAK = JAK + 1
                        ISP = ISPA(JAK)
                        ISQ = ISQA(JAK)
                        IFP = ( ISP - 1  ) / 2
                        IFQ = (-ISQ - 1  ) / 2
                        IPQ = ( ISP - ISQ) / 2
!                    
                        IX = I - MY_IS_GLB + 1
                        JX = J - MY_JS_GLB + 1
!                    
                        IF (HTM(IX+IHE(JX)+IFP,JX+ISP    ,K) *                                    &
    &                       HTM(IX+IHE(JX)+IFQ,JX    +ISQ,K) *                                    &
    &                       HTM(IX        +IPQ,JX+ISP+ISQ,K) > 0.1) GOTO 172
!                    
                        IF (HTM(IX+IHE(JX)+IFP,JX+ISP    ,K) +                                    &
    &                       HTM(IX+IHE(JX)+IFQ,JX    +ISQ,K) +                                    &
    &                       HTM(IX        +IPQ,JX+ISP+ISQ,K) < 0.1) THEN
!
                            TST(IX+IHE(JX)+IFP,JX+ISP    ) = TST(IX,JX)
                            TST(IX+IHE(JX)+IFQ,JX    +ISQ) = TST(IX,JX)
                            TST(IX        +IPQ,JX+ISP+ISQ) = TST(IX,JX)
!
                        ELSE IF (HTM(IX+IHE(JX)+IFP,JX+ISP    ,K) +                               &
    &                            HTM(IX        +IPQ,JX+ISP+ISQ,K) < 0.99) THEN
!
                            TST(IX+IHE(JX)+IFP,JX+ISP    ) = TST(IX            ,JX    )
                            TST(IX        +IPQ,JX+ISP+ISQ) = TST(IX+IHE(JX)+IFQ,JX+ISQ)
!
                        ELSE IF (HTM(IX+IHE(JX)+IFQ,JX    +ISQ,K) +                               &
    &                            HTM(IX+IPQ        ,JX+ISP+ISQ,K) < 0.99) THEN

                            TST(IX+IHE(JX)+IFQ,JX    +ISQ) = TST(IX            ,JX    )
                            TST(IX        +IPQ,JX+ISP+ISQ) = TST(IX+IHE(JX)+IFP,JX+ISP)
!
                        ELSE IF (HTM(IX+IHE(JX)+IFP,JX+ISP    ,K) +                               &
    &                            HTM(IX+IHE(JX)+IFQ,JX    +ISQ,K) < 0.99) THEN
!
                            TST(IX+IHE(JX)+IFP,JX+ISP    ) = 0.5                                  &
    &                                                      * (TST(IX    ,JX        )              &
    &                                                      +  TST(IX+IPQ,JX+ISP+ISQ))
!
                            TST(IX+IHE(JX)+IFQ,JX    +ISQ) = TST(IX+IHE(JX)+IFP,JX+ISP)
!
                        ELSE IF (HTM(IX+IHE(JX)+IFP,JX+ISP,K) < 0.99) THEN
!
                            TST(IX+IHE(JX)+IFP,JX+ISP) = TST(IX            ,JX        )           &
    &                                                  + TST(IX+IPQ        ,JX+ISP+ISQ)           &
    &                                                  - TST(IX+IHE(JX)+IFQ,JX    +ISQ)
!
                        ELSE IF (HTM(IX+IHE(JX)+IFQ,JX+ISQ,K) < 0.99) THEN
!
                            TST(IX+IHE(JX)+IFQ,JX+ISQ) = TST(IX            ,JX        )           &
    &                                                  + TST(IX+IPQ        ,JX+ISP+ISQ)           &
    &                                                  - TST(IX+IHE(JX)+IFP,JX+ISP    )
                        ELSE
                            TST(IX+IPQ,JX+ISP+ISQ) = TST(IX+IHE(JX)+IFP,JX+ISP    )               &
    &                                              + TST(IX+IHE(JX)+IFQ,JX    +ISQ)               &
    &                                              - TST(IX            ,JX        )
                        END IF
!                    
                    172 CONTINUE
!                    
                        F0 = ARRAY0(JAK)
                        F1 = ARRAY1(JAK)
                        F2 = ARRAY2(JAK)
                        F3 = ARRAY3(JAK)
!
                        ADT(IX,JX) = F0 * TST(IX,JX) + F1 * TST(IX+IHE(JX)+IFP,JX+ISP    )        &
    &                              +                   F2 * TST(IX+IHE(JX)+IFQ,JX    +ISQ)        &
    &                              +                   F3 * TST(IX+IPQ        ,JX+ISP+ISQ)
                    END IF
!
                END DO
        173 END DO
!        
            DO 175 JAK=1,JAKONE
                I = IHLAS(JAK)
                J = JHLAS(JAK)
!            
                IX = I - MY_IS_GLB + 1
                JX = J - MY_JS_GLB + 1
!            
                ISP = ISPA(JAK)
                ISQ = ISQA(JAK)
!
                IFP = ( ISP-1  ) / 2
                IFQ = (-ISQ-1  ) / 2
                IPQ = ( ISP-ISQ) / 2
!
                F0 = ARRAY0(JAK)
                F1 = ARRAY1(JAK)
                F2 = ARRAY2(JAK)
                F3 = ARRAY3(JAK)
!
                ADQ2M(IX,JX) = F0 * Q2M(IX,JX) + F1 *  Q2M(IX+IHE(JX)+IFP,JX+ISP    )             &
    &                        +                   F2 *  Q2M(IX+IHE(JX)+IFQ,JX    +ISQ)             &
    &                        +                   F3 *  Q2M(IX+IPQ        ,JX+ISP+ISQ)
!
                ADQ2L(IX,JX) = F0 * Q2L(IX,JX) + F1 * Q2L(IX+IHE(JX)+IFP,JX+ISP    )              &
    &                        +                   F2 * Q2L(IX+IHE(JX)+IFQ,JX    +ISQ)              &
    &                        +                   F3 * Q2L(IX+IPQ        ,JX+ISP+ISQ)
        175 END DO
!        
        END IF
!----------------------------------------------- 
! END OF THIS UPSTREAM REGION
!
! CALCULATION OF MOMENTUM ADVECTION COMPONENTS 
!
! THE FOLLOWING EW AND NS ARRAYS ARE ON H POINTS
!-----------------------------------------------
        DO 180 J=MYJS4_P4,MYJE4_P4
            DO 180 I=MYIS_P4,MYIE_P4
               UEW(I,J) = (FEW(I+IHW(J),J  ) + FEW(I+IHE(J),J  ))                                 &
    &                   * (UST(I+IHE(J),J  ) - UST(I+IHW(J),J  ))
!
               UNS(I,J) = (FNS(I+IHW(J),J  ) + FNS(I+IHE(J),J  ))                                 &
    &                   * (UST(I       ,J+1) - UST(I       ,J-1))
!
               VEW(I,J) = (FEW(I       ,J-1) + FEW(I       ,J+1))                                 &
    &                   * (VST(I+IHE(J),J  ) - VST(I+IHW(J),J  ))
!
               VNS(I,J) = (FNS(I       ,J-1) + FNS(I       ,J+1))                                 &
    &                   * (VST(I       ,J+1) - VST(I       ,J-1))
!---------------------------------------------------- 
! THE FOLLOWING NE AND SE ARRAYS ARE TIED TO V POINTS
!---------------------------------------------------- 
               UNE(I,J) = (FNE(I+IVW(J),J  ) + FNE(I+IVE(J),J  ))                                 &
    &                   * (UST(I+IVE(J),J+1) - UST(I      ,J   ))
!
               USE(I,J) = (FSE(I+IVW(J),J  ) + FSE(I+IVE(J),J  ))                                 &
    &                   * (UST(I+IVE(J),J-1) - UST(I       ,J  ))
!
               VNE(I,J) = (FNE(I       ,J-1) + FNE(I       ,J+1))                                 &
    &                   * (VST(I+IVE(J),J+1) - VST(I       ,J  ))
!
               VSE(I,J) = (FSE(I       ,J-1) + FSE(I       ,J+1))                                 &
    &                   * (VST(I+IVE(J),J-1) - VST(I       ,  J))
    180 END DO
!----------------------------  
! EQUATION OF MOTION 
!
! ADU AND ADV ARE ON V POINTS
!----------------------------  
        DO 200 J=MYJS5_P2,MYJE5_P2
            DO 200 I=MYIS_P2,MYIE_P2
                ADU(I,J) =  (UEW(I+IVW(J),J  ) + UEW(I+IVE(J),J  )                                &
    &                    +   UNS(I       ,J-1) + UNS(I       ,J+1)                                &
    &                    +   UNE(I+IVW(J),J-1) + UNE(I       ,J  )                                &
    &                    +   USE(I       ,J  ) + USE(I+IVW(J),J+1))                               &
    &                    * RDPDX(I       ,J  ) * FAD(I+IVW(J),J  )
!
                ADV(I,J) =  (VEW(I+IVW(J),J  ) + VEW(I+IVE(J),J  )                                &
    &                    +   VNS(I       ,J-1) + VNS(I       ,J+1)                                &
    &                    +   VNE(I+IVW(J),J-1) + VNE(I       ,J  )                                &
    &                    +   VSE(I       ,J  ) + VSE(I+IVW(J),J+1))                               &
    &                    * RDPDY(I       ,J  ) * FAD(I+IVW(J),J  )
    200 END DO
!------------------------------------------     
! UPSTREAM ADVECTION OF VELOCITY COMPONENTS 
!------------------------------------------    
        IF (UPSTRM) THEN
            DO 205 JAK=1,JAKTWO
                JA = KVHAS(JAK)
                I  = IVLAS(JAK)
                J  = JVLAS(JAK)
!            
                IX = I - MY_IS_GLB + 1
                JX = J - MY_JS_GLB + 1
!            
                TTA = EM(JA) * UST(IX,JX)
                TTB = EN     * VST(IX,JX)
!
                PP = -TTA - TTB
                QP =  TTA - TTB
!            
                IF (PP < 0.) THEN
                    ISP = -1
                ELSE
                    ISP =  1
                END IF
!            
                IF (QP < 0.) THEN
                    ISQ = -1
                ELSE
                    ISQ =  1
                END IF
!            
                IFP = ( ISP - 1  ) / 2
                IFQ = (-ISQ - 1  ) / 2
                IPQ = ( ISP - ISQ) / 2
!
                PP = ABS(PP)
                QP = ABS(QP)
!
                F3 = PP * QP
                F0 = F3 - PP - QP
                F1 = PP - F3
                F2 = QP - F3
!
                ADU(IX,JX) = F0 * UST(IX,JX) + F1 * UST(IX+IVE(JX)+IFP,JX+ISP    )                &
    &                      +                   F2 * UST(IX+IVE(JX)+IFQ,JX    +ISQ)                &
    &                      +                   F3 * UST(IX        +IPQ,JX+ISP+ISQ)
!
                ADV(IX,JX) = F0 * VST(IX,JX) + F1 * VST(IX+IVE(JX)+IFP,JX+ISP    )                &
    &                                        + F2 * VST(IX+IVE(JX)+IFQ,JX    +ISQ)                &
    &                                        + F3 * VST(IX        +IPQ,JX+ISP+ISQ)
        205 END DO
        END IF
!---------------------------- 
! END OF THIS UPSTREAM REGION
!----------------------------
        IF (ITER2) GOTO 235
!
        DO 220 J=MYJS2_P2,MYJE2_P2
            DO 220 I=MYIS1_P2,MYIE1_P2
                TST(I,J) = ADT  (I,J) * (HM(I,J) * TLC) + TST(I,J)
                Q2M(I,J) = ADQ2M(I,J) * (HM(I,J) * TLC) + Q2M(I,J)
                Q2L(I,J) = ADQ2L(I,J) * (HM(I,J) * TLC) + Q2L(I,J)
    220 END DO
!    
        DO 230 J=MYJS2_P2,MYJE2_P2
            DO 230 I=MYIS1_P2,MYIE1_P2
                UST(I,J) = ADU(I,J) * VM(I,J) * TLC + UST(I,J)
                VST(I,J) = ADV(I,J) * VM(I,J) * TLC + VST(I,J)
    230 END DO
!
        ITER2 = .TRUE.
! 
        GOTO 125
!
    235 DO 240 J=MYJS2,MYJE2
        DO 240 I=MYIS1,MYIE1
            T(I,J,K) = ADT(I,J) * (2.0*HM(I,J)) + T(I,J,K)
    240 END DO
!    
        DO 250 J=MYJS2,MYJE2
            DO 250 I=MYIS1,MYIE1
                U(I,J,K) = ADU(I,J) * (2.0*VM(I,J)) + U(I,J,K)
                V(I,J,K) = ADV(I,J) * (2.0*VM(I,J)) + V(I,J,K)
    250 END DO
!
        IF (K == LM) THEN
            DO 260 J=MYJS2,MYJE2
                DO 260 I=MYIS1,MYIE1
                    ADQ2HL(I,J,LM-1) = ADQ2M(I,J)
        260 END DO
!
        ELSE IF (K > 1 .AND. K < LM) THEN
            DO 270 J=MYJS2,MYJE2
                DO 270 I=MYIS1,MYIE1
                    ADQ2HL(I,J,K-1) = ADQ2M(I,J)
                        Q2(I,J,K  ) = ADQ2L(I,J) * HM(I,J) + Q2(I,J,K)
        270 END DO
!
        ELSE
            DO 280 J=MYJS2,MYJE2
                DO 280 I=MYIS1,MYIE1
                    Q2(I,J,K) = ADQ2L(I,J) * HM(I,J) + Q2(I,J,K)
        280 END DO
        END IF
!
500 END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (HM)
!
    DO 600 K=2,LM-1
        DO J=MYJS2,MYJE2
            DO I=MYIS1,MYIE1
                HM(I,J)   =    HTM(I,J,K) * HBM2(I,J)
                Q2(I,J,K) = ADQ2HL(I,J,K) *   HM(I,J) + Q2(I,J,K)
            END DO
        END DO
600 END DO
!
    RETURN
!
    END SUBROUTINE HZADV_LM1
