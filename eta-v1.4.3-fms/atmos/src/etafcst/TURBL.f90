    SUBROUTINE TURBL
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE TURBL
!>
!> SUBROUTINE: TURBL - VERTICAL TURBULENT EXCHANGE
!> PROGRAMMER: JANJIC
!> ORG: W/NP2
!> DATE: 95-03-20
!>
!> ABSTRACT:
!> TURBL UPDATES THE TURBULENT KINETIC ENERGY WITH THE PRODUCTION/DISSIPATION TERM AND THE VERTICAL 
!> DIFFUSION TERM DIFFUSION TERM (USING AN IMPLICIT FORMULATION). EXCHANGE COEFFICIENTS FOR THE 
!> SURFACE AND FOR ALL LAYER INTERFACES ARE THEN COMPUTED AND THE EXCHANGE IS EXECUTED.
!>
!> PROGRAM HISTORY LOG:
!> 95-03-15  JANJIC     - ORIGINATOR
!> 95-03-28  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 96-03-29  BLACK      - ADDED EXTERNAL EDGE; REMOVED SCRCH COMMON
!> 96-07-19  MESINGER   - ADDED Z0 EFFECTIVE
!> 98-??-??  TUCCILLO   - MODIFIED FOR CLASS VIII PARALLELISM
!> 98-10-27  BLACK      - PARALLEL CHANGES INTO MOST RECENT CODE
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
!> USE MODULES: CTLBLK
!>              DYNAM
!>              EXCHM
!>              F77KINDS
!>              GLB_TABLE
!>              INDX
!>              LOOPS
!>              MAPPINGS
!>              MASKS
!>              MOMENTO
!>              MPPCOM
!>              PARMETA
!>              PARM_TBL
!>              PHYS
!>              PVRBLS
!>              TEMPCOM
!>              TIMMING
!>              TOPO
!>              VRBLS
!>              Z0EFFT         
!>
!> DRIVER     : EBU
!>              NEWFLT
!>
!> CALLS      : DIFCOF
!>              EXCH
!>              MIXLEN
!>              PRODQ2
!>              SGETMO
!>              SFCDIF
!>              SURFCE
!>              VDIFH
!>              VDIFQ
!>              VDIFV
!>              ZERO2
!>              ZERO3
!>              ZERO3_T
!>--------------------------------------------------------------------------------------------------
    USE CTLBLK
    USE DYNAM
    USE EXCHM
    USE F77KINDS
    USE GLB_TABLE
    USE INDX
    USE LOOPS
    USE MAPPINGS
    USE MASKS
    USE MOMENTO
    USE MPPCOM
    USE PARMETA
    USE PARM_TBL, ONLY : ITB, JTB , ITBQ , JTBQ
    USE PHYS
    USE PVRBLS
    USE SOIL, ONLY : IVGTYP
    USE TEMPCOM
    USE TIMMING
    USE TOPO
    USE VRBLS
    USE Z0EFFT
!
    IMPLICIT NONE
!
    INCLUDE "EXCHM.h"
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: KTMQ2 =   1
!
    REAL   (KIND=R4KIND), PARAMETER :: CAPA  =   0.28589641
    REAL   (KIND=R4KIND), PARAMETER :: G     =   9.8 
    REAL   (KIND=R4KIND), PARAMETER :: RG    =   1.   /  G 
    REAL   (KIND=R4KIND), PARAMETER :: ROG   = 287.04 /  G
    REAL   (KIND=R4KIND), PARAMETER :: EPSZ  =   1.E-4
    REAL   (KIND=R4KIND), PARAMETER :: EPSQ2 =   0.12
!
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM  = IM * JM - JM / 2
    INTEGER(KIND=I4KIND), PARAMETER :: NHRZ  = (IDIM2 - IDIM1 + 1) * (JDIM2 - JDIM1 + 1)
! 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & CKLQ    , CT      , UZ0H    ,  VZ0H   , AKMSV
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & APE     , UCOL    , VCOL
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM1)                              ::&
    & AKH     , AKM     , AKMCOL  , AKHCOL
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LP1)                              ::&
    & ZINT    , ZCOL
!
    REAL   (KIND=R4KIND), DIMENSION(LM1, IDIM1:IDIM2, JDIM1:JDIM2)                              ::&
    & AKH_T   , AKM_T
!
    REAL   (KIND=R4KIND), DIMENSION(LM, IDIM1:IDIM2, JDIM1:JDIM2)                               ::&
    & APECOL_T, UCOL_T  , VCOL_T  , TCOL_T  , QCOL_T  , Q2COL_T
!
    REAL   (KIND=R4KIND), DIMENSION(LP1, IDIM1:IDIM2, JDIM1:JDIM2)                              ::&
    & ZCOL_T  , ZCOL_T2
!
    REAL   (KIND=R4KIND), DIMENSION(LM1)                                                        ::&
    & GM      , GH      , EL
! 
    REAL   (KIND=R4KIND), DIMENSION(4)                                                          ::&
    & ZEFF
!------------------------------------------------ 
! THE FOLLOWING ARE USED FOR TIMIMG PURPOSES ONLY
!------------------------------------------------ 
    REAL   (KIND=R8KIND)                                                                        ::&
    & TIMEF
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K       , LMHK    , LMHP    , LPBL    , IIM     , JJM     , JJ      ,   &
    & II      , LMVK    , LMHM    , IVGTYPIJ
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & PDSL    , APESTS  , WMSK    , RWMSK   , HPBL    , ULM     , VLM     , WSTAR   , BTIM    ,   &
    & EVPR
!
    CALL ZERO3(AKM, LM1)
    CALL ZERO3(ZINT, LP1)
    CALL ZERO3_T(AKH_T, LM1)
    CALL ZERO3_T(AKM_T, LM1)
    CALL ZERO2(UZ0H)
    CALL ZERO2(VZ0H)
    CALL ZERO2(AKM10)
!--------------------------------------------------------------------- 
! AKMS WAS ZEROED OUT BY MISTAKE, SEE THE CORRESP. WITH SM, FEB.22, 07
!---------------------------------------------------------------------
!
!-------------------------------------------------------------------     
! COMPUTE THE HEIGHTS OF THE LAYER INTERFACES AND THE EXNER FUNCTION
!------------------------------------------------------------------- 
!
!$omp parallel do 
!
    DO J=MYJS_P1,MYJE_P1
        DO I=MYIS_P1,MYIE_P1
            ZINT(I,J,LP1) = EPSZ
            IF (SIGMA) ZINT(I,J,LP1) = RG * FIS(I,J)
        END DO
    END DO
!
    DO 10 K=LM,1,-1
!
!$omp parallel do private (APESTS, PDSL)
!
        DO J=MYJS1_P1,MYJE1_P1 
            DO I=MYIS_P1,MYIE_P1
                PDSL   = PD(I,J) * RES(I,J)
                APESTS = PDSL    * AETA(K) + PT
!            
                ZINT(I,J,K) =  ZINT(I,J,K+1) + T(I,J,K)  / APESTS     * PDSL * DETA(K) * ROG      &
    &                       *    (Q(I,J,K)   * 0.608     + 1.)
!
                ZINT(I,J,K) = (ZINT(I,J,K)   - DFRLG(K)) * HTM(I,J,K) + DFRLG(K)
!            
                 APE(I,J,K) = (1.E5/APESTS)  ** CAPA
            END DO
        END DO
!
 10 END DO
!------------------- 
! REMOVE NEGATIVE Q2
!------------------- 
!
!$omp parallel do 
!
    DO 40 K=1,LM
        DO J=MYJS_P1,MYJE_P1
            DO I=MYIS_P1,MYIE_P1
                Q2(I,J,K) = AMAX1(Q2(I,J,K) * HBM2(I,J), EPSQ2)
            END DO
        END DO
 40 END DO
!
!$omp parallel do 
!
    DO J=MYJS2_P1,MYJE2_P1
        DO I=MYIS_P1,MYIE_P1
            UZ0H(I,J) = (UZ0(I+IHE(J),J) + UZ0(I+IHW(J),J) + UZ0(I,J+1) + UZ0(I,J-1))             &
    &                 * HBM2(I       ,J) * 0.25
!
            VZ0H(I,J) = (VZ0(I+IHE(J),J) + VZ0(I+IHW(J),J) + VZ0(I,J+1) + VZ0(I,J-1))             &
    &                 * HBM2(I       ,J) * 0.25
        END DO
    END DO
!---------------------------------- 
! PREPARE THE EXCHANGE COEFFICIENTS
!---------------------------------- 
!
!---------------------------------------- 
! COMPUTE VELOCITY COMPONENTS AT H POINTS
!---------------------------------------- 
!
!$omp parallel do private (RWMSK, WMSK)
!
    DO 60 K=1,LM
!    
        DO J=MYJS2_P1,MYJE2_P1
            DO I=MYIS_P1,MYIE_P1
                WMSK = VTM(I+IHE(J),J,K) + VTM(I+IHW(J),J,K) + VTM(I,J+1,K) + VTM(I,J-1,K)
                IF (WMSK > 0.) THEN
                    RWMSK       = 0.25
                    UCOL(I,J,K) = (U(I+IHE(J),J  ,K) * VTM(I+IHE(J),J  ,K)                         &
    &                           +  U(I+IHW(J),J  ,K) * VTM(I+IHW(J),J  ,K)                         &
    &                           +  U(I       ,J+1,K) * VTM(I       ,J+1,K)                         &
    &                           +  U(I       ,J-1,K) * VTM(I       ,J-1,K))                        &
    &                           * RWMSK
!                    
                    VCOL(I,J,K) = (V(I+IHE(J),J  ,K) * VTM(I+IHE(J),J  ,K)                         &
    &                           +  V(I+IHW(J),J  ,K) * VTM(I+IHW(J),J  ,K)                         &
    &                           +  V(I       ,J+1,K) * VTM(I       ,J+1,K)                         &
    &                           +  V(I       ,J-1,K) * VTM(I       ,J-1,K))                        &
    &                           *  RWMSK
                ELSE
                    UCOL(I,J,K) = 0.
                    VCOL(I,J,K) = 0.
                END IF
            END DO
        END DO
!
 60 END DO
!----------------------- 
! FILL TRANSPOSED ARRAYS
!----------------------- 
! 
!$omp parallel sections
!$omp section
!
    CALL SGETMO(T   , NHRZ, NHRZ, LM ,   TCOL_T, LM)
!
!$omp section
!
    CALL SGETMO(Q   , NHRZ, NHRZ, LM ,   QCOL_T, LM)
!
!$omp section
!
    CALL SGETMO(APE , NHRZ, NHRZ, LM , APECOL_T, LM)
!
!$omp section
!
    CALL SGETMO(Q2  , NHRZ, NHRZ, LM ,  Q2COL_T, LM)
!
!$omp section
!
    CALL SGETMO(ZINT, NHRZ, NHRZ, LP1,   ZCOL_T, LP1)
!
!$omp section
!
    CALL SGETMO(UCOL, NHRZ, NHRZ, LM ,   UCOL_T, LM)
!
!$omp section
!
    CALL SGETMO(VCOL, NHRZ, NHRZ, LM ,   VCOL_T, LM)
!
!$omp end parallel sections
!
!----------------------- 
! FIND THE MIXING LENGTH
!----------------------- 
!
!$omp parallel do private (EL, GH, GM, HPBL, LMHK, LMHM, LMHP, LPBL)
!$omp private (ULM, VLM, WSTAR, ZEFF)
!
    DO 100 J=MYJS2_P1,MYJE2_P1
        DO 100 I=MYIS_P1,MYIE1_P1
!
            LMHK = LMH(I,J)
            LMHP = LMHK + 1
            LMHM = LMHK - 1
	    IVGTYPIJ = IVGTYP(I,J)
!        
            CALL MIXLEN(LMHK, LPBL, HPBL,                                                         &
    &                      UCOL_T(1,I,J),  VCOL_T(1,I,J),   TCOL_T(1,I,J),                        &
    &                      QCOL_T(1,I,J), Q2COL_T(1,I,J), APECOL_T(1,I,J),                        &
    &                      ZCOL_T(1,I,J),                                                         &
    &                   GM  , GH  , EL)
!--------------------------------------------------------------------- 
! SOLVE FOR THE PRODUCTION/DISSIPATION OF THE TURBULENT KINETIC ENERGY
!---------------------------------------------------------------------      
            CALL PRODQ2(LMHK, DTQ2,                                                               &
    &                   USTAR(I,J),                                                               &
    &                   GM, GH, EL,                                                               &
    &                   Q2COL_T(1,I,J))
!------------------------------------------------------  
! FIND THE EXCHANGE COEFFICIENTS IN THE FREE ATMOSPHERE
!------------------------------------------------------ 
            CALL DIFCOF(LMHK, GM, GH, EL,                                                         &
    &                   Q2COL_T(1,I,J), ZCOL_T(1,I,J), AKM_T(1,I,J), AKH_T(1,I,J))
!------------------------------------------------------------- 
! CARRY OUT THE VERTICAL DIFFUSION OF TURBULENT KINETIC ENERGY
!-------------------------------------------------------------      
            CALL VDIFQ(LMHK, KTMQ2, DTQ2,                                                         &
    &                  Q2COL_T(1,I,J),                                                            &
    &                  EL,                                                                        &
    &                  ZCOL_T(1,I,J))
!---------------------- 
! FIND THE Z0 EFFECTIVE
!----------------------  
            ZEFF(1) = ZEFFIJ(I,J,1)
            ZEFF(2) = ZEFFIJ(I,J,2)
            ZEFF(3) = ZEFFIJ(I,J,3)
            ZEFF(4) = ZEFFIJ(I,J,4)
!---------------------------------------  
! FIND THE SURFACE EXCHANGE COEFFICIENTS
!--------------------------------------- 
            ULM = UCOL(I,J,LMHK)
            VLM = VCOL(I,J,LMHK)
!        
            CALL SFCDIF(LMHK         ,        SM(I,J),        THS(I,J),       QS(I,J),            &
    &                       UZ0H(I,J),      VZ0H(I,J),       THZ0(I,J),      QZ0(I,J),            &
    &                      USTAR(I,J),  WSTAR        ,         Z0(I,J),  ZEFF        ,            &
    &                       AKMS(I,J),      AKHS(I,J), HPBL           ,       CT(I,J),            &
    &                        U10(I,J),       V10(I,J),     TSHLTR(I,J),     TH10(I,J),            &    
    &                     QSHLTR(I,J),       Q10(I,J),                                            &
!---------- 
! GSM V100M
!----------
    &                    TH100(I,J)  ,      Q100(I,J),       U100(I,J),     V100(I,J),            &
    &                   ULM          , VLM           ,                                            & 
    &                   TCOL_T(1,I,J),  QCOL_T(1,I,J), APECOL_T(1,I,J), ZCOL_T(1,I,J),            &
    &                         PD(I,J), PT            ,        T(I,J,LMHK),  IVGTYPIJ)
!-------- 
! MORELLI
!--------
            XMOMFLUX(I,J) = UMFLX
            YMOMFLUX(I,J) = VMFLX
               AKM10(I,J) = AKMS10
!       
100 END DO
!-------------------------------------------- 
! FILL STANDARD ARRAYS FROM TRANSPOSED ARRAYS
!--------------------------------------------
!
!$omp parallel sections
!$omp section
!
    CALL SGETMO(Q2COL_T, LM , LM , NHRZ, Q2 , NHRZ)
!
!$omp section
!
    CALL SGETMO(  AKH_T, LM1, LM1, NHRZ, AKH, NHRZ)
!
!$omp section
!
    CALL SGETMO(  AKM_T, LM1, LM1, NHRZ, AKM, NHRZ)
!
!$omp end parallel sections
!
!-------------------------------------------------------------- 
! UNCOMPUTED LOCATIONS MUST BE FILLED IN FOR THE POST-PROCESSOR
!-------------------------------------------------------------- 
    IIM = IM - MY_IS_GLB + 1
    JJM = JM - MY_JS_GLB + 1
!------------------------  
! EASTERN GLOBAL BOUNDARY
!------------------------
    IF (MY_IE_GLB == IM) THEN
        DO J=1,JM
            IF (J >= MY_JS_GLB .AND. J <= MY_JE_GLB) THEN
                           JJ  = J - MY_JS_GLB + 1
!
                  TH10(IIM,JJ) =   TH10(IIM-1,JJ)
                   Q10(IIM,JJ) =    Q10(IIM-1,JJ)
                   U10(IIM,JJ) =    U10(IIM-1,JJ)
                   V10(IIM,JJ) =    V10(IIM-1,JJ)
                TSHLTR(IIM,JJ) = TSHLTR(IIM-1,JJ)
                QSHLTR(IIM,JJ) = QSHLTR(IIM-1,JJ)
!---------- 
! GSM V100M
!----------
                 TH100(IIM,JJ) =  TH100(IIM-1,JJ)
                  Q100(IIM,JJ) =   Q100(IIM-1,JJ)
                  U100(IIM,JJ) =   U100(IIM-1,JJ)
                  V100(IIM,JJ) =   V100(IIM-1,JJ)
            END IF
        END DO
    END IF
!-------------------------
! SOUTHERN GLOBAL BOUNDARY
!-------------------------
    IF (MY_JS_GLB == 1) THEN
        DO J=1,2
            DO I=1,IM
                IF (I >= MY_IS_GLB .AND. I <= MY_IE_GLB) THEN
                           II    = I - MY_IS_GLB + 1
!
                      TH10(II,J) =   TH10(II,3)
                       Q10(II,J) =    Q10(II,3)
                       U10(II,J) =    U10(II,3)
                       V10(II,J) =    V10(II,3)
                    TSHLTR(II,J) = TSHLTR(II,3)
                    QSHLTR(II,J) = QSHLTR(II,3)
!---------- 
! GSM V100M
!----------
                     TH100(II,J) = TH100(II,3)
                      Q100(II,J) =  Q100(II,3)
                      U100(II,J) =  U100(II,3)
                      V100(II,J) =  V100(II,3)
                END IF
            END DO
        END DO
    END IF
!-------------------------
! NORTHERN GLOBAL BOUNDARY
!-------------------------
    IF (MY_JE_GLB == JM) THEN
        DO J=JM-1,JM
            JJ = J - MY_JS_GLB + 1
            DO I=1,IM
                IF (I >= MY_IS_GLB .AND. I <= MY_IE_GLB) THEN
                          II      = I - MY_IS_GLB + 1
!
                      TH10(II,JJ) =   TH10(II,JJM-2)
                       Q10(II,JJ) =    Q10(II,JJM-2)
                       U10(II,JJ) =    U10(II,JJM-2)
                       V10(II,JJ) =    V10(II,JJM-2)
                    TSHLTR(II,JJ) = TSHLTR(II,JJM-2)
                    QSHLTR(II,JJ) = QSHLTR(II,JJM-2)
!---------- 
! GSM V100M
!----------
                     TH100(II,JJ) =  TH100(II,JJM-2)
                      Q100(II,JJ) =   Q100(II,JJM-2)
                      U100(II,JJ) =   U100(II,JJM-2)
                      V100(II,JJ) =   V100(II,JJM-2)
                END IF
            END DO
        END DO
    END IF
!
    BTIM = TIMEF()
!
    CALL EXCH(UZ0H, 1, VZ0H, 1, 1, 1)
!
    EXCH_TIM = EXCH_TIM + TIMEF() - BTIM
!-------------------------------------
! AVERAGE UZ0 AND VZ0 BACK TO V POINTS
!-------------------------------------
!
!$omp parallel do
!
    DO 125 J=MYJS2,MYJE2
        DO 125 I=MYIS,MYIE
            UZ0(I,J) = (UZ0H(I+IVE(J),J  ) * HBM2(I+IVE(J),J  )                                   &
    &                +  UZ0H(I+IVW(J),J  ) * HBM2(I+IVW(J),J  )                                   &
    &                +  UZ0H(I       ,J+1) * HBM2(I       ,J+1)                                   &
    &                +  UZ0H(I       ,J-1) * HBM2(I       ,J-1)) * 0.25
 !           
            VZ0(I,J) = (VZ0H(I+IVE(J),J  ) * HBM2(I+IVE(J),J  )                                   & 
    &                +  VZ0H(I+IVW(J),J  ) * HBM2(I+IVW(J),J  )                                   &
    &                +  VZ0H(I       ,J+1) * HBM2(I       ,J+1)                                   &
    &                +  VZ0H(I       ,J-1) * HBM2(I       ,J-1)) * 0.25
125 END DO
!----------------------------- 
! EXECUTE THE GROUND PROCESSES
!----------------------------- 
    CALL SURFCE(APE(IDIM1, JDIM1, 1), ZINT(IDIM1, JDIM1, 1), CKLQ(IDIM1, JDIM1))
!------------------------------ 
! EXECUTE THE VERTICAL EXCHANGE
!------------------------------ 
    BTIM=TIMEF()
!
    CALL EXCH(AKMS, 1, AKM, LM1, ZINT, LP1, 1, 1)
!
    EXCH_TIM = EXCH_TIM + TIMEF() - BTIM
!
!$omp parallel do 
!
    DO K=1,LM1
        DO J=MYJS2,MYJE2
            DO I=MYIS,MYIE
                AKMCOL(I,J,K) = (AKM(I+IVE(J),J  ,K) * HBM2(I+IVE(J),J  )                         &
    &                         +  AKM(I+IVW(J),J  ,K) * HBM2(I+IVW(J),J  )                         &
    &                         +  AKM(I       ,J+1,K) * HBM2(I       ,J+1)                         &
                              +  AKM(I       ,J-1,K) * HBM2(I       ,J-1))                        & 
    &                         *  VTM(I       ,J  ,K) * VBM2(I       ,J  ) * 0.25
!                
                AKHCOL(I,J,K) =  AKH(I       ,J  ,K) *  HTM(I       ,J,K) * HBM2(I,J)
            END DO
        END DO
    END DO
!
!$omp parallel do
! 
    DO J=MYJS2,MYJE2
        DO I=MYIS,MYIE
             THZ0(I,J) =   (1.-SM(I       ,J))  *  THS(I       ,J)                                &
    &                  +       SM(I       ,J)   * THZ0(I       ,J)
!             
             QZ0 (I,J) =   (1.-SM(I       ,J))  *  QS (I       ,J)                                &
    &                  +       SM(I       ,J)   * QZ0 (I       ,J)
!           
            AKMSV(I,J) =    (AKMS(I+IVE(J),J)   * HBM2(I+IVE(J),J)                                &
    &                  +     AKMS(I+IVW(J),J)   * HBM2(I+IVW(J),J)                                &
    &                  +     AKMS(I       ,J+1) * HBM2(I       ,J+1)                              &
    &                  +     AKMS(I       ,J-1) * HBM2(I       ,J-1)) * VBM2(I,J) * 0.25
        END DO
    END DO
!
!$omp parallel do 
!
    DO K=1,LP1
        DO J=MYJS2,MYJE2
            DO I=MYIS,MYIE
                ZCOL(I,J,K) = 0.25 * (ZINT(I+IVE(J),J  ,K) + ZINT(I+IVW(J),J  ,K)                 &
    &                       +         ZINT(I       ,J+1,K) + ZINT(I       ,J-1,K))
            END DO
        END DO
    END DO
!---------
! MORELLI
!---------
    BTIM=TIMEF()
!
    CALL EXCH(AKM10, 1, 1, 1)
!
    EXCH_TIM = EXCH_TIM + TIMEF() - BTIM
!    
    DO J=MYJS2,MYJE2
        DO I=MYIS,MYIE
                   LMVK =    LMV(I       ,J)
!
            AKM10V(I,J) = (AKM10(I+IVE(J),J  ) * HBM2(I+IVE(J),J  )                               & 
    &                   +  AKM10(I+IVW(J),J  ) * HBM2(I+IVW(J),J  )                               &
    &                   +  AKM10(I       ,J+1) * HBM2(I       ,J+1)                               &
    &                   +  AKM10(I       ,J-1) * HBM2(I       ,J-1)) * VBM2(I,J) * 0.25
!
            IF (AKM10V(I,J) == 0) THEN
                   U10(I,J) =  0.
                   V10(I,J) =  0.
            ELSE
                   U10(I,J) = AKMSV(I,J) * (U(I,J,LMVK) - UZ0(I,J)) / AKM10V(I,J) + UZ0(I,J)
                   V10(I,J) = AKMSV(I,J) * (V(I,J,LMVK) - VZ0(I,J)) / AKM10V(I,J) + VZ0(I,J)
            END IF
        END DO
    END DO
!-----------------
! TRANSPOSE ARRAYS
!-----------------
!
!$omp parallel sections
!$omp section
!
    CALL SGETMO(ZCOL  , NHRZ, NHRZ, LP1, ZCOL_T2, LP1)
!
!$omp section
!
    CALL SGETMO(U     , NHRZ, NHRZ, LM , UCOL_T , LM)
!
!$omp section
!
    CALL SGETMO(V     , NHRZ, NHRZ, LM , VCOL_T , LM)
!
!$omp section
!
    CALL SGETMO(AKHCOL, NHRZ, NHRZ, LM1,  AKH_T , LM1)
!
!$omp section
!
    CALL SGETMO(AKMCOL, NHRZ, NHRZ, LM1,  AKM_T , LM1)
!
!$omp end parallel sections
!
!$omp parallel do private (LMHK,LMVK)
!
    DO 200 J=MYJS2,MYJE2
        DO 200 I=MYIS,MYIE1
!        
            LMHK = LMH(I,J)
            LMVK = LMV(I,J)
!----------------------------------------------------------------
! CARRY OUT THE VERTICAL DIFFUSION OF TEMPERATURE AND WATER VAPOR
!----------------------------------------------------------------
            CALL VDIFH(LMHK, KTMQ2, DTQ2,                                                         &
    &                    THZ0(I,J),      QZ0(I,J),    AKHS(I,J),                                  &
    &                      CT(I,J),     CKLQ(I,J),                                                &
    &                TCOL_T(1,I,J), QCOL_T(1,I,J), AKH_T(1,I,J),                                  &
    &              APECOL_T(1,I,J), ZCOL_T(1,I,J), EVPR)
!            
            PD(I,J) = PD(I,J) + PD(I,J) * RES(I,J) * EVPR
!-------------------------------------------------------- 
! CARRY OUT THE VERTICAL DIFFUSION OF VELOCITY COMPONENTS
!-------------------------------------------------------- 
            CALL VDIFV(LMVK         , KTMQ2        , DTQ2        ,                                &
    &                       UZ0(I,J),      VZ0(I,J),   AKMSV(I,J),                                &
    &                  UCOL_T(1,I,J), VCOL_T(1,I,J), AKM_T(1,I,J), ZCOL_T2(1,I,J))
!
200 END DO
!--------------------------------------------  
! FILL STANDARD ARRAYS FROM TRANSPOSED ARRAYS
!-------------------------------------------- 
!
!$omp parallel sections
!$omp section
!
    CALL SGETMO(QCOL_T, LM, LM, NHRZ, Q, NHRZ)
!
!$omp section
!
    CALL SGETMO(TCOL_T, LM, LM, NHRZ, T, NHRZ)
!
!$omp section
!
    CALL SGETMO(UCOL_T, LM, LM, NHRZ, U, NHRZ)
!
!$omp section
!
    CALL SGETMO(VCOL_T, LM, LM, NHRZ, V, NHRZ)
!
!$omp end parallel sections
!
!------------------- 
! REMOVE NEGATIVE Q2
!------------------- 
!
!$omp parallel do 
!
    DO K=1,LM
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                Q2(I,J,K) = AMAX1(Q2(I,J,K) * HBM2(I,J), EPSQ2)
            END DO
        END DO
    END DO
!
    RETURN
!
    END SUBROUTINE TURBL
