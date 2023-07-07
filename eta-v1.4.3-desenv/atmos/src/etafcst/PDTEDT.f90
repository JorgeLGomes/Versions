    SUBROUTINE PDTEDT
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE PDTEDT
!>
!> SUBROUTINE: PDTEDT - SURFACE PRESSURE TENDENCY CALC
!> PROGRAMMER: JANJIC
!> ORG: W/NMC2
!> DATE: 96-07-??
!>
!> ABSTRACT:
!> PDTEDT VERTICALLY INTEGRATES THE MASS FLUX DIVERGENCE TO OBTAIN THE SURFACE PRESSURE TENDENCY AND
!> ETADOT ON THE LAYER INTERFACES.
!> THEN IT UPDATES THE HYDROSTATIC SURFACE PRESSURE, THE NONHYDROSTATIC PRESSURE, AND ADDS THE LOCAL
!> TIME DERIVATIVE AND VERTICAL ADVECTION OF NONHYDROSTATIC PRESSURE CONTRIBUTION TO THE OMEGA-ALPHA 
!> TERM OF THE THERMODYNAMIC EQUATION. 
!> ALSO, THE OMEGA-ALPHA TERM IS COMPUTED FOR DIAGNOSTICS.
!>
!> PROGRAM HISTORY LOG:
!> 87-06-??  JANJIC     - ORIGINATOR
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 96-05-??  JANJIC     - ADDED NONHYDROSTATIC EFFECTS & MERGED THE PREVIOUS SUBROUTINES PDTE AND 
!>                        PDNEW
!> 00-01-03  BLACK      - DISTRIBUTED MEMORY AND THREADS
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
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
!> USE MODULES: CONTIN
!>              CTLBLK
!>              DYNAM
!>              EXCHM
!>              F77KINDS
!>              GLB_TABLE
!>              INDX
!>              LOOPS
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              NHYDRO
!>              PARMETA
!>              TEMPCOM
!>              TIMMING
!>              TOPO
!>              VRBLS     
!>
!> DRIVER     : DIGFLT
!>              EBU
!>
!> CALLS      : EXCH
!>              ZERO2
!>--------------------------------------------------------------------------------------------------
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
    USE TEMPCOM
    USE TIMMING
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    INCLUDE "EXCHM.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM  = IM * JM - JM / 2
    INTEGER(KIND=I4KIND), PARAMETER :: ITRMX = 0
!
    REAL   (KIND=R4KIND), PARAMETER :: WC    =  2. / 3.
    REAL   (KIND=R4KIND), PARAMETER :: RWCQ  = (1. - WC) * 0.25
    REAL   (KIND=R4KIND), PARAMETER :: RWC   =  1. / WC
!
    INTEGER(KIND=I4KIND), PARAMETER :: KSMUD = 7
    INTEGER(KIND=I4KIND), PARAMETER :: LNSDT = 7
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & PRET    , RPSL    ,                                                                         &
    & FNE     , FSE     ,                                                                         &
    & HBMS    ,                                                                                   &
    & TTB     ,                                                                                   &
    & APDT    , PPDT    ,                                                                         &
    & TPM
!------------------------------------------------ 
! THE FOLLOWING ARE USED FOR TIMIMG PURPOSES ONLY
!------------------------------------------------ 
    REAL   (KIND=R8KIND)                                                                        ::&
    & TIMEF
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K       ,                                                               &
    & IX      , JX      ,                                                                         &
    & KS      ,                                                                                   &
    & NSMUD   ,                                                                                   &
    & JHL     , JHH     ,                                                                         &
    & IHL     , IHH
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & DWDTP   , TPMP    , TTAL    , RHS     , ETADTL  , BTIM
!
    CALL ZERO2(PDSLO)
!------------------------------------------------
! COMPUTATION OF PRESSURE TENDENCY & PREPARATIONS 
!------------------------------------------------
    DO 100 K=2,LM
!
!$omp parallel do  
!  
        DO J=MYJS_P2,MYJE_P2
            DO I=MYIS_P2,MYIE_P2
                DIV(I,J,K) = DIV(I,J,K-1) + DIV(I,J,K)
            END DO
        END DO
!    
100 END DO
!
!$omp parallel do  
!  
    DO J=MYJS_P2,MYJE_P2
        DO I=MYIS_P2,MYIE_P2
             PSDT(I,J) = -DIV(I,J,LM)
             APDT(I,J) = PSDT(I,J)
             PPDT(I,J) = PSDT(I,J)
            PDSLO(I,J) = PDSL(I,J)
             RPSL(I,J) = 1. / PDSL(I,J)
        END DO
    END DO
!
!$omp parallel do  
! 
    DO J=MYJS_P2,MYJE_P2
        DO I=MYIS_P2,MYIE_P2
            PRET(I,J) = PSDT(I,J) * RES(I,J)
            PDSL(I,J) =   PD(I,J) * RES(I,J)
!        
            PINT(I,J,1) = PT
!        
            TPM(I,J) = PT + PINT(I,J,2)
            TTB(I,J) = 0.
        END DO
    END DO
!--------------------- 
! COMPUTATION OF ETADT 
!--------------------- 
!
!$omp parallel do  
!  
    DO 300 K=1,LM1
!   
        DO J=MYJS_P2,MYJE_P2
            DO I=MYIS_P2,MYIE_P2
                ETADT(I,J,K) = - (PRET(I,J)     *  ETA(K+1) +  DIV(I,J,K))                        &
    &                        *     HTM(I,J,K+1) * HBM2(I,J) * RPSL(I,J)
            END DO
        END DO
300 END DO
!---------------------------------------------- 
! KINETIC ENERGY GENERATION TERMS IN T EQUATION 
!---------------------------------------------- 
!
!$omp parallel do private (DWDTP  , RHS     , TPMP    , TTAL)
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            DWDTP = DWDT(I,J,1)
            TPMP  = PINT(I,J,2) + PINT(I,J,3)
!
            TTAL = 0.
!        
            RHS = -DIV(I,J,1) * RTOP(I,J,1) * HTM(I,J,1) * DWDTP * EF4T
!
            OMGALF(I,J,1) = OMGALF(I,J,1) + RHS
                 T(I,J,1) = (TTAL * RDETA(1) + RHS) * HBM2(I,J) + T(I,J,1)
              PINT(I,J,2) = PRET(I,J) * (ETA(1) + ETA(2)) * DWDTP * DT + TPM(I,J) - PINT(I,J,1)
!        
            TPM(I,J) = TPMP
            TTB(I,J) = TTAL
        END DO
    END DO
!
    DO 410 K=2,LM1
!
!$omp parallel do private (DWDTP  , RHS     , TPMP    , TTAL) 
! 
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                DWDTP = DWDT(I,J,K  )
                TPMP  = PINT(I,J,K+1) + PINT(I,J,K+2)
!
                TTAL = 0.
!           
                RHS = -(DIV(I,J,K-1) + DIV(I,J,K)) * RTOP(I,J,K) * HTM(I,J,K) * DWDTP * EF4T
!
                OMGALF(I,J,K  ) =  OMGALF(I,J,K) + RHS
                     T(I,J,K  ) = ((TTAL + TTB(I,J)) * RDETA(L) + RHS) * HBM2(I,J) + T(I,J,K)
                  PINT(I,J,K+1) = PRET(I,J) * (ETA(K) + ETA(K+1)) * DWDTP * DT + TPM(I,J)         &
    &                           - PINT(I,J,K)
            
                TPM(I,J) = TPMP
                TTB(I,J) = TTAL
            END DO
        END DO
    
    410 END DO
!
!$omp parallel do private (DWDTP  , RHS)
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            DWDTP = DWDT(I,J,LM)
!        
            RHS = -(DIV(I,J,LM1) + DIV(I,J,LM)) * RTOP(I,J,LM) * HTM(I,J,LM) * DWDTP * EF4T
!
            OMGALF(I,J,LM  ) = OMGALF(I,J,LM) + RHS
                 T(I,J,LM  ) = (TTB(I,J) * RDETA(LM) + RHS) * HBM2(I,J) + T(I,J,LM)
              PINT(I,J,LM+1) = PRET(I,J) * (ETA(LM) + ETA(LM+1)) * DWDTP * DT + TPM(I,J)          &
    &                        - PINT(I,J,LM)
        END DO
    END DO
!--------------------------------------- 
! REGENERATE THE UNINTEGRATED DIVERGENCE 
!--------------------------------------- 
    DO 425 K=LM,2,-1
!
!$omp parallel do
!
        DO J=MYJS,MYJE2
            DO I=MYIS,MYIE
                DIV(I,J,K) = DIV(I,J,K) - DIV(I,J,K-1)
            END DO
        END DO
!   
425 END DO
!--------------------------------------------- 
! SMOOTHING VERTICAL VELOCITY ALONG BOUNDARIES 
!--------------------------------------------- 
    IF (.NOT. HYDRO .AND. KSMUD > 0) THEN
!
        NSMUD = KSMUD
!
!$omp parallel do
!
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                HBMS(I,J) = HBM2(I,J)
            END DO
        END DO
!    
        JHL = LNSDT
        JHH = JM - JHL + 1
!    
        DO J=JHL,JHH
            IF (J >= MY_JS_GLB .AND. J <= MY_JE_GLB) THEN
                IHL = JHL / 2   + 1
                IHH = IM  - IHL + MOD(J,2)
!            
                DO I=IHL,IHH
                    IF (I >= MY_IS_GLB .AND. I <= MY_IE_GLB) THEN
                        IX = I - MY_IS_GLB + 1
                        JX = J - MY_JS_GLB + 1
                        HBMS(IX,JX) = 0.
                    END IF
                END DO
!            
            END IF
        END DO
!
        DO KS=1,NSMUD
!
!$omp parallel do (ETADTL         , FNE     , FSE)
!
            DO 450 K=1,LM-1
!            
                DO J=MYJS_P1,MYJE1_P1
                    DO I=MYIS_P1,MYIE1_P1
                        FNE(I,J) = (ETADT(I+IHE(J),J+1,K  ) - ETADT(I       ,J  ,K  ))            &
    &                            *    HTM(I       ,J  ,K+1) *   HTM(I+IHE(J),J+1,K+1)
                    END DO
                END DO
!            
                DO J=MYJS1_P1,MYJE_P1
                    DO I=MYIS_P1,MYIE1_P1
                        FSE(I,J) = (ETADT(I+IHE(J),J-1,K  ) - ETADT(I       ,J  ,K  ))            &
    &                            *    HTM(I+IHE(J),J-1,K+1) *   HTM(I       ,J  ,K+1)
                    END DO
                END DO
!            
                DO J=MYJS2,MYJE2
                    DO I=MYIS1,MYIE1
                        ETADTL = (FNE(I,J) - FNE(I+IHW(J),J-1)  +  FSE(I,J)                       &
    &                          -             FSE(I+IHW(J),J+1)) * HBM2(I,J)
!
                        ETADT(I,J,K) = ETADTL * HBMS(I,J) * 0.125 + ETADT(I,J,K)
                    END DO
                END DO
!            
        450 END DO
!        
            BTIM = TIMEF()
!
            CALL EXCH(ETADT, LM-1, 2, 2)
!
            EXCH_TIM = EXCH_TIM + TIMEF() - BTIM
!
        END DO
!
    END IF
! 
    RETURN
!
    END SUBROUTINE PDTEDT
