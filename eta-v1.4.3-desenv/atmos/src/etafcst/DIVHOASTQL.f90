    SUBROUTINE DIVHOASTQL
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE DIVHOASTQL
!>
!> SUBPROGRAM: DIVHOA - DIVERGENCE / HORIZONTAL OMEGA-ALPHA
!> PROGRAMMER: JANJIC
!> ORG: W/NP22
!> DATE: 93-10-28
!>
!> ABSTRACT:
!> DIVHOA COMPUTES THE DIVERGENCE INCLUDING THE MODIFICATION PREVENTING GRAVITY WAVE GRID SEPARATION
!> AND CALCULATES THE HORIZONTAL PART OF THE OMEGA-ALPHA TERM (THE PART PROPORTIONAL TO THE 
!> ADVECTION OF MASS ALONG ETA / SIGMA SURFACES).
!>
!> MODIFIED TO INCLUDE DIVERGENCE RESULTING FROM SLOPING STEPS ("S": DIVHOAST)
!> EXPANDED TO ALSO SLANTWISE T ADVECTION ("T": DIVHOAST), ALONG WITH APPROPRIATE OMEGA-ALPHA
!> CHANGES TO THE ADVECTED T.
!>
!> WARNING: THIS MAKES SUBROUTINE SLADVT REDUNDANT (NEEDS TO BE COMENTED OUT OR REMOVED).
!>
!> PROGRAM HISTORY LOG:
!> 87-06-??  JANJIC             - ORIGINATOR
!> 95-03-25  BLACK              - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 96-03-29  BLACK              - ADDED EXTERNAL EDGE
!> 97-03-17  MESINGER           - SPLIT FROM PFDHT
!> 98-10-30  BLACK              - MODIFIED FOR DISTRIBUTED MEMORY
!> 00-10-20  BLACK              - INCORPORATED PRESSURE GRADIENT METHOD FROM MESO MODEL
!> 06-04-04  MESINGER AND JOVIC - ADDED DIVERGENCE CALCULATION DUE TO SLOPES
!> 06-05-??  MESINGER, AT CPTEC - EXPANDED TO INCLUDE SLANTWISE T ADV.
!> 06-06-??  MESINGER AND JOVIC - REVISED OMEGA-ALPHA DUE TO SLOPES
!> 18-01-15  LUCCI              - MODERNIZATION OF THE CODE, INCLUDING:
!>                                * F77 TO F90/F95
!>                                * INDENTATION & UNIFORMIZATION CODE
!>                                * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                                * DOCUMENTATION WITH DOXYGEN
!>                                * OPENMP FUNCTIONALITY
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
!>              DYNAM0
!>              F77KINDS
!>              GLB_TABLE
!>              INDX
!>              LOOPS
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              NHYDRO
!>              PARMETA
!>              SLOPES
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
    USE NHYDRO
    USE PARMETA
    USE SLOPES
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND), PARAMETER :: RFCP = 1.E0 / (4.E0 * 1004.6E0) 
!
    REAL   (KIND=R8KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM+1)                             ::&
    & PINTLG
!
    REAL   (KIND=R8KIND)                                                                        ::&
    & ETA_DP , PDSL_DP  , PT_DP     , PINT_DP!
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    &  FIM    ,                                                                                   &
    &  FILO   , RDPD    ,                                                                         &
    &  ADPDX  , RDPDX   ,                                                                         &
    &  ADPDY  , RDPDY   ,                                                                         &
    &  ADPDNE , ADPDSE  ,                                                                         &
    &  PEW    , PNS     ,                                                                         &
    &  PCEW   , PCNS    ,                                                                         &
    &  DPFEW  , DPFNS   ,                                                                         &
    &  FNS    , TNS     ,                                                                         &
    &  HM     , VM      ,                                                                         &
    &  EDIV   , DIVL  
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    &   DPDE  ,                                                                                   &
    &   APEL  , PCXC    ,                                                                         &
    &   ALP1  ,                                                                                   &
    &   DFDZ  ,                                                                                   &
    &   UDY   , VDX     ,                                                                         &
    &   TEW   , FEW     ,                                                                         &
    &   TNE   , TSE     ,                                                                         &
    &   FNE   , FSE     ,                                                                         &
    &   PNE   , PSE     ,                                                                         &
    &   CNE   , CSE     ,                                                                         &
    &   PPNE  , PPSE    ,                                                                         &
    &   PCNE  , PCSE    ,                                                                         &
    &   DIVS  , DPDEP1  ,                                                                         &
    &   ATC   , ATCP1   ,                                                                         &
!-----------------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER CONTRIBUTIONS
!-----------------------------------------------------------------
    & AQC     , AQCP1   ,                                                                         &
    & ACWMC   , ACWMCP1  
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K 
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ALP1P   , ALP1X   , DFI     , RDPDS   , FIUPK   , ALP1PL  , ALP2P   , ALP2PL  , DPNEK   ,   &
    & DPSEK   , DCNEK   , DCSEK   , PVNEK   , PVSEK   , TOASL   , OAADP1  , TFSTR   , FSTR    ,   &
    & STSEP1  , STNE    , STSE    , STNEM1  , STEW    , STNS    , WDEP    , DPT     , WDUTD   ,   &
    & WARR    , WDUCWMD , WDUQD
!
    CALL ZERO2(ALP1)
    CALL ZERO2(DPDE)
    CALL ZERO2(APEL)
    CALL ZERO2(PCXC)
    CALL ZERO2(DFDZ)
    CALL ZERO2(UDY)
    CALL ZERO2(VDX)
    CALL ZERO2(TEW)
    CALL ZERO2(FEW)
    CALL ZERO2(TNE)
    CALL ZERO2(TSE)
    CALL ZERO2(FNE)
    CALL ZERO2(FSE)
    CALL ZERO2(PNE)
    CALL ZERO2(PSE)
    CALL ZERO2(CNE)
    CALL ZERO2(CSE)
    CALL ZERO2(PPNE)
    CALL ZERO2(PPSE)
    CALL ZERO2(PCNE)
    CALL ZERO2(PCSE)
!------------------------- 
! PREPARATORY CALCULATIONS 
!------------------------- 
    IF (SIGMA) THEN
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO J=MYJS_P5,MYJE_P5
            DO I=MYIS_P5,MYIE_P5
                FILO(I,J) = FIS(I,J)
                PDSL(I,J) =  PD(I,J)
            END DO
        END DO
!
    ELSE
!------- 
! OPENMP
!-------
! 
!$omp parallel do
!
        DO J=MYJS_P5,MYJE_P5
            DO I=MYIS_P5,MYIE_P5
                FILO(I,J) = 0.0
                PDSL(I,J) = RES(I,J) * PD(I,J)
            END DO
        END DO
!
    END IF
!
    IF (HYDRO) THEN
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO K=1,LM+1
            DO J=MYJS_P5,MYJE_P5
                DO I=MYIS_P5,MYIE_P5
                    PINTLG(I,J,K) = ALOG(ETA(K) * PDSL(I,J) + PT)
                END DO
            END DO
        END DO
!
!
    ELSE
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO K=1,LM+1
            DO J=MYJS_P5,MYJE_P5
                DO I=MYIS_P5,MYIE_P5
                    PINTLG(I,J,K) = ALOG(PINT(I,J,K))
                END DO
            END DO
        END DO
!
    END IF
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
    DO K=1,LM
        DO J=MYJS_P4,MYJE_P4
            DO I=MYIS_P4,MYIE_P4
                OMGALF(I,J,K) = 0.
                   DIV(I,J,K) = 0.
            END DO
        END DO
    END DO
!
    DO J=MYJS_P4,MYJE_P4
        DO I=MYIS_P4,MYIE_P4
               DIVS(I,J) = 0.
!--------------------------------------------------------------------------------------------------
! DIVERGENCE TO SAVE IS NEEDED TO COLLECT CONTRIBUTIONS TO DIV DUE TO SLOPES AT H POINTS OF THE 
! LAYER BEING PROCESSED (L).
! THEY CAN BE ADDED ONLY ONCE THE CODE MOVES ONE LAYER UP.
! THE REASON IS THAT THE COLLECTION IS DONE FOR POINTS ALSO SIDEWAYS OF THE H POINT THAT THE CODE
! IS PROCESSING, SO THAT IF ADDED IMMEDIATELLY THESE WOULD BE LOST WHEN THE CODE VISITS THESE 
! SIDEWAYS LOCATED POINTS, AND PUTS THE DIVERGENCE "CORRECTION"INTO DIV (LOOP 270)
!
! ACCUMULATED T CHANGES AT L AND AT L+1
!--------------------------------------------------------------------------------------------------
                ATC(I,J) = 0.0
              ATCP1(I,J) = 0.0
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                AQC(I,J) = 0.0
              AQCP1(I,J) = 0.0
              ACWMC(I,J) = 0.0
            ACWMCP1(I,J) = 0.0
        END DO
    END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (ALP1X)
!
    DO J=MYJS_P5,MYJE_P5
        DO I=MYIS_P5,MYIE_P5
            ALP1X     = PINTLG(I,J,LM+1)
            ALP1(I,J) = ALP1X
        END DO
    END DO
!------------------------------- 
! MAIN VERTICAL INTEGRATION LOOP  
!------------------------------- 
    DO K=LM,1,-1
!--------------------------- 
! INTEGRATE THE GEOPOTENTIAL
!--------------------------- 
        FIM = 0.  
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (ALP1P  , DFI     , FIUPK   , RDPDS)
!
        DO J=MYJS_P5,MYJE_P5
            DO I=MYIS_P5,MYIE_P5
!            
                ALP1P = PINTLG(I,J,K)
!            
!CHOU CWM in Tv 20210401
                DFI =  (Q(I,J,K) * 0.608 + 1.-CWM(I,J,K)) * T(I,J,K) * R * (ALP1(I,J) - ALP1P)    &
    &                 / DWDT(I,J,K)
!           
                RDPDS       = 1. / (DETA(K) * PDSL(I,J))
                RTOP(I,J,K) = RDPDS * DFI
                FIUPK       = FILO(I,J) + DFI
                 FIM(I,J)   = FILO(I,J) + FIUPK
                IF (ABS(FIM(I,J)) <= 5.E+10) THEN
!
                ELSE
!
                    WRITE(6,*) 'BAD FIM ', I, J, FIM(I,J), FILO(I,J), DFI
                    WRITE(6,*) 'Q,T,ALP1,ALP1P,DWDT: ',   Q(I,J,K), T(I,J,K), ALP1(I,J), ALP1P,   &
    &                                                  DWDT(I,J,K)
                    STOP
!
                END IF
!            
                FILO(I,J) = (FIUPK - DFL(K)) * HTM(I,J,K) + DFL(K)
                ALP1(I,J) = ALP1P
!
            END DO
        END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (ALP1P  , ALP1PL  , ALP2P   , ALP2PL  , DFI)
!
        DO J=MYJS_P5,MYJE_P5
            DO I=MYIS_P5,MYIE_P5
                HM(I,J) = HTM(I,J,K) * HBM2(I,J)
                VM(I,J) = VTM(I,J,K) * VBM2(I,J)
!            
                ALP1P  = PINTLG(I,J,K)
                ALP1PL = PINTLG(I,J,K+1)
                ALP2P  = ALP1P  * ALP1P
                ALP2PL = ALP1PL * ALP1PL
!            
!CHOU CWM in Tv, but no need (1 + CWM(I,J,K))
!CHOU                 DFI       = ((Q(I,J,K) * 0.608 + 1. - CWM(I,J,K)) * T(I,J,K)                      &
!CHOU     &                     * R * (ALP1PL - ALP1P) / (1 + CWM(I,J,K))) / DWDT(I,J,K)
                DFI       = ((Q(I,J,K) * 0.608 + 1. - CWM(I,J,K)) * T(I,J,K)                      &
    &                     * R * (ALP1PL - ALP1P) ) / DWDT(I,J,K)

                DFDZ(I,J) = DFI * DWDT(I,J,K) / (ALP2PL - ALP2P)
                APEL(I,J) = (ALP2PL + ALP2P) * 0.5
            END DO
        END DO
!------- 
! OPENMP
!-------
! 
!$omp parallel do
!
        DO J=MYJS_P4,MYJE_P4
            DO I=MYIS_P4,MYIE_P4
                DPDE(I,J) = DETA(K) * PDSL(I,J)
                DIVL(I,J) = 0.
                EDIV(I,J) = 0.
            END DO
        END DO
!    
        IF (K < LM) THEN
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
            DO J=MYJS_P4,MYJE_P4
                DO I=MYIS_P4,MYIE_P4
                    DPDEP1(I,J) = DETA(K+1) * PDSL(I,J)
                END DO
            END DO
!
        END IF
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO J=MYJS_P1,MYJE_P1
            DO I=MYIS_P1,MYIE_P1
                RDPD(I,J) = 1. / DPDE(I,J)
            END DO
        END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO J=MYJS1_P3,MYJE1_P3
            DO I=MYIS_P3,MYIE_P3
                ADPDX(I,J) = DPDE(I+IVW(J),J  ) + DPDE(I+IVE(J),J  )
                ADPDY(I,J) = DPDE(I       ,J-1) + DPDE(I       ,J+1)
                RDPDX(I,J) = 1. / ADPDX(I,J)
                RDPDY(I,J) = 1. / ADPDY(I,J)
            END DO
        END DO
!--------------------------------------------------   
! DIAGONAL CONTRIBUTIONS TO PRESSURE GRADIENT FORCE 
!--------------------------------------------------    
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO J=MYJS_P4,MYJE_P4
            DO I=MYIS_P4,MYIE_P4
                ADPDNE(I,J) =  DPDE(I+IHE(J),J+1)   + DPDE(I,J)
                   PNE(I,J) =  (FIM(I+IHE(J),J+1)   -  FIM(I,J))                                  &
    &                       * (DWDT(I+IHE(J),J+1,K) + DWDT(I,J,K))
!JLG Utilizado para debug  na ocorrencia de NaN 
!            IF (ieee_is_nan(PNE(I,J))) THEN
!               WRITE(6,*)MYPE,I,J,K,I+IHE(J), LMH(I,J), FIM(I+IHE(J),J+1), FIM(I,J), DWDT(I+IHE(J),J+1,K), DWDT(I,J,K)
!               STOP
!            END IF
!
                IF (ABS(PNE(I,J)) <= 5.E10) THEN
!
                ELSE
                    WRITE(6,*) 'CRAZY PNE ', I, J, PNE(I,J)
                    WRITE(6,*) 'PIECES', I+IHE(J), J+1, FIM(I+IHE(J),J+1)
                END IF
!
                PPNE(I,J) = PNE(I,J) * ADPDNE(I,J)
                 CNE(I,J) =      (DFDZ(I+IHE(J),J+1) + DFDZ(I,J))                                 &
    &                     * 2. * (APEL(I+IHE(J),J+1) - APEL(I,J))
                PCNE(I,J) = CNE(I,J) * ADPDNE(I,J)
            END DO
        END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO J=MYJS1_P4,MYJE_P4
            DO I=MYIS_P4,MYIE1_P4
                ADPDSE(I,J) =  DPDE(I+IHE(J),J-1)   + DPDE(I,J)
                   PSE(I,J) =  (FIM(I+IHE(J),J-1)   -  FIM(I,J))                                  &
    &                       * (DWDT(I+IHE(J),J-1,K) + DWDT(I,J,K))
!
                  PPSE(I,J) = PSE(I,J) * ADPDSE(I,J)
                   CSE(I,J) =      (DFDZ(I+IHE(J),J-1) + DFDZ(I,J))                               &
    &                       * 2. * (APEL(I+IHE(J),J-1) - APEL(I,J))
                  PCSE(I,J) = CSE(I,J) * ADPDSE(I,J)
            END DO
        END DO
!---------------------------------    
! CONTINUITY EQUATION MODIFICATION
!-------------------------------
!------- 
! OPENMP
!-------
! 
!$omp parallel do
!
        DO J=MYJS1_P1,MYJE1_P1
            DO I=MYIS_P1,MYIE_P1
                PCXC(I,J) = VBM3(I,J) * VTM(I,J,K)                                                &
    &                     * (PNE(I+IVW(J),J) + CNE(I+IVW(J),J  ) +  PSE(I+IVW(J),J  )             &
    &                     +  CSE(I+IVW(J),J) - PNE(I       ,J-1) -  CNE(I       ,J-1)             &
    &                     -                    PSE(I       ,J+1) -  CSE(I       ,J+1))
            END DO
        END DO
!
        DO J=MYJS2,MYJE2
            DO I=MYIS1,MYIE1
                DIV(I,J,K) = DETA(K) * WPDAR(I,J) * (PCXC(I+IHE(J),J) - PCXC(I,J+1)               &
    &                      +                         PCXC(I+IHW(J),J) - PCXC(I,J-1))
            END DO
        END DO
!---------------------------------------    
! LAT AND LONG PRESSURE FORCE COMPONENTS
!---------------------------------------
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (DCNEK  , DCSEK   , DPNEK   , DPSEK)
!
        DO J=MYJS1_P3,MYJE1_P3
            DO I=MYIS_P3,MYIE_P3
                DPNEK     = PNE(I+IVW(J),J) +   PNE(I,J-1)
                DPSEK     = PSE(I+IVW(J),J) +   PSE(I,J+1)
                PEW(I,J)  = DPNEK + DPSEK
                PNS(I,J)  = DPNEK - DPSEK
                DCNEK     = CNE(I+IVW(J),J) +   CNE(I,J-1)
                DCSEK     = CSE(I+IVW(J),J) +   CSE(I,J+1)
                PCEW(I,J) = (DCNEK + DCSEK) * ADPDX(I,J)
                PCNS(I,J) = (DCNEK - DCSEK) * ADPDY(I,J)
            END DO
        END DO
!-----------------------------------------------    
! LAT AND  LON FLUXES AND OMEGA-ALPHA COMPONENTS
!-----------------------------------------------
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO J=MYJS1_P3,MYJE1_P3
            DO I=MYIS_P3,MYIE_P3
                UDY(I,J) = DY       *     U(I,J,K)
                FEW(I,J) = UDY(I,J) * ADPDX(I,J)
                TEW(I,J) = UDY(I,J) *  PCEW(I,J)
                VDX(I,J) =  DX(I,J) *     V(I,J,K)
                FNS(I,J) = VDX(I,J) * ADPDY(I,J)
                TNS(I,J) = VDX(I,J) *  PCNS(I,J)
            END DO
        END DO
!---------------------------------------------    
! DIAGONAL FLUXES AND DIAGONALLY AVERAGED WIND 
!---------------------------------------------    
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (PVNEK)
!
        DO J=MYJS1_P2,MYJE2_P2
            DO I=MYIS_P2,MYIE1_P2
                PVNEK    = (UDY(I+IHE(J),J) + VDX(I+IHE(J),J)) + (UDY(I,J+1) + VDX(I,J+1))
                FNE(I,J) = PVNEK * ADPDNE(I,J)
                TNE(I,J) = PVNEK *   PCNE(I,J) * 2.
            END DO
        END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (PVSEK)
!
        DO J=MYJS2_P2,MYJE1_P2
            DO I=MYIS_P2,MYIE1_P2
                PVSEK    = (UDY(I+IHE(J),J) - VDX(I+IHE(J),J)) + (UDY(I,J-1) - VDX(I,J-1))
                FSE(I,J) = PVSEK * ADPDSE(I,J)
                TSE(I,J) = PVSEK *   PCSE(I,J) * 2.
            END DO
        END DO
!
        IF (K < LM-1) THEN
            DO J=MYJS2_P2,MYJE2_P2
                DO I=MYIS1_P2,MYIE1_P2
                      DIV(I,J,K+1) = DIV(I,J,K+1) +    DIVS(I,J)
                     DIVS(I,J)     = 0.0
                        T(I,J,K+2) =   T(I,J,K+2) +   ATCP1(I,J)
                        Q(I,J,K+2) =   Q(I,J,K+2) +   AQCP1(I,J)   ! CHOU
                      CWM(I,J,K+2) = CWM(I,J,K+2) + ACWMCP1(I,J)   ! CHOU
                    ATCP1(I,J)     = ATC(I,J)
                      ATC(I,J)     = 0.0
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                      AQCP1(I,J) =   AQC(I,J) ! CHOU
                        AQC(I,J) = 0.0        ! CHOU
                    ACWMCP1(I,J) = ACWMC(I,J) ! CHOU
                      ACWMC(I,J) = 0.0        ! CHOU
                END DO
            END DO
        END IF
!-----------------------------------------------    
! HORIZONTAL PART OF OMEGA-ALPHA AND DIVERGENCE 
!-----------------------------------------------
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO J=MYJS2_P1,MYJE2_P1
            DO I=MYIS1_P1,MYIE1_P1
                OMGALF(I,J,K) = (TEW(I+IHE(J),J  ) + TEW(I+IHW(J),J  )                            &
    &                         +  TNS(I       ,J+1) + TNS(I       ,J-1)                            &
    &                         +  TNE(I       ,J  ) + TNE(I+IHW(J),J-1)                            &
    &                         +  TSE(I       ,J  ) + TSE(I+IHW(J),J+1))                           & 
    &                         * RDPD(I,J) * FCP(I,J) * HM(I,J)
!
                IF (K < 3) THEN
                      T(I,J,K) = OMGALF(I,J,K) +      T(I,J,K)
                ELSE IF (K < LM) THEN
                    ATC(I,J)   =    ATC(I,J)   + OMGALF(I,J,K)
                END IF
!            
                IF (K < LM .AND. (HTM(I,J,K+1) * HBM2(I,J)) > 0) THEN
!--------------------------------------------------------------------------------------------------
! THE T POINT BELOW IS IN THE ATMOSPHERE. SHOULD IT RECEIVE OMEGA-ALPHA CONTRIBUTIONS FROM FLUXES 
! OF LAYER L ?
!
! CHECK EACH OF THE FOUR SURROUNDING V POINTS IF THEY ARE THE POINTS JUST ABOVE THE SLOPE
!--------------------------------------------------------------------------------------------------
                    TOASL = 0.0
!
                    IF (VTMS(I+IHE(J),J,K+1) == 1) THEN
                        IF (ISLD(I+IHE(J),J) == 1)                                                &
    &                       TOASL = TOASL + 0.50 * TEW(I+IHE(J),J)
!
                        IF (ISLD(I+IHE(J),J) == 8)                                                &
    &                       TOASL = TOASL + 0.50 * TEW(I+IHE(J),J) + 0.25 * TSE(I,J)
!
                        IF (ISLD(I+IHE(J),J) == 2)                                                &
    &                       TOASL = TOASL + 0.50 * TEW(I+IHE(J),J) + 0.25 * TNE(I,J)
                    END IF
!
                    IF (VTMS(I,J+1,K+1) == 1) THEN
                        IF (ISLD(I,J+1) == 3)                                                     &
    &                       TOASL = TOASL + 0.50 * TNS(I,J+1)
!
                        IF (ISLD(I,J+1) == 2)                                                     &
    &                       TOASL = TOASL + 0.50 * TNS(I,J+1) + 0.25 * TNE(I       ,J  )
!
                        IF (ISLD(I,J+1) == 4)                                                     &
    &                       TOASL = TOASL + 0.50 * TNS(I,J+1) + 0.25 * TSE(I+IHW(J),J+1)
                    END IF
!
                    IF (VTMS(I+IHW(J),J,K+1) == 1) THEN
                        IF (ISLD(I+IHW(J),J) == 5)                                                &
    &                       TOASL = TOASL + 0.50 * TEW(I+IHW(J),J)
!
                        IF (ISLD(I+IHW(J),J) == 4)                                                &
    &                       TOASL = TOASL + 0.50 * TEW(I+IHW(J),J) + 0.25 * TSE(I+IHW(J),J+1)
                        IF (ISLD(I+IHW(J),J) == 6)                                                &
    &                       TOASL = TOASL + 0.50 * TEW(I+IHW(J),J) + 0.25 * TNE(I+IHW(J),J-1)
                    END IF
!
                    IF (VTMS(I,J-1,K+1) == 1) THEN
                        IF (ISLD(I,J-1) == 7)                                                     &
    &                       TOASL = TOASL + 0.50 * TNS(I,J-1)
!
                        IF (ISLD(I,J-1) == 6)                                                     &
    &                       TOASL = TOASL + 0.50 * TNS(I,J-1) + 0.25 * TNE(I+IHW(J),J-1)
!
                        IF (ISLD(I,J-1) == 8)                                                     &
    &                       TOASL = TOASL + 0.50 * TNS(I,J-1) + 0.25 * TSE(I       ,J  )
                    END IF
!
                    OAADP1          = TOASL * DETA(K+1) / DETA(K) * RDPD(I,J) * FCP(I,J)
                    OMGALF(I,J,K+1) = OMGALF(I,J,K+1) + OAADP1
!
                    IF (K < 3) THEN
                            T(I,J,K+1) =     T(I,J,K+1) + OAADP1
                    ELSE IF (K == LM-1) THEN
                        ATCP1(I,J)     = ATCP1(I,J)     + OMGALF(I,J,K+1)
                    ELSE
                        ATCP1(I,J)     = ATCP1(I,J)     + OAADP1
                    END IF
!
                END IF
!            
                EDIV(I,J) = ((FEW(I+IHE(J),J) + FNS(I,J+1) + FNE(I,J)                             &
    &                     +                                  FSE(I,J))                            &
    &                     -  (FEW(I+IHW(J),J) + FNS(I,J-1) + FNE(I+IHW(J),J-1)                    &
    &                     +                                  FSE(I+IHW(J),J+1)))                  &
    &                     *   FDIV(I,J)
!
                DIVL(I,J) =   EDIV(I,J) * HBM2(I,J)
!--------------------------------------------------------------------------------------------------
! SLANTWISE MASS AND T ADVECTION, WITH T INCREMENT TO ACCOUNT FOR THE OMEGA-ALPHA IMPACT, T CHANGE 
! DUE TO MOVEMENT TO A DIFFERENT P
!--------------------------------------------------------------------------------------------------
                IF (K < LM .AND. K == LMV(I+IHE(J),J) .AND. ISLD(I+IHE(J),J) > 0) THEN
!
                    TFSTR = 0.5 * DETA(K+1) / DETA(K) * FDIV(I,J) * HM(I,J)
                    FSTR  = 0.5 * TFSTR
!
                    STSEP1 = FSTR  * FSE(I+IHE(J),J+1)
                    STNE   = FSTR  * FNE(I       ,J  )
                    STSE   = FSTR  * FSE(I       ,J  )
                    STNEM1 = FSTR  * FNE(I+IHE(J),J-1)
                    STEW   = TFSTR * FEW(I+IHE(J),J  )
                    STNS   = TFSTR * FNS(I+IHE(J),J  )
!
                    IF  (ISLD(I+IHE(J),J) == 1) THEN
                         DIV(I+IHE(J),J+1,K+1) =  DIV(I+IHE(J),J+1,K+1) + STSEP1
                        DIVS(I+1,J)            = DIVS(I+1,J)            - STSEP1
!
                        WDEP = ABS(STSEP1) * DT
!
                        DPT =   (RTOP(I+IHE(J),J+1,K+1) + RTOP(I+1,J,K))                          &
    &                       * (DPDEP1(I+IHE(J),J+1)     + DPDE(I+1,J))                            &
    &                       *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,STSEP1) * (T(I+IHE(J),J+1,K+1) - DPT              &
    &                         -                            T(I+1,J,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0, STSEP1) *   (Q(I+IHE(J),J+1,K+1) -   Q(I+1,J,K))
                        WDUCWMD = WDEP * SIGN(1.0, STSEP1) * (CWM(I+IHE(J),J+1,K+1) - CWM(I+1,J,K))
!
                        IF ( STSEP1 > 0.0) THEN
                            WARR = DPDE(I+1,J)
                              ATC(I+1,J)        =   ATC(I+1,J) + WDUTD   / (WARR+WDEP)
                              AQC(I+1,J)        =   AQC(I+1,J) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+1,J)        = ACWMC(I+1,J) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR = DPDEP1(I+IHE(J),J+1)
                              ATCP1(I+IHE(J),J+1) =   ATCP1(I+IHE(J),J+1) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+IHE(J),J+1) =   AQCP1(I+IHE(J),J+1) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+IHE(J),J+1) = ACWMCP1(I+IHE(J),J+1) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                         DIV(I+IHE(J),J-1,K+1) =  DIV(I+IHE(J),J-1,K+1) + STNEM1
                        DIVS(I+1,J)            = DIVS(I+1,J)            - STNEM1
!
                        WDEP  = ABS(STNEM1) * DT
                        DPT   =   (RTOP(I+IHE(J),J-1,K+1) + RTOP(I+1,J,K))                        &
    &                         * (DPDEP1(I+IHE(J),J-1)     + DPDE(I+1,J))                          &
    &                         *  RFCP
                        WDUTD = WDEP * SIGN(1.0,STNEM1) * (T(I+IHE(J),J-1,K+1) - DPT              &
    &                         -                            T(I+1,J,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0, STNEM1) * (  Q(I+IHE(J),J-1,K+1) -   Q(I+1,J,K))
                        WDUCWMD = WDEP * SIGN(1.0, STNEM1) * (CWM(I+IHE(J),J-1,K+1) - CWM(I+1,J,K))

                        IF (STNEM1 > 0.0) THEN
                            WARR         =  DPDE(I+1,J)
                              ATC(I+1,J) =   ATC(I+1,J) + WDUTD   / (WARR+WDEP)
                              AQC(I+1,J) =   AQC(I+1,J) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+1,J) = ACWMC(I+1,J) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR                  =  DPDEP1(I+IHE(J),J-1)
                              ATCP1(I+IHE(J),J-1) =   ATCP1(I+IHE(J),J-1) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+IHE(J),J-1) =   AQCP1(I+IHE(J),J-1) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+IHE(J),J-1) = ACWMCP1(I+IHE(J),J-1) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                         DIV(I,J,K+1) =  DIV(I,J,K+1) + STEW !BLUE
                        DIVS(I+1,J)   = DIVS(I+1,J)   - STEW !BLUE
!
                        WDEP  = ABS(STEW) * DT
                        DPT   =   (RTOP(I,J,K+1) + RTOP(I+1,J,K))                                 & 
    &                         * (DPDEP1(I,J)     + DPDE(I+1,J))                                   &
    &                         *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,STEW) * (T(I,J,K+1) - DPT                         &
    &                         -                          T(I+1,J,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0, STEW) * (  Q(I,J,K+1) -   Q(I+1,J,K))
                        WDUCWMD = WDEP * SIGN(1.0, STEW) * (CWM(I,J,K+1) - CWM(I+1,J,K))
!
                        IF (STEW > 0.0) THEN
                            WARR         =  DPDE(I+1,J)
                              ATC(I+1,J) =   ATC(I+1,J) + WDUTD   / (WARR+WDEP)
                              AQC(I+1,J) =   AQC(I+1,J) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+1,J) = ACWMC(I+1,J) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR         =  DPDEP1(I,J)
                              ATCP1(I,J) =   ATCP1(I,J) + WDUTD   / (WARR+WDEP)
                              AQCP1(I,J) =   AQCP1(I,J) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I,J) = ACWMCP1(I,J) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                    ELSE IF (ISLD(I+IHE(J),J) == 2) THEN
                         DIV(I,J,K+1)      =   DIV(I,J,K+1)      + STNE
                        DIVS(I+IHE(J),J+1) =  DIVS(I+IHE(J),J+1) - STNE
!
                        WDEP = ABS(STNE) * DT
                        DPT  =   (RTOP(I,J,K+1) + RTOP(I+IHE(J),J+1,K))                           &
    &                        * (DPDEP1(I,J)     + DPDE(I+IHE(J),J+1))                             &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,STNE) * (T(I,J,K+1) - DPT                         &
    &                         -                          T(I+IHE(J),J+1,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,STNE) * (  Q(I,J,K+1) -   Q(I+IHE(J),J+1,K))
                        WDUCWMD = WDEP * SIGN(1.0,STNE) * (CWM(I,J,K+1) - CWM(I+IHE(J),J+1,K))
!
                        IF (STNE > 0.0) THEN
                            WARR                =  DPDE(I+IHE(J),J+1)
                              ATC(I+IHE(J),J+1) =   ATC(I+IHE(J),J+1) + WDUTD   / (WARR+WDEP)
                              AQC(I+IHE(J),J+1) =   AQC(I+IHE(J),J+1) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+IHE(J),J+1) = ACWMC(I+IHE(J),J+1) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR         =  DPDEP1(I,J) 
                              ATCP1(I,J) =   ATCP1(I,J) + WDUTD   / (WARR+WDEP)
                              AQCP1(I,J) =   AQCP1(I,J) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I,J) = ACWMCP1(I,J) + WDUCWMD / (WARR+WDEP)
                        END IF
!                       		   
                         DIV(I+IHE(J),J-1,K+1) =  DIV(I+IHE(J),J-1,K+1) + STNEM1
                        DIVS(I+1,J)            = DIVS(I+1,J)            - STNEM1
!
                        WDEP = ABS(STNEM1) * DT
                        DPT  =   (RTOP(I+IHE(J),J-1,K+1) + RTOP(I+1,J,K))                         &
    &                        * (DPDEP1(I+IHE(J),J-1)     + DPDE(I+1,J))                           &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,STNEM1) * (T(I+IHE(J),J-1,K+1) - DPT              &
    &                         -                            T(I+1,J,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,STNEM1) * (  Q(I+IHE(J),J-1,K+1) -  Q(I+1,J,K))
                        WDUCWMD = WDEP * SIGN(1.0,STNEM1) * (CWM(I+IHE(J),J-1,K+1) -CWM(I+1,J,K))
!
                        IF (STNEM1 > 0.0) THEN
                            WARR         =  DPDE(I+1,J)
                              ATC(I+1,J) =   ATC(I+1,J) + WDUTD   / (WARR+WDEP)
                              AQC(I+1,J) =   AQC(I+1,J) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+1,J) = ACWMC(I+1,J) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR                  =  DPDEP1(I+IHE(J),J-1)
                              ATCP1(I+IHE(J),J-1) =   ATCP1(I+IHE(J),J-1) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+IHE(J),J-1) =   AQCP1(I+IHE(J),J-1) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+IHE(J),J-1) = ACWMCP1(I+IHE(J),J-1) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                         DIV(I,J,K+1)   =  DIV(I,J,K+1)   + STEW !BLUE         
                        DIVS(I+1,J)     = DIVS(I+1,J)     - STEW !BLUE
!
                        WDEP = ABS(STEW) * DT
                        DPT  =   (RTOP(I,J,K+1) + RTOP(I+1,J,K))                                  &
    &                        * (DPDEP1(I,J)     + DPDE(I+1,J))                                    &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,STEW) * (T(I,J,K+1) - DPT                         &
    &                         -                          T(I+1,J,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,STEW) * (  Q(I,J,K+1) -   Q(I+1,J,K))
                        WDUCWMD = WDEP * SIGN(1.0,STEW) * (CWM(I,J,K+1) - CWM(I+1,J,K))

!
                        IF (STEW > 0.0) THEN
                            WARR         =  DPDE(I+1,J)
                              ATC(I+1,J) =   ATC(I+1,J) + WDUTD   / (WARR+WDEP)
                              ATC(I+1,J) =   ATC(I+1,J) + WDUTD   / (WARR+WDEP)
                              AQC(I+1,J) =   AQC(I+1,J) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+1,J) = ACWMC(I+1,J) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR         =  DPDEP1(I,J)
                              ATCP1(I,J) =   ATCP1(I,J) + WDUTD   / (WARR+WDEP)
                              AQCP1(I,J) =   AQCP1(I,J) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I,J) = ACWMCP1(I,J) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                         DIV(I+IHE(J),J-1,K+1) =  DIV(I+IHE(J),J-1,K+1) + STNS !YELLOW
                        DIVS(I+IHE(J),J+1)     = DIVS(I+IHE(J),J+1)     - STNS !YELLOW
!
                        WDEP = ABS(STNS) * DT
                        DPT  =   (RTOP(I+IHE(J),J-1,K+1) + RTOP(I+IHE(J),J+1,K))                  &
    &                        * (DPDEP1(I+IHE(J),J-1)     + DPDE(I+IHE(J),J+1))                    &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,STNS) * (T(I+IHE(J),J-1,K+1) - DPT                &
    &                         -                          T(I+IHE(J),J+1,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,STNS) * (  Q(I+IHE(J),J-1,K+1)                  &
    &                           -                            Q(I+IHE(J),J+1,K  ))
!
                        WDUCWMD = WDEP * SIGN(1.0,STNS) * (CWM(I+IHE(J),J-1,K+1)                  &
    &                           -                          CWM(I+IHE(J),J+1,K  ))
!
                        IF (STNS > 0.0) THEN
                            WARR                =  DPDE(I+IHE(J),J+1)
                              ATC(I+IHE(J),J+1) =   ATC(I+IHE(J),J+1) + WDUTD   / (WARR+WDEP)
                              AQC(I+IHE(J),J+1) =   AQC(I+IHE(J),J+1) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+IHE(J),J+1) = ACWMC(I+IHE(J),J+1) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR                  =  DPDEP1(I+IHE(J),J-1)
                              ATCP1(I+IHE(J),J-1) =   ATCP1(I+IHE(J),J-1) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+IHE(J),J-1) =   AQCP1(I+IHE(J),J-1) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+IHE(J),J-1) = ACWMCP1(I+IHE(J),J-1) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                    ELSE IF (ISLD(I+IHE(J),J) == 3) THEN
                         DIV(I+1,J,K+1)    =  DIV(I+1,J,K+1)    - STSEP1
                        DIVS(I+IHE(J),J+1) = DIVS(I+IHE(J),J+1) + STSEP1
!
                        WDEP = ABS(-STSEP1) * DT
                        DPT  =   (RTOP(I+1,J,K+1) + RTOP(I+IHE(J),J+1,K))                         &
    &                        * (DPDEP1(I+1,J)     + DPDE(I+IHE(J),J+1))                           &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,-STSEP1) * (T(I+1,J,K+1)       - DPT              &
    &                         -                             T(I+IHE(J),J+1,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,-STSEP1) * (  Q(I+1,J,K+1) -   Q(I+IHE(J),J+1,K))
                        WDUCWMD = WDEP * SIGN(1.0,-STSEP1) * (CWM(I+1,J,K+1) - CWM(I+IHE(J),J+1,K))

                        IF (-STSEP1 > 0.0) THEN
                            WARR                = DPDE(I+IHE(J),J+1)
                              ATC(I+IHE(J),J+1) =  ATC(I+IHE(J),J+1) + WDUTD   / (WARR+WDEP)
                              AQC(I+IHE(J),J+1) =  AQC(I+IHE(J),J+1) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+IHE(J),J+1) =ACWMC(I+IHE(J),J+1) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR           =  DPDEP1(I+1,J)
                              ATCP1(I+1,J) =   ATCP1(I+1,J) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+1,J) =   AQCP1(I+1,J) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+1,J) = ACWMCP1(I+1,J) + WDUCWMD / (WARR+WDEP)
                        END IF
!                        		   
                        DIV(I,J,K+1)       = DIV(I,J,K+1)       + STNE
                        DIVS(I+IHE(J),J+1) = DIVS(I+IHE(J),J+1) - STNE
!
                        WDEP = ABS(STNE) * DT
                        DPT  = (RTOP(I,J,K+1) + RTOP(I+IHE(J),J+1,K))                             &
    &                        * (DPDEP1(I,J)   + DPDE(I+IHE(J),J+1))                               &
    &                        * RFCP
!
                        WDUTD = WDEP * SIGN(1.0,STNE) * (T(I,J,K+1)       - DPT                   &
    &                         -                        T(I+IHE(J),J+1,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        IF (STNE > 0.0) THEN
                            WARR                =  DPDE(I+IHE(J),J+1)
                              ATC(I+IHE(J),J+1) =   ATC(I+IHE(J),J+1) + WDUTD   / (WARR+WDEP)
                              AQC(I+IHE(J),J+1) =   AQC(I+IHE(J),J+1) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+IHE(J),J+1) = ACWMC(I+IHE(J),J+1) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR         =  DPDEP1(I,J)
                            ATCP1(I,J)   =   ATCP1(I,J) + WDUTD   / (WARR+WDEP)
                            AQCP1(I,J)   =   AQCP1(I,J) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I,J) = ACWMCP1(I,J) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                         DIV(I+IHE(J),J-1,K+1) =  DIV(I+IHE(J),J-1,K+1) + STNS !YELLOW
                        DIVS(I+IHE(J),J+1)     = DIVS(I+IHE(J),J+1)     - STNS !YELLOW
!
                        WDEP = ABS(STNS) * DT
                        DPT  =   (RTOP(I+IHE(J),J-1,K+1) + RTOP(I+IHE(J),J+1,K))                  &
    &                        * (DPDEP1(I+IHE(J),J-1)     + DPDE(I+IHE(J),J+1))                    &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,STNS) * (T(I+IHE(J),J-1,K+1) - DPT                &
    &                         -                          T(I+IHE(J),J+1,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,STNS) * (  Q(I+IHE(J),J-1,K+1)                  &
    &                           -                            Q(I+IHE(J),J+1,K  ))
!
                        WDUCWMD = WDEP * SIGN(1.0,STNS) * (CWM(I+IHE(J),J-1,K+1)                  &
    &                           -                          CWM(I+IHE(J),J+1,K  ))
!
                        IF (STNS > 0.0) THEN
                            WARR                =  DPDE(I+IHE(J),J+1)
                              ATC(I+IHE(J),J+1) =   ATC(I+IHE(J),J+1) + WDUTD   / (WARR+WDEP)
                              AQC(I+IHE(J),J+1) =   AQC(I+IHE(J),J+1) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+IHE(J),J+1) = ACWMC(I+IHE(J),J+1) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR                  =  DPDEP1(I+IHE(J),J-1)
                              ATCP1(I+IHE(J),J-1) =   ATCP1(I+IHE(J),J-1) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+IHE(J),J-1) =   AQCP1(I+IHE(J),J-1) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+IHE(J),J-1) = ACWMCP1(I+IHE(J),J-1) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                    ELSE IF (ISLD(I+IHE(J),J) == 4) THEN
                         DIV(I+1,J,K+1)    =  DIV(I+1,J,K+1)    - STSEP1
                        DIVS(I+IHE(J),J+1) = DIVS(I+IHE(J),J+1) + STSEP1
!
                        WDEP = ABS(-STSEP1) * DT
                        DPT  =   (RTOP(I+1,J,K+1) + RTOP(I+IHE(J),J+1,K))                         &
    &                        * (DPDEP1(I+1,J)     + DPDE(I+IHE(J),J+1))                           &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,-STSEP1) * (T(I+1,J,K+1)       - DPT              &
    &                         -                             T(I+IHE(J),J+1,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,-STSEP1) * (  Q(I+1,J,K+1) -   Q(I+IHE(J),J+1,K))
                        WDUCWMD = WDEP * SIGN(1.0,-STSEP1) * (CWM(I+1,J,K+1) - CWM(I+IHE(J),J+1,K))
!
                        IF (-STSEP1 > 0.0) THEN
                            WARR                =  DPDE(I+IHE(J),J+1)
                              ATC(I+IHE(J),J+1) =   ATC(I+IHE(J),J+1) + WDUTD   / (WARR+WDEP)
                              AQC(I+IHE(J),J+1) =   AQC(I+IHE(J),J+1) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+IHE(J),J+1) = ACWMC(I+IHE(J),J+1) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR           =  DPDEP1(I+1,J)
                              ATCP1(I+1,J) =   ATCP1(I+1,J) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+1,J) =   AQCP1(I+1,J) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+1,J) = ACWMCP1(I+1,J) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                         DIV(I+IHE(J),J-1,K+1) =  DIV(I+IHE(J),J-1,K+1) - STSE
                        DIVS(I,J)              = DIVS(I,J)              + STSE
                        WDEP = ABS(-STSE) * DT
                        DPT  =   (RTOP(I+IHE(J),J-1,K+1) + RTOP(I,J,K))                           &
    &                        * (DPDEP1(I+IHE(J),J-1)     + DPDE(I,J))                             &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,-STSE) * (T(I+IHE(J),J-1,K+1) - DPT               &
    &                         -                           T(I,J,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,-STSE) * (  Q(I+IHE(J),J-1,K+1) -   Q(I,J,K))
                        WDUCWMD = WDEP * SIGN(1.0,-STSE) * (CWM(I+IHE(J),J-1,K+1) - CWM(I,J,K))
!
                        IF (-STSE > 0.0) THEN
                            WARR       =  DPDE(I,J)
                              ATC(I,J) =   ATC(I,J) + WDUTD   / (WARR+WDEP)
                              AQC(I,J) =   AQC(I,J) + WDUQD   / (WARR+WDEP)
                            ACWMC(I,J) = ACWMC(I,J) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR                  =  DPDEP1(I+IHE(J),J-1)
                              ATCP1(I+IHE(J),J-1) =   ATCP1(I+IHE(J),J-1) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+IHE(J),J-1) =   AQCP1(I+IHE(J),J-1) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+IHE(J),J-1) = ACWMCP1(I+IHE(J),J-1) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                        DIV(I+IHE(J),J-1,K+1) =  DIV(I+IHE(J),J-1,K+1) + STNS !YELLOW
                        DIVS(I+IHE(J),J+1)    = DIVS(I+IHE(J),J+1)     - STNS !YELLOW
!
                        WDEP = ABS(STNS) * DT
                        DPT  =   (RTOP(I+IHE(J),J-1,K+1) + RTOP(I+IHE(J),J+1,K))                  &
    &                        * (DPDEP1(I+IHE(J),J-1)     + DPDE(I+IHE(J),J+1))                    &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,STNS) * (T(I+IHE(J),J-1,K+1) - DPT                &
    &                         -                          T(I+IHE(J),J+1,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,STNS) * (  Q(I+IHE(J),J-1,K+1)                  &
    &                           -                            Q(I+IHE(J),J+1,K  ))
!
                        WDUCWMD = WDEP * SIGN(1.0,STNS) * (CWM(I+IHE(J),J-1,K+1)                  &
    &                           -                          CWM(I+IHE(J),J+1,K  ))
!
                        IF (STNS > 0.0) THEN
                            WARR                =  DPDE(I+IHE(J),J+1)
                              ATC(I+IHE(J),J+1) =   ATC(I+IHE(J),J+1) + WDUTD   / (WARR+WDEP)
                              AQC(I+IHE(J),J+1) =   AQC(I+IHE(J),J+1) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+IHE(J),J+1) = ACWMC(I+IHE(J),J+1) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR                  =  DPDEP1(I+IHE(J),J-1)
                              ATCP1(I+IHE(J),J-1) =   ATCP1(I+IHE(J),J-1) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+IHE(J),J-1) =   AQCP1(I+IHE(J),J-1) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+IHE(J),J-1) = ACWMCP1(I+IHE(J),J-1) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                        DIV(I+1,J,K+1) =  DIV(I+1,J,K+1) - STEW !GREEN
                        DIVS(I,J)      = DIVS(I,J)       + STEW !GREEN
!
                        WDEP  = ABS(-STEW) * DT
                        DPT   =   (RTOP(I+1,J,K+1) + RTOP(I,J,K))                                 &
    &                         * (DPDEP1(I+1,J)     + DPDE(I,J))                                   &
    &                         *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,-STEW) * (T(I+1,J,K+1) - DPT                      &
    &                         -                           T(I  ,J,K  ))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,-STEW) * (  Q(I+1,J,K+1) -   Q(I,J,K))
                        WDUCWMD = WDEP * SIGN(1.0,-STEW) * (CWM(I+1,J,K+1) - CWM(I,J,K))
!
                        IF (-STEW > 0.0) THEN
                            WARR       =  DPDE(I,J)
                              ATC(I,J) =   ATC(I,J) + WDUTD   / (WARR+WDEP)
                              AQC(I,J) =   AQC(I,J) + WDUQD   / (WARR+WDEP)
                            ACWMC(I,J) = ACWMC(I,J) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR           =  DPDEP1(I+1,J)
                              ATCP1(I+1,J) =   ATCP1(I+1,J) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+1,J) =   AQCP1(I+1,J) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+1,J) = ACWMCP1(I+1,J) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                    ELSE IF (ISLD(I+IHE(J),J) == 5) THEN
                         DIV(I+IHE(J),J+1,K+1) =  DIV(I+IHE(J),J+1,K+1) - STNE
                        DIVS(I,J)              = DIVS(I,J)              + STNE
!
                        WDEP = ABS(-STNE) * DT
                        DPT  =   (RTOP(I+IHE(J),J+1,K+1) + RTOP(I,J,K))                           &
    &                        * (DPDEP1(I+IHE(J),J+1)     + DPDE(I,J))                             &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,-STNE) * (T(I+IHE(J),J+1,K+1) - DPT               &
    &                         -                           T(I       ,J  ,K  ))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,-STNE) * (  Q(I+IHE(J),J+1,K+1) -   Q(I,J,K))
                        WDUCWMD = WDEP * SIGN(1.0,-STNE) * (CWM(I+IHE(J),J+1,K+1) - CWM(I,J,K))
!
                        IF (-STNE > 0.0) THEN
                            WARR       =  DPDE(I,J)
                              ATC(I,J) =   ATC(I,J) + WDUTD   / (WARR+WDEP)
                              AQC(I,J) =   AQC(I,J) + WDUQD   / (WARR+WDEP)
                            ACWMC(I,J) = ACWMC(I,J) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR                  =  DPDEP1(I+IHE(J),J+1)
                              ATCP1(I+IHE(J),J+1) =   ATCP1(I+IHE(J),J+1) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+IHE(J),J+1) =   AQCP1(I+IHE(J),J+1) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+IHE(J),J+1) = ACWMCP1(I+IHE(J),J+1) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                         DIV(I+IHE(J),J-1,K+1) =  DIV(I+IHE(J),J-1,K+1) - STSE
                        DIVS(I,J)              = DIVS(I,J)              + STSE
!
                        WDEP = ABS(-STSE) * DT
                        DPT  =   (RTOP(I+IHE(J),J-1,K+1) + RTOP(I,J,K))                           &
    &                        * (DPDEP1(I+IHE(J),J-1)     + DPDE(I,J))                             &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,-STSE) * (T(I+IHE(J),J-1,K+1) - DPT               &
    &                         -                           T(I       ,J  ,K  ))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,-STSE) * (  Q(I+IHE(J),J-1,K+1) -   Q(I,J,K))
                        WDUCWMD = WDEP * SIGN(1.0,-STSE) * (CWM(I+IHE(J),J-1,K+1) - CWM(I,J,K))
!
                        IF (-STSE > 0.0) THEN
                            WARR       =  DPDE(I,J)
                              ATC(I,J) =   ATC(I,J) + WDUTD   / (WARR+WDEP)
                              AQC(I,J) =   AQC(I,J) + WDUQD   / (WARR+WDEP)
                            ACWMC(I,J) = ACWMC(I,J) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR                  =  DPDEP1(I+IHE(J),J-1)
                              ATCP1(I+IHE(J),J-1) =   ATCP1(I+IHE(J),J-1) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+IHE(J),J-1) =   AQCP1(I+IHE(J),J-1) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+IHE(J),J-1) = ACWMCP1(I+IHE(J),J-1) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                         DIV(I+1,J,K+1) =  DIV(I+1,J,K+1) - STEW !GREEN
                        DIVS(I,J)       = DIVS(I,J)       + STEW !GREEN
!
                        WDEP = ABS(-STEW) * DT
                        DPT  =   (RTOP(I+1,J,K+1) + RTOP(I,J,K))                                  &
    &                        * (DPDEP1(I+1,J)     + DPDE(I,J))                                    &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,-STEW) * (T(I+1,J,K+1) - DPT                      &
    &                         -                           T(I,J,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,-STEW) * (  Q(I+1,J,K+1) -   Q(I,J,K))
                        WDUCWMD = WDEP * SIGN(1.0,-STEW) * (CWM(I+1,J,K+1) - CWM(I,J,K))
!
                        IF (-STEW > 0.0) THEN
                            WARR       =  DPDE(I,J)
                              ATC(I,J) =   ATC(I,J) + WDUTD   / (WARR+WDEP)
                              AQC(I,J) =   AQC(I,J) + WDUQD   / (WARR+WDEP)
                            ACWMC(I,J) = ACWMC(I,J) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR           =  DPDEP1(I+1,J)
                              ATCP1(I+1,J) =   ATCP1(I+1,J) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+1,J) =   AQCP1(I+1,J) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+1,J) = ACWMCP1(I+1,J) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                    ELSE IF (ISLD(I+IHE(J),J) == 6) THEN
                         DIV(I+IHE(J),J+1,K+1) =  DIV(I+IHE(J),J+1,K+1) - STNE
                        DIVS(I,J)              = DIVS(I,J)              + STNE
!
                        WDEP = ABS(-STNE) * DT
                        DPT  =   (RTOP(I+IHE(J),J+1,K+1) + RTOP(I,J,K))                           &
    &                        * (DPDEP1(I+IHE(J),J+1)     + DPDE(I,J))                             &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,-STNE) * (T(I+IHE(J),J+1,K+1) - DPT               &
    &                         -                           T(I,J,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,-STNE) * (  Q(I+IHE(J),J+1,K+1) -   Q(I,J,K))
                        WDUCWMD = WDEP * SIGN(1.0,-STNE) * (CWM(I+IHE(J),J+1,K+1) - CWM(I,J,K))
!
                        IF (-STNE > 0.0) THEN
                            WARR       =  DPDE(I,J)
                              ATC(I,J) =   ATC(I,J) + WDUTD   / (WARR+WDEP)
                              AQC(I,J) =   AQC(I,J) + WDUQD   / (WARR+WDEP)
                            ACWMC(I,J) = ACWMC(I,J) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR                  =  DPDEP1(I+IHE(J),J+1)
                              ATCP1(I+IHE(J),J+1) =   ATCP1(I+IHE(J),J+1) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+IHE(J),J+1) =   AQCP1(I+IHE(J),J+1) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+IHE(J),J+1) = ACWMCP1(I+IHE(J),J+1) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                         DIV(I+1,J,K+1)    =  DIV(I+1,J,K+1)    - STNEM1
                        DIVS(I+IHE(J),J-1) = DIVS(I+IHE(J),J-1) + STNEM1
!
                        WDEP = ABS(-STNEM1) * DT
!
                        DPT =   (RTOP(I+1,J,K+1) + RTOP(I+IHE(J),J-1,K))                          &
    &                       * (DPDEP1(I+1,J)     + DPDE(I+IHE(J),J-1))                            &
    &                       *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,-STNEM1) * (T(I+1,J,K+1)       - DPT              &
    &                         -                             T(I+IHE(J),J-1,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,-STNEM1) * (  Q(I+1,J,K+1) -   Q(I+IHE(J),J-1,K))
                        WDUCWMD = WDEP * SIGN(1.0,-STNEM1) * (CWM(I+1,J,K+1) - CWM(I+IHE(J),J-1,K))
!
                        IF (-STNEM1 > 0.0) THEN
                            WARR                =  DPDE(I+IHE(J),J-1)
                              ATC(I+IHE(J),J-1) =   ATC(I+IHE(J),J-1) + WDUTD   / (WARR+WDEP)
                              AQC(I+IHE(J),J-1) =   AQC(I+IHE(J),J-1) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+IHE(J),J-1) = ACWMC(I+IHE(J),J-1) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR           =  DPDEP1(I+1,J)
                            ATCP1(I+1,J)   =   ATCP1(I+1,J) + WDUTD   / (WARR+WDEP)
                            AQCP1(I+1,J)   =   AQCP1(I+1,J) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+1,J) = ACWMCP1(I+1,J) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                         DIV(I+1,J,K+1) =  DIV(I+1,J,K+1) - STEW !GREEN
                        DIVS(I,J)       = DIVS(I,J)       + STEW !GREEN
!
                        WDEP = ABS(-STEW) * DT
                        DPT =   (RTOP(I+1,J,K+1) + RTOP(I,J,K))                                   &
    &                       * (DPDEP1(I+1,J)     + DPDE(I,J))                                     &
    &                       *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,-STEW) * (T(I+1,J,K+1) - DPT                      &
    &                         -                           T(I,J,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,-STEW) * (  Q(I+1,J,K+1) -   Q(I,J,K))
                        WDUCWMD = WDEP * SIGN(1.0,-STEW) * (CWM(I+1,J,K+1) - CWM(I,J,K))
!
                        IF (-STEW > 0.0) THEN
                            WARR       =  DPDE(I,J)
                              ATC(I,J) =   ATC(I,J) + WDUTD   / (WARR+WDEP)
                              AQC(I,J) =   AQC(I,J) + WDUQD   / (WARR+WDEP)
                            ACWMC(I,J) = ACWMC(I,J) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR           =  DPDEP1(I+1,J)
                              ATCP1(I+1,J) =   ATCP1(I+1,J) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+1,J) =   AQCP1(I+1,J) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+1,J) = ACWMCP1(I+1,J) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                         DIV(I+IHE(J),J+1,K+1) =  DIV(I+IHE(J),J+1,K+1) - STNS !PINK
                        DIVS(I+IHE(J),J-1)     = DIVS(I+IHE(J),J-1)     + STNS !PINK
!
                        WDEP = ABS(-STNS) * DT
                        DPT  =   (RTOP(I+IHE(J),J+1,K+1) + RTOP(I+IHE(J),J-1,K))                  &
    &                        * (DPDEP1(I+IHE(J),J+1)     + DPDE(I+IHE(J),J-1))                    &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,-STNS) * (T(I+IHE(J),J+1,K+1) - DPT               &
    &                         -                           T(I+IHE(J),J-1,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,-STNS) * (  Q(I+IHE(J),J+1,K+1)                 &
    &                           -                             Q(I+IHE(J),J-1,K  ))
!
                        WDUCWMD = WDEP * SIGN(1.0,-STNS) * (CWM(I+IHE(J),J+1,K+1)                 &
    &                           -                           CWM(I+IHE(J),J-1,K  ))
!
                        IF (-STNS > 0.0) THEN
                            WARR                =  DPDE(I+IHE(J),J-1)
                              ATC(I+IHE(J),J-1) =   ATC(I+IHE(J),J-1) + WDUTD   / (WARR+WDEP)
                              AQC(I+IHE(J),J-1) =   AQC(I+IHE(J),J-1) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+IHE(J),J-1) = ACWMC(I+IHE(J),J-1) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR                  =  DPDEP1(I+IHE(J),J+1)
                              ATCP1(I+IHE(J),J+1) =   ATCP1(I+IHE(J),J+1) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+IHE(J),J+1) =   AQCP1(I+IHE(J),J+1) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+IHE(J),J+1) = ACWMCP1(I+IHE(J),J+1) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                    ELSE IF (ISLD(I+IHE(J),J) == 7) THEN
                         DIV(I,J,K+1)      =  DIV(I,J,K+1)      + STSE
                        DIVS(I+IHE(J),J-1) = DIVS(I+IHE(J),J-1) - STSE
!
                        WDEP = ABS(STSE) * DT
                        DPT  =   (RTOP(I,J,K+1) + RTOP(I+IHE(J),J-1,K))                           &
    &                        * (DPDEP1(I,J)     + DPDE(I+IHE(J),J-1))                             &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,STSE) * (T(I,J,K+1) - DPT                         &
    &                         -                          T(I+IHE(J),J-1,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0, STSE) * (  Q(I,J,K+1) -   Q(I+IHE(J),J-1,K))
                        WDUCWMD = WDEP * SIGN(1.0, STSE) * (CWM(I,J,K+1) - CWM(I+IHE(J),J-1,K))
!
                        IF (STSE > 0.0) THEN
                            WARR                =  DPDE(I+IHE(J),J-1)
                              ATC(I+IHE(J),J-1) =   ATC(I+IHE(J),J-1) + WDUTD   / (WARR+WDEP)
                              AQC(I+IHE(J),J-1) =   AQC(I+IHE(J),J-1) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+IHE(J),J-1) = ACWMC(I+IHE(J),J-1) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR         =  DPDEP1(I,J)
                              ATCP1(I,J) =   ATCP1(I,J) + WDUTD   / (WARR+WDEP)
                              AQCP1(I,J) =   AQCP1(I,J) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I,J) = ACWMCP1(I,J) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                        DIV(I+1,J,K+1)     = DIV(I+1,J,K+1)     - STNEM1
                        DIVS(I+IHE(J),J-1) = DIVS(I+IHE(J),J-1) + STNEM1
!
                        WDEP = ABS(-STNEM1) * DT
                        DPT  =   (RTOP(I+1,J,K+1) + RTOP(I+IHE(J),J-1,K))                         &
    &                        * (DPDEP1(I+1,J)     + DPDE(I+IHE(J),J-1))                           &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,-STNEM1) * (T(I+1     ,J  ,K+1)       - DPT       &
    &                         -                             T(I+IHE(J),J-1,K  ))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,-STNEM1) * (  Q(I+1,J,K+1) -   Q(I+IHE(J),J-1,K))
                        WDUCWMD = WDEP * SIGN(1.0,-STNEM1) * (CWM(I+1,J,K+1) - CWM(I+IHE(J),J-1,K))
!
                        IF (-STNEM1 > 0.0) THEN
                            WARR                =  DPDE(I+IHE(J),J-1)
                              ATC(I+IHE(J),J-1) =   ATC(I+IHE(J),J-1) + WDUTD   / (WARR+WDEP)
                              AQC(I+IHE(J),J-1) =   AQC(I+IHE(J),J-1) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+IHE(J),J-1) = ACWMC(I+IHE(J),J-1) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR           =  DPDEP1(I+1,J)
                              ATCP1(I+1,J) =   ATCP1(I+1,J) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+1,J) =   AQCP1(I+1,J) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+1,J) = ACWMCP1(I+1,J) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                         DIV(I+IHE(J),J+1,K+1) =  DIV(I+IHE(J),J+1,K+1) - STNS !PINK
                        DIVS(I+IHE(J),J-1)     = DIVS(I+IHE(J),J-1)     + STNS !PINK
!
                        WDEP = ABS(-STNS) * DT
                        DPT  =   (RTOP(I+IHE(J),J+1,K+1) + RTOP(I+IHE(J),J-1,K))                  &
    &                        * (DPDEP1(I+IHE(J),J+1)     + DPDE(I+IHE(J),J-1))                    &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,-STNS) * (T(I+IHE(J),J+1,K+1) - DPT               &
    &                         -                           T(I+IHE(J),J-1,K  ))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,-STNS) * (  Q(I+IHE(J),J+1,K+1)                 &
    &                           -                             Q(I+IHE(J),J-1,K  )) 
!
                        WDUCWMD = WDEP * SIGN(1.0,-STNS) * (CWM(I+IHE(J),J+1,K+1)                 &
    &                           -                           CWM(I+IHE(J),J-1,K  ))
!
                        IF (-STNS > 0.0) THEN
                            WARR                =  DPDE(I+IHE(J),J-1)
                              ATC(I+IHE(J),J-1) =   ATC(I+IHE(J),J-1) + WDUTD   / (WARR+WDEP)
                              AQC(I+IHE(J),J-1) =   AQC(I+IHE(J),J-1) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+IHE(J),J-1) = ACWMC(I+IHE(J),J-1) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR                  =  DPDEP1(I+IHE(J),J+1)
                              ATCP1(I+IHE(J),J+1) =   ATCP1(I+IHE(J),J+1) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+IHE(J),J+1) =   AQCP1(I+IHE(J),J+1) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+IHE(J),J+1) = ACWMCP1(I+IHE(J),J+1) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                    ELSE IF (ISLD(I+IHE(J),J) == 8) THEN
                         DIV(I+IHE(J),J+1,K+1) =  DIV(I+IHE(J),J+1,K+1) + STSEP1
                        DIVS(I+1,J)            = DIVS(I+1,J)            - STSEP1
!
                        WDEP = ABS( STSEP1) * DT
                        DPT  =   (RTOP(I+IHE(J),J+1,K+1) + RTOP(I+1,J,K))                         &
    &                        * (DPDEP1(I+IHE(J),J+1)     + DPDE(I+1,J))                           &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,STSEP1) * (T(I+IHE(J),J+1,K+1) - DPT              &
    &                         -                            T(I+1     ,J  ,K  ))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,STSEP1) * (  Q(I+IHE(J),J+1,K+1) -   Q(I+1,J,K))
                        WDUCWMD = WDEP * SIGN(1.0,STSEP1) * (CWM(I+IHE(J),J+1,K+1) - CWM(I+1,J,K))
!
                        IF ( STSEP1 > 0.0) THEN
                            WARR         =  DPDE(I+1,J)
                              ATC(I+1,J) =   ATC(I+1,J) + WDUTD   / (WARR+WDEP)
                              AQC(I+1,J) =   AQC(I+1,J) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+1,J) = ACWMC(I+1,J) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR                  =  DPDEP1(I+IHE(J),J+1)
                              ATCP1(I+IHE(J),J+1) =   ATCP1(I+IHE(J),J+1) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+IHE(J),J+1) =   AQCP1(I+IHE(J),J+1) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+IHE(J),J+1) = ACWMCP1(I+IHE(J),J+1) + WDUCWMD / (WARR+WDEP)
                        END IF
!                        		   
                         DIV(I,J,K+1)      =  DIV(I,J,K+1)      + STSE
                        DIVS(I+IHE(J),J-1) = DIVS(I+IHE(J),J-1) - STSE
!
                        WDEP = ABS(STSE) * DT
                        DPT  =   (RTOP(I,J,K+1) + RTOP(I+IHE(J),J-1,K))                           &
    &                        * (DPDEP1(I,J)     + DPDE(I+IHE(J),J-1))                             &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,STSE) * (T(I,J,K+1) - DPT                         &
    &                         -                          T(I+IHE(J),J-1,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0, STSE) * (  Q(I,J,K+1) -   Q(I+IHE(J),J-1,K))
                        WDUCWMD = WDEP * SIGN(1.0, STSE) * (CWM(I,J,K+1) - CWM(I+IHE(J),J-1,K))
!
                        IF (STSE > 0.0) THEN
                            WARR                =  DPDE(I+IHE(J),J-1)
                              ATC(I+IHE(J),J-1) =   ATC(I+IHE(J),J-1) + WDUTD   / (WARR+WDEP)
                              AQC(I+IHE(J),J-1) =   AQC(I+IHE(J),J-1) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+IHE(J),J-1) = ACWMC(I+IHE(J),J-1) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR         =  DPDEP1(I,J)
                              ATCP1(I,J) =   ATCP1(I,J) +WDUTD    / (WARR+WDEP)
                              AQCP1(I,J) =   AQCP1(I,J) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I,J) = ACWMCP1(I,J) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                         DIV(I,J,K+1) =  DIV(I,J,K+1) + STEW !BLUE
                        DIVS(I+1,J)   = DIVS(I+1,J)   - STEW !BLUE
!
                        WDEP = ABS(STEW) * DT
                        DPT  =   (RTOP(I,J,K+1) + RTOP(I+1,J,K))                                  &
    &                        * (DPDEP1(I,J)     + DPDE(I+1,J))                                    &
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,STEW) * (T(I  ,J,K+1) - DPT                       &
    &                         -                          T(I+1,J,K  ))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,STEW) * (  Q(I,J,K+1) -   Q(I+1,J,K))
                        WDUCWMD = WDEP * SIGN(1.0,STEW) * (CWM(I,J,K+1) - CWM(I+1,J,K))
!
                        IF (STEW > 0.0) THEN
                            WARR         =  DPDE(I+1,J)
                              ATC(I+1,J) =   ATC(I+1,J) + WDUTD   / (WARR+WDEP)
                              AQC(I+1,J) =   AQC(I+1,J) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+1,J) = ACWMC(I+1,J) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR         =  DPDEP1(I,J)
                              ATCP1(I,J) =   ATCP1(I,J) + WDUTD   / (WARR+WDEP)
                              AQCP1(I,J) =   AQCP1(I,J) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I,J) = ACWMCP1(I,J) + WDUCWMD / (WARR+WDEP)
                        END IF
!
                         DIV(I+IHE(J),J+1,K+1) =  DIV(I+IHE(J),J+1,K+1) - STNS !PINK
                        DIVS(I+IHE(J),J-1)     = DIVS(I+IHE(J),J-1)     + STNS !PINK
!
                        WDEP = ABS(-STNS) * DT
                        DPT  =   (RTOP(I+IHE(J),J+1,K+1) + RTOP(I+IHE(J),J-1,K))                  &
    &                        * (DPDEP1(I+IHE(J),J+1)     + DPDE(I+IHE(J),J-1))                    & 
    &                        *  RFCP
!
                        WDUTD = WDEP * SIGN(1.0,-STNS) * (T(I+IHE(J),J+1,K+1) - DPT               &
    &                         -                           T(I+IHE(J),J-1,K))
!----------------------------------------------------------
! CHOU: INCLUDING SPECIFIC HUMIDITY AND LIQUID WATER FLUXES
!----------------------------------------------------------
                        WDUQD   = WDEP * SIGN(1.0,-STNS) * (  Q(I+IHE(J),J+1,K+1)                 &
    &                           -                             Q(I+IHE(J),J-1,K  ))
!
                        WDUCWMD = WDEP * SIGN(1.0,-STNS) * (CWM(I+IHE(J),J+1,K+1)                 &
    &                           -                           CWM(I+IHE(J),J-1,K  ))
!
                        IF (-STNS > 0.0) THEN
                            WARR                =  DPDE(I+IHE(J),J-1)
                              ATC(I+IHE(J),J-1) =   ATC(I+IHE(J),J-1) + WDUTD   / (WARR+WDEP)
                              AQC(I+IHE(J),J-1) =   AQC(I+IHE(J),J-1) + WDUQD   / (WARR+WDEP)
                            ACWMC(I+IHE(J),J-1) = ACWMC(I+IHE(J),J-1) + WDUCWMD / (WARR+WDEP)
                        ELSE
                            WARR                  =  DPDEP1(I+IHE(J),J+1)
                              ATCP1(I+IHE(J),J+1) =   ATCP1(I+IHE(J),J+1) + WDUTD   / (WARR+WDEP)
                              AQCP1(I+IHE(J),J+1) =   AQCP1(I+IHE(J),J+1) + WDUQD   / (WARR+WDEP)
                            ACWMCP1(I+IHE(J),J+1) = ACWMCP1(I+IHE(J),J+1) + WDUCWMD / (WARR+WDEP)
                        END IF
                    END IF
!-------
! K < LM
!-------
                END IF 
!            
            END DO
        END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                DIV(I,J,K) = (DIV(I,J,K) + DIVL(I,J)) * HM(I,J)
            END DO
        END DO
!
    END DO
!
    RETURN
!
    END SUBROUTINE DIVHOASTQL
