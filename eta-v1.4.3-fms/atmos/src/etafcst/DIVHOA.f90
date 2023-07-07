    SUBROUTINE DIVHOA
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE DIVHOA
!>
!> SUBPROGRAM: DIVHOA - DIVERGENCE/HORIZONTAL OMEGA-ALPHA
!> PROGRAMMER: JANJIC
!> ORG: W/NP22
!> DATE: 93-10-28
!>
!> ABSTRACT:
!> DIVHOA COMPUTES THE DIVERGENCE INCLUDING THE MODIFICATION PREVENTING GRAVITY WAVE GRID SEPARATION
!> AND CALCULATES THE HORIZONTAL PART OF THE OMEGA-ALPHA TERM (THE PART PROPORTIONAL TO THE 
!> ADVECTION OF MASS ALONG ETA/SIGMA SURFACES).
!>
!> PROGRAM HISTORY LOG:
!> 87-06-??  JANJIC     - ORIGINATOR
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 96-03-29  BLACK      - ADDED EXTERNAL EDGE
!> 97-03-17  MESINGER   - SPLIT FROM PFDHT
!> 98-10-30  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
!> 00-10-20  BLACK      - INCORPORATED PRESSURE GRADIENT METHOD FROM MESO MODEL
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
!>              NHYDRO
!>              PARMETA
!>              TEMPCOM 
!>              TOPO
!>              VRBLS
!>  
!> DRIVER     : DIGFLT
!>              EBU
!>              NEWFLT
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
    USE MAPOT
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE NHYDRO
    USE PARMETA
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM+1)                             ::&
    & PINTLG
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::& 
    & FIM     ,                                                                                   &
    & FILO    , RDPD    ,                                                                         &
    & ADPDX   , RDPDX   ,                                                                         &
    & ADPDY   , RDPDY   ,                                                                         &
    & ADPDNE  , ADPDSE  ,                                                                         &
    & PEW     , PNS     ,                                                                         &
    & PCEW    , PCNS    ,                                                                         &
    & DPFEW   , DPFNS   ,                                                                         &
    & FNS     , TNS     ,                                                                         &
    & HM      , VM      ,                                                                         &
    & EDIV    , DIVL  
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & DPDE    ,                                                                                   &
    & APEL    , PCXC    ,                                                                         &
    & ALP1    ,                                                                                   &
    & DFDZ    ,                                                                                   &
    & UDY     , VDX     ,                                                                         &
    & TEW     , FEW     ,                                                                         &
    & TNE     , TSE     ,                                                                         &
    & FNE     , FSE     ,                                                                         &
    & PNE     , PSE     ,                                                                         &
    & CNE     , CSE     ,                                                                         &
    & PPNE    , PPSE    ,                                                                         & 
    & PCNE    , PCSE  
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K 
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ALP1P   , ALP1X   , DFI     , RDPDS   , FIUPK   , ALP1PL  , ALP2P   , ALP2PL  , DPNEK   ,   &
    & DPSEK   , DCNEK   , DCSEK   , PVNEK   , PVSEK
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
        DO 50 J=MYJS_P5,MYJE_P5
            DO 50 I=MYIS_P5,MYIE_P5
                FILO(I,J) = FIS(I,J)
                PDSL(I,J) = PD (I,J)
     50 END DO
!
    ELSE
!------- 
! OPENMP
!-------
! 
!$omp parallel do
!
        DO 100 J=MYJS_P5,MYJE_P5
            DO 100 I=MYIS_P5,MYIE_P5
                FILO(I,J) = 0.0
                PDSL(I,J) = RES(I,J) * PD(I,J)
    100 END DO
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
                DIV   (I,J,K) = 0.
            END DO
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
            ALP1X = PINTLG(I,J,LM+1)
            ALP1(I,J) = ALP1X
        END DO
    END DO
!-------------------------------
! MAIN VERTICAL INTEGRATION LOOP  
!------------------------------- 
    DO 400 K=LM,1,-1
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
        DO 125 J=MYJS_P5,MYJE_P5
            DO 125 I=MYIS_P5,MYIE_P5
!            
                ALP1P = PINTLG(I,J,K)
!
                DFI = (Q(I,J,K) * 0.608 + 1.) * T(I,J,K) * R * (ALP1(I,J) - ALP1P)                &
    &               / (1 + CWM(I,J,K)) / DWDT(I,J,K)
!            
                RDPDS       = 1. / (DETA(K) * PDSL(I,J))
                RTOP(I,J,L) = RDPDS * DFI
                FIUPK       = FILO(I,J) + DFI
                FIM(I,J)    = FILO(I,J) + FIUPK
!
                IF (ABS(FIM(I,J)) <= 5.E+10) THEN
!
                ELSE
                    WRITE(6,*) 'BAD FIM ', I, J, FIM(I,J), FILO(I,J), DFI
                    WRITE(6,*) 'Q,T,ALP1,ALP1P,DWDT: ', Q(I,J,K),    T(I,J,K), ALP1(I,J),         &
    &                                                   ALP1P   , DWDT(I,J,K)
                    STOP
!
                END IF
!            
                FILO(I,J) = (FIUPK - DFL(K)) * HTM(I,J,K) + DFL(K)
                ALP1(I,J) = ALP1P
    125 END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (ALP1P  , ALP1PL  , ALP2P   , ALP2PL  , DFI)
!
        DO 205 J=MYJS_P5,MYJE_P5
            DO 205 I=MYIS_P5,MYIE_P5
                HM(I,J) = HTM(I,J,K) * HBM2(I,J)
                VM(I,J) = VTM(I,J,K) * VBM2(I,J)
!            
                ALP1P  = PINTLG(I,J,K  )
                ALP1PL = PINTLG(I,J,K+1)
                ALP2P  = ALP1P  * ALP1P
                ALP2PL = ALP1PL * ALP1PL
!
                DFI = (Q(I,J,K) * 0.608 + 1.) * T(I,J,K) * R * (ALP1PL-ALP1P)                     &
    &               / (1 + CWM(I,J,K)) / DWDT(I,J,K)
                DFDZ(I,J) = DFI * DWDT(I,J,K) / (ALP2PL-ALP2P)
                APEL(I,J) = (ALP2PL + ALP2P) * 0.5
    205 END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO 210 J=MYJS_P4,MYJE_P4
            DO 210 I=MYIS_P4,MYIE_P4
                DPDE(I,J) = DETA(K) * PDSL(I,J)
                DIVL(I,J) = 0.
                EDIV(I,J) = 0.
    210 END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO 215 J=MYJS_P1,MYJE_P1
            DO 215 I=MYIS_P1,MYIE_P1
                RDPD(I,J) = 1. / DPDE(I,J)
    215 END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO 220 J=MYJS1_P3,MYJE1_P3
            DO 220 I=MYIS_P3,MYIE_P3
                ADPDX(I,J) = DPDE(I+IVW(J),J) + DPDE(I+IVE(J),J)
                ADPDY(I,J) = DPDE(I,J-1) + DPDE(I,J+1)
                RDPDX(I,J) = 1. / ADPDX(I,J)
                RDPDY(I,J) = 1. / ADPDY(I,J)
    220 END DO
!--------------------------------------------------    
! DIAGONAL CONTRIBUTIONS TO PRESSURE GRADIENT FORCE
!--------------------------------------------------     
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO 240 J=MYJS_P4,MYJE_P4
            DO 240 I=MYIS_P4,MYIE_P4
                ADPDNE(I,J) = DPDE(I+IHE(J),J+1) + DPDE(I,J)
                PNE   (I,J) = (FIM(I+IHE(J),J+1) - FIM (I,J)) * (DWDT(I+IHE(J),J+1,K)+DWDT(I,J,K))
!
                IF ( ABS(PNE(I,J)) <= 5.E10) THEN
!
                ELSE
!
                    WRITE(6,*) 'CRAZY PNE ', I, J, PNE(I,J)
                    WRITE(6,*) 'PIECES', I+IHE(J), J+1, FIM(I+IHE(J),J+1)
!
                END IF
!
                PPNE(I,J) = PNE(I,J) * ADPDNE(I,J)
                CNE(I,J)  = (DFDZ(I+IHE(J),J+1)+DFDZ(I,J)) * 2. * (APEL(I+IHE(J),J+1)-APEL(I,J))
                PCNE(I,J) = CNE(I,J) * ADPDNE(I,J)
    240 END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO 250 J=MYJS1_P4,MYJE_P4
            DO 250 I=MYIS_P4,MYIE1_P4
                ADPDSE(I,J) = DPDE(I+IHE(J),J-1) + DPDE(I,J)
                PSE( I,J)   = (FIM(I+IHE(J),J-1)-FIM(I,J)) * (DWDT(I+IHE(J),J-1,K)+DWDT(I,J,K))
                PPSE(I,J)   = PSE(I,J) * ADPDSE(I,J)
                CSE (I,J)   = (DFDZ(I+IHE(J),J-1)+DFDZ(I,J)) * 2. * (APEL(I+IHE(J),J-1)-APEL(I,J))
                PCSE(I,J)   = CSE(I,J) * ADPDSE(I,J)
    250 END DO
!---------------------------------     
! CONTINUITY EQUATION MODIFICATION 
!---------------------------------
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO 260 J=MYJS1_P1,MYJE1_P1
            DO 260 I=MYIS_P1,MYIE_P1
                PCXC(I,J) = VBM3(I,J) * VTM(I,J,K) * (PNE(I+IVW(J),J) + CNE(I+IVW(J),J)           &
    &                     +                           PSE(I+IVW(J),J) + CSE(I+IVW(J),J)           &
    &                     -                           PNE(I,J-1) - CNE(I,J-1)                     &
    &                     -                           PSE(I,J+1) - CSE(I,J+1))
    260 END DO
!
        DO 270 J=MYJS2,MYJE2
            DO 270 I=MYIS1,MYIE1
                DIV(I,J,K) = DETA(K) * WPDAR(I,J) * (PCXC(I+IHE(J),J) - PCXC(I,J+1)               &
    &                      +                         PCXC(I+IHW(J),J) - PCXC(I,J-1))
    270 END DO
!---------------------------------------      
! LAT AND LONG PRESSURE FORCE COMPONENTS
!---------------------------------------      
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (DCNEK  , DCSEK   , DPNEK   , DPSEK)
!
        DO 280 J=MYJS1_P3,MYJE1_P3
            DO 280 I=MYIS_P3,MYIE_P3
                DPNEK     = PNE(I+IVW(J),J) + PNE(I,J-1)
                DPSEK     = PSE(I+IVW(J),J) + PSE(I,J+1)
                PEW(I,J)  = DPNEK           + DPSEK
                PNS(I,J)  = DPNEK           - DPSEK
                DCNEK     = CNE(I+IVW(J),J) + CNE(I,J-1)
                DCSEK     = CSE(I+IVW(J),J) + CSE(I,J+1)
                PCEW(I,J) = (DCNEK+DCSEK)   * ADPDX(I,J)
                PCNS(I,J) = (DCNEK-DCSEK)   * ADPDY(I,J)
    280 END DO
!----------------------------------------------      
! LAT AND LON FLUXES AND OMEGA-ALPHA COMPONENTS
!----------------------------------------------    
!------- 
! OPENMP
!------- 
!$omp parallel do
!
        DO 310 J=MYJS1_P3,MYJE1_P3
            DO 310 I=MYIS_P3,MYIE_P3
                UDY(I,J) = DY       * U    (I,J,K)
                FEW(I,J) = UDY(I,J) * ADPDX(I,J)
                TEW(I,J) = UDY(I,J) * PCEW (I,J)
                VDX(I,J) = DX (I,J) * V    (I,J,K)
                FNS(I,J) = VDX(I,J) * ADPDY(I,J)
                TNS(I,J) = VDX(I,J) * PCNS (I,J)
    310 END DO
!---------------------------------------------      
! DIAGONAL FLUXES AND DIAGONALLY AVERAGED WIND
!---------------------------------------------   
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (PVNEK)
!
        DO 320 J=MYJS1_P2,MYJE2_P2
            DO 320 I=MYIS_P2,MYIE1_P2
                PVNEK    = (UDY(I+IHE(J),J) + VDX(I+IHE(J),J)) + (UDY(I,J+1) + VDX(I,J+1))
                FNE(I,J) = PVNEK * ADPDNE(I,J)
                TNE(I,J) = PVNEK * PCNE(I,J) * 2.
    320 END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (PVSEK)
!
        DO 330 J=MYJS2_P2,MYJE1_P2
            DO 330 I=MYIS_P2,MYIE1_P2
                PVSEK    = (UDY(I+IHE(J),J) - VDX(I+IHE(J),J)) + (UDY(I,J-1) - VDX(I,J-1))
                FSE(I,J) = PVSEK * ADPDSE(I,J)
                TSE(I,J) = PVSEK * PCSE(I,J) * 2.
    330 END DO
!----------------------------------------------     
! HORIZONTAL PART OF OMEGA-ALPHA AND DIVERGENCE
!----------------------------------------------  
!------- 
! OPENMP
!-------
! 
!$omp parallel do
!
        DO 340 J=MYJS2_P1,MYJE2_P1
            DO 340 I=MYIS1_P1,MYIE1_P1
                OMGALF(I,J,K) = (TEW(I+IHE(J),J) + TEW(I+IHW(J),J) + TNS(I,J+1) + TNS(I,J-1)      &
    &                         +  TNE(I,J) + TNE(I+IHW(J),J-1)                                     &
    &                         +  TSE(I,J) + TSE(I+IHW(J),J+1))                                    &
    &                         * RDPD(I,J) * FCP(I,J) * HM(I,J)
!
                T(I,J,K)      = OMGALF(I,J,K) + T(I,J,K)
!
                EDIV(I,J) = ((FEW(I+IHE(J),J  ) + FNS(I       , J+1)                              &
    &                     +   FNE(I       ,J  ) + FSE(I       , J  ))                             &
    &                     -  (FEW(I+IHW(J),J  ) + FNS(I       , J-1)                              &
    &                     +   FNE(I+IHW(J),J-1) + FSE(I+IHW(J), J+1)))                            &
    &                     *  FDIV(I       ,J  )
!
                DIVL(I,J) = EDIV(I,J) * HBM2(I,J)
    340 END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO 390 J=MYJS,MYJE
            DO 390 I=MYIS,MYIE
                DIV(I,J,K) = (DIV(I,J,K) + DIVL(I,J)) * HM(I,J)
    390 END DO
!
400 END DO
!
    RETURN 
!
    END SUBROUTINE DIVHOA
