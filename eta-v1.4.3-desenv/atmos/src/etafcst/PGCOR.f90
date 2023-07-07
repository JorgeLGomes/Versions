    SUBROUTINE PGCOR
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE PGCOR
!>
!> SUBROUTINE: PGCOR - PRESSURE GRADIENT/CORIOLIS CALC
!> PROGRAMMER: JANJIC
!> ORG: W/NP22
!> DATE: 93-10-28
!>
!> ABSTRACT:
!> PGCOR CALCULATES THE PRESSURE GRADIENT FORCE, UPDATES THE VELOCITY COMPONENTS DUE TO THE EFFECT
!> OF THE PRESSURE GRADIENT AND CORIOLIS FORCES.
!>
!> PROGRAM HISTORY LOG:
!> 87-06-??  JANJIC     - ORIGINATOR
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 96-03-29  BLACK      - ADDED EXTERNAL EDGE
!> 97-03-17  MESINGER   - SPLIT FROM PFDHT
!> 98-10-28  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
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
!> OUTPUT FILES:
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
    REAL   (KIND=R8KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM+1)                             ::&
    & PINTLG

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
    & HM      , VM    
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & DPDE    ,                                                                                   &
    & PNE     , PSE     ,                                                                         &
    & CNE     , CSE     ,                                                                         &
    & PPNE    , PPSE    ,                                                                         &
    & PCNE    , PCSE  
!
    REAL   (KIND=R8KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & APEL    ,                                                                                   &
    & ALP1    ,                                                                                   &
    & DFDZ    
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & RDPDS   , FIUPK   , F0K     ,                                                               &
    & DPNEK   , DPSEK   ,                                                                         &
    & DCNEK   , DCSEK   ,                                                                         &
    & DPFNEK  , DPFSEK  ,                                                                         &
    & UPK     , VPK     ,                                                                         &
    & UTK     , VTK 
!
    REAL   (KIND=R8KIND)                                                                        ::&
    & ALP1P   , ALP1PL  ,                                                                         &
    & ALP2P   , ALP2PL  ,                                                                         &
    & DFI
!
    CALL ZERO2(DPDE)
    CALL ZERO2(APEL)
    CALL ZERO2(ADPDX)
    CALL ZERO2(ADPDY)
    CALL ZERO2(DFDZ)
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
!
!$omp parallel do 
!
        DO 50 J=MYJS_P5,MYJE_P5
            DO 50 I=MYIS_P5,MYIE_P5
                FILO(I,J) = FIS(I,J)
                PDSL(I,J) =  PD(I,J)
     50 END DO
!
    ELSE
!
!$omp parallel do  
!
        DO 100 J=MYJS_P5,MYJE_P5
            DO 100 I=MYIS_P5,MYIE_P5
                FILO(I,J) = 0.0
                PDSL(I,J) = RES(I,J) * PD(I,J)
    100 END DO
    END IF
!
    IF (HYDRO) THEN
!
!$omp parallel do 
!
        DO K=1,LM+1
            DO J=MYJS_P5,MYJE_P5
                DO I=MYIS_P5,MYIE_P5
                    PINTLG(I,J,K) = ALOG(ETA(K) * PDSL(I,J) + PT)
                    IF (K > 1) THEN
                        IF (PINTLG(I,J,K) == PINTLG(I,J,K-1)) THEN
                            WRITE(6,*) 'SAME PINTLG AT DIFFERENT LEVELS: ', MYPE,I,J,K,LMH(I,J), PINT(I,J,K  ), PINT(I,J,K-1)
!
                            PINT(I,J,K)   = PDSL(I,J) * (DETA(K-1)/2.) + PINT(I,J,K-1)
                            PINTLG(I,J,K) = ALOG(PINT(I,J,K))
!
                            WRITE(6,*) 'NEW PINT: ', I,J,K, PINT(I,J,K), PDSL(I,J), DETA(K), PINT(I,J,K-1), ALOG(PINT(I,J,K)), ALOG(PINT(I,J,K-1))
!
                        END IF
                    END IF

                END DO
            END DO
        END DO
    ELSE
!
!$omp parallel do  
!
        DO K=1,LM+1
            DO J=MYJS_P5,MYJE_P5
                DO I=MYIS_P5,MYIE_P5
                    PINTLG(I,J,K) = ALOG(PINT(I,J,K))
                    IF (K > 1) THEN
                        IF (PINTLG(I,J,K) == PINTLG(I,J,K-1)) THEN
                            WRITE(6,*) 'SAME PINTLG AT DIFFERENT LEVELS: ', MYPE,I,J,K,LMH(I,J), PINT(I,J,K  ), PINT(I,J,K-1)
!
                            PINT(I,J,K)   = PDSL(I,J) * (DETA(K-1)/2.) + PINT(I,J,K-1)
                            PINTLG(I,J,K) = ALOG(PINT(I,J,K))
!
                            WRITE(6,*) 'NEW PINT: ', I,J,K, PINT(I,J,K), PDSL(I,J), DETA(K), PINT(I,J,K-1), ALOG(PINT(I,J,K)), ALOG(PINT(I,J,K-1))
!
                        END IF
                    END IF
                END DO
            END DO
        END DO
    END IF
!
!$omp parallel do private (ALP1P)
!
    DO J=MYJS_P5,MYJE_P5
        DO I=MYIS_P5,MYIE_P5
            ALP1P     = PINTLG(I,J,LM+1)
            ALP1(I,J) = ALP1P
        END DO
    END DO
!------------------------------- 
! MAIN VERTICAL INTEGRATION LOOP 
!------------------------------- 
    FIM = 0.
!
    do400: DO K=LM,1,-1
!--------------------------- 
! INTEGRATE THE GEOPOTENTIAL
!--------------------------- 
!
!$omp parallel do private (ALP1P  , DFI     , FIUPK   , RDPDS)
!
        doout125: DO J=MYJS_P5,MYJE_P5
            doin125: DO I=MYIS_P5,MYIE_P5
!            
                ALP1P = PINTLG(I,J,K)
!            
                DFI = (Q(I,J,K) * 0.608 + 1.) * T(I,J,K) * R * (ALP1(I,J) - ALP1P) / DWDT(I,J,K)
!
                IF (ABS(DFI) < 2.E13) THEN
                ELSE
                    WRITE(6,*) 'BAD DFI: '  , DFI
                    WRITE(6,*) 'MYPE,I,J,K,LMH,Q,T: ', MYPE,I,J,K,LMH(I,J)  , Q(I,J,K)   , T(I,J,K)
                    WRITE(6,*) 'MYPE,I,J,K,PDSL,PD,RES: ', MYPE,I,J,K,PDSL(I,J) , PD(I,J)   , RES(I,J)
                    WRITE(6,*) 'ALP VALS: ' , ALP1(I,J)  , ALP1P
                    WRITE(6,*) 'DWDT= '     , DWDT(I,J,K)
                END IF
!           
                RDPDS = 1. / (DETA(K) * PDSL(I,J))
                FIUPK = FILO(I,J) + DFI
!
                IF (ABS(FIUPK) < 2.E13) THEN
                ELSE
                    WRITE(6,*) 'BAD FIUPK.  FILO, DFI ', FILO(I,J), DFI
                END IF
!
                 FIM(I,J) = FILO(I,J) + FIUPK
            
                FILO(I,J) = (FIUPK - DFL(K)) * HTM(I,J,K) + DFL(K)
!
                IF (ABS(FILO(I,J)) < 20000000.) THEN
                ELSE
                    WRITE(6,*) 'BAD FILO VALUE ', FILO(I,J),' ON PE: ' ,MYPE, 'AT ', I, J
                    WRITE(6,*) 'FIUPK,DFL: ', FIUPK, DFL(K)
                    STOP 999
                END IF
!
                ALP1(I,J) = ALP1P
            END DO doin125
         END DO doout125
!
!$omp parallel do private (ALP1P  , ALP1PL  , ALP2P   , ALP2PL  , DFI)
!
         doout205: DO J=MYJS_P5,MYJE_P5
            doin205:DO I=MYIS_P5,MYIE_P5
                HM(I,J) = HTM(I,J,K) * HBM2(I,J)
                VM(I,J) = VTM(I,J,K) * VBM2(I,J)
!            
                ALP1P  = PINTLG(I,J,K  )
                ALP1PL = PINTLG(I,J,K+1)
                ALP2P  = ALP1P  * ALP1P
                ALP2PL = ALP1PL * ALP1PL
!            
!CHOU CWM in Tv                DFI = (Q(I,J,K) * 0.608 + 1.) * T(I,J,K) * R * (ALP1PL - ALP1P) / DWDT(I,J,K)
                DFI = (Q(I,J,K) * 0.608 + 1. - CWM(I,J,K)) * T(I,J,K) * R * (ALP1PL - ALP1P)      &
    &               / DWDT(I,J,K)
!
                IF (ABS(DFI) <= 2.E13) THEN
                ELSE
                    WRITE(6,*) 'BAD DFI     '   , DFI
                END IF
!
                IF (ABS(DWDT(I,J,K)) <= 2.E13) THEN
                ELSE
                    WRITE(6,*) 'BAD DWDTI     ' , DWDT(I,J,K)
                END IF

                IF (ABS(ALP2PL) <= 2.E13) THEN
                ELSE
                    WRITE(6,*) 'BAD ALP2PL     ', ALP2PL
                END IF

                IF (ABS(ALP2P) <= 2.E13) THEN
                ELSE
                    WRITE(6,*) 'BAD ALP2P     ' , ALP2P
                END IF
!            	
                DFDZ(I,J) = DFI * DWDT(I,J,K) / (ALP2PL - ALP2P)
!
                IF (ABS(DFDZ(I,J)) <= 2.E13) THEN
                ELSE
                    WRITE(6,*) 'ON PE: '     , MYPE
                    WRITE(6,*) 'AT = '       , I,J,K
                    WRITE(6,*) 'DFDZ= '      , DFDZ(I,J)
                    WRITE(6,*) 'DFI= '       , DFI
                    WRITE(6,*) 'DWDT= '      , DWDT(I,J,K)
                    WRITE(6,*) 'DENOM= '     , ALP2PL - ALP2P
                    WRITE(6,*) 'PINTLG(K) '  , PINTLG(I,J,K)
                    WRITE(6,*) 'PINTLG(K+1) ', PINTLG(I,J,K+1)
                END IF
!
                APEL(I,J) = (ALP2PL + ALP2P) * 0.5
                DPDE(I,J) = DETA(K) * PDSL(I,J)
            END DO doin205
         END DO doout205
!
!$omp parallel do
! 
         doout215: DO J=MYJS_P1,MYJE_P1
            doin215: DO I=MYIS_P1,MYIE_P1
                RDPD(I,J) = 1. / DPDE(I,J)
            END DO doin215
         END DO doout215
!
!$omp parallel do 
!
         doout220: DO J=MYJS1_P3,MYJE1_P3
            doin220: DO I=MYIS_P3,MYIE_P3
                ADPDX(I,J) = DPDE(I+IVW(J),J  ) + DPDE(I+IVE(J),J  )
                ADPDY(I,J) = DPDE(I       ,J-1) + DPDE(I       ,J+1)
                RDPDX(I,J) = 1./ ADPDX(I,J)
                RDPDY(I,J) = 1./ ADPDY(I,J)
            END DO doin220
         END DO doout220
!--------------------------------------------------     
! DIAGONAL CONTRIBUTIONS TO PRESSURE GRADIENT FORCE 
!--------------------------------------------------  
!  
!$omp parallel do 
!
        DO 240 J=MYJS_P4,MYJE_P4
            DO 240 I=MYIS_P4,MYIE_P4
                ADPDNE(I,J) = DPDE(I+IHE(J),J+1) + DPDE(I,J)
!
                IF (ABS(FIM (I+IHE(J),J+1)) < 2000000.) THEN
                ELSE
                    WRITE(6,*) 'USING FIM VAL: ', FIM (I+IHE(J),J+1), 'AT POINT',I+IHE(J),J+1,    &
    &                          'ON PE: ', MYPE
                END IF
!
                 PNE(I,J) = (FIM (I+IHE(J),J+1)   - FIM (I,J))                                    &
    &                     * (DWDT(I+IHE(J),J+1,K) + DWDT(I,J,K))
!
                PPNE(I,J) = PNE(I,J) * ADPDNE(I,J)
!
                 CNE(I,J) = (DFDZ(I+IHE(J),J+1) + DFDZ(I,J)) * 2.                                 &
    &                     * (APEL(I+IHE(J),J+1) - APEL(I,J))
!
                PCNE(I,J) = CNE(I,J) * ADPDNE(I,J)
    240 END DO
!
!$omp parallel do 
!
        DO 250 J=MYJS1_P4,MYJE_P4
            DO 250 I=MYIS_P4,MYIE1_P4
                ADPDSE(I,J) = DPDE(I+IHE(J),J-1) + DPDE(I,J)
!
                   PSE(I,J) = (FIM (I+IHE(J),J-1)   - FIM (I,J))                                  &
    &                       * (DWDT(I+IHE(J),J-1,K) + DWDT(I,J,K))
!
                  PPSE(I,J) = PSE(I,J)*ADPDSE(I,J)
!
                   CSE(I,J) = (DFDZ(I+IHE(J),J-1) + DFDZ(I,J)) * 2.                               &
    &                       * (APEL(I+IHE(J),J-1) - APEL(I,J))
!
                  PCSE(I,J) = CSE(I,J) * ADPDSE(I,J)
    250 END DO
!---------------------------------------    
! LAT AND LONG PRESSURE FORCE COMPONENTS 
!---------------------------------------    
! 
!$omp parallel do private (DCNEK  , DCSEK   , DPNEK   , DPSEK)
!
        DO 280 J=MYJS1_P3,MYJE1_P3
            DO 280 I=MYIS_P3,MYIE_P3
                DPNEK = PNE(I+IVW(J),J) + PNE(I,J-1)
                DPSEK = PSE(I+IVW(J),J) + PSE(I,J+1)
!
                PEW(I,J) = DPNEK + DPSEK
                PNS(I,J) = DPNEK - DPSEK
!
                DCNEK = CNE(I+IVW(J),J) + CNE(I,J-1)
                DCSEK = CSE(I+IVW(J),J) + CSE(I,J+1)
!
                PCEW(I,J) = (DCNEK + DCSEK) * ADPDX(I,J)
                PCNS(I,J) = (DCNEK - DCSEK) * ADPDY(I,J)
    280 END DO
!----------------------------------      
! UPDATE U AND V (CORIOLIS AND PGF) 
!----------------------------------  
!   
!$omp parallel do private (DPFNEK , DPFSEK)
!
        DO 290 J=MYJS2_P3,MYJE2_P3
            DO 290 I=MYIS_P3,MYIE1_P3
                DPFNEK = ((PPNE(I+IVW(J),J)+PPNE(I,J-1))                                          &
    &                  +  (PCNE(I+IVW(J),J)+PCNE(I,J-1))) * 2.
!
                DPFSEK = ((PPSE(I+IVW(J),J)+PPSE(I,J+1))                                          &
    &                  +  (PCSE(I+IVW(J),J)+PCSE(I,J+1))) * 2.
!
                DPFEW(I,J) = DPFNEK + DPFSEK
                DPFNS(I,J) = DPFNEK - DPFSEK
    290 END DO
!
!$omp parallel do private (F0K    , UPK     , UTK     , VPK     , VTK)
!
        DO 300 J=MYJS2_P2,MYJE2_P2
            DO 300 I=MYIS_P2,MYIE1_P2
                F0K     =   U(I,J,K) * CURV(I,J) + F(I,J)
                VM(I,J) = VTM(I,J,K) * VBM2(I,J)
!
                UPK = ((DPFEW(I,J) + PCEW(I,J)) * RDPDX(I,J)   + PEW(I,J))                        &
    &               *   CPGFU(I,J) + F0K        *     V(I,J,K) +   U(I,J,K)
!
                VPK = ((DPFNS(I,J) + PCNS(I,J)) * RDPDY(I,J)   + PNS(I,J))                        &
    &               *   CPGFV      - F0K        *     U(I,J,K) +   V(I,J,K)
!
                UTK = U(I,J,K)
                VTK = V(I,J,K)
!
                U(I,J,K) = ((F0K * VPK + UPK)     / (F0K * F0K + 1.) - UTK) * VM(I,J) + UTK
                V(I,J,K) =  (VPK - F0K * U(I,J,K)                    - VTK) * VM(I,J) + VTK
    300 END DO
!
    END DO do400
!
    RETURN
!
    END SUBROUTINE PGCOR
