    SUBROUTINE HADZ
!>-------------------------------------------------------------------------------------------------- 
!> SUBROUTINE HADZ
!> 
!> SUBPROGRAM: HADZ - HORIZONTAL ADVECTION OF HEIGHT
!> PROGRAMMER: JANJIC 
!> ORG: W/NP22
!> DATE: 96-05-??
!>
!> ABSTRACT:
!> HADZ CALCULATES DIAGNOSTICALLY THE CONTRIBUTION OF THE HORIZONTAL ADVECTION OF HEIGHT
!>
!> PROGRAM HISTORY LOG:
!> 96-05-??  JANJIC  - ORIGINATOR
!> 00-01-04  BLACK   - DISTRIBUTED MEMORY AND THREADS
!> 18-01-15  LUCCI   - MODERNIZATION OF THE CODE, INCLUDING:
!>                     * F77 TO F90/F95
!>                     * INDENTATION & UNIFORMIZATION CODE
!>                     * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                     * DOCUMENTATION WITH DOXYGEN
!>                     * OPENMP FUNCTIONALITY
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
!>              TEMPCOM
!>              TOPO
!>              VRBLS 
!>
!> DRIVER     : EBU
!>
!> CALLS      : -----           
!>--------------------------------------------------------------------------------------------------
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
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM  = IM * JM - JM / 2
    INTEGER(KIND=I4KIND), PARAMETER :: NTSHY = 2
!
    REAL   (KIND=R4KIND), PARAMETER :: G     = 9.8

!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & HBMS    , DPDE    ,                                                                         &
    & UDY     , VDX     ,                                                                         &
    & UNED    , USED    ,                                                                         &
    & ZEW     , ZNS     ,                                                                         &
    & ZNE     , ZSE
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K       , IX      , JX      , IVH     , IVL     , IHH     , IHL                          
!
    IF (NTSD <= NTSHY .OR. HYDRO) THEN
!------- 
! OPENMP
!-------
! 
!$omp parallel do
!
        DO K=1,LM
            DO J=MYJS,MYJE
                DO I=MYIS,MYIE
                    W(I,J,K) = 0.
                END DO
            END DO
        END DO
!
        RETURN
!
    END IF
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
!$omp private (DPDE     , IHH     , IHL     , IVH     , IVL     , IX      , JX      , UDY     ,   &
!$omp          UNED     , USED    , VDX     , ZEW     , ZNE     , ZNS     , ZSE     )
!
    DO 200 K=1,LM
!
        DO J=MYJS_P3,MYJE_P3
            DO I=MYIS_P3,MYIE_P3
                DPDE(I,J) = PDSL(I,J) * DETA(K)
            END DO
        END DO
!------------------------------------------------- 
! MASS FLUXES AND MASS POINTS ADVECTION COMPONENTS 
!------------------------------------------------- 
        DO 125 J=2,JM-1
            IF (J >= MY_JS_GLB-JBPAD2 .AND. J <= MY_JE_GLB+JTPAD2) THEN
                JX  = J  - MY_JS_GLB + 1
                IVL = 2  - MOD(J,2)
                IVH = IM - 1 
!            
                DO 120 I=IVL,IVH
                    IF (I >= MY_IS_GLB-ILPAD2 .AND. I <= MY_IE_GLB+IRPAD2) THEN
                        IX = I - MY_IS_GLB + 1
                        UDY(IX,JX)  =   U(IX,JX,K) * DY
                        ZEW(IX,JX)  = UDY(IX,JX)   * (DPDE(IX + IVW(JX), JX)                      &
    &                               +                DPDE(IX + IVE(JX), JX))                      &
    &                               *                  (Z(IX + IVE(JX), JX,K)                     &
    &                               -                   Z(IX + IVW(JX), JX,K))
!
                        VDX(IX,JX)  =   V(IX,JX,K) *    DX(IX,JX)
                        ZNS(IX,JX)  = VDX(IX,JX)   * (DPDE(IX,JX-1)   + DPDE(IX,JX+1))            &
    &                                              *    (Z(IX,JX+1,K) -    Z(IX,JX-1,K))
                        UNED(IX,JX) = UDY(IX,JX)   +   VDX(IX,JX)
                        USED(IX,JX) = UDY(IX,JX)   -   VDX(IX,JX)
                    END IF
            120 END DO
!
            END IF
    125 END DO
!--------------------------------------------- 
! DIAGONAL FLUXES AND DIAGONALLY AVERAGED WIND 
!--------------------------------------------- 
        DO 145 J=2,JM-2
            IF (J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
                JX  = J  - MY_JS_GLB + 1
                IHL = 2  -             MOD(J+1,2)
                IHH = IM - 2         + MOD(J  ,2)
!            
                DO 140 I=IHL,IHH
                    IF (I >= MY_IS_GLB-ILPAD1 .AND. I <= MY_IE_GLB+IRPAD1) THEN
                        IX = I - MY_IS_GLB + 1
!
                        ZNE(IX,JX) = (UNED(IX + IHE(JX), JX  )                                    &
    &                              +  UNED(IX          , JX+1))                                   &
    &                              * (DPDE(IX          , JX  )                                    &
    &                              +  DPDE(IX + IHE(JX), JX+1))                                   &
    &                              *    (Z(IX + IHE(JX), JX+1, K)                                 &
    &                              -     Z(IX          , JX  , K))
                    END IF
            140 END DO
!
            END IF
    145 END DO
!    
        DO 165 J=3,JM-1
            IF (J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
                JX  = J  - MY_JS_GLB + 1
                IHL = 2  -             MOD(J+1, 2)
                IHH = IM - 2         + MOD(J  , 2)
!            
                DO 160 I=IHL,IHH
                    IF (I >= MY_IS_GLB-ILPAD1 .AND. I <= MY_IE_GLB+IRPAD1) THEN
                        IX = I - MY_IS_GLB + 1
!
                        ZSE(IX,JX) = (USED(IX + IHE(JX), JX  )                                    &
    &                              +  USED(IX          , JX-1))                                   &
    &                              * (DPDE(IX          , JX  )                                    &
    &                              +  DPDE(IX + IHE(JX), JX-1))                                   &
    &                              *    (Z(IX + IHE(JX), JX-1, K)                                 &
    &                              -     Z(IX          , JX  , K))
                    END IF 
            160 END DO
!
            END IF
    165 END DO
!--------------- 
! ADVECTION OF Z 
!--------------- 
        DO 175 J=3,JM-2
            IF (J >= MY_JS_GLB .AND. J <= MY_JE_GLB) THEN
                JX  = J - MY_JS_GLB + 1
                IHL = 2
                IHH = IM - 2 + MOD(J,2)
!            
                DO 170 I=IHL,IHH
                    IF (I >= MY_IS_GLB .AND. I <= MY_IE_GLB) THEN
                        IX = I - MY_IS_GLB + 1
!
                        W(IX,JX,K) = -(ZEW(IX + IHW(JX), JX)   + ZEW(IX + IHE(JX), JX  )          &
    &                              +   ZNS(IX          , JX-1) + ZNS(IX          , JX+1)          &
    &                              +   ZNE(IX + IHW(JX), JX-1) + ZNE(IX          , JX  )          &
    &                              +   ZSE(IX          , JX  ) + ZSE(IX + IHW(JX), JX+1))         &
    &                              *   FAD(IX          , JX  )                                    &
    &                              *   HTM(IX,JX,K) * HBM2(IX,JX) / (DPDE(IX,JX) * DT)            &
    &                              +     W(IX,JX,K)
                    END IF
            170 END DO
!
            END IF
    175 END DO
!
200 END DO
!
    RETURN
!
    END SUBROUTINE HADZ
