    SUBROUTINE SLPSIGSPLINE(PD, FIS, TSIG, QSIG, SPL, LSL, DETA, PT, PSLP)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SLPSIGSPLINE
!>
!> SUBROUTINE: SLPSIGSPLINE - SLP REDUCTION
!> PROGRAMMER: JANJIC, D. JOVIC, S. NICKOVIC
!> ORG: ?????
!> DATE: 81-02-??
!>
!> PROGRAMMER: JANJIC, D. JOVIC, S. NICKOVIC
!> ORG: ?????
!> DATE: 81-02-??
!>
!> ABSTRACT: THIS ROUTINE COMPUTES THE SEA LEVEL PRESSURE REDUCTION USING CUBIC SPLINE FITTING
!>           METHOD OR THE STANDARD NCEP REDUCTION FOR SIGMA COORDINATES.
!>
!> PROGRAM HISTORY LOG:
!> 81-02-??  JANJIC       - ORIGINATOR
!> 00-09-??  H CHUANG     - INCLUDED IN OPERATIONAL MPI QUILT 
!> 18-01-15  LUCCI        - MODERNIZATION OF THE CODE, INCLUDING:
!>                          * F77 TO F90/F95
!>                          * INDENTATION & UNIFORMIZATION CODE
!>                          * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                          * DOCUMENTATION WITH DOXYGEN
!>                          * OPENMP FUNCTIONALITY
!>
!> NOTE: CURRENT POST ONLY PROCESS UP TO 1000 MB PRESSURE LEVEL 2 ADDITIONAL LEVELS (1025 1050) ARE
!>       ADDED DURING SPLINE FITTING COMPUTATION BUT THE FIELDS ON THESE TWO LEVELS ARE NOT OUTPUT 
!>       TO GRID FILES.
!>
!> INPUT ARGUMENT LIST:
!> PD   - SFC PRESSURE MINUS PTOP
!> FIS  - SURFACE GEOPOTENTIAL
!> TSIG - TEMPERATURE
!> QSIG - SPECIFIC HUMIDITY
!> FI   - GEOPOTENTIAL
!> PT   - TOP PRESSURE OF DOMAIN
!>
!> OUTPUT ARGUMENT LIST:
!> PSLP - THE FINAL REDUCED SEA LEVEL PRESSURE ARRAY
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MPPCOM
!>              PARA
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : QUILT
!>
!> CALLS      : SPLINEF
!>              UPDATE 
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARA
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE	
!
    INCLUDE "mpif.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: LMP1   = LM  + 1
    INTEGER(KIND=I4KIND), PARAMETER :: NSMUD  =  25 
    INTEGER(KIND=I4KIND), PARAMETER :: LP2    = LM  + 2
    INTEGER(KIND=I4KIND), PARAMETER :: LSMP2  = LSM + 2
    INTEGER(KIND=I4KIND), PARAMETER :: III    =  22
    INTEGER(KIND=I4KIND), PARAMETER :: JJJ    = 293
!
    REAL   (KIND=R4KIND), PARAMETER :: RD     = 287.04
    REAL   (KIND=R4KIND), PARAMETER :: ROG    = RD / 9.8
    REAL   (KIND=R4KIND), PARAMETER :: PQ0    = 379.90516
    REAL   (KIND=R4KIND), PARAMETER :: A2     =  17.2693882 
    REAL   (KIND=R4KIND), PARAMETER :: A3     = 273.16
    REAL   (KIND=R4KIND), PARAMETER :: A4     =  35.86
    REAL   (KIND=R4KIND), PARAMETER :: GAMMA  =   6.5E-3
    REAL   (KIND=R4KIND), PARAMETER :: RGAMOG = GAMMA * ROG
    REAL   (KIND=R4KIND), PARAMETER :: H1M12  =   1.E-12
    REAL   (KIND=R4KIND), PARAMETER :: G      =   9.8
!
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED)                    , INTENT(IN)          ::&
    & PD      , FIS     
!
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED)                                          ::&
    & TMASK   , FSLL
!
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED, LSMP2)                                   ::&
    & FI
!
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED, LM)                , INTENT(IN)          ::&
    & TSIG    , QSIG
!
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED)                    , INTENT(INOUT)       ::&
    & PSLP
!
    REAL   (KIND=R4KIND), DIMENSION(LSM)                                  , INTENT(IN)          ::&
    & SPL
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                   , INTENT(IN)          ::&
    & DETA 
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                                         ::&
    & AETA
!
    REAL   (KIND=R4KIND), DIMENSION(LMP1)                                                       ::&
    & ETA
!
    REAL   (KIND=R4KIND), DIMENSION(LM+2)                                                       ::&
    & ZTH     , HCOL    , Y2      , PHLD    , QHLD
!
    REAL   (KIND=R4KIND), DIMENSION(LSMP2)                                                      ::&
    & ZTSL    , OVRLX   , HCOLSL
!
    REAL   (KIND=R4KIND), DIMENSION(3)                                                          ::&
    & HCOL3   , PCOL3   , Y23
!
    DATA Y2 /LP2 * 0./
    DATA OVRLX /LSMP2 * 0.175/
!
    INTEGER(KIND=I4KIND), DIMENSION(JM)                                                         ::&
    & IHE     , IHW
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & SIGMA, STDRD
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LSL
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & PT
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , L       , IHL     , IHH     , IVI     , LMD     , N       , LLL     ,   &
    & LMA
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ZTBOT   , PDP     , ALP1L   , ALP2L   , DH      , X       , DUM     , D2      , RLX     ,   &
    & HREF    , GAR     , SLP1    , ALP1U   , ALP2U   , ALPP1   , SLOP    , SLPP    , TTT
!
    STDRD = .FALSE. 
!------------------------------------------ 
!CALCULATE THE I-INDEX EAST-WEST INCREMENTS
!------------------------------------------ 
    DO J=1,JM
        IHE(J) = MOD(J+1,2)
        IHW(J) = IHE(J)-1
    END DO
!---------------------------------------------------------- 
! INITIALIZE ARRAYS.  LOAD SLP ARRAY WITH SURFACE PRESSURE.
!----------------------------------------------------------
!
!$omp parallel do
!
    DO J=JSTA_I,JEND_I
        DO I=1,IM
            PSLP(I,J) = PD(I,J) + PT
        END DO
    END DO
!
    DO 100 L=1,LSM
        ZTSL(L) = ALOG(SPL(L) ) ** 2
100 END DO
!
    ZTSL(LSM+1) = ALOG(102500.) ** 2
    ZTSL(LSMP2) = ALOG(105000.) ** 2
!
    ZTBOT = ZTSL(LSMP2)
!---------------------  
! COMPUTE ETA AND AETA
!--------------------- 
    ETA(1) = 0.0
!
    DO L=2,LM+1
         ETA(L) =        ETA(L-1) + DETA(L-1)
    END DO
!
    DO L=1,LM
        AETA(L) = 0.5 * (ETA(L  ) +  ETA(L+1))
    END DO
!
    DO 200 J=JSTA_I,JEND_I
        IHL = 1
        IHH = IM-1 + MOD(J,2)
!
        DO 201 I=IHL,IHH
            PDP        =  PD(I,J)
            HCOL(LM+1) = FIS(I,J) / G
!
            ALP1L      =  ALOG(PT  + PDP)
            ZTH(LM+1)  = (ALOG(PDP + PT )) ** 2
!
            DO 220 IVI=1,LM
                L       = LM + 1 - IVI
                ALP1U   =  ALOG(ETA(L) * PDP + PT)
                ALP2U   = (ALOG(ETA(L) * PDP + PT)) ** 2
                DH      = (QSIG(I,J,L) * 0.608 + 1.) * (ALP1L - ALP1U) * TSIG(I,J,L) * ROG
!
                 ZTH(L) = ALP2U
                HCOL(L) = HCOL(L+1) + DH
!
                ALP1L   = ALP1U
        220 END DO
!
            IF (ZTH(LM+1) >= ZTBOT) THEN
                LMD = LM + 1
            ELSE
                LMD        = LM + 2
                ZTH(LM+2)  = ZTBOT
                X          = ZTBOT - ZTH(LM+1)
                DUM        = (HCOL(LM+1) - HCOL(LM-5)) / (ZTH(LM+1) - ZTH(LM-5))
                D2         = 0.
                HCOL(LM+2) = D2 * X * X + DUM * X + HCOL(LM+1)
            END IF
!
            Y2(LMD) = 0.
!
            CALL SPLINEF(LMD, ZTH, HCOL, Y2, LSMP2, ZTSL, HCOLSL, PHLD, QHLD)
!
            DO 240 L=1,LSMP2
                FI(I,J,L) = HCOLSL(L)
        240 END DO
!
    201 END DO
!
 200 END DO
!-------------------------------------------------- 
! FILLING REMAINING UNDEFINED MASS POINT BOUNDARIES 
!-------------------------------------------------- 
    DO 300 L=1,LSMP2
!
        IF (MOD(JSTA_I,2) < 1) THEN !EVEN JSTA_I
!
            DO 351 J=JSTA_I,JEND_I,2
                FI(IM,J,L) = FI(IM-1,J,L)
        351 END DO
!
        ELSE
!
            DO 352 J=JSTA_I+1,JEND_I,2
                FI(IM,J,L) = FI(IM-1,J,L)
            352 END DO
!
        END IF
!
        RLX = OVRLX(L)
!
        DO 310 J=JSTA_I,JEND_I
            IHL = 1
            IHH = IM - 1 + MOD(J,2)
!
            DO 311 I=IHL,IHH
                HREF = FIS(I,J) / G
!
                IF (HREF > 300.) THEN
                    HREF = HREF + 1500.
                END IF
!
                IF (HREF >= FI(I,J,L)) THEN
                    TMASK(I,J) = RLX
                ELSE
                    TMASK(I,J) = 0.
                END IF
!
        311 END DO
!
    310 END DO
!
        DO 320 N=1,NSMUD
!
            CALL UPDATE(FI(1,MY_JSD,L))
!
            DO 330 J=JSTA_IM,JEND_IM
                IHL = 1 + MOD(J,2)
                IHH = IM-1
!
                DO 331 I=IHL,IHH
                    IF (TMASK(I,J) > 0.05) THEN
                        FSLL(I,J) = FI(I + IHW(J),J-1,L) + FI(I + IHE(J),J-1,L)                   &
    &                             + FI(I + IHW(J),J+1,L) + FI(I + IHE(J),J+1,L) - FI(I,J,L) * 4.
                    END IF
            331 END DO
!
        330 END DO
!------------------------- 
! PRINT*,'AFTER SMOOTHING'
!------------------------- 
            DO 340 J=JSTA_IM,JEND_IM
                IHL = 1 + MOD(J,2)
                IHH = IM-1
!
                DO 341 I=IHL,IHH
                    IF (TMASK(I,J) > 0.05) THEN
                        FI(I,J,L) = FSLL(I,J) * TMASK(I,J) + FI(I,J,L)
                    END IF
            341 END DO
!
        340 END DO
!
    320 END DO
!
300 END DO
!-------------------  
! SEA LEVEL PRESSURE 
!-------------------
    DO 600 J=JSTA_I,JEND_I
        IHL = 1
        IHH = IM - 1 + MOD(J,2)
!
        DO 601 I=IHL,IHH
            IF (FIS(I,J) > -1. .AND. FIS(I,J) < 1.) THEN
                PSLP(I,J) = PT + PD(I,J)
            ELSE
                IF (FI(I,J,LSMP2) > 0) THEN
!---------------------- 
! ASSUME PSLP < 1100 MB
!---------------------- 
                    GAR = EXP(ZTSL(LSMP2) ** 0.5) + 5000. 
                    HCOL3(1) = FI(I,J,LSMP2) + (FI(I,J,LSMP2) - FI(I,J,LSMP2-5))                  &
    &                        * (ALOG(GAR) ** 2. - ZTSL(LSMP2)) / (ZTSL(LSMP2) - ZTSL(LSMP2-5))
                    HCOL3(2) = FI(I,J,LSMP2  )
                    HCOL3(3) = FI(I,J,LSMP2-1)
                    PCOL3(1) = GAR
                    PCOL3(2) = EXP(ZTSL(LSMP2  ) ** 0.5)
                    PCOL3(3) = EXP(ZTSL(LSMP2-1) ** 0.5)
                ELSE
                    DO L=LSMP2,2,-1
                        IF (FI(I,J,L) <= 0. .AND. FI(I,J,L-1) > 0.) THEN
                            LLL = L
                            GOTO 603
                        END IF
                    END DO
!
                603 HCOL3(1) = FI(I,J,LLL  )
                    HCOL3(2) = FI(I,J,LLL-1)
                    HCOL3(3) = FI(I,J,LLL-2)
!
                    PCOL3(1) = EXP(ZTSL(LLL  ) ** 0.5)
                    PCOL3(2) = EXP(ZTSL(LLL-1) ** 0.5)
                    PCOL3(3) = EXP(ZTSL(LLL-2) ** 0.5)
                END IF
!
                DO L=1,3
                    Y23(L) = 0.0
                END DO
!
                IF (HCOL3(1) == HCOL3(2) .OR. HCOL3(2) == HCOL3(3))                               &
    &               PRINT*, 'PROBLEM SLP', I, J, HCOL3(1), HCOL3(2), HCOL3(3)
!
                CALL SPLINEF(3, HCOL3, PCOL3, Y23, 1, 0.0, SLP1, PHLD, QHLD)
                PSLP(I,J) = SLP1
            END IF
!
    601 END DO
!
600 END DO
!-------------------------- 
! SKIP THE STANDARD SCHEME.
!-------------------------- 
    GOTO 430
!-------------------------------------------------------------------------
! IF YOU WANT THE "STANDARD" ETA/SIGMA REDUCTION THIS IS WHERE IT IS DONE.
!-------------------------------------------------------------------------
    DO 410 J=JSTA_I,JEND_I
        DO 410 I=1,IM
            IF (FIS(I,J) >= 1.) THEN
                LMA   = LM
                ALPP1 = ALOG(PD(I,J) + PT)
                SLOP  = 0.0065 * ROG * TSIG(I,J,LM)
!
                IF (SLOP < 0.50) THEN
                    SLPP = ALPP1 + FIS(I,J) / (RD * TSIG(I,J,LMA))
                ELSE
                    TTT = -(ALOG(PD(I,J) + PT) + ALPP1) *   SLOP * 0.50 + TSIG(I,J,LMA)
                    SLPP = (-TTT + SQRT(TTT * TTT + 2. * SLOP * (FIS(I,J) / RD                    &
    &                    +  (TTT + 0.50 * SLOP * ALPP1) * ALPP1))) / SLOP
                END IF
!
                PSLP(I,J) = EXP(SLPP)
!
            END IF
410 END DO
!--------------------------------------------------------------------------------------------------
! AT THIS POINT WE HAVE A SEA LEVEL PRESSURE FIELD BY EITHER METHOD. 5-POINT AVERAGE THE FIELD ON
! THE E-GRID.
!--------------------------------------------------------------------------------------------------
!
430 CONTINUE
!
    RETURN
!
    END SUBROUTINE SLPSIGSPLINE
!
!
!
    SUBROUTINE SPLINEF(NOLD, XOLD, YOLD, Y2, NNEW, XNEW, YNEW, P, Q)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SPLINEF
!>
!> SUBROUTINE: SLPSIGSPLINE - SLP REDUCTION
!> PROGRAMMER: JZ. JANJIC, YUGOSLAV FED. HYDROMET.
!> ORG: INST. BEOGRAD
!> DATE: 81-02-??
!>
!> ABSTRACT:
!> THIS IS A ONE-DIMENSIONAL CUBIC SPLINE FITTING ROUTINE PROGRAMED FOR A SMALL SCALAR MACHINE.
!>
!> PROGRAM HISTORY LOG:
!> 81-02-??  JANJIC       - ORIGINATOR
!> 18-01-15  LUCCI        - MODERNIZATION OF THE CODE, INCLUDING:
!>                          * F77 TO F90/F95
!>                          * INDENTATION & UNIFORMIZATION CODE
!>                          * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                          * DOCUMENTATION WITH DOXYGEN
!>                          * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NOLD - 
!> XOLD -
!> YOLD -
!> NNEW -
!> XNEW - 
!> YNEW -
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> Y2   -
!> P    -
!> Q    - 
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : SLPSIGSPLINE
!>
!> CALLS      : -----
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
! NOLD - NUMBER OF GIVEN VALUES OF THE FUNCTION. MUST BE GE 3.
! XOLD - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE FUNCTION ARE GIVEN. MUST BE IN 
!        ASCENDING ORDER.
! YOLD - THE GIVEN VALUES OF THE FUNCTION AT THE POINTS XOLD. 
! Y2   - THE SECOND DERIVATIVES AT THE POINTS XOLD. IF NATURAL SPLINE IS FITTED Y2(1)=0. 
!        AND Y2(NOLD)=0. MUST BE SPECIFIED.
! NNEW - NUMBER OF VALUES OF THE FUNCTION TO BE CALCULATED.
! XNEW - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE FUNCTION ARE CALCULATED. XNEW(K) 
!        MUST BE GE XOLD(1) AND LE XOLD(NOLD).
! YNEW - THE VALUES OF THE FUNCTION TO BE CALCULATED. P, Q - AUXILIARY VECTORS OF THE LENGTH NOLD-2
!--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & NOLD    , NNEW
!
    REAL   (KIND=R4KIND), DIMENSION(NOLD)                                 , INTENT(IN)          ::&
    & XOLD    , YOLD    
!
    REAL   (KIND=R4KIND), DIMENSION(NOLD)                                 , INTENT(INOUT)       ::&
    & Y2      , P       , Q
!
    REAL   (KIND=R4KIND), DIMENSION(NNEW)                                 , INTENT(IN)          ::&
    & XNEW     
!
    REAL   (KIND=R4KIND), DIMENSION(NNEW)                                 , INTENT(INOUT)       ::&
    & YNEW
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NOLDM1  , K       , K1      , K2      , KOLD     
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & DXL     , DXR     , DYDXL   , DYDXR   , RTDXC   , DXC     , DEN     , XK      , Y2K     ,   &
    & Y2KP1   , DX      , RDX     , AK      , BK      , CK      , X       , XSQ
!
    NOLDM1 = NOLD - 1
!
    DXL   =  XOLD(2)  - XOLD(1)
    DXR   =  XOLD(3)  - XOLD(2)
    DYDXL = (YOLD(2)  - YOLD(1)) / DXL
    DYDXR = (YOLD(3)  - YOLD(2)) / DXR
    RTDXC = .5 / (DXL + DXR)
!
    P(1) =  RTDXC * (6. * (DYDXR - DYDXL) - DXL * Y2(1))
    Q(1) = -RTDXC * DXR
!
    IF (NOLD == 3) GOTO 700
!
    K = 3
!
100 DXL   = DXR
    DYDXL = DYDXR
    DXR   = XOLD(K+1) - XOLD(K)
!
    DYDXR = (YOLD(K+1) - YOLD(K)) / DXR
    DXC   = DXL + DXR
    DEN   = 1. / (DXL * Q(K-2) + DXC + DXC)
!
    P(K-1) =  DEN * (6. * (DYDXR - DYDXL) - DXL * P(K-2))
    Q(K-1) = -DEN * DXR
!
    K = K + 1
!
    IF (K < NOLD) GOTO 100
!
700 K = NOLDM1
!
200 Y2(K) = P(K-1) + Q(K-1) * Y2(K+1)
!
    K = K - 1
!
    IF (K > 1) GOTO 200
!
    K1 = 1
!
300 XK = XNEW(K1)
!
    DO 400 K2=2,NOLD
        IF (XOLD(K2) <= XK) GOTO 400
        KOLD = K2 - 1
        GOTO 450
400 END DO
!
    YNEW(K1) = YOLD(NOLD)
!
    GOTO 600
!
450 IF (K1 == 1)   GOTO 500
!
    IF (K == KOLD) GOTO 550
!
500 K = KOLD
!
    Y2K   = Y2(K  )
    Y2KP1 = Y2(K+1)
!
    DX    = XOLD(K+1) - XOLD(K)
    RDX   = 1. / DX
!
    AK = .1666667 * RDX * (Y2KP1 - Y2K)
    BK = .5 * Y2K
    CK = RDX * (YOLD(K+1) - YOLD(K)) - .1666667 * DX * (Y2KP1 + Y2K + Y2K)
!
550 X   = XK - XOLD(K)
    XSQ = X * X
!
    YNEW(K1) = AK * XSQ * X + BK * XSQ + CK * X + YOLD(K)
!
600 K1 = K1 + 1
!
    IF (K1 <= NNEW) GOTO 300
!
    RETURN
!
    END SUBROUTINE SPLINEF
