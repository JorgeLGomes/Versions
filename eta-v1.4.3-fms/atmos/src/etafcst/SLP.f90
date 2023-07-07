    SUBROUTINE SLP(NHB, PD, RES, FIS, T, Q, NTSD, NEST, PSLP)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SLP
!> 
!> SUBROUTINE: SLP - SLP REDUCTION
!> PROGRAMMER: BLACK
!> ORG: W/NP22
!> DATE: 99-04-22
!>
!> ABSTRACT:
!> THIS ROUTINE COMPUTES THE SEA LEVEL PRESSURE REDUCTION USING EITHER THE MESINGER RELAXATION 
!> METHOD OR THE STANDARD NCEP REDUCTION.
!>
!> PROGRAM HISTORY LOG:
!> 99-04-22  T BLACK      - ORIGINATOR
!> 00-01-10  JIM TUCCILLO - MPI VERSION
!> 00-11-02  T BLACK      - SLP FOR NEST BOUNDARIES
!> 18-01-15  LUCCI        - MODERNIZATION OF THE CODE, INCLUDING:
!>                          * F77 TO F90/F95
!>                          * INDENTATION & UNIFORMIZATION CODE
!>                          * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                          * DOCUMENTATION WITH DOXYGEN
!>                          * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NHB  - UNIT NUMBER FOR READING THE NHB FILE
!> PD   - SFC PRESSURE MINUS PTOP
!> RES  - RECIPROCAL OF ETA AT THE GROUND
!> FIS  - SURFACE GEOPOTENTIAL
!> NTSD - THE TIMESTEP
!> NEST - IF A NESTED RUN THEN NEST IS .TRUE.
!>
!> OUTPUT ARGUMENT LIST:
!> PSLP - THE FINAL REDUCED SEA LEVEL PRESSURE ARRAY
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> T    - TEMPERATURE
!> Q    - SPECIFIC HUMIDITY
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
!> CALLS      : MPI_ALLREDUCE
!>              MPI_BCAST
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
    INTEGER(KIND=I4KIND), PARAMETER :: LP1    = LM  + 1
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM   = IM  * JM - JM / 2
    INTEGER(KIND=I4KIND), PARAMETER :: JM2    = JM  - 2
    INTEGER(KIND=I4KIND), PARAMETER :: KBI    = 2   * IM + JM - 3
    INTEGER(KIND=I4KIND), PARAMETER :: KBI2   = KBI - 4
    INTEGER(KIND=I4KIND), PARAMETER :: NRLX1  = 500
    INTEGER(KIND=I4KIND), PARAMETER :: NRLX2  = 100
    INTEGER(KIND=I4KIND), PARAMETER :: KSLPD  = 1
!
    REAL   (KIND=R4KIND), PARAMETER :: OVERRC = 1.50
    REAL   (KIND=R4KIND), PARAMETER :: AD05   = OVERRC * 0.05
    REAL   (KIND=R4KIND), PARAMETER :: CFT0   = OVERRC - 1.
    REAL   (KIND=R4KIND), PARAMETER :: ROG    = 287.04 / 9.8
!
    REAL   (KIND=R4KIND), ALLOCATABLE       , DIMENSION(:,:,:)                                  ::&
    & HTM

    REAL   (KIND=R4KIND)                    , DIMENSION(IM, MY_JSD:MY_JED)                      ::&
    & TTV     , PDSL1   , PBI     , SLPX
!
    REAL   (KIND=R4KIND)                    , DIMENSION(LM)                                     ::&
    & DETA    , RDETA   , AETA    , F4Q2    
!
    REAL   (KIND=R4KIND)                    , DIMENSION(LP1)                                    ::&
    & ETA     , DFL
!
    REAL   (KIND=R4KIND)                    , DIMENSION(KBI,LM)                                 ::&
    & TSLPB
!
    REAL   (KIND=R4KIND)                    , DIMENSION(KBI2,LM)                                ::& 
    & TSLPB2  
!
    REAL   (KIND=R4KIND)                    , DIMENSION(KBI)                                    ::&
    & PSLPB
!
    REAL   (KIND=R4KIND)                    , DIMENSION(KBI2)                                   ::&
    & PSLPB2
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & KMNTM
!
    INTEGER(KIND=I4KIND)                    , DIMENSION(IM, JM)                                 ::&
    & LMH
!
    INTEGER(KIND=I4KIND)                    , DIMENSION(IMJM)                                   ::&
    & IMNT    , JMNT
!
    INTEGER(KIND=I4KIND)                    , DIMENSION(JM)                                     ::&
    & IHE     , IHW     , IVE     , IVW
!
    LOGICAL(KIND=L4KIND)                                                       , INTENT(IN)     ::&
    & NEST
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & SIGMA   , STDRD   ,NO_REDUCE, EXBC
!
    REAL   (KIND=R4KIND)                    , DIMENSION(IM, MY_JSD:MY_JED)     , INTENT(IN)     ::&
    & PD      , RES     , FIS
!
    REAL   (KIND=R4KIND)                    , DIMENSION(IM, MY_JSD:MY_JED)     , INTENT(INOUT)    ::&
    & PSLP
!
    REAL   (KIND=R4KIND)                    , DIMENSION(IM, MY_JSD:MY_JED,LM)  , INTENT(INOUT)  ::&
    & T       , Q
!
    REAL   (KIND=R4KIND)                    , DIMENSION(IM, JM)                                 ::&
    & DUM
!
    REAL   (KIND=R4KIND)                    , DIMENSION(2 * KBI * (6 * LM + 1))                 ::&
    & BCDUM
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & LFRST
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                       , INTENT(IN)     ::&
    & NHB     , NTSD    
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NFCST   , NBC     , I       , J       , L       , LREC    , LXXX    , LHMNT   , IERR    ,   &
    & LMV     , IDTAD   , LIST    , LL      , LMAP1   , LRECBC  , NREC    , IRTN    , NRLX    ,   &
    & KMN     , LMST    , LI      , N       , KMM     , KM      , LMA     , KS      , IHH2    
! 
    REAL   (KIND=R4KIND)                                                                        ::&
    & DT      , DY      , CPGFV   , EN      , ENT     , R       , PT      , TDDAMP  , F4D     ,   &
    & F4Q     , EF4T    , TGSS    , BCHR    , PBIN    , TRTV    , PHTI    , DPOSP   , PRBIN   ,   &
    & PRTIN   , ALPP1   , SLOP    , SLPP    , TTT     , PTIN    , PHBI 
!
    DATA LFRST /.TRUE./
!
    SAVE
!
    IF (LFRST) THEN
        LFRST = .FALSE. 
!------------------------------------------------- 
! READ IN THE ARRAYS AND CONSTANTS THAT ARE NEEDED
!-------------------------------------------------
        REWIND NHB
!    
        READ(NHB) NFCST, NBC, LIST, DT, IDTAD, SIGMA
        READ(NHB) LMH
        READ(NHB) LMV
        READ(NHB)
        READ(NHB)
        READ(NHB)
        READ(NHB)
        READ(NHB)
!    
        STDRD     = .FALSE. 
        NO_REDUCE = .FALSE. 
!
        IF (SIGMA) STDRD = .TRUE. 
!-------------------------------------------------------- 
! FIND THE MOST ELEVATED GLOBAL LAYER WHERE THERE IS LAND
!--------------------------------------------------------    
        DO L=1,LM
            READ(NHB) ((DUM(I,J),I=1,IM),J=1,JM)
!        
            DO J=JSTA_I,JEND_I
                DO I=1,IM
                    IF (DUM(I,J) < 0.5) GOTO 666
                END DO
            END DO
!------------------------------------------         
! IF WE GET TO HERE, WE HAVE ALL ATMOSPHERE
!------------------------------------------        
            GOTO 667
!
        666 CONTINUE
!----------------------------------------------        
! IF WE GET HERE, WE HAVE FOUND A NON_ATM POINT
!----------------------------------------------        
            LHMNT = L
!
            GOTO 668
!
        667 CONTINUE
!
        END DO
!    
        LHMNT = LM+1
!
    668 CONTINUE
!    
        CALL MPI_ALLREDUCE(LHMNT, LXXX, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_COMP, IERR)
!
        LREC = MIN(LHMNT,LM)
        DO I=1,LREC-LXXX+1
            BACKSPACE NHB
        END DO
!    
        LHMNT = LXXX
!
        ALLOCATE(HTM(IM, MY_JSD:MY_JED, LHMNT:LM))
!    
        DO L=LHMNT,LM
            READ(NHB) ((DUM(I,J),I=1,IM),J=1,JM)
            DO J=MY_JSD,MY_JED
                DO I=1,IM
                    HTM(I,J,L) = DUM(I,J)
                END DO
            END DO
        END DO
!    
    669 CONTINUE
    
        DO L=1,LM
            READ(NHB)
        END DO
!    
        READ(NHB) DY , CPGFV, EN, ENT, R, PT, TDDAMP, F4D, F4Q, EF4T, DETA, RDETA, AETA, F4Q2,    &
    &             ETA, DFL
!---------------------- 
! END OF LFRST IF BLOCK
!---------------------- 
    END IF 
!------------------------------------------- 
! CALCULATE THE I-INDEX EAST-WEST INCREMENTS
!------------------------------------------- 
    DO J=1,JM
        IHE(J) = MOD(J+1,2)
        IHW(J) = IHE(J) - 1
        IVE(J) = MOD(J  ,2)
        IVW(J) = IVE(J) - 1
    END DO
!---------------------------------------------------------- 
! INITIALIZE ARRAYS.  LOAD SLP ARRAY WITH SURFACE PRESSURE.
!---------------------------------------------------------- 
    DO J=JSTA_I,JEND_I
        DO I=1,IM
            PSLP(I,J) = 0.
             TTV(I,J) = 0.
        END DO
    END DO
!
    DO 110 J=JSTA_I,JEND_I
        DO 110 I=1,IM
            PDSL1(I,J) =  RES(I,J) * PD(I,J)
             PSLP(I,J) =   PD(I,J) + PT
             PBI (I,J) = PSLP(I,J)
110 END DO
!---------------------------------------------------------------------------------------- 
! CALCULATE SEA LEVEL PRESSURE FOR PROFILES (AND POSSIBLY FOR POSTING BY POST PROCESSOR).
!
! "STDRD" REFERS TO THE "STANDARD" SLP REDUCTION SCHEME.
! THIS IS THE ONLY SCHEME AT PRESENT AVAILABLE FOR A SIGMA=.TRUE. ETA MODEL RUN.
!---------------------------------------------------------------------------------------- 
220 IF (LHMNT == LP1) THEN
        NO_REDUCE = .TRUE. 
    END IF
!
    IF (NO_REDUCE) GOTO 430
!
    IF (STDRD)     GOTO 400
!
    LL = LM
!
    DO J=JSTA_IM2,JEND_IM2
        DO I=2,IM-1
            IF (HTM(I,J,LL) <= 0.5) THEN
                TGSS  = FIS(I,J) / (R * ALOG((PDSL1(I,J) + PT) / (PD(I,J) + PT)))
                LMAP1 = LMH(I,J) + 1
!            
                DO 260 L=LMAP1,LM
                    T(I,J,L) = TGSS
            260 END DO
!           
            END IF
        END DO
    END DO
!--------------------------------------------------------------------------------------------------
! IF THIS IS A NESTED RUN, READ IN THE PARENT TEMPERATURE VALUES (INTERPOLATED THE NEST VERTICAL
! DISTRIBUTION) FOR THE SLP RELAXATION PROCEDURE.
!--------------------------------------------------------------------------------------------------
    IF (NEST) THEN
        IF (ME == 0) THEN
            LRECBC = 4 * (1 + (1 + 6 * LM) * KBI * 2 + (KBI + KBI2) * (LM + 1))
!
            OPEN(UNIT=NBC, ACCESS='DIRECT', RECL=LRECBC)
        
            NREC = NINT((NTSD-1) * DT / 3600.) + 2
!
            READ(NBC,REC=NREC) BCHR, BCDUM, TSLPB, TSLPB2, PSLPB, PSLPB2
!
            CLOSE(NBC)
        END IF
!    
        CALL MPI_BCAST(TSLPB , KBI*LM , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
        CALL MPI_BCAST(TSLPB2, KBI2*LM, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
        CALL MPI_BCAST(PSLPB , KBI    , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
        CALL MPI_BCAST(PSLPB2, KBI2   , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
!    
    END IF
!-----------------------------------------------------------------------------------------
! CREATE A TEMPORARY TV ARRAY, AND FOLLOW BY SEQUENTIAL OVERRELAXATION, DOING NRLX PASSES.
!-----------------------------------------------------------------------------------------
    NRLX = NRLX1
!
    DO 300 L=LHMNT,LM
!
          KMN = 0
        KMNTM = 0
!    
        DO 240 J=JSTA_IM2,JEND_IM2
            DO 240 I=2,IM-1
                IF (HTM(I,J,L) > 0.5) GOTO 240
                KMN = KMN+1
                IMNT(KMN) = I
                JMNT(KMN) = J
    240 END DO
!    
        KMNTM = KMN
!    
        DO 270 J=JSTA_I,JEND_I
            DO 270 I=1,IM
                TTV(I,J) = T(I,J,L)
    270 END DO  
!--------------------------------------------------------------------------------------------------
! FOR GRID BOXES NEXT TO MOUNTAINS REPLACE TTV BY AN "EQUIVALENT" TV, ONE WHICH CORRESPONDS TO THE
! CHANGE IN P BETWEEN REFERENCE INTERFACE GEOPOTENTIALS, INSTEAD OF BETWEEN LAYER INTERFACES
!--------------------------------------------------------------------------------------------------   
        DO J=JSTA_IM2,JEND_IM2
            DO I=2,IM-1
!
                IF (HTM(I,J,L) > 0.5 .AND.                                                        &
    &               HTM(I+IHW(J),J-1,L) * HTM(I+IHE(J),J-1,L) *                                   &
    &               HTM(I+IHW(J),J+1,L) * HTM(I+IHE(J),J+1,L) *                                   &
    &               HTM(I-1     ,J  ,L) * HTM(I+1     ,J  ,L) *                                   &
    &               HTM(I       ,J-2,L) * HTM(I       ,J+2,L) < 0.5) THEN
                    LMST = LMH(I,J)
!-------------------------------------------------------------
! FIND P AT THE REFERENCE INTERFACE GEOPOTENTIAL AT THE BOTTOM
!-------------------------------------------------------------
                    PBIN = PT + PD(I,J)
                    PHBI = DFL(LMST+1)
!                
                    DO LI=LMST,1,-1
                        PTIN = PBIN - DETA(LI) * PD(I,J) * RES(I,J)
                        TRTV = 2. * R * T(I,J,LI) * (1. + 0.608 * Q(I,J,LI))
                        PHTI = PHBI + TRTV * (PBIN - PTIN) / (PBIN + PTIN)
!
                        IF (PHTI >= DFL(L+1)) GO TO 273
!
                        PBIN = PTIN
                        PHBI = PHTI
                    END DO
!                
                273 DPOSP = (PHTI - DFL(L+1)) / TRTV
                    PRBIN = (1. + DPOSP) / (1. - DPOSP) * PTIN
!---------------------------------------------------------- 
! FIND P AT THE REFERENCE INTERFACE GEOPOTENTIAL AT THE TOP
!---------------------------------------------------------- 
                    PBIN = PT + PD(I,J)
                    PHBI = DFL(LMST+1)
!                
                    DO LI=LMST,1,-1
                        PTIN = PBIN - DETA(LI) * PD(I,J) * RES(I,J)
                        TRTV = 2. * R * T(I,J,LI) * (1. + 0.608 * Q(I,J,LI))
                        PHTI = PHBI + TRTV * (PBIN - PTIN) / (PBIN + PTIN)
!
                        IF (PHTI >= DFL(L)) GO TO 275
!
                        PBIN = PTIN
                        PHBI = PHTI
                    END DO
!                
                275 DPOSP = (PHTI - DFL(L)) / TRTV
                    PRTIN = (1. + DPOSP) / (1. - DPOSP) * PTIN
!                
                    TTV(I,J) = (DFL(L) - DFL(L+1)) / (2. * R) * (PRBIN + PRTIN) / (PRBIN - PRTIN)
                END IF
!
            END DO
        END DO
!--------------------------------------------------------------------------------------------------
! FOR POINTS IN THE OUTER TWO BOUNDARY ROWS THAT ARE NEXT TO MOUNTAINS, SET TTVS EQUAL TO THE 
! TEMPERATURES DERIVED EARLIER FROM THE PARENT GRIDS SLP AND SET PSLP EQUAL TO THE PARENTS SLP.
!--------------------------------------------------------------------------------------------------
        IF (NEST) THEN
!--------------------------------             
! FIRST DO THE OUTER BOUNDARY ROW
!--------------------------------        
            N = 1
!        
            DO I=1,IM
!            
                IF (JSTA_I == 1) THEN ! SOUTHERN EDGE
                    IF (HTM(I,3,L) < 0.5) THEN
                        TTV(I,1) = TSLPB(N,L)
                    END IF
!                
                    IF (L == LM) PSLP(I,1)=PSLPB(N)
                END IF
                N = N + 1
            END DO
!        
            DO I=1,IM
!            
                IF (JEND_I == JM) THEN ! NORTHERN EDGE
                    IF (HTM(I,JM-2,L) < 0.5) THEN
                        TTV(I,JM) = TSLPB(N,L)
                    END IF
!                
                    IF (L == LM) PSLP(I,JM) = PSLPB(N)
                END IF
                N = N + 1
            END DO
!        
            DO J=3,JM-2,2 ! WESTERN EDGE
!            
                IF (J >= JSTA_I .AND. J <= JEND_I) THEN
                    IF (HTM(3,J,L) < 0.5) THEN
                        TTV(1,J) = TSLPB(N,L)
                    END IF
!                
                    IF (L == LM) PSLP(1,J) = PSLPB(N)
                END IF
                N = N + 1
            END DO
!        
            DO J=3,JM-2,2 ! EASTERN EDGE
!            
                IF (J >= JSTA_I .AND. J <= JEND_I) THEN
                    IF (HTM(IM-2,J,L) < 0.5) THEN
                        TTV(IM,J) = TSLPB(N,L)
                    END IF
!                
                    IF (L == LM) PSLP(IM,J) = PSLPB(N)
                END IF
                N = N + 1
            END DO
!------------------------------------------        
! NOW DO THE INNER 2ND (INNER) BOUNDARY ROW
!------------------------------------------        
            N = 1
            DO I=1,IM-1 ! 1ST ROW FROM SOUTHERN EDGE
!            
                IF (JSTA_I == 1) THEN
                    IF (HTM(I,3,L) < 0.5 .OR. HTM(I+1,3,L) < 0.5) THEN
                        TTV(I,2) = TSLPB2(N,L)
                    END IF
!                
                    IF (L == LM) PSLP(I,2) = PSLPB2(N)
                END IF
                N = N + 1
            END DO
!--------------------------- 
! 1ST ROW FROM NORTHERN EDGE
!--------------------------- 
            DO I=1,IM-1 
!            
                IF (JEND_I == JM) THEN
                    IF (HTM(I,JM-2,L) < 0.5 .OR. HTM(I+1,JM-2,L) < 0.5) THEN
                        TTV(I,JM-1) = TSLPB2(N,L)
                    END IF
!                
                    IF (L == LM) PSLP(I,JM-1) = PSLPB2(N)
                END IF
                N = N + 1
            END DO
!-------------------------- 
! 1ST ROW FROM WESTERN EDGE
!-------------------------- 
            DO J=4,JM-3,2 
!    
                IF (J >= JSTA_I .AND. J <= JEND_I) THEN
                    IF (HTM(2,J-1,L) < 0.5 .OR. HTM(2,J+1,L) < 0.5) THEN
                        TTV(1,J) = TSLPB2(N,L)
                    END IF
!                
                    IF (L == LM) PSLP(1,J) = PSLPB2(N)
                END IF
                N = N + 1
            END DO
!--------------------------
! 1ST ROW FROM EASTERN EDGE
!--------------------------
            DO J=4,JM-3,2 
!            
                IF (J >= JSTA_I .AND. J <= JEND_I) THEN
                    IF (HTM(IM-1,J-1,L) < 0.5 .OR. HTM(IM-1,J+1,L) < 0.5) THEN
                        TTV(IM-1,J) = TSLPB2(N,L)
                    END IF
!                
                    IF (L == LM) PSLP(IM-1,J) = PSLPB2(N)
                END IF
                N = N + 1
            END DO
!------------------ 
! END OF NEST BLOCK
!------------------ 
        END IF 
!
        KMM = KMNTM
!---------------------------- 
! HERE IS THE RELAXATION LOOP
!---------------------------- 
        DO 285 N=1,NRLX
!        
            CALL UPDATE(TTV) ! EXCHANGE HALOES
!        
            DO 280 KM=1,KMM
                I = IMNT(KM)
                J = JMNT(KM)
                TTV(I,J) = AD05 * (4. * (TTV(I+IHW(J),J-1) + TTV(I+IHE(J),J-1)                    &
    &                    +               TTV(I+IHW(J),J+1) + TTV(I+IHE(J),J+1))                   &
    &                    +               TTV(I-1     ,J  ) + TTV(I+1     ,J  )                    &
    &                    +               TTV(I       ,J-2) + TTV(I       ,J+2))                   &
    &                    -               CFT0              * TTV(I       ,J  )
        280 END DO
!        
    285 END DO
!
        DO 290 KM=1,KMM
            I = IMNT(KM)
            J = JMNT(KM)
            T(I,J,L) = TTV(I,J)
    290 END DO
!
300 END DO
!-----------------------------------------------------------------
! CALCULATE THE SEA LEVEL PRESSURE AS PER THE NEW SCHEME.
!
! VALUES FOR IMNT AND JMNT ARE FOR LAYER LM - THIS IS WHAT WE WANT
!-----------------------------------------------------------------
    KMM = KMNTM
!
    DO 320 KM=1,KMM
        I = IMNT(KM)
        J = JMNT(KM)
        LMAP1 = LMH(I,J) + 1
        PBIN = PT + PD(I,J)
!    
        DO L=LMAP1,LM
            PTIN = PBIN
            DPOSP = (DFL(L) - DFL(L+1)) / (2. * R * T(I,J,L))
            PBIN = (1. + DPOSP) / (1. - DPOSP) * PTIN
        END DO
!    
        PSLP(I,J) = PBIN
320 END DO
!-------------------------- 
! SKIP THE STANDARD SCHEME.
!--------------------------
    GOTO 430
!-------------------------------------------------------------------------
! IF YOU WANT THE "STANDARD" ETA/SIGMA REDUCTION THIS IS WHERE IT IS DONE.
!-------------------------------------------------------------------------
400 CONTINUE
!
    DO 410 J=JSTA_I,JEND_I
        DO 410 I=1,IM
            IF (FIS(I,J) >= 1.) THEN
                LMA   = LMH(I,J)
                ALPP1 = ALOG(PDSL1(I,J) * ETA(LMA+1) + PT)
                SLOP  = 0.0065 * ROG * T(I,J,LMA)
!
                IF (SLOP < 0.50) THEN
                    SLPP = ALPP1 + FIS(I,J) / (R * T(I,J,LMA))
                ELSE
                    TTT = -(ALOG(PDSL1(I,J) * ETA(LMA) + PT) + ALPP1) * SLOP * 0.50 + T(I,J,LMA)
                    SLPP = (-TTT+SQRT(TTT*TTT + 2. * SLOP * (FIS(I,J) / R                         &
    &                    + (TTT + 0.50 * SLOP * ALPP1) * ALPP1))) / SLOP
                END IF
                PSLP(I,J) = EXP(SLPP)
            END IF
410 END DO
!--------------------------------------------------------------------------------------------------
! AT THIS POINT WE HAVE A SEA LEVEL PRESSURE FIELD BY EITHER METHOD. 5-POINT AVERAGE THE FIELD ON
! THE E-GRID.
!--------------------------------------------------------------------------------------------------
430 CONTINUE
!
    DO 440 J=JSTA_I,JEND_I
        DO 440 I=1,IM
            SLPX(I,J) = PSLP(I,J)
440 END DO
!
    DO 480 KS=1,KSLPD
!---------------- 
! EXCHANGE HALOES
!---------------- 
        CALL UPDATE(PSLP)
!    
        DO 460 J=JSTA_IM2,JEND_IM2
            IHH2 = IM-1 - MOD(J+1,2)
            DO 460 I=2,IHH2
!--------------------------------------------------------            
! EXTRA AVERAGING UNDER MOUNTAINS TAKEN OUT, FM, MARCH 96
!--------------------------------------------------------              
                SLPX(I,J) = 0.125 * (PSLP(I+IHW(J),J-1) + PSLP(I+IHE(J),J-1)                      &
    &                     +          PSLP(I+IHW(J),J+1) + PSLP(I+IHE(J),J+1)                      &
    &                     + 4.    *  PSLP(I       , J ))
    460 END DO
!    
        DO J=JSTA_I,JEND_I
            DO I=1,IM
                PSLP(I,J) = SLPX(I,J)
            END DO
        END DO
!    
480 END DO
!
    END SUBROUTINE SLP
