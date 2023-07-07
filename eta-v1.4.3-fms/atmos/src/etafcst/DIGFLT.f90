    SUBROUTINE DIGFLT
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE DIGFLT
!> 
!> SUBPROGRAM: DIGFLT - ?????
!> PROGRAMMER: ?????
!> ORG: ??????
!> DATE: ??-??-??       
!>     
!> ABSTRACT:
!> ?????
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????      - ORIGINATOR
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                       * F77 TO F90/F95
!>                       * INDENTATION & UNIFORMIZATION CODE
!>                       * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                       * DOCUMENTATION WITH DOXYGEN
!>                       * OPENMP FUNCTIONALITY
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
!>              EXCHM
!>              F77KINDS
!>              GLB_TABLE
!>              LOOPS
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              PARMETA
!>              PARM_TBL
!>              PVRBLS
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>
!> DRIVER     : EBU
!>
!> CALLS      : BOCOHF
!>              BOCOV
!>              DIVHOA
!>              DSTRB
!>              EXCH
!>              HZADV
!>              LOC2GLB
!>              MPI_BARRIER
!>              MPI_BCAST
!>              PDTEDT
!>              PGCOR
!>              VTADVF
!>              ZERO2            
!>--------------------------------------------------------------------------------------------------   
    USE CONTIN
    USE CTLBLK
    USE DYNAM
    USE EXCHM
    USE F77KINDS
    USE GLB_TABLE
    USE LOOPS    , ONLY : JAM
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE PARMETA
    USE PARM_TBL
    USE PVRBLS
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    INCLUDE "EXCHM.h"
    INCLUDE "mpif.h"
!
#include "sp.h"
!
    REAL   (KIND=R4KIND), PARAMETER :: PIQ  =    3.141592654
!--------------------------------------- 
! NTIM IS A HALF-SPAN IN TIME STEP UNITS
! ACTUAL TIME IS NTIM*DT
!--------------------------------------- 
    INTEGER(KIND=I4KIND), PARAMETER :: NTIM =  110
!
    REAL   (KIND=R4KIND), PARAMETER :: CM1  = 2937.4
    REAL   (KIND=R4KIND), PARAMETER :: CM2  =    4.9283 
    REAL   (KIND=R4KIND), PARAMETER :: CM3  =   23.5518
    REAL   (KIND=R4KIND), PARAMETER :: EPS  =    0.622
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & PDQ     , RELHUM  , PDSL0
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & TQ      , UQ      , VQ
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & KNT     , IUNRH   , IUNDF   , I       , J       , K       , IRTN    , NT      , JK
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & CLOGES  , ESE     , PRS     , QIJ     , E       , HSPAN   , TIME    , QNT     , TETC    ,   &
    & CSTT    , SQ      , FNT     , AS1     , AS2     , ASS1    , ASS2    , BNT     , BS1     ,   &
    & BS2     , BSS1    , BSS2    , HNTSD   , TSPH    , BCHR
!
    DATA KNT/0/, IUNRH/51/, IUNDF/52/
!
    NBOCO = 0
!
    REWIND IUNRH
    REWIND IUNDF
!-------------------------------------------- 
! COMPUTE AND WRITE OUT THE RELATIVE HUMIDITY
!--------------------------------------- ----
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            PDSL0(I,J) = RES(I,J) * PD(I,J)
        END DO
    END DO
!
    DO K=1,LM
!    
        CALL ZERO2(RELHUM)
!    
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                IF (HTM(I,J,K) > 0.5) THEN
                    CLOGES = - CM1 / T(I,J,K) - CM2 * ALOG10(T(I,J,K)) + CM3
                    ESE = 10. ** (CLOGES + 2.)
                    PRS = AETA(K) * PDSL0(I,J) + PT
                    QIJ = Q(I,J,K)
                    E   = PRS * QIJ / (EPS * (1.-QIJ) + QIJ)
                    RELHUM(I,J) = AMIN1(E/ESE,0.98)
                ELSE
                    RELHUM(I,J) = 0.
                END IF
            END DO
        END DO
!-----------------------------------------
! SMOOTH THE INITIAL RH FIELD THEN SAVE IT
!--------------------------------------- --
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                RELHUM(I,J) = AMIN1(RELHUM(I,J),0.98)
                RELHUM(I,J) = AMAX1(RELHUM(I,J),0.02)
            END DO
        END DO
!-----------------------  
! ASSEMBLE GLOBAL RELHUM
!-----------------------  
        CALL LOC2GLB(RELHUM,TEMP1)
!    
        IF (MYPE == 0) WRITE(IUNRH) TEMP1
!-------------------------------------------------------------------------     
! SMOOTH THE INITIAL TEMPERATURE FIELD BEFORE EXECUTING THE DIGITAL FILTER
!-------------------------------------------------------------------------    
    END DO
!------------------------------------------------------- 
! SAVE CURRENT PROG VARIABLES IN WORK FILE FOR LATER USE
!--------------------------------------- ---------------
    CALL LOC2GLB(PD, TEMP1)
!
    IF (MYPE == 0) WRITE(IUNDF) TEMP1
!
    DO K=1,LM
!
        CALL LOC2GLB(T(IDIM1,JDIM1,K), TEMP1)
!
        IF(MYPE == 0)  WRITE(IUNDF) TEMP1
!    
        CALL LOC2GLB(U(IDIM1,JDIM1,K), TEMP1)
!
        IF (MYPE == 0) WRITE(IUNDF) TEMP1
!    
        CALL LOC2GLB(V(IDIM1,JDIM1,K), TEMP1)
!
        IF (MYPE == 0) WRITE(IUNDF) TEMP1
!
    END DO
!---------------------------------------------- 
! SET UP ARRAYS TO HOLD SUMS FOR DIGITAL FILTER
!----------------------------------------------
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            PDQ(I,J) = PD(I,J)
        END DO
    END DO
!
    DO K=1,LM
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                TQ(I,J,K) = T(I,J,K)
                UQ(I,J,K) = U(I,J,K)
                VQ(I,J,K) = V(I,J,K)
            END DO
        END DO
    END DO
!
    HSPAN = FLOAT(NTIM) * DT / 3600.
!
    IF (MYPE == 0) WRITE(6,100) HSPAN
    100 FORMAT(' ','INITIALIZATION CALLED WITH HALF-SPAN',F5.1,' HOURS')
!
    CALL MPI_BARRIER(MPI_COMM_COMP, IRTN)
!-------------------------------------------------------------
! RUN THE FORECAST MODEL FORWARD AND BACKWARD (DIGITAL FILTER)
!
! FIRST, FORWARD INTEGRATION
!
!--------------------------------------------------------------------------------------------------
! CALCULATE NORM NEEDED FOR LANCZOS WINDOW
!
! 'TETC' DEFINES CUT-OFF FREQUENCY FOR THE WINDOW (FACTOR 2 APPEARS ONLY TO SHOW GENERAL FORMULA) 
! (IT SHOULD BE TETC=PIQ/FLOAT(NTIM))
! 'NTIM' IS A NUMBER OF TIME STEPS IN A HALF-SPAN
! 'TIME' CORRESPONDS TO NUMBER OF TIME STEPS OF A SPAN
!--------------------------------------------------------------------------------------------------
    TIME = 2. * FLOAT(NTIM)
    QNT  =      FLOAT(NTIM+1)
    TETC = 2. * PIQ / TIME
    CSTT =     TETC / PIQ
    SQ   =     CSTT

    DO NT=1,NTIM
        FNT  = FLOAT(NT)
        AS1  = PIQ * FNT / QNT
        AS2  = FNT * TETC
        ASS1 = SIN(AS1) / AS1
        ASS2 = SIN(AS2) / AS2
        SQ   = SQ + 2. * CSTT * ASS1 * ASS2
    END DO

!-------------------------------------
! NORMALIZATION OF THE WINDOW FUNCTION
!-------------------------------------
    CSTT = CSTT / SQ
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            PDQ(I,J) = PDQ(I,J) * CSTT
        END DO
    END DO
!
    DO K=1,LM
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                TQ(I,J,K) = TQ(I,J,K) * CSTT
                UQ(I,J,K) = UQ(I,J,K) * CSTT
                VQ(I,J,K) = VQ(I,J,K) * CSTT
            END DO
        END DO
    END DO
!
    NTSD = 0
!------------------------- 
! ENTRY INTO THE TIME LOOP
!------------------------- 
    2010 NTSD = NTSD + 1
          KNT = KNT  + 1
    IF (MYPE == 0) WRITE(6,2015) NTSD, NTIM
    2015 FORMAT(' NTSD=',I5,'  NTSTM=',I4)
!-------------------------------------------------------
! DIVERGENCE AND HORIZONTAL PART OF THE OMEGA-ALPHA TERM
!-------------------------------------------------------
    IF (NTSD > 1) CALL EXCH(T, LM, U, LM, V, LM, 2, 2) !EXCHANGE T, U AND V
    CALL DIVHOA
!----------------------------------------------
! PRESS. TEND.,ETA DOT AND VERTICAL OMEGA-ALPHA
!----------------------------------------------
    CALL PDTEDT
!------------------- 
! VERTICAL ADVECTION
!-------------------
    IF (MOD(NTSD-1,IDTAD) == 0) THEN
        CALL EXCH(ETADT, LM-1, 1, 1)
!
        CALL VTADVF
!    
        CALL EXCH(T, LM, U, LM, V, LM, Q2, LM, 1, 1)
    END IF
!------------------------------------- 
! UPDATE SURFACE PRESSURE (MINUS PTOP)
!------------------------------------- 
    CALL PDNEW
!------------------------- 
! UPDATE H BOUNDARY POINTS
!------------------------- 
    IF (MOD(NTSD,IDTAD) == 0) THEN
        CALL EXCH(T, LM, Q2, LM, 1, 1)
    END IF
!
    CALL EXCH(PD, 1, 1, 1)
!
    CALL BOCOHF
!-------------------------------------------
! PRESSURE GRADIENT AND CORIOLIS FORCE TERMS
!-------------------------------------------
    CALL EXCH(PD, 1, T, LM, 2, 2) !EXCHANGE PD AND T
!
    CALL PGCOR
!
    CALL EXCH(PDSL, 1, 5, 5)
!------------------------- 
! UPDATE V BOUNDARY POINTS
!------------------------- 
    CALL EXCH(U, LM, V, LM, 1, 1) !EXCHANGE U AND V
!
    CALL BOCOV
!---------------------
! HORIZONTAL ADVECTION
!---------------------
    IF (MOD(NTSD,IDTAD) == 0) THEN
        CALL EXCH(T, LM, U, LM, V, LM, 4, 4) !EXCHANGE T, U, AND V
!        
        CALL EXCH(Q2, LM, 5, 5)
!    
        CALL HZADV
 !   
        CALL EXCH(U, LM, V, LM, 2, 2) !EXCHANGE U AND V
    END IF
!--------------------------- 
! CALCULATE WINDOW PARAMETER
!---------------------------
    BNT   = FLOAT(NTSD)
    BS1   = BNT * PIQ / QNT
    BS2   = BNT * TETC
    BSS1  = SIN(BS1) / BS1
    BSS2  = SIN(BS2) / BS2
    HNTSD = CSTT * BSS1 * BSS2
!
    DO K=1,LM
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                TQ(I,J,K) = TQ(I,J,K) + HNTSD * T(I,J,K)
                UQ(I,J,K) = UQ(I,J,K) + HNTSD * U(I,J,K)
                VQ(I,J,K) = VQ(I,J,K) + HNTSD * V(I,J,K)
            END DO
        END DO
    END DO
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            PDQ(I,J) = PDQ(I,J) + HNTSD * PD(I,J)
        END DO
    END DO
!
    IF (NTSD == NTIM) GOTO 2013
!
    GOTO 2010
!------------------------ 
! EXIT FROM THE TIME LOOP 
!------------------------
    2013 CONTINUE
!
    CALL MPI_BARRIER(MPI_COMM_COMP, IRTN)
!---------------------------------------------------------
! NOW BACKWARD INTEGRATION, STARTING FROM THE INITIAL TIME
!---------------------------------------------------------
!
!-----------------------------------------------
! CHANGE (SIGN ONLY OF) IMPORTANT TIME CONSTANTS
!-----------------------------------------------
    DT    = -DT
    CPGFV = -CPGFV
    EN    = -EN
    ENT   = -ENT
    F4D   = -F4D
    F4Q   = -F4Q
    EF4T  = -EF4T
!
    DO JK=1,JAM
        EM (JK) = - EM(JK)
        EMT(JK) = -EMT(JK)
    END DO
!
    DO K=1,LM
        F4Q2(K) = -F4Q2(K)
    END DO
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            WPDAR(I,J) = -WPDAR(I,J)
            CPGFU(I,J) = -CPGFU(I,J)
             CURV(I,J) =  -CURV(I,J)
              FCP(I,J) =   -FCP(I,J)
              FAD(I,J) =   -FAD(I,J)
                F(I,J) =     -F(I,J)
        END DO
    END DO
!
    CALL EXCH(WPDAR, 1, CPGFU, 1, 2, 2)
    CALL EXCH(CURV , 1, FCP  , 1, 2, 2)
    CALL EXCH(FAD  , 1, F    , 1, 2, 2)
!-------------- 
! END OF CHANGE 
!--------------  
!
!-----------------------------------------------------------------------------------------------
! DEFINE INITIAL DATA FOR BACKWARD INTEGRATION (PDQ,TQ,UQ,VQ ARE SUMS NEEDED FOR DIGITAL FILTER)
!-----------------------------------------------------------------------------------------------
    IF (MYPE == 0) THEN
        REWIND IUNDF
        READ(IUNDF) TEMP1
    END IF
!
    CALL DSTRB(TEMP1, PD, 1, 1, 1)
!
    DO K=1,LM
        IF (MYPE == 0) THEN
            READ(IUNDF) TEMP1
            READ(IUNDF) TEMP2
            READ(IUNDF) TEMP3
        END IF
!
        CALL DSTRB(TEMP1, T, 1, LM, K)
        CALL DSTRB(TEMP2, U, 1, LM, K)
        CALL DSTRB(TEMP3, V, 1, LM, K)
    END DO
!
    CALL EXCH(T , LM, U, LM, V, LM, 5, 5) !EXCHANGE T, U, AND V
    CALL EXCH(PD,  1, 5, 5)
!
    NTSD = 0
    FIRST = .TRUE. 
    TSPH = -3600. / DT
    IF (MYPE == 0) THEN
        IF (.NOT. NEST) THEN
            REWIND NBC
            READ(NBC)
            READ(NBC) BCHR
        ELSE
            READ(NBC,REC=2) BCHR
        END IF
    END IF
!
    CALL MPI_BCAST(BCHR, 1, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    NBOCO = INT(BCHR * TSPH + 0.5)
!------------------------- 
! ENTRY INTO THE TIME LOOP
!------------------------- 
    2020 NTSD = NTSD + 1
          KNT = KNT  + 1
!
    IF (MYPE == 0) WRITE(6,2015) NTSD, NTIM
!-------------------------------------------------------
! DIVERGENCE AND HORIZONTAL PART OF THE OMEGA-ALPHA TERM
!-------------------------------------------------------
    IF (NTSD > 1) CALL EXCH(T, LM, U, LM, V, LM, 2, 2) !EXCHANGE T, U, AND V
!
    CALL DIVHOA
!
!--------------------------------------------
! PRESS. TEND.,ETA DOT & VERTICAL OMEGA-ALPHA
!--------------------------------------------
    CALL PDTEDT
!------------------- 
! VERTICAL ADVECTION
!------------------- 
    IF (MOD(NTSD-1,IDTAD) == 0) THEN
        CALL EXCH(ETADT, LM-1, 1, 1)
!    
        CALL VTADVF
!    
        CALL EXCH(T, LM, U, LM, V, LM, Q2, LM, 1, 1) !EXCHANGE T, U, AND V
    END IF
!------------------------------------- 
! UPDATE SURFACE PRESSURE (MINUS PTOP)
!-------------------------------------
    CALL PDNEW
!------------------------- 
! UPDATE H BOUNDARY POINTS
!------------------------- 
    IF (MOD(NTSD,IDTAD) == 0) THEN
        CALL EXCH(T, LM, Q2, LM, 1, 1)
    END IF
!
    CALL EXCH(PD, 1, 1, 1)
!
    CALL BOCOHF
!-------------------------------------------
! PRESSURE GRADIENT AND CORIOLIS FORCE TERMS
!-------------------------------------------
    CALL EXCH(PD, 1, T, LM, 2, 2) !EXCHANGE PD AND T
!
    CALL PGCOR
!
    CALL EXCH(PDSL, 1, 5, 5)
!------------------------- 
! UPDATE V BOUNDARY POINTS
!-------------------------
    CALL EXCH(U, LM, V, LM, 1, 1) !EXCHANGE U AND V
!
    CALL BOCOV
!--------------------- 
! HORIZONTAL ADVECTION
!--------------------- 
    IF (MOD(NTSD,IDTAD) == 0) THEN
!    
        CALL EXCH(T, LM, U, LM, V, LM, 4, 4) !EXCHANGE T, U, AND V
!
        CALL EXCH(Q2, LM, 5, 5)
!    
        CALL HZADV
!    
        CALL EXCH(U, LM, V, LM, 2, 2) !EXCHANGE U AND V
    END IF
!---------------------------
! CALCULATE WINDOW PARAMETER
!---------------------------
    BNT   = FLOAT(NTSD)
    BS1   = BNT * PIQ / QNT
    BS2   = BNT * TETC
    BSS1  = SIN(BS1) / BS1
    BSS2  = SIN(BS2) / BS2
    HNTSD = CSTT * BSS1 * BSS2
!
    DO K=1,LM
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                TQ(I,J,K) = TQ(I,J,K) + HNTSD * T(I,J,K)
                UQ(I,J,K) = UQ(I,J,K) + HNTSD * U(I,J,K)
                VQ(I,J,K) = VQ(I,J,K) + HNTSD * V(I,J,K)
            END DO
        END DO
    END DO
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            PDQ(I,J) = PDQ(I,J) + HNTSD * PD(I,J)
        END DO
    END DO
!
    IF (NTSD == NTIM) GOTO 2022
!
    GOTO 2020
!------------------------
! EXIT FROM THE TIME LOOP
!------------------------
    2022 CONTINUE
!-------------------------------------------------
! CHANGE BACK (SIGN ONLY) IMPORTANT TIME CONSTANTS
!-------------------------------------------------
    DT    = -DT
    CPGFV = -CPGFV
    EN    = -EN
    ENT   = -ENT
    F4D   = -F4D
    F4Q   = -F4Q
    EF4T  = -EF4T
!
    DO JK=1,JAM
         EM(JK) = - EM(JK)
        EMT(JK) = -EMT(JK)
    END DO
!
    DO K=1,LM
        F4Q2(K) = -F4Q2(K)
    END DO
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            WPDAR(I,J) = -WPDAR(I,J)
            CPGFU(I,J) = -CPGFU(I,J)
             CURV(I,J) =  -CURV(I,J)
              FCP(I,J) =   -FCP(I,J)
              FAD(I,J) =   -FAD(I,J)
                F(I,J) =     -F(I,J)
        END DO
    END DO
!
    CALL EXCH(WPDAR, 1, CPGFU, 1, 2, 2)
!
    CALL EXCH(CURV , 1, FCP  , 1, 2, 2)
!
    CALL EXCH(FAD  , 1, F    , 1, 2, 2)
!----------------------------------------
! UPDATE INITIALIZED PROGNOSTIC VARIALBES
!----------------------------------------
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            PD(I,J) = PDQ(I,J)
        END DO
    END DO
!
    DO K=1,LM
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                T(I,J,K) = TQ(I,J,K)
                U(I,J,K) = UQ(I,J,K)
                V(I,J,K) = VQ(I,J,K)
            END DO
        END DO
    END DO
!---------------------------------------------------------------------------------------------------------
! RETRIEVE THE INITIAL RELATIVE HUMIDITY AND COMPUTE Q SO AS TO MAINTAIN THE RH GIVEN THE NEW TEMPERATURES
!---------------------------------------------------------------------------------------------------------
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            PDSL0(I,J) = RES(I,J) * PD(I,J)
        END DO
    END DO
!
    IF (MYPE == 0) REWIND IUNRH
!
    DO K=1,LM
        IF (MYPE == 0) READ(IUNRH) TEMP1
!
        CALL DSTRB(TEMP1, RELHUM, 1, 1, 1)
!    
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                IF (HTM(I,J,K) > 0.5) THEN
                    CLOGES   = -CM1 / T(I,J,K) - CM2 * ALOG10(T(I,J,K)) + CM3
                    ESE      = 10. ** (CLOGES + 2.)
                    PRS      = AETA(K) * PDSL0(I,J) + PT
                    E        = RELHUM(I,J) * ESE
                    Q(I,J,K) = EPS * E / (PRS + E * (EPS - 1.))
                END IF
            END DO
        END DO
    END DO
!
    CALL EXCH(T, LM, Q, LM, U, LM, V, LM, 4, 4) !EXCHANGE T, Q, U, AND V
!
    CALL EXCH(PD, 1, 5, 5)
!---------------------------------
! RETURN BC FILE TO START FORECAST
!---------------------------------
    NTSD = 0
    IF (MYPE == 0) THEN
        IF (.NOT. NEST) THEN
            REWIND NBC
            READ(NBC)
            READ(NBC) BCHR
        ELSE
            READ(NBC,REC=2) BCHR
        END IF
    END IF
!
    CALL MPI_BCAST(BCHR, 1, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
!
    NBOCO = INT(BCHR * TSPH + 0.5)
!-------------- 
! END OF CHANGE  
!-------------- 
    RETURN
!
    END SUBROUTINE DIGFLT
