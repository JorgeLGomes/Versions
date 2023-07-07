    SUBROUTINE NEWFLT
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE NEWFLT
!>
!> SUBROUTINE: NEWFLT - ?????
!> PROGRAMMER: ?????
!> ORG: ?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> ?????
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????    - ORIGINATOR
!> 18-01-15  LUCCI    - MODERNIZATION OF THE CODE, INCLUDING:
!>                      * F77 TO F90/F95
!>                      * INDENTATION & UNIFORMIZATION CODE
!>                      * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                      * DOCUMENTATION WITH DOXYGEN
!>                      * OPENMP FUNCTIONALITY
!>
!> INPUT  ARGUMENT LIST:
!> NONE
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!> 
!> USE MODULES: BOCO
!>              CLDWTR
!>              CONTIN
!>              CTLBLK
!>              DYNAM
!>              EXCHM
!>              F77KINDS
!>              GLB_TABLE
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
!> DRIVER     : -----
!>
!> CALLS      : BOCOHF
!>              BOCOV
!>              DDAMP
!>              DIVHOA
!>              EXCH
!>              EXIT
!>              HDIFF
!>              HZADV
!>              HZADV2
!>              MPI_BARRIER
!>              MPI_BCAST
!>              PDNEW
!>              PDTE
!>              PGCOR
!>              PRECPD
!>              RADTN
!>              RDTEMP
!>              TURBL
!>              VTADV
!>              VTADVF             
!>-------------------------------------------------------------------------------------------------- 
    USE BOCO
    USE CLDWTR
    USE CONST
    USE CONTIN
    USE CTLBLK
    USE DYNAM
    USE EXCHM
    USE F77KINDS
    USE GLB_TABLE
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
    REAL   (KIND=R4KIND), PARAMETER :: PIQ = 3.141592654
!
    INTEGER(KIND=I4KIND), PARAMETER :: LP1  = LM +  1
    INTEGER(KIND=I4KIND), PARAMETER :: JAM  =  6 +  2 * (JM - 10)  
    INTEGER(KIND=I4KIND), PARAMETER :: LB   =  2 * IM +  JM -  3 
!---------------------------------------------------------------
! NTIM IS A SPAN IN TIME STEP UNITS FOR THE BACKWARD INTEGRATION
!---------------------------------------------------------------
    INTEGER(KIND=I4KIND), PARAMETER :: NTIM =  270
!
    REAL   (KIND=R4KIND), PARAMETER :: CM1  = 2937.4
    REAL   (KIND=R4KIND), PARAMETER :: CM2  =    4.9283
    REAL   (KIND=R4KIND), PARAMETER :: CM3  =   23.5518
    REAL   (KIND=R4KIND), PARAMETER :: EPS  =    0.622
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & DIVHOA_STATUS     ,                                                                         &
    &   PDTE_STATUS     ,                                                                         &
    & VTADVF_STATUS     ,                                                                         &
    &  PGCOR_STATUS     ,                                                                         &
    &  HZADV_STATUS 
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & RUN     , FIRST   , RESTRT  , SIGMA
!
    DATA KNT/0/, IUNRH/51/, IUNDF/52/
!
    NBOCO = 0
!-----------------------------------------------------------------
! SMOOTH THE INITIAL TEMPERATURE FIELD BEFORE EXECUTING THE FILTER
!-----------------------------------------------------------------
!
!    CALL FILT25(T(IDIM1,JDIM1,L),HTM(IDIM1,JDIM1,L),5)
!
    SPAN = FLOAT(NTIM) * DT / 3600.
!
    IF (MYPE == 0) WRITE(6,100) SPAN
    100 FORMAT(' ','INITIALIZATION CALLED WITH SPAN',F5.1,' HOURS')
!------------------------------------ 
! RUN THE MODEL BACKWARD THEN FORWARD
!------------------------------------ 
!
!--------------------------------------------------------------- 
! ADIABATIC BACKWARD INTEGRATION, STARTING FROM THE INITIAL TIME
!--------------------------------------------------------------- 
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
        EM (JK) = -EM (JK)
        EMT(JK) = -EMT(JK)
    END DO
!
    DO L=1,LM
        F4Q2(L) = -F4Q2(L)
    END DO
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            WPDAR(I,J) = -WPDAR(I,J)
            CPGFU(I,J) = -CPGFU(I,J)
            CURV (I,J) = -CURV (I,J)
            FCP  (I,J) = -FCP  (I,J)
            FAD  (I,J) = -FAD  (I,J)
            F    (I,J) = -F    (I,J)
        END DO
    END DO
!
    NTSD  = 0
    FIRST = .TRUE. 
    TSPH  = -3600. / DT
!------------------------- 
! ENTRY INTO THE TIME LOOP 
!------------------------- 
101 NTSD = NTSD + 1
    KNT  = KNT  + 1
!
    IF (MYPE == 0) WRITE(6,2015) NTSD, NTIM
!------------------------------------------------------- 
! DIVERGENCE AND HORIZONTAL PART OF THE OMEGA-ALPHA TERM
!-------------------------------------------------------
!
!--------------------- 
! EXCHANGE T, U, AND V
!--------------------- 
    IF (NTSD > 1) CALL EXCH(T, LM, U, LM, V, LM, 2, 2) 
!
    CALL DIVHOA
!
    CALL EXIT(DIVHOA_STATUS)
!----------------------------------------------- 
! PRESS. TEND., ETA DOT AND VERTICAL OMEGA-ALPHA
!-----------------------------------------------
    CALL PDTE
!
    CALL EXIT(PDTE_STATUS)
!-------------------  
! VERTICAL ADVECTION
!-------------------
    IF (MOD(NTSD-1,IDTAD) == 0) THEN
        CALL EXCH(ETADT, LM-1, 1, 1)
!
        CALL VTADVF
!
        CALL EXIT(VTADVF_STATUS)
!--------------------- 
! EXCHANGE T, U, AND V
!--------------------- 
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
!
!------------------  
! EXCHANGE PD AND T
!------------------ 
    CALL EXCH(PD, 1, T, LM, 2, 2) 
!
    CALL PGCOR
!
    CALL EXIT(PGCOR_STATUS)
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
!--------------------- 
! EXCHANGE T, U, AND V
!--------------------- 
        CALL EXCH(T , LM, U, LM, V, LM, 4, 4) 
!
        CALL EXCH(Q2, LM, 5, 5)
!
        CALL HZADV
!
        CALL EXIT(HZADV_STATUS)
!-----------------  
! EXCHANGE U AND V
!-----------------    
        CALL EXCH(U, LM, V, LM, 2, 2) 
    END IF
!
    IF (NTSD == NTIM) GOTO 200
!
    GOTO 101
!------------------------ 
! EXIT FROM THE TIME LOOP
!------------------------ 
200 CONTINUE
!------------------------------ 
! READY FOR FORWARD INTEGRATION
!------------------------------ 
    NTSD  = 0
    FIRST = .TRUE. 
    TSPH  = -3600. / DT
!
    IF (MYPE == 0) THEN
!
        REWIND NBC
!
        READ(NBC)
        READ(NBC) BCHR
        READ(NBC) PDB
        READ(NBC) TB
        READ(NBC) QB
        READ(NBC) UB
        READ(NBC) VB
        READ(NBC) Q2B
        READ(NBC) CWMB
    END IF
!
    CALL MPI_BCAST(BCHR       , 1      , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST( PDB(1,1)  , 2*LB   , MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST(  TB(1,1,1), 2*LB*LM, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST(  QB(1,1,1), 2*LB*LM, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST(  UB(1,1,1), 2*LB*LM, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST(  VB(1,1,1), 2*LB*LM, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST( Q2B(1,1,1), 2*LB*LM, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
    CALL MPI_BCAST(CWMB(1,1,1), 2*LB*LM, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
!
    CALL MPI_BARRIER(MPI_COMM_COMP, IRTN)
!
    IF (MYPE == 0) WRITE(LIST,*)'  READ UNIT NBC=', NBC
!------------------------------------------------- 
! COMPUTE THE 1ST TIME FOR BOUNDARY CONDITION READ
!------------------------------------------------- 
    NBOCO = INT(BCHR * TSPH + 0.5)

    IF (MYPE == 0) WRITE(6,*)' NBOCO=', NBOCO
!
    IF (NTSD == NBOCO) THEN
        IF (MYPE == 0) THEN
            BACKSPACE NBC
            BACKSPACE NBC
            BACKSPACE NBC
            BACKSPACE NBC
            BACKSPACE NBC
            BACKSPACE NBC
            BACKSPACE NBC
            WRITE(LIST,*)'  BACKSPACE UNIT NBC=', NBC
        END IF
    END IF
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
        EM (JK) = -EM (JK)
        EMT(JK) = -EMT(JK)
    END DO
!
    DO L=1,LM
        F4Q2(L) = -F4Q2(L)
    END DO
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            WPDAR(I,J) = -WPDAR(I,J)
            CPGFU(I,J) = -CPGFU(I,J)
            CURV (I,J) = -CURV (I,J)
            FCP  (I,J) = -FCP  (I,J)
            FAD  (I,J) = -FAD  (I,J)
            F    (I,J) = -F    (I,J)
        END DO
    END DO
!------------------------------------ 
! INTEGRATE FORWARD WITH FULL PHYSICS
!------------------------------------ 
    NTSD = 0
!------------------------- 
! ENTRY INTO THE TIME LOOP 
!------------------------- 
 300 NTSD = NTSD + 1
      KNT = KNT  + 1
!
    IF (MYPE == 0) WRITE(6,2015) NTSD, NTIM
    2015 FORMAT(' NTSD=',I5,'  NTSTM=',I4)
!---------- 
! RADIATION 
!---------- 
    IF (MOD(NTSD-1,NRADS) == 0 .OR. MOD(NTSD-1,NRADL) == 0) THEN
        CALL RADTN
    END IF
!------------------------------------------------------- 
! DIVERGENCE AND HORIZONTAL PART OF THE OMEGA-ALPHA TERM
!------------------------------------------------------- 
!--------------------- 
! EXCHANGE T, U, AND V
!--------------------- 
    IF (NTSD > 1) CALL EXCH(T, LM, U, LM, V, LM, 2, 2)
!
    CALL DIVHOA
!---------------------------------------------- 
! PRESS. TEND.,ETA DOT AND VERTICAL OMEGA-ALPHA
!----------------------------------------------
    CALL PDTE
!------------------- 
! VERTICAL ADVECTION
!-------------------
    IF (MOD(NTSD-1,IDTAD) == 0) THEN
        CALL EXCH(ETADT, LM-1, 1, 1)
!
        CALL VTADV
!--------------------- 
! EXCHANGE T, U, AND V
!---------------------     
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
        CALL EXCH(T, LM, Q, LM, Q2, LM, 1, 1)
    END IF
!
    CALL EXCH(PD, 1, 1, 1)
!
    CALL MPI_BARRIER(MPI_COMM_COMP, ISTAT)
!
    CALL EXCH(CWM, LM, 1, 1)
!
    CALL BOCOH
!-------------------------------------------
! PRESSURE GRADIENT AND CORIOLIS FORCE TERMS
!-------------------------------------------
!------------------
! EXCHANGE PD AND T
!------------------
    CALL EXCH(PD, 1, T, LM, Q, LM, 2, 2)
!
    CALL PGCOR
!
    CALL EXCH(PDSL, 1, 5, 5)
!------------------- 
! DIVERGENCE DAMPING 
!------------------- 
    IF (MOD(NTSD,NTDDMP) == 0) THEN
!--------------------- 
! EXCHANGE T, U, AND V
!--------------------- 
        CALL EXCH(T, LM, U, LM, V, LM, DIV, LM, 1, 1)
!    
        CALL DDAMP
    END IF
!------------------------- 
! UPDATE V BOUNDARY POINTS
!------------------------- 
!--------------------- 
! EXCHANGE U, AND V
!--------------------- 
    CALL EXCH(U, LM, V, LM, 1, 1)
!
    CALL BOCOV
!-------------------------------------------- 
! APPLY TEMPERATURE TENDENCY DUE TO RADIATION 
!-------------------------------------------- 
    CALL RDTEMP
!------------------ 
! LATERAL DIFFUSION 
!------------------ 
!----------------------- 
! EXCHANGE T, U, V AND Q
!----------------------- 
    CALL EXCH(T , LM, U, LM, V, LM, Q, LM, 2, 2)
!
    CALL EXCH(Q2, LM, 1, 1)
!
    CALL HDIFF
!---------------------  
! HORIZONTAL ADVECTION
!--------------------- 
    IF (MOD(NTSD,IDTAD) == 0) THEN
!--------------------- 
! EXCHANGE T, U, AND V
!--------------------- 
        CALL EXCH(T , LM, U, LM, V, LM, 4, 4)  
!
        CALL EXCH(Q2, LM, 5, 5)
!    
        CALL HZADV
!-----------------  
! EXCHANGE U AND V
!-----------------    
        CALL EXCH(U, LM, V, LM, CWM, LM, 2, 2) 
!    
        CALL HZADV2
    END IF
!-------------------------------------- 
! TURBULENT PROCESSES AND PRECIPITATION 
!-------------------------------------- 
    IF (MOD(NTSD-NPHS/2,NPHS) == 0) THEN
!--------------------------- 
! EXCHANGE PD, T, U, V AND Q
!--------------------------- 
        CALL EXCH(PD, 1, UZ0, 1, VZ0, 1, T, LM, U, LM, V, LM, Q, LM, 1, 1) 
!----------------------- 
! CONTAINS CALLS TO EXCH
!-----------------------     
        CALL TURBL
    END IF
!---------------------------------------- 
! CONDENSATION/EVAPORATION OF CLOUD WATER 
!---------------------------------------- 
    IF (MOD(NTSD-NPHS/2,NPHS) == 0) THEN
        CALL GSCOND
    END IF
!------------------------- 
! CONVECTIVE PRECIPITATION 
!------------------------- 
    IF (MOD(NTSD-NCNVC/2,NCNVC) == 0) THEN
        CALL CUCNVC
    END IF
!------------------------ 
! GRIDSCALE PRECIPITATION
!------------------------ 
    IF (MOD(NTSD-NPHS/2,NPHS) == 0) THEN
        CALL PRECPD
    END IF
!
    IF (NTSD == NTIM) GOTO 400
!
    GOTO 300
!------------------------
! EXIT FROM THE TIME LOOP 
!------------------------
400 CONTINUE
!
    CALL MPI_BARRIER(MPI_COMM_COMP, IRTN)
!---------------------------------
! RETURN BC FILE TO START FORECAST
!---------------------------------
    NTSD = 0
    IF (MYPE == 0) THEN
        REWIND NBC
!    
        READ(NBC)
        READ(NBC) BCHR
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
    END SUBROUTINE NEWFLT
