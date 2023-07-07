    SUBROUTINE MPI_FIRST
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE MPI_FIRST
!>
!>
!> SUBROUTINE: MPI_FIRST - INTIALIZES MPI STUFF
!> PROGRAMMER: TUCCILLO
!> ORG: IBM
!> DATE: 00-01-20
!>
!> ABSTRACT:
!> INTIALIZES MPI STUFF
!>
!> PROGRAM HISTORY LOG:
!> 00-01-20  TUCCILLO - ORIGINATOR
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
!> EXIT STATES:
!> COND =   0 - NORMAL EXIT
!>
!> USE MODULES: GLB_TABLE
!>              F77KINDS
!>              MAPPINGS
!>              MPPCOM
!>              PARA
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : QUILT
!>
!> CALLS      : MPI_ABORT
!>              MPI_COMM_RANK
!>              MPI_COMM_SIZE           
!>              PARA_RANGE
!>--------------------------------------------------------------------------------------------------
    USE GLB_TABLE
    USE F77KINDS
    USE MAPPINGS
    USE MPPCOM
    USE PARA
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
    INCLUDE 'mpif.h'
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IERR
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & LNIP    , LNJP    ,                                                                         &
    & JS_X    , JE_X    ,                                                                         &
    & IRMND   , JRMND   ,                                                                         &
    & IPOSN   , JPOSN   ,                                                                         &
    & I       ,                                                                                   &
    & IPE
!
    LNIP = IM / INPES
    LNJP = JM / JNPES
!-------------------------------------------------------------------------
! NUM_PROCS IS THE NUMBER OF TASKS DOING THE QUILTING IN THIS SERVER GROUP
!-------------------------------------------------------------------------
    CALL MPI_COMM_SIZE(MPI_COMM_COMP, NUM_PROCS, IERR)
    CALL MPI_COMM_RANK(MPI_COMM_COMP, ME       , IERR)
!
    PRINT *, ' NUM_PROCS = ', NUM_PROCS
!
    IF ( NUM_PROCS > JNPES ) THEN
        PRINT *, ' TOO MANY MPI TASKS, MAX IS ',JNPES, ', STOPPING'
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    END IF
!
    IF ( NUM_PROCS > 2800 ) THEN
        PRINT *, ' TOO MANY MPI TASKS, MAX IS 1024, STOPPING'
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    END IF
!--------------------------------------------------------------------------------------------------
! JS_X AND JS_Y ARE THE STARTING AND ENDING ROWS OF TASKS IN THE MODEL FORECAST DECOMPOSITION THAT 
! WILL BE SENDING TO EACH QUILT TASK
!
! JSTA IS THE FIRST FORECAST TASK AND JEND IS THE LAST FORECAST TASK IN THE ENTIRE RANGE OF 
! FORECAST TASKS THAT WILL BE SENDING TO EACH QUILT TASK. 
! REMEMBER THAT AN INTEGER NUMBER OF FORECAST TASK ROWS IS SENT TO EACH QUILT TASK.
!--------------------------------------------------------------------------------------------------
    DO I = 0, NUM_PROCS-1
        CALL PARA_RANGE(0, JNPES-1, NUM_PROCS, I, JS_X, JE_X)
!
        JSTA(I) = JS_X * INPES
        JEND(I) = JSTA(I) + (JE_X-JS_X+1) * INPES-1
!
        IF (ME == 0) THEN
            PRINT *, ' TASK ID, JSTA, END = ', I, JSTA(I), JEND(I)
        END IF
    END DO
!--------------------------------------------------------------------------------------------------
! LOCATIONS
!
! PARAMETER LNIP (LNJP) IS SMALLEST THAT THE I (J) EXTENT OF EACH SUBDOMAIN CAN BE.
! IRMND (JRMND) IS THE NUMBER OF "REMAINDER" I (J) POINTS THAT WILL BE GIVEN TO THE LEFTMOST 
! (LOWERMOST) PES.
!--------------------------------------------------------------------------------------------------
    IRMND = MOD(IM, INPES)
    JRMND = MOD(JM, JNPES)
!
    DO IPE=0,NPES-1
!    
        IPOSN = MOD(IPE, INPES) + 1
        JPOSN = IPE / INPES + 1
!------------------------------------    
! GLOBAL LIMITS OF THIS PES SUBDOMAIN
!------------------------------------     
        MY_IS_GLB_A(IPE) = (IPOSN-1) * LNIP + MIN(IRMND, IPOSN-1) + 1
        MY_IE_GLB_A(IPE) = MY_IS_GLB_A(IPE) + LNIP - 1
!
        IF (IPOSN <= IRMND) MY_IE_GLB_A(IPE) = MY_IE_GLB_A(IPE) + 1
!    
        MY_JS_GLB_A(IPE) = (JPOSN-1) * LNJP + MIN(JRMND, JPOSN-1) + 1
        MY_JE_GLB_A(IPE) = MY_JS_GLB_A(IPE) + LNJP - 1
!
        IF (JPOSN <= JRMND) MY_JE_GLB_A(IPE) = MY_JE_GLB_A(IPE) + 1
!
    END DO
!------------------------- 
! DIMENSIONING INFORMATION
!------------------------- 
    MY_ISD = 1
    MY_IED = IM
    MY_JSD = MY_JS_GLB_A(JSTA(ME)) - 2
    MY_JED = MY_JE_GLB_A(JEND(ME)) + 2
!
    IF ( MY_JSD < 1  ) MY_JSD = 1
    IF ( MY_JED > JM ) MY_JED = JM
!
    PRINT *, ' ME, MY_ISD,MY_IED,MY_JSD,MY_JED = ', ME, MY_ISD, MY_IED, MY_JSD, MY_JED
!
    JSTA_I = MY_JS_GLB_A(JSTA(ME))
    JEND_I = MY_JE_GLB_A(JEND(ME))
!
    JSTA_IM  = JSTA_I
    JSTA_IM2 = JSTA_I
    JEND_IM  = JEND_I
    JEND_IM2 = JEND_I
!
    IF (ME == 0) THEN
        JSTA_IM  = 2
        JSTA_IM2 = 3
    END IF
!
    IF (ME == NUM_PROCS-1) THEN
        JEND_IM=JM-1
        JEND_IM2=JM-2
    END IF
!
    PRINT *, ' JSTA_I, JEND_I, JSTA_IM, JEND_IM, JSTA_IM2, JEND_IM2= ',                           &
    &          JSTA_I, JEND_I, JSTA_IM, JEND_IM, JSTA_IM2, JEND_IM2
!----------  
! NEIGHBORS
!----------
    IUP = ME + 1
    IDN = ME - 1
!
    IF (ME == 0) THEN
        IDN = MPI_PROC_NULL
    END IF
!
    IF (ME == NUM_PROCS-1) THEN
        IUP = MPI_PROC_NULL
    END IF
!
    END SUBROUTINE MPI_FIRST

