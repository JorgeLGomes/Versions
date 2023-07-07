    SUBROUTINE SETUP_SERVERS(NPES_MOD     , MYPE          , NPES                ,                 &
    &                        IQUILT_GROUP , INUMQ         ,                                       &
    &                        MPI_COMM_COMP, MPI_COMM_INTER, MPI_COMM_INTER_ARRAY)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SETUP_SERVERS
!> 
!> SUBROUTINE: SETUP_SERVERS - SETUP I/O SERVERS
!> PROGRAMMER: TUCCILLO
!> ORG: IBM
!> DATE: 00-03-20
!>
!> ABSTRACT:
!> SETUP I/O SERVERS
!>
!> PROGRAM HISTORY LOG:
!> 00-03-11  TUCCILLO - ORIGINATOR
!> 18-01-15  LUCCI    - MODERNIZATION OF THE CODE, INCLUDING:
!>                      * F77 TO F90/F95
!>                      * INDENTATION & UNIFORMIZATION CODE
!>                      * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                      * DOCUMENTATION WITH DOXYGEN
!>                      * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NPES_MOD - NUMBER OF MPI TASKS FOR MODEL INTEGRATION FROM INPES AND JNPES THIS IS THE NUMBER
!>            OF MPI TASKS THE EXECUTABLE HAS BEEN BUILT FOR.
!>            NPES, RETURNED FROM MPI_COMM_SIZE, MUST BE AT LEAST THIS SIZE OTHERWISE THE
!>            INTEGRATION CANNOT PROCEED. THE DIFFERENCE BETWEEN NPES_MOD AND NPES IS THE NUMBER
!>            OF MPI TASKS THAT ARE AVAILABLE FOR I/O SERVING. THIS CAN BE ZERO, IN WHICH CASE 
!>            CHKOUT WILL WRITE A DIRECT ACCESS FILE THAT CAN BE SEPARTELY "QUILTED".
!>            IN ORDER TO SKIP THE SEPARATE QUILTING STEP, MAKE SURE THAT THE NUMBER OF MPI TASKS 
!>            THAT THE CODE IS INITIATED WITH IS AT LEAST ONE GREATER THAN NPES_MOD.
!>
!> OUTPUT ARGUMENT LIST:
!> MYPE                 - MY RANK
!> IQUILT_GROUP         - NUMBER OF I/O SERVER GROUPS
!> INUMQ                - ARRAY THAT HOLDS THE NUMBER OF SERVERS IN EACH GROUP
!> NPES                 - NUMBER OF MPI TASKS FOR MODEL INTEGRATION
!> MPI_COMM_COMP        - THE NEW INTRACOMMUNICATOR FOR ALL TASKS
!> MPI_COMM_INTER       - THE INTERCOMMUNICATOR FOR THE I/O SERVERS
!> MPI_COMM_INTER_ARRAY - THE ARRAY OF INTERCOMMUNICATORS FOR THE INTEGRATION TASKS
!>
!> INPUT FILES:  NONE
!>
!> OUTPUT FILES:
!> NONE BUT THE CODE DOES ATTEMPT TO READ THE ENVIRONMENT VARIABLE "SERVER_GROUPS".
!> THIS IS THE NUMBER OF INDEPENDENT GROUPS OF SERVER TASKS.
!> THE DEFAULT IS ONE AND SHOULD BE OK FOR MOST APPLICATIONS OF THE ETA CODE.
!> IF ONE SET OF I/O SERVERS CAN NOT COMPLETE BEFORE THE NEXT OUPUT TIME THEN ADDITIONAL I/O SERVER
!> GROUPS WOULD BE USEFUL.
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : EBU
!> 
!> CALLS      : GETENV
!>              MPI_ABORT
!>              MPI_BARRIER
!>              MPI_COMM_CREATE
!>              MPI_COMM_DUP
!>              MPI_COMM_GROUP
!>              MPI_COMM_RANK
!>              MPI_COMM_SIZE
!>              MPI_COMM_SPLIT
!>              MPI_GROUP_CREATE
!>              MPI_GROUP_EXCL
!>              MPI_GROUP_FREE
!>              MPI_INIT   
!>              MPI_INTERCOMM_CREATE    
!>              PARA_RANGE
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    INCLUDE 'mpif.h'
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & COMDUP
!
    INTEGER(KIND=I4KIND), DIMENSION(:), ALLOCATABLE                                             ::&
    & IRANK
!
    INTEGER(KIND=I4KIND), DIMENSION(*)                                    , INTENT(INOUT)       ::& 
    & MPI_COMM_INTER_ARRAY        , INUMQ
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & YES
!
    CHARACTER(4)                                                                                ::&
    & GET
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & NPES_MOD
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(INOUT)       ::&
    & MYPE    , NPES    , IQUILT_GROUP      , MPI_COMM_COMP     , MPI_COMM_INTER     
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IERR    , IQSERVER, I       , ISTAQ   , IENDQ   , ICOLOR  , ISTAXX  , IENDXX  , IRLR    ,   &
    & ICC     , ISS     , JJ      , ISSL    , KK      , IWORLD  , IGROUP  , IGROUP_X,             &
    & IWORLD_MINUS      , IXX
!---------------------------------------------------
! INITIALIZE MPI
! RETRIEVE THE NUMBER OF TOTAL MPI TASKS AND MY RANK
!---------------------------------------------------
    CALL MPI_INIT(IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYPE, IERR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPES, IERR)
!--------------------------------------------------------------------------------------------------
! AT THIS POINT NPES IS THE TOTAL NUMBER OF MPI TASKS. 
! WE WILL RESET THIS AT THE END OF THE SUBROUTINE TO THE NUMBER OF MPI TASKS
! THAT ARE WORKING ON THE MODEL INTGRATION.
!
! FIRST, HOWEVER, WE NEED TO MAKE SURE THAT A SUFFICIENT NUMBER OF MPI TASKS HAVE BEEN INITIATED. 
! IF NOT, WE WILL STOP.
!--------------------------------------------------------------------------------------------------
    IF (NPES < NPES_MOD) THEN
        PRINT *, ' *******************************************************'
        PRINT *, ' *******************************************************'
        PRINT *, ' ******************** MAJOR PROBLEM ********************'
        PRINT *, ' ******************** MAJOR PROBLEM ********************'
        PRINT *, ' ******************** MAJOR PROBLEM ********************'
        PRINT *, ' ******************** MAJOR PROBLEM ********************'
        PRINT *
        PRINT *, '                      STOPPING NOW                      '
        PRINT *, '      THERE ARE INSUFFICIENT MPI TASKS TO CONTINUE      '
        PRINT *, '      YOU MUST SPECIFY AT LEAST ',NPES_MOD,' TASKS      '
        PRINT *, '                      STOPPING NOW                      '
        PRINT *
        PRINT *, ' ******************** MAJOR PROBLEM ********************'
        PRINT *, ' ******************** MAJOR PROBLEM ********************'
        PRINT *, ' ******************** MAJOR PROBLEM ********************'
        PRINT *, ' ******************** MAJOR PROBLEM ********************'
        PRINT *, ' *******************************************************'
        PRINT *, ' *******************************************************'
!
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
!
    END IF
!-----------------------------------------------------------------------------------
! OK, WE HAVE A SUFFICIENT NUMBER OF MPI TASKS TO CONTINUE
!
! HOW MANY GROUPS OF SERVERS ? THE DEFAULT IS 1 GROUP
! THE ENVIRONMENT VARIABLE, SERVER_GROUPS, CAN BE USED TO SPECIFY MORE SERVER GROUPS
!-----------------------------------------------------------------------------------
    GET = '1'
!
    CALL GETENV('SERVER_GROUPS', GET)
    READ(GET,FMT='(I4)') IQUILT_GROUP
!
    IQUILT_GROUP = MAX(IQUILT_GROUP, 1)
!------------------------------------------------------------------ 
! ERROR CHECK NUMBER OF GROUPS - THE MAXIMUM IS 100 - THIS IS A LOT
!------------------------------------------------------------------ 
    IF (IQUILT_GROUP > 100) THEN
        PRINT *, ' ***** IQUILT_GROUP IS GREATER THAN 100'
        PRINT *, ' ***** DO YOU REALLY WANT THIS ?'
        PRINT *, ' ***** IF SO THEN INCREASE SIZE IN MPP.H'
        PRINT *, ' ***** ALSO, CHANGE IF CHECK IN SETUP_SERVERS'
        PRINT *, ' ***** RESETTING THE NUMBER OF SERVER GROUPS TO 100'
        PRINT *, ' ***** WE ARE CONTINUING ....   '
        IQUILT_GROUP = 100
    END IF
!
    IF (MYPE == 0) THEN
        PRINT *, ' WE WILL TRY TO RUN WITH ',IQUILT_GROUP,' SERVER GROUPS'
    END IF
!--------------------------------------------------------------------------------------------------  
! COMPUTE THE NUMBER OF SERVERS PER GROUP
! ALL MPI TASKS BEYOND NPES_MOD WILL BE SERVERS
! IF THE NUMBER OF SERVERS IS NOT EQUALLY DIVISIBLE BY THE NUMBER OF GROUPS OF SERVERS THEN SOME
! GROUPS MAY HAVE MORE SERVERS THEN 
! OTHERS - THIS IS FINE
! NOTE THAT WE REQUIRE AT LEAST ONE SERVER PER GROUP
! WE MAY NEED TO REDUCE THE NUMBER OF SERVER GROUPS IF IT EXCEEDS THE NUMBER OF SERVERS
!--------------------------------------------------------------------------------------------------
    IQSERVER = NPES - NPES_MOD
!
    IF (IQSERVER == 0) THEN
        IF (MYPE == 0) THEN
            PRINT *, ' YOU SPECIFIED 0 I/O SERVERS '
            PRINT *, ' CHKOUT WILL WRITE A FILE'
        END IF
!
        IQUILT_GROUP = 0
    END IF
!
    IF (IQUILT_GROUP > IQSERVER)  THEN
        IQUILT_GROUP = IQSERVER
        PRINT *, ' NOT ENOUGH SERVERS'
        PRINT *, ' WE NEED TO REDUCE THE NUMB OF SERVER GROUPS'
        PRINT *, ' NUMB OF SERVER GROUPS IS ', IQUILT_GROUP
    END IF
!
    DO I = 0, IQUILT_GROUP - 1
        CALL PARA_RANGE(1, IQSERVER, IQUILT_GROUP, I, ISTAQ, IENDQ)
!
        INUMQ(I+1) = IENDQ - ISTAQ + 1
        IF (MYPE == 0) PRINT *, ' I, INUMQ = ',I+1, INUMQ(I+1)
    END DO
!------------------------------------------------------------------------------- 
! SETUP THE "COLOR" FOR MPI_COMM_SPLIT
! THOSE TASKS WHICH WILL DO MODEL INTEGRATION WILL BE COLOR 0
! THE SERVER TASKS WILL HAVE THE COLOR OF THE GROUP NUMBER THAT THEY WILL BELONG
!-------------------------------------------------------------------------------
    IF (MYPE < NPES_MOD) THEN
        ICOLOR = 0
    ELSE
        ISTAXX = NPES_MOD
        DO I=1,IQUILT_GROUP
            IENDXX = ISTAXX + INUMQ(I) - 1
            IF (MYPE >= ISTAXX .AND. MYPE <= IENDXX) THEN
                ICOLOR = I
            END IF
            ISTAXX = IENDXX + 1
        END DO
    END IF
!--------------------------------------------------------------------------------------------------
! SPLIT THE COMMUNICATOR - THE NEW INTRACOMMUNICATOR FOR ALL TASKS IS MPI_COMM_COMP. 
! MPI_COMM_WORLD IS STILL AVAILABLE BUT IT DOES REFER TO ALL THE MPI TASKS ( MODEL INTEGRATION AND
! I/O SERVING )
!--------------------------------------------------------------------------------------------------
    CALL MPI_COMM_DUP(MPI_COMM_WORLD, COMDUP, IERR)
    CALL MPI_COMM_SPLIT(COMDUP, ICOLOR, MYPE, MPI_COMM_COMP, IERR)
!--------------------------------------------------------------------------------------------------
! AT THIS POINT WE HAVE A NEW COMMUNICATOR, MPI_COMM_COMP, THAT CAN BE USED BY THE FORECASTS TASKS 
! AND THE I/O SERVER TASKS FOR THEIR INTERNAL COMMUNICATIONS. ONTO THE INTERCOMMUNICATORS ...
!
! NOW WE MUST CREATE THE INTERCOMMUNICATORS FOR USE BETWEEN THE MPI TASKS DOING THE MODEL
! INTEGRATION AND THE MPI TASKS FOR EACH SERVER GROUP. THE FIRST STEP IS TO EXCLUDE THE TASKS THAT
! DONT BELONG. WE WILL DO THIS FOR EACH SERVER GROUP BY EXCLUDING THE TASKS FROM ALL OF THE OTHER
! SERVER GROUPS.
!--------------------------------------------------------------------------------------------------
    ALLOCATE (IRANK(IQSERVER))
    IXX = NPES_MOD
    DO I=1,IQUILT_GROUP
        YES = .TRUE. 
        IF (MYPE < NPES_MOD) THEN
            IRLR = IXX
        ELSE
            IRLR = 0
        END IF
!
        ICC = 0
        ISS = NPES_MOD
!---------------------------------------------------------- 
! THIS IS THE FIRST POSSIBLE TASK ID THAT COULD BE EXCLUDED
!----------------------------------------------------------
        DO JJ=1,IQUILT_GROUP
            IF (JJ /= I) THEN
                ISSL = ISS
                DO KK=1,INUMQ(JJ)
                    ICC        = ICC + 1
                    IRANK(ICC) = ISSL
                    IF (MYPE == ISSL) YES = .FALSE. 
                    ISSL = ISSL + 1 
                END DO
            END IF
            ISS = ISS + INUMQ(JJ)
        END DO
!--------------------------------------------------------------------------------------------------    
! AT THIS POINT WE HAVE AN ARRAY, IRANK, WITH TASK IDS TO EXCLUDE THERE ARE ICC OF THEM.
! CREATE A NEW GROUP WITH THE TASKS FROM THE OTHER SERVER GROUPS EXCLUDED AND THEN CREATE A NEW 
! COMMUNICATOR ( IWORLD_MINUS ) THAT CONTAINS ONLY THE MPI TASKS DOING THE MODEL INTEGRATION AND 
! THE TASKS THAT BLONG TO THE SERVER GROUP WE ARE CONSIDERING.
!--------------------------------------------------------------------------------------------------     
        IWORLD = MPI_COMM_WORLD
!
        CALL MPI_COMM_GROUP (IWORLD, IGROUP, IERR)
        CALL MPI_GROUP_EXCL (IGROUP, ICC     ,IRANK        , IGROUP_X, IERR)
        CALL MPI_COMM_CREATE(IWORLD, IGROUP_X, IWORLD_MINUS, IERR)
        CALL MPI_GROUP_FREE (IGROUP  , IERR)
        CALL MPI_GROUP_FREE (IGROUP_X, IERR)
!--------------------------------------------------------------------------------------------------     
! AT THIS POINT WE HAVE A COMMUNICATOR THAT EXCLUDES THE TASKS WE DONT WANT.
! CREATE AN INTERCOMMUNICATOR FOR USE BETWEEN THE MPI TASKS DOING THE MODEL INTEGRATION AND THE I/O
! SERVER GROUP WE ARE CONSIDERING. THIS PROCESS IS A COLLECTIVE ROUTINE SO IT CAN ONLY BE DONE BY
! THE TASKS THAT HAVE NOT BEEN EXCLUDED. SAVE THIS NEW COMMUNICATOR IN MPI_COMM_INTER FOR USE BY 
! THE TASKS THAT BELONG TO THE SERVER GROUP THAT WE ARE CONSIDERING. 
! THE TASKS THAT ARE PERFORMING THE MODEL INTEGRATION WILL REFERENCE MPI_COMM_INTER_ARRAY() SINCE 
! WE WILL NEED TO SELECT WHICH 
! SERVER GROUP WE WISH TO COMMUNICATE WITH.
!--------------------------------------------------------------------------------------------------   
        IF (YES) THEN
            CALL MPI_INTERCOMM_CREATE                                                             &
    &                     (MPI_COMM_COMP, 0, IWORLD_MINUS, IRLR, 0, MPI_COMM_INTER_ARRAY(I), IERR)
!
            MPI_COMM_INTER = MPI_COMM_INTER_ARRAY(I)
        END IF
!    
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
!------------------------------------------------- 
! END DO FOR LOOP OVER THE NUMBER OF SERVER GROUPS
!-------------------------------------------------
    END DO 
!--------------------------------------------------------------------
! NPES IS REALLY THE NUMBER OF TASKS WORKING ON THE MODEL INTEGRATION
!--------------------------------------------------------------------
    NPES = NPES - IQSERVER
!
    IF (MYPE == 0) THEN
        PRINT *, ' THE MODEL INTEGRATION IS USING ',NPES,' MPI TASK'
        PRINT *, ' THERE ARE ',IQSERVER,' I/O SERVERS'
    END IF
!
    DEALLOCATE (IRANK)
!
    END SUBROUTINE SETUP_SERVERS
