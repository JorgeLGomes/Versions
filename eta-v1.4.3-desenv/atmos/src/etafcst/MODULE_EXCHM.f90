    MODULE EXCHM
!--------------------------------------------------------------------------------------------------
    CONTAINS
!
    SUBROUTINE EXCH0(ARR1, LL1, IHALO, JHALO)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE EXCH0
!>
!> SUBPROGRAM: EXCH0 - FUNDAMENTAL EXCHANGE ROUTINES
!> PROGRAMMER: BLACK
!> ORG: W/NP2
!> DATE: 97-06-25
!>
!> ABSTRACT:
!> SUBROUTINE EXCH0 IS USED TO EXCHANGE HALOS BETWEEN PROCESSORS
!>
!>  CURRENTLY SUPPORTED INTERFACES
!>
!> 0, 00, 01, 011, 0001111, 1, 11, 111, 1111, 11111, 111111
!>
!> WHERE 0 REFERS TO A 2-D ARRAY AND 1 REFERS TO A 3-D ARRAY
!>
!> PROGRAM HISTORY LOG:
!> 97-05-??  MEYS       - ORIGINATOR
!> 97-06-25  BLACK      - CONVERTED FROM SHMEM TO MPI
!> 98-??-??  TUCCILLO   - REMOVED EXPLICIT EXCHANGES OF CORNERS
!> 98-??-??  TUCCILLO   - REWROTE TO USE NON_BLOCKING MPI ROUTINES
!> 99-??-??  BLACK      - ADDED VARIABLE HALO SIZES
!> 00-03-10  TUCCILLO   - CHANGED TO USE MODULE PROCDURES FOR INCREASED MESSAGE SIZES AND A UNIFORM
!>                        INTERFACE FOR ALL CALLS
!> 01-02-25  TUCCILLO   - SOME PERFORMANCE IMPROVEMENTS
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> LL1  - THE VERTICAL DIMENSION OF ARR
!> IHALO - THE NUMBER OF POINTS IN THE X DIRECTION TO EXCHANGE IN THE HALO
!> JHALO - THE NUMBER OF POINTS IN THE Y DIRECTION TO EXCHANGE IN THE HALO
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> ARR1 - THE ARRAY TO BE EXCHANGED
!>
!> USE MODULES: EXCH_BUF_REAL
!>              F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : EBU
!>              NEWFLT
!>
!> CALLS      : MPI_IRECV
!>              MPI_ISEND
!>              MPI_WAIT
!>--------------------------------------------------------------------------------------------------
    USE EXCH_BUF_REAL
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
    SAVE
!
    INCLUDE "mpif.h"
!---------------------------------------
! LOCAL VARIABLE CAUSED BY IMPLICIT NONE
!---------------------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LL1     , IHALO   , JHALO
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ITYPE   , IRECV   , IBEG    , IEND    , IC      , J       , I       , ISEND   , IERR
!------------------------
! STANDARD LOCAL VARIABLE 
!------------------------
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!
    INTEGER(KIND=I4KIND), DIMENSION(4)                                                          ::&
    & IHANDLE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(INOUT)       ::&
    & ARR1
!
    ITYPE = MPI_REAL
!------------ 
! NORTH/SOUTH
!------------ 
!
!------------------- 
! RECEIVE FROM NORTH
!------------------- 

    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(1),                                          &
    &                                         MY_NEB(1), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------- 
! RECEIVE FROM SOUTH
!-------------------
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(3),                                          &
    &                                         MY_NEB(3), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!--------------  
! SEND TO NORTH
!-------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR1(I,MYJE-J)
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(1), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!-------------- 
! SEND TO SOUTH
!-------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR1(I,MYJS+J)
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(3), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!------------------------- 
! STORE RESULTS FROM SOUTH
!------------------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,MYJS-J-1) = BUF1(IC)
            END DO
        END DO
    END IF
!----------------- 
! STORE FROM NORTH
!----------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1),ISTAT,IERR)
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,MYJE+J+1) = BUF0(IC)
            END DO
        END DO
    END IF
!
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!--------- 
!EAST/WEST
!--------- 
!
!------------------ 
! RECEIVE FROM WEST
!------------------ 

    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(4),                                          &
    &                                         MY_NEB(4), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------ 
! RECEIVE FROM EAST
!------------------ 
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(2),                                          &
    &                                         MY_NEB(2), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!------------- 
! SEND TO EAST
!------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE - IHALO + 1
        IEND = MYIE
        IC   = 0
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR1(I,J)
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(2), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!------------- 
! SEND TO WEST
!------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS
        IEND = MYIS + IHALO - 1
        IC   = 0
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR1(I,J)
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(4), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!---------------- 
! STORE FROM WEST
!---------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIS - 1
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,J) = BUF0(IC)
            END DO
        END DO
    END IF
!---------------- 
! STORE FROM EAST
!---------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE+1
        IEND = MYIE+IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,J) = BUF1(IC)
            END DO
        END DO
    END IF
!
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    END SUBROUTINE
!
!
    SUBROUTINE EXCH1(ARR1, LL1, IHALO, JHALO)
!--------------------------------------------------------------------------------------------------
    USE EXCH_BUF_REAL
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!---------------------------------------
! LOCAL VARIABLE CAUSED BY IMPLICIT NONE
!---------------------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LL1     , IHALO   , JHALO
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ITYPE   , IRECV   , IBEG    , IEND    , IC      , J       , I       , K       , ISEND   ,   &
    & IERR
!------------------------
! STANDARD LOCAL VARIABLE 
!------------------------
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!
    INTEGER(KIND=I4KIND), DIMENSION(4)                                                          ::&
    & IHANDLE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR1
!
    ITYPE = MPI_REAL
!------------ 
! NORTH/SOUTH
!------------ 
!
!------------------- 
! RECEIVE FROM NORTH
!------------------- 

    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(1),                                          &
    &                                         MY_NEB(1), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------- 
! RECEIVE FROM SOUTH
!------------------- 
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(3),                                          &
    &                                         MY_NEB(3), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!--------------
! SEND TO NORTH
!-------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR1(I,MYJE-J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(1), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!-------------- 
! SEND TO SOUTH
!-------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR1(I,MYJS+J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(3), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!------------------------- 
! STORE RESULTS FROM SOUTH
!------------------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!----------------- 
! STORE FROM NORTH
!----------------- 

    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!---------- 
! EAST/WEST
!---------- 
!
!------------------ 
! RECEIVE FROM WEST
!------------------ 

    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(4),                                          &
    &                                         MY_NEB(4), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------ 
! RECEIVE FROM EAST
!------------------ 
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(2),                                          &
    &                                         MY_NEB(2), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!------------- 
! SEND TO EAST
!------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE - IHALO + 1
        IEND = MYIE
        IC   = 0
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR1(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(2), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!------------- 
! SEND TO WEST
!------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS
        IEND = MYIS + IHALO - 1
        IC   = 0
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR1(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(4), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!---------------- 
! STORE FROM WEST
!---------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIS - 1
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!---------------- 
! STORE FROM EAST
!---------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE + 1
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    END SUBROUTINE
!
!
    SUBROUTINE EXCH01(ARR1, LL1, ARR2, LL2, IHALO, JHALO)
!--------------------------------------------------------------------------------------------------
    USE EXCH_BUF_REAL
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!---------------------------------------
! LOCAL VARIABLE CAUSED BY IMPLICIT NONE
!---------------------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LL1     , LL2     , IHALO   , JHALO
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ITYPE   , IRECV   , IBEG    , IEND    , IC      , J       , I       , K       , ISEND   ,   &
    & IERR
!------------------------
! STANDARD LOCAL VARIABLE 
!------------------------
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!
    INTEGER(KIND=I4KIND), DIMENSION(4)                                                          ::&
    & IHANDLE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(INOUT)       ::&
    & ARR1
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR2
!
    ITYPE = MPI_REAL
!------------ 
! NORTH/SOUTH
!------------ 
!
!------------------- 
! RECEIVE FROM NORTH
!------------------- 
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(1),                                          &
    &                                         MY_NEB(1), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------- 
! RECEIVE FROM SOUTH
!------------------- 
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(3),                                          &
    &                                         MY_NEB(3), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!-------------- 
! SEND TO NORTH
!-------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR1(I,MYJE-J)
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR2(I,MYJE-J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(1), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!-------------- 
! SEND TO SOUTH
!-------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR1(I,MYJS+J)
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR2(I,MYJS+J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(3), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!------------------------- 
! STORE RESULTS FROM SOUTH
!------------------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,MYJS-J-1) = BUF1(IC)
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!----------------- 
! STORE FROM NORTH
!----------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,MYJE+J+1) = BUF0(IC)
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!----------  
! EAST/WEST
!---------- 
!
!------------------ 
! RECEIVE FROM WEST
!------------------ 
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(4),                                          &
    &                                         MY_NEB(4), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------ 
! RECEIVE FROM EAST
!------------------  
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(2),                                          &
    &                                         MY_NEB(2), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!------------- 
! SEND TO EAST
!------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE - IHALO + 1
        IEND = MYIE
        IC   = 0
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR1(I,J)
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR2(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(2), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!------------- 
! SEND TO WEST
!------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS
        IEND = MYIS + IHALO -1 
        IC   = 0
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR1(I,J)
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR2(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(4), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!---------------- 
! STORE FROM WEST
!---------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS-IHALO
        IEND = MYIS-1
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1),ISTAT,IERR)
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,J) = BUF0(IC)
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!---------------- 
! STORE FROM EAST
!---------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE + 1
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,J) = BUF1(IC)
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    END SUBROUTINE
!
!
    SUBROUTINE EXCH00(ARR1, LL1, ARR2, LL2, IHALO, JHALO)
!--------------------------------------------------------------------------------------------------
    USE EXCH_BUF_REAL
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
! 
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!---------------------------------------
! LOCAL VARIABLE CAUSED BY IMPLICIT NONE
!---------------------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LL1     , LL2     , IHALO   , JHALO
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ITYPE   , IRECV   , IBEG    , IEND    , IC      , J       , I       , K       , ISEND   ,   &
    & IERR
!------------------------
! STANDARD LOCAL VARIABLE 
!------------------------
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!
    INTEGER(KIND=I4KIND), DIMENSION(4)                                                          ::&
    & IHANDLE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(INOUT)       ::&
    & ARR1
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(INOUT)       ::&
    & ARR2
!
    ITYPE = MPI_REAL
!------------  
! NORTH/SOUTH
!------------ 
!
!------------------- 
! RECEIVE FROM NORTH
!------------------- 
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(1),                                          &
    &                                         MY_NEB(1), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------- 
! RECEIVE FROM SOUTH
!------------------- 
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(3),                                          &
    &                                         MY_NEB(3), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!-------------- 
! SEND TO NORTH
!-------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR1(I,MYJE-J)
            END DO
        END DO
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR2(I,MYJE-J)
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(1), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!-------------- 
! SEND TO SOUTH
!-------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR1(I,MYJS+J)
            END DO
        END DO
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR2(I,MYJS+J)
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(3), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!------------------------- 
! STORE RESULTS FROM SOUTH
!------------------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,MYJS-J-1) = BUF1(IC)
            END DO
        END DO
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR2(I,MYJS-J-1) = BUF1(IC)
            END DO
        END DO
    END IF

!-----------------
! STORE FROM NORTH
!----------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,MYJE+J+1) = BUF0(IC)
            END DO
        END DO
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR2(I,MYJE+J+1) = BUF0(IC)
            END DO
        END DO
    END IF
!
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!----------  
! EAST/WEST
!---------- 
!
!------------------ 
! RECEIVE FROM WEST
!------------------ 
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(4),                                          &
    &                                         MY_NEB(4), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------ 
! RECEIVE FROM EAST
!------------------ 
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(2),                                          &
    &                                         MY_NEB(2), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!------------- 
! SEND TO EAST
!------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE - IHALO + 1
        IEND = MYIE
        IC   = 0
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR1(I,J)
            END DO
        END DO
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR2(I,J)
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(2), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!------------- 
! SEND TO WEST
!------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS
        IEND = MYIS + IHALO - 1
        IC   =  0
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR1(I,J)
            END DO
        END DO
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR2(I,J)
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(4), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!---------------- 
! STORE FROM WEST
!---------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIS-1
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,J) = BUF0(IC)
            END DO
        END DO
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR2(I,J) = BUF0(IC)
            END DO
        END DO
    END IF
!---------------- 
! STORE FROM EAST
!---------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE + 1
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,J) = BUF1(IC)
            END DO
        END DO
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR2(I,J) = BUF1(IC)
            END DO
        END DO
    END IF
!
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    END SUBROUTINE
!
!
    SUBROUTINE EXCH11(ARR1, LL1, ARR2, LL2, IHALO, JHALO)
!--------------------------------------------------------------------------------------------------
    USE EXCH_BUF_REAL
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!---------------------------------------
! LOCAL VARIABLE CAUSED BY IMPLICIT NONE
!---------------------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LL1     , LL2     , IHALO   , JHALO
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ITYPE   , IRECV   , IBEG    , IEND    , IC      , J       , I       , K       , ISEND   ,   &
    & IERR
!------------------------
! STANDARD LOCAL VARIABLE 
!------------------------
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!
    INTEGER(KIND=I4KIND), DIMENSION(4)                                                          ::&
    & IHANDLE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR1
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR2
!
    ITYPE = MPI_REAL
!------------ 
! NORTH/SOUTH
!------------
!
!-------------------
! RECEIVE FROM NORTH
!------------------- 
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(1),                                          &
    &                                         MY_NEB(1), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------- 
! RECEIVE FROM SOUTH
!------------------- 
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(3),                                          &
    &                                         MY_NEB(3), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!-------------- 
! SEND TO NORTH
!-------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR1(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR2(I,MYJE-J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(1), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!-------------- 
! SEND TO SOUTH
!-------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR1(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR2(I,MYJS+J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(3), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!------------------------- 
! STORE RESULTS FROM SOUTH
!------------------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!----------------- 
! STORE FROM NORTH
!----------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF

    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!---------- 
! EAST/WEST
!----------
!
!------------------ 
! RECEIVE FROM WEST
!------------------ 
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(4),                                          &
    &                                         MY_NEB(4), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!--------------------------------------------------------------------
! RECEIVE FROM EAST
!--------------------------------------------------------------------
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(2),                                          &
    &                                         MY_NEB(2), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!------------- 
! SEND TO EAST
!------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE - IHALO + 1
        IEND = MYIE
        IC   = 0
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR1(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR2(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(2), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!------------- 
! SEND TO WEST
!------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS
        IEND = MYIS + IHALO - 1
        IC   = 0
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR1(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR2(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(4), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!---------------- 
! STORE FROM WEST
!---------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIS - 1
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!----------------
! STORE FROM EAST
!----------------
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE + 1
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    END SUBROUTINE
!
!
    SUBROUTINE EXCH111(ARR1, LL1, ARR2, LL2, ARR3, LL3, IHALO, JHALO)
!--------------------------------------------------------------------------------------------------
    USE EXCH_BUF_REAL
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!---------------------------------------
! LOCAL VARIABLE CAUSED BY IMPLICIT NONE
!---------------------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LL1     , LL2     , LL3     , IHALO   , JHALO
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ITYPE   , IRECV   , IBEG    , IEND    , IC      , J       , I       , K       , ISEND   ,   &
    & IERR
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NEBPE
!------------------------
! STANDARD LOCAL VARIABLE 
!------------------------
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!
    INTEGER(KIND=I4KIND), DIMENSION(4)                                                          ::&
    & IHANDLE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR1
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR2
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR3
!
    ITYPE = MPI_REAL
!------------ ARR1
! NORTH/SOUTH
!------------ 
!
!------------------- 
! RECEIVE FROM NORTH
!------------------- 
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(1),                                          &
    &                                         MY_NEB(1), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------- 
! RECEIVE FROM SOUTH
!------------------- 
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(3),                                          &
    &                                         MY_NEB(3), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!-------------- 
! SEND TO NORTH
!-------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR1(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR2(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR3(I,MYJE-J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(1), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!-------------- 
! SEND TO SOUTH
!-------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR1(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR2(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR3(I,MYJS+J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(3), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!------------------------- 
! STORE RESULTS FROM SOUTH
!------------------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,MYJS-J-1,K) = BUF1(IC)
                END DO 
            END DO
        END DO
    END IF
!----------------- 
! STORE FROM NORTH
!----------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF

    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!---------- 
! EAST/WEST
!---------- 
!
!------------------ 
! RECEIVE FROM WEST
!------------------ 
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(4),                                          &
    &                                         MY_NEB(4), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------ 
! RECEIVE FROM EAST
!------------------ 
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(2),                                          &
    &                                         MY_NEB(2), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!------------- 
! SEND TO EAST
!------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE - IHALO + 1
        IEND = MYIE
        IC   = 0
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR1(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR2(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR3(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(2), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!------------- 
! SEND TO WEST
!------------- 
    IF (MY_NEB(4) >= 0) THEN
        NEBPE = MY_NEB(4)
        IBEG  = MYIS
        IEND  = MYIS + IHALO - 1
        IC    = 0
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR1(I,J,K)
                END DO 
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR2(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR3(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(4), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!---------------- 
! STORE FROM WEST
!---------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIS-1
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!---------------- 
! STORE FROM EAST
!---------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE + 1
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    END SUBROUTINE
!
!
    SUBROUTINE EXCH1111(ARR1, LL1, ARR2, LL2, ARR3, LL3, ARR4, LL4, IHALO, JHALO)
!--------------------------------------------------------------------------------------------------
    USE EXCH_BUF_REAL
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
! 
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!---------------------------------------
! LOCAL VARIABLE CAUSED BY IMPLICIT NONE
!---------------------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LL1     , LL2     , LL3     , LL4     , IHALO   , JHALO
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ITYPE   , IRECV   , IBEG    , IEND    , IC      , J       , I       , K       , ISEND   ,   &
    & IERR
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NEBPE
!------------------------
! STANDARD LOCAL VARIABLE 
!------------------------
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!
    INTEGER(KIND=I4KIND), DIMENSION(4)                                                          ::&
    & IHANDLE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR1
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR2
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR3
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR4
!
    ITYPE = MPI_REAL
!------------ 
! NORTH/SOUTH
!------------ 
!
!------------------- 
! RECEIVE FROM NORTH
!------------------- 
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(1),                                          &
    &                                         MY_NEB(1), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------- 
! RECEIVE FROM SOUTH
!------------------- 
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(3),                                          &
    &                                         MY_NEB(3), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!-------------- 
! SEND TO NORTH
!-------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR1(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR2(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR3(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR4(I,MYJE-J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(1), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!-------------- 
! SEND TO SOUTH
!-------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR1(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR2(I,MYJS+J,K)
                END DO 
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR3(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR4(I,MYJS+J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(3), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!-------------------------
! STORE RESULTS FROM SOUTH
!-------------------------
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!----------------- 
! STORE FROM NORTH
!----------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!---------- 
! EAST/WEST
!---------- 
!
!------------------ 
! RECEIVE FROM WEST
!------------------ 
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(4),                                          &
    &                                         MY_NEB(4), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------ 
! RECEIVE FROM EAST
!------------------ 
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(2),                                          &
    &                                         MY_NEB(2), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!------------- 
! SEND TO EAST
!------------- 
    IF (MY_NEB(2) >= 0) THEN
        NEBPE = MY_NEB(2)
        IBEG  = MYIE - IHALO + 1
        IEND  = MYIE
        IC    = 0
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR1(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR2(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR3(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR4(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(2), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!------------- 
! SEND TO WEST
!------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS
        IEND = MYIS + IHALO - 1
        IC   = 0
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR1(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR2(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR3(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR4(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(4), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!---------------- 
! STORE FROM WEST
!---------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIS - 1
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!---------------- 
! STORE FROM EAST
!---------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE + 1
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    END SUBROUTINE
!
!
    SUBROUTINE EXCH11111(ARR1, LL1, ARR2, LL2, ARR3, LL3, ARR4, LL4, ARR5, LL5, IHALO, JHALO)
!--------------------------------------------------------------------------------------------------
    USE EXCH_BUF_REAL
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!---------------------------------------
! LOCAL VARIABLE CAUSED BY IMPLICIT NONE
!---------------------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LL1     , LL2     , LL3     , LL4     , LL5     , IHALO   , JHALO
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ITYPE   , IRECV   , IBEG    , IEND    , IC      , J       , I       , K       , ISEND   ,   &
    & IERR
!------------------------
! STANDARD LOCAL VARIABLE 
!------------------------
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!
    INTEGER(KIND=I4KIND), DIMENSION(4)                                                          ::&
    & IHANDLE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR1
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR2
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR3
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR4
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR5
!
    ITYPE = MPI_REAL
!------------  
! NORTH/SOUTH
!------------
!
!------------------- 
! RECEIVE FROM NORTH
!------------------- 
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(1),                                          &
    &                                         MY_NEB(1), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!-------------------
! RECEIVE FROM SOUTH
!------------------- 
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(3),                                          &
    &                                         MY_NEB(3), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!-------------- 
! SEND TO NORTH
!-------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR1(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR2(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR3(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR4(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR5(I,MYJE-J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(1), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!-------------- 
! SEND TO SOUTH
!-------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR1(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR2(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR3(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR4(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR5(I,MYJS+J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(3), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!------------------------- 
! STORE RESULTS FROM SOUTH
!------------------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2),ISTAT,IERR)
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR5(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!----------------- 
! STORE FROM NORTH
!----------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR5(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!---------- 
! EAST/WEST
!---------- 
!
!------------------ 
! RECEIVE FROM WEST
!------------------ 
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(4),                                          &
    &                                         MY_NEB(4), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------ 
! RECEIVE FROM EAST
!------------------ 
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(2),                                          &
    &                                         MY_NEB(2), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!------------- 
! SEND TO EAST
!------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE - IHALO + 1
        IEND = MYIE
        IC   = 0
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR1(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR2(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR3(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR4(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR5(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(2), MYPE,MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!------------- 
! SEND TO WEST
!------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS
        IEND = MYIS + IHALO - 1
        IC   = 0
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR1(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR2(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR3(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR4(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR5(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(4), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!---------------- 
! STORE FROM WEST
!---------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIS - 1
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR5(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!---------------- 
! STORE FROM EAST
!---------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE + 1
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR5(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    END SUBROUTINE
!
!
    SUBROUTINE EXCH111111(ARR1, LL1, ARR2, LL2, ARR3, LL3,                                        &
    &                     ARR4, LL4, ARR5, LL5, ARR6, LL6, IHALO, JHALO)
!--------------------------------------------------------------------------------------------------
    USE EXCH_BUF_REAL
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!---------------------------------------
! LOCAL VARIABLE CAUSED BY IMPLICIT NONE
!---------------------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LL1     , LL2     , LL3     , LL4     , LL5     , LL6     , IHALO   , JHALO
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ITYPE   , IRECV   , IBEG    , IEND    , IC      , J       , I       , K       , ISEND   ,   &
    & IERR
!------------------------
! STANDARD LOCAL VARIABLE 
!------------------------
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!
    INTEGER(KIND=I4KIND), DIMENSION(4)                                                          ::&
    & IHANDLE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR1
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR2
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR3
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR4
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR5
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR6
!
    ITYPE = MPI_REAL
!------------
! NORTH/SOUTH
!------------
!
!------------------- 
! RECEIVE FROM NORTH
!------------------- 
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(1),                                          &
    &                                         MY_NEB(1), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------- 
! RECEIVE FROM SOUTH
!------------------- 
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(3),                                          &
    &                                         MY_NEB(3), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!-------------- 
! SEND TO NORTH
!-------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR1(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR2(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR3(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR4(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR5(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR6(I,MYJE-J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(1), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!-------------- 
! SEND TO SOUTH
!-------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR1(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR2(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR3(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR4(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR5(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR6(I,MYJS+J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(3), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!------------------------- 
! STORE RESULTS FROM SOUTH
!------------------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR5(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR6(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!----------------- 
! STORE FROM NORTH
!----------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO K=1,LL1
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR5(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR6(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!---------- 
! EAST/WEST
!---------- 
!
!------------------ 
! RECEIVE FROM WEST
!------------------ 
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(4),                                          &
    &                                         MY_NEB(4), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------ 
! RECEIVE FROM EAST
!------------------ 
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(2),                                          &
    &                                         MY_NEB(2), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!------------- 
! SEND TO EAST
!------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE - IHALO + 1
        IEND = MYIE
        IC   = 0
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR1(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR2(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR3(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR4(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR5(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR6(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(2), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!------------- 
! SEND TO WEST
!------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS
        IEND = MYIS + IHALO - 1 
        IC   = 0
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR1(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR2(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR3(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR4(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR5(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR6(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(4), MYPE ,MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!---------------- 
! STORE FROM WEST
!---------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIS - 1
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR5(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR6(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!---------------- 
! STORE FROM EAST
!---------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE + 1
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO K=1,LL1
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR1(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR5(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR6(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    END SUBROUTINE
!
!
    SUBROUTINE EXCH011(ARR1, LL1, ARR2, LL2, ARR3, LL3, IHALO, JHALO)
!--------------------------------------------------------------------------------------------------
    USE EXCH_BUF_REAL
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!---------------------------------------
! LOCAL VARIABLE CAUSED BY IMPLICIT NONE
!---------------------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LL1     , LL2     , LL3     , IHALO   , JHALO
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ITYPE   , IRECV   , IBEG    , IEND    , IC      , J       , I       , K       , ISEND   ,   &
    & IERR
!------------------------
! STANDARD LOCAL VARIABLE 
!------------------------
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!
    INTEGER(KIND=I4KIND), DIMENSION(4)                                                          ::&
    & IHANDLE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(INOUT)       ::&
    & ARR1
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR2
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR3
!
    ITYPE = MPI_REAL
!------------ 
! NORTH/SOUTH
!------------ 
!
!------------------- 
! RECEIVE FROM NORTH
!------------------- 
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(1),                                          &
    &                                         MY_NEB(1), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------- 
! RECEIVE FROM SOUTH
!------------------- 
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(3),                                          &
    &                                         MY_NEB(3), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!-------------- 
! SEND TO NORTH
!-------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR1(I,MYJE-J)
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR2(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR3(I,MYJE-J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(1), MYPE,MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!-------------- 
! SEND TO SOUTH
!-------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR1(I,MYJS+J)
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR2(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR3(I,MYJS+J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(3), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!------------------------- 
! STORE RESULTS FROM SOUTH
!------------------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,MYJS-J-1) = BUF1(IC)
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!----------------- 
! STORE FROM NORTH
!----------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,MYJE+J+1) = BUF0(IC)
            END DO
        END DO
        DO K=1,LL2
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!---------- 
! EAST/WEST
!---------- 
!
!------------------ 
! RECEIVE FROM WEST
!------------------ 
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(4),                                          &
    &                                         MY_NEB(4), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------ 
! RECEIVE FROM EAST
!------------------ 
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(2),                                          &
    &                                         MY_NEB(2), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!-------------
! SEND TO EAST
!-------------
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE - IHALO + 1
        IEND = MYIE
        IC   = 0
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR1(I,J)
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR2(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR3(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(2), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!------------- 
! SEND TO WEST
!------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS
        IEND = MYIS + IHALO - 1
        IC   = 0
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR1(I,J)
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR2(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR3(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(4), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!---------------- 
! STORE FROM WEST
!---------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIS - 1
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,J) = BUF0(IC)
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!---------------- 
! STORE FROM EAST
!---------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE + 1
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2),ISTAT,IERR)
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,J) = BUF1(IC)
            END DO
        END DO
        DO K=1,LL2
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR2(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL3
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR3(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    END SUBROUTINE
!
!
    SUBROUTINE IEXCH(ARR1, LL1, IHALO, JHALO)
!--------------------------------------------------------------------------------------------------
    USE EXCH_BUF_INTEGER
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!---------------------------------------
! LOCAL VARIABLE CAUSED BY IMPLICIT NONE
!---------------------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LL1     , IHALO   , JHALO
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ITYPE   , IRECV   , IBEG    , IEND    , IC      , J       , I       , ISEND   , IERR
!------------------------
! STANDARD LOCAL VARIABLE 
!------------------------
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!
    INTEGER(KIND=I4KIND), DIMENSION(4)                                                          ::&
    & IHANDLE
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(INOUT)       ::&
    & ARR1
!
    ITYPE = MPI_INTEGER
!------------ 
! NORTH/SOUTH
!------------ 
!
!------------------- 
! RECEIVE FROM NORTH
!------------------- 
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(1),                                          &
    &                                         MY_NEB(1), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------- 
! RECEIVE FROM SOUTH
!-------------------
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(3),                                          &
    &                                         MY_NEB(3), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!-------------- 
! SEND TO NORTH
!--------------
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC)=ARR1(I,MYJE-J)
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(1), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!-------------- 
! SEND TO SOUTH
!-------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC)=ARR1(I,MYJS+J)
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(3), MYPE,MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!------------------------- 
! STORE RESULTS FROM SOUTH
!------------------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,MYJS-J-1) = BUF1(IC)
            END DO
        END DO
    END IF
!----------------- 
! STORE FROM NORTH
!----------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,MYJE+J+1) = BUF0(IC)
            END DO
        END DO
    END IF
!
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!---------- 
! EAST/WEST
!----------
!
!------------------ 
! RECEIVE FROM WEST
!------------------
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(4),                                          &
    &                                         MY_NEB(4), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------ 
! RECEIVE FROM EAST
!------------------ 
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(2),                                          &
    &                                         MY_NEB(2), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!------------- 
! SEND TO EAST
!------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE - IHALO + 1
        IEND = MYIE
        IC   = 0
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR1(I,J)
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(2), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!------------- 
! SEND TO WEST
!-------------
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS
        IEND = MYIS + IHALO - 1
        IC   = 0
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC)=ARR1(I,J)
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(4), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!---------------- 
! STORE FROM WEST
!---------------- 

    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIS - 1
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1),ISTAT,IERR)
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,J) = BUF0(IC)
            END DO
        END DO
    END IF
!----------------  
! STORE FROM EAST
!---------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE + 1
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,J) = BUF1(IC)
            END DO
        END DO
    END IF
!
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    END SUBROUTINE
!
!
    SUBROUTINE EXCH0001111(ARR1, LL1, ARR2, LL2, ARR3, LL3,                                       &
    &                      ARR4, LL4, ARR5, LL5, ARR6, LL6, ARR7, LL7, IHALO, JHALO)
!--------------------------------------------------------------------------------------------------
    USE EXCH_BUF_REAL
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!---------------------------------------
! LOCAL VARIABLE CAUSED BY IMPLICIT NONE
!---------------------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LL1     , LL2     , LL3     , LL4     , LL5     , LL6     , LL7
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & IHALO   , JHALO
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ITYPE   , IRECV   , IBEG    , IEND    , IC      , J       ,  I      , K       , ISEND   ,   &
    & IERR
!------------------------
! STANDARD LOCAL VARIABLE 
!------------------------
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!
    INTEGER(KIND=I4KIND), DIMENSION(4)                                                          ::&
    & IHANDLE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(INOUT)       ::&
    & ARR1
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(INOUT)       ::&
    & ARR2
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(INOUT)       ::&
    & ARR3
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR4
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR5
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR6
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, *)          , INTENT(INOUT)       ::&
    & ARR7
!
    ITYPE = MPI_REAL
!------------ 
! NORTH/SOUTH
!------------ 
!
!------------------- 
! RECEIVE FROM NORTH
!------------------- 

    IF(MY_NEB(1) >= 0)THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(1),                                          &
    &                                         MY_NEB(1), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------ 
! RECEIVE FROM SOUTH
!------------------- 
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(3),                                          &
    &                                         MY_NEB(3), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!-------------- 
! SEND TO NORTH
!-------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR1(I,MYJE-J)
            END DO
        END DO
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR2(I,MYJE-J)
            END DO
        END DO
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR3(I,MYJE-J)
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR4(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR5(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR6(I,MYJE-J,K)
                END DO
            END DO
        END DO
        DO K=1,LL7
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR7(I,MYJE-J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(1), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!-------------- 
! SEND TO SOUTH
!-------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR1(I,MYJS+J)
            END DO
        END DO
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR2(I,MYJS+J)
            END DO
        END DO
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR3(I,MYJS+J)
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR4(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR5(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR6(I,MYJS+J,K)
                END DO
            END DO
        END DO
        DO K=1,LL7
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR7(I,MYJS+J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(3), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!------------------------- 
! STORE RESULTS FROM SOUTH
!------------------------- 
    IF (MY_NEB(3) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,MYJS-J-1) = BUF1(IC)
            END DO
        END DO
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR2(I,MYJS-J-1) = BUF1(IC)
            END DO
        END DO
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR3(I,MYJS-J-1) = BUF1(IC)
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR5(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR6(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL7
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR7(I,MYJS-J-1,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!----------------- 
! STORE FROM NORTH
!----------------- 
    IF (MY_NEB(1) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,MYJE+J+1) = BUF0(IC)
            END DO
        ENDDO
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR2(I,MYJE+J+1) = BUF0(IC)
            END DO
        END DO
        DO J=0,JHALO-1
            DO I=IBEG,IEND
                IC = IC + 1
                ARR3(I,MYJE+J+1) = BUF0(IC)
            END DO
        END DO
        DO K=1,LL4
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR5(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR6(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL7
            DO J=0,JHALO-1
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR7(I,MYJE+J+1,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(1) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(3) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!---------- 
! EAST/WEST
!---------- 
!
!------------------ 
! RECEIVE FROM WEST
!------------------ 

    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_IRECV(BUF0, IBUFEXCH, ITYPE, MY_NEB(4),                                          &
    &                                         MY_NEB(4), MPI_COMM_COMP, IHANDLE(1), IRECV)
    END IF
!------------------ 
! RECEIVE FROM EAST
!------------------ 
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_IRECV(BUF1, IBUFEXCH, ITYPE, MY_NEB(2),                                          &
    &                                         MY_NEB(2), MPI_COMM_COMP, IHANDLE(2), IRECV)
    END IF
!------------- 
! SEND TO EAST
!------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE - IHALO + 1
        IEND = MYIE
        IC   = 0
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR1(I,J)
            END DO
        END DO
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR2(I,J)
            END DO
        END DO
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF2(IC) = ARR3(I,J)
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR4(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR5(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR6(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL7
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF2(IC) = ARR7(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF2, IC, ITYPE, MY_NEB(2), MYPE, MPI_COMM_COMP, IHANDLE(3), ISEND)
    END IF
!------------- 
! SEND TO WEST
!------------- 

    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS
        IEND = MYIS + IHALO - 1
        IC   = 0
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR1(I,J)
            END DO
        END DO
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR2(I,J)
            END DO
        END DO
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                BUF3(IC) = ARR3(I,J)
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR4(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR5(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR6(I,J,K)
                END DO
            END DO
        END DO
        DO K=1,LL7
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    BUF3(IC) = ARR7(I,J,K)
                END DO
            END DO
        END DO
        CALL MPI_ISEND(BUF3, IC, ITYPE, MY_NEB(4), MYPE, MPI_COMM_COMP, IHANDLE(4), ISEND)
    END IF
!---------------- 
! STORE FROM WEST
!---------------- 
    IF (MY_NEB(4) >= 0) THEN
        IBEG = MYIS - IHALO
        IEND = MYIS - 1
        IC   = 0
        CALL MPI_WAIT(IHANDLE(1), ISTAT, IERR)
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,J) = BUF0(IC)
            END DO
        END DO
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR2(I,J) = BUF0(IC)
            END DO
        END DO
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR3(I,J) = BUF0(IC)
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR5(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR6(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
        DO K=1,LL7
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR7(I,J,K) = BUF0(IC)
                END DO
            END DO
        END DO
    END IF
!---------------- 
! STORE FROM EAST
!---------------- 
    IF (MY_NEB(2) >= 0) THEN
        IBEG = MYIE + 1
        IEND = MYIE + IHALO
        IC   = 0
        CALL MPI_WAIT(IHANDLE(2), ISTAT, IERR)
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR1(I,J) = BUF1(IC)
            END DO
        END DO
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR2(I,J) = BUF1(IC)
            END DO
        END DO
        DO J=MYJS-JHALO,MYJE+JHALO
            DO I=IBEG,IEND
                IC = IC + 1
                ARR3(I,J) = BUF1(IC)
            END DO
        END DO
        DO K=1,LL4
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR4(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL5
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR5(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL6
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR6(I,J,K)=BUF1(IC)
                END DO
            END DO
        END DO
        DO K=1,LL7
            DO J=MYJS-JHALO,MYJE+JHALO
                DO I=IBEG,IEND
                    IC = IC + 1
                    ARR7(I,J,K) = BUF1(IC)
                END DO
            END DO
        END DO
    END IF
!
    IF (MY_NEB(4) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(4), ISTAT, IERR)
    END IF
!
    IF (MY_NEB(2) >= 0) THEN
        CALL MPI_WAIT(IHANDLE(3), ISTAT, IERR)
    END IF
!
    END SUBROUTINE
!
    END MODULE
