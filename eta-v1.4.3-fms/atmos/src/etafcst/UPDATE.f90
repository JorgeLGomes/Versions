    SUBROUTINE UPDATE(A)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE UPDATE
!>
!> SUBROUTINE: UPDATE - EXCHANGE TWO HALO ROWS
!> PROGRAMMER: TUCCILLO
!> ORG: IBM
!> DATE: 00-01-06
!>
!> ABSTRACT:
!> EXCHANGE TWO HALO ROWS
!>
!> PROGRAM HISTORY LOG:
!> 00-01-06  TUCCILLO - ORIGINAL
!> 18-01-15  LUCCI    - MODERNIZATION OF THE CODE, INCLUDING:
!>                      * F77 TO F90/F95
!>                      * INDENTATION & UNIFORMIZATION CODE
!>                      * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                      * DOCUMENTATION WITH DOXYGEN
!>                      * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> A - ARRAY TO HAVE HALOS EXCHANGED
!> NONE
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> A - ARRAY TO HAVE HALOS EXCHANGE
!>
!> OUTPUT FILES:
!> STDOUT  - RUN TIME STANDARD OUT.
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
!> DRIVER     : SLP
!>              SLPSIG
!>              SLPSIGSPLINE
!>
!> CALLS      : MPI_SENDRECV
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
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED)                    , INTENT(INOUT)       ::&
    & A
!
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & STATUS
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IERR
!
    IF (NUM_PROCS == 1) RETURN
!
    CALL MPI_SENDRECV(A(1,JEND_I-1), 2*IM                    , MPI_REAL                         , & 
    &    IUP                       , 1                       , A(1,JSTA_I-2)                    , &
    &    2*IM                      , MPI_REAL                , IDN                              , &
    &    1                         , MPI_COMM_COMP           , STATUS                           , &
    &    IERR)
!
    CALL MPI_SENDRECV(A(1,JSTA_I  ), 2*IM                    , MPI_REAL                         , &
    &    IDN                       , 1                       , A(1,JEND_I+1)                    , &
    &    2*IM                      , MPI_REAL                , IUP                              , &
    &    1                         , MPI_COMM_COMP           , STATUS                           , &
    &    IERR)
!
    END SUBROUTINE UPDATE
