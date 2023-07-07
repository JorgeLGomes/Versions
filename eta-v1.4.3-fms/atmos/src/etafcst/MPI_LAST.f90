    SUBROUTINE MPI_LAST
!>-------------------------------------------------------------------------------------------------- 
!> SUBROUTINE MPI_LAST
!>
!> SUBROUTINE: MPI_LAST - SHUTS DOWN MPI
!> PROGRAMMER: TUCCILLO
!> ORG: IBM
!> DATE: 00-01-20
!>
!> ABSTRACT:
!> SHUTS DOWN MPI
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
!> USE MODULES: -----
!>
!> DRIVER     : -----
!>
!> CALLS      : ----- 
!-------------------------------------------------------------------------------------------------- 
    IMPLICIT NONE
!
    CALL MPI_FINALIZE(IERR)
!
    END SUBROUTINE MPI_LAST
