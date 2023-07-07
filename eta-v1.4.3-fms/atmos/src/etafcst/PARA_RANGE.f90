    SUBROUTINE PARA_RANGE(N1, N2, NPROCS, IRANK, ISTA, IEND)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE PARA_RANGE
!>
!> SUBROUTINE: PARA_RANGE - SET UP DECOMPOSITION VALUES
!> PROGRAMMER: TUCCILLO
!> ORG: IBM
!> DATE: 00-01-06
!>
!> ABSTRACT:
!> SETS UP DECOMOSITION VALUES
!>
!> PROGRAM HISTORY LOG:
!> 00-01-06  TUCCILLO - ORIGINAL
!> 18-01-15  LUCCI    - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> N1     - FIRST INTERATE VALUE
!> N2     - LAST INTERATE VALUE
!> NPROCS - NUMBER OF MPI TASKS
!> IRANK  - MY TAKS ID
!>
!> OUTPUT ARGUMENT LIST:
!> ISTA    - FIRST LOOP VALUE
!> IEND    - LAST LOOP VALUE
!>
!> OUTPUT FILES:
!> STDOUT  - RUN TIME STANDARD OUT.
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : CHKOUT
!>              SETUP_SERVERS
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & N1      , N2      , NPROCS  , IRANK 
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(INOUT)       ::&
    & ISTA    , IEND
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IWORK1  , IWORK2
!
    IWORK1 =    (N2 - N1 + 1) / NPROCS
    IWORK2 = MOD(N2 - N1 + 1,   NPROCS)
!
    ISTA = IRANK * IWORK1 + N1 + MIN(IRANK, IWORK2)
!
    IEND = ISTA + IWORK1 - 1
!
    IF (IWORK2 > IRANK) IEND = IEND + 1
!
    END SUBROUTINE PARA_RANGE

