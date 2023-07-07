    MODULE M_PRINTVEG
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE M_PRINTVEG
!>
!> SUBROUTINE: M_PRINTVEG - ?????
!> PROGRAMMER: ?????
!> ORG: ?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> MODULE RESPONSIBLE FOR THE WRITING OF MATRICES IN FILES IN AN ENVIRONMENT WITH SUPPORT TO MPI. 
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
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : ----- 
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    PRIVATE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & MPI_MYRANK
!
    PUBLIC INICIA, WRITEFILE2D
!
    CONTAINS
!
!
!
    SUBROUTINE INICIA(MRANK)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE INICIA
!>
!> SUBROUTINE: INICIA - ?????
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
!> MRANK - 
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!> 
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : -----
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & MRANK
!
    MPI_MYRANK = MRANK
!
    END SUBROUTINE INICIA
!
!
!
    SUBROUTINE WRITEFILE2D(FILE, MATRIX, DIM1, DIM2)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE WRITEFILE2D
!>
!> SUBROUTINE: WRITEFILE2D - ?????
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
!> FILE   - 
!> MATRIX - 
!> DIM1   -
!> DIM2   -
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!> 
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : -----
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    CHARACTER(*)                                                          , INTENT(IN)          ::&
    & FILE
!
    REAL   (KIND=R4KIND), DIMENSION(:,:)                                  , INTENT(IN)          ::&
    & MATRIX
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NTSD
!
    REAL   (KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & DIM1    , DIM2
!
    CHARACTER(LEN=1000)                                                                         ::&
    & FINALFILE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J           
!
    WRITE(FINALFILE, "(A,I6.6, I4.4, A)") TRIM(FILE), NTSD, MPI_MYRANK, ".OUT"
!
    OPEN(UNIT=1000+MPI_MYRANK, FILE=FINALFILE, STATUS='REPLACE')
!
    WRITE(UNIT=1000+MPI_MYRANK, FMT='(A, I2, A)') "PROCESSADOR ", MPI_MYRANK, ":"
!
    DO J=1,DIM2
       WRITE(UNIT=1000+MPI_MYRANK, FMT='(800(F10.5))') (MATRIX(I,J), I=1,DIM1)
    END DO
!
    CLOSE(UNIT=1000+MPI_MYRANK)
!
    END SUBROUTINE WRITEFILE2D
!
!
!
    END MODULE M_PRINTVEG

