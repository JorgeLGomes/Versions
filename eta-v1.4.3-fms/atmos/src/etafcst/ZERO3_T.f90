    SUBROUTINE ZERO3_T(ARRAY, LL)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE ZERO3_T
!> 
!> SUBROUTINE: ZERO3_T - ZERO OUT TRANSPOSED 3-D ARRAYS
!> PROGRAMMER: BLACK
!> ORG: W/NP2
!> DATE: 99-07-06
!>
!> ABSTRACT:
!> SET THE VALUES OF THE ARTIFICIAL EXTERNAL (OUT-OF-BOUNDS) EDGES TO ZERO
!>
!> PROGRAM HISTORY LOG:
!> 99-07-06  BLACK      - ORIGINATOR
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> LL    -
!>
!> OUTPUT ARGUMENT LIST:
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> ARRAY - 
!>
!> OUTPUT FILES:
!> NONE
!>
!> USE MODULES: F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              TEMPCOM
!>              TOP
!> 
!> DRIVER     : TURBL
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
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
#include "sp.h"
!
    REAL   (KIND=R4KIND), DIMENSION(LL, IDIM1:IDIM2, JDIM1:JDIM2)         , INTENT(INOUT)       ::&
    & ARRAY
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LL
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , L
!
    DO J=JDIM1,JDIM2
        DO I=IDIM1,IDIM2
            DO L=1,LL
                ARRAY(L,I,J) = 0.
            END DO
        END DO
    END DO
!
    RETURN
!
    END SUBROUTINE ZERO3_T
