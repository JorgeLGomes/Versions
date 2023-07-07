    SUBROUTINE ZERO3(ARRAY, LL)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE ZERO3
!>
!> SUBROUTINE: ZERO3 - ZERO OUT 3-D ARRAYS
!> PROGRAMMER: BLACK
!> ORG: W/NP2
!> DATE: 96-03-28
!>
!> ABSTRACT:
!> SET THE VALUES OF THE ARTIFICIAL EXTERNAL (OUT-OF-BOUNDS) EDGES TO ZERO
!>
!> PROGRAM HISTORY LOG:
!> 96-03-28  BLACK      - ORIGINATOR
!> 97-06-??  MEYS       - MODIFIED FOR DISTRIBUTED MEMORY
!> 99-07-06  BLACK      - FULL ARRAY AND NOT JUST EDGES
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
!>              TOPO
!>
!> DRIVER     : INIT
!>              INITS
!>              TURBL
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
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LL)         , INTENT(INOUT)       ::&
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
    DO L=1,LL
        DO J=JDIM1,JDIM2
            DO I=IDIM1,IDIM2
                ARRAY(I,J,L) = 0.
            END DO
        END DO
    END DO
!
    RETURN
!
    END SUBROUTINE ZERO3
