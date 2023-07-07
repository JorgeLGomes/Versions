    SUBROUTINE ZERO2(ARRAY)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE ZERO2
!>
!> SUBROUTINE: ZERO2 - ZERO OUT 2-D ARRAYS
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
!> 99-07-06  BLACK      - FULL ARRAY RATHER THAN JUST EDGES
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NONE
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> ARRAY - THE DUMMY ARRAY NAME
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
!> DRIVER     : CUCNVC
!>              DDAMP
!>              DIGFLT
!>              DIVHOA
!>              DIVHOAST
!>              HDIFF
!>              HDIFFS
!>              HZADAV
!>              HZADAV_LM1
!>              HZADVS
!>              INIT
!>              INITS
!>              PDTEDT
!>              PGCOR
!>              RADTN
!>              SURFCE
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
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(INOUT)       ::&
    & ARRAY
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J 
!
    DO J=JDIM1,JDIM2
        DO I=IDIM1,IDIM2
            ARRAY(I,J) = 0.
        END DO
    END DO
!
    RETURN
!
    END SUBROUTINE ZERO2
