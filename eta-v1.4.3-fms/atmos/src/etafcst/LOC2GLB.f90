    SUBROUTINE LOC2GLB(ARRL, ARRG)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE LOC2GLB
!>
!> SUBPROGRAM: LOC2GLB - REATE GLOBAL ARRAYS
!> PROGRAMMER: BLACK
!> ORG: W/NP22
!> DATE: 97-10-28
!>
!> ABSTRACT:
!> LOC2GLB CREATES A SINGLE GLOBAL ARRAY FROM MANY LOCAL ONES
!>
!> PROGRAM HISTORY LOG:
!> 97-10-28  BLACK   - ORIGINATOR
!> 18-01-15  LUCCI   - MODERNIZATION OF THE CODE, INCLUDING:
!>                     * F77 TO F90/F95
!>                     * INDENTATION & UNIFORMIZATION CODE
!>                     * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                     * DOCUMENTATION WITH DOXYGEN
!>                     * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> ARRL - THE LOCAL ARRAY
!>
!> OUTPUT ARGUMENT LIST:
!> ARRG - THE GLOBAL ARRAYS
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT FILES:
!> NONE
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
!> DRIVER     : DIGFLT
!>
!> CALLS      : MPI_RECV
!>              MPI_SEND
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
    INCLUDE "mpif.h"
!
#include "sp.h"
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(IN)          ::&
    & ARRL
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & ARRX
!
    REAL   (KIND=R4KIND), DIMENSION(IM, JM)                               , INTENT(INOUT)       ::&
    & ARRG 
!
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NUMVAL  , I       , J       , IPE     , IRECV   , JKNT    , JGLB    , IKNT    , IGLB    ,   &
    & ISEND 
!
    NUMVAL = (IDIM2 - IDIM1 + 1) * (JDIM2 - JDIM1 + 1)
!
    IF (MYPE /= 0) THEN
        CALL MPI_SEND(ARRL, NUMVAL, MPI_REAL, 0, MYPE, MPI_COMM_COMP, ISEND)
    ELSE
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                ARRG(I + MY_IS_GLB - 1, J + MY_JS_GLB - 1) = ARRL(I,J)
            END DO
        END DO
!    
        DO IPE=1,NPES-1
!
            CALL MPI_RECV(ARRX, NUMVAL, MPI_REAL, IPE, IPE, MPI_COMM_COMP, ISTAT, IRECV)
!        
            JKNT = 0
!
            DO J=JS_LOC_TABLE(IPE),JE_LOC_TABLE(IPE)
                JGLB = JS_GLB_TABLE(IPE) + JKNT
                IKNT = 0
!
                DO I=IS_LOC_TABLE(IPE),IE_LOC_TABLE(IPE)
                    IGLB            = IS_GLB_TABLE(IPE) + IKNT
                    ARRG(IGLB,JGLB) = ARRX(I,J)
                    IKNT            = IKNT + 1
                END DO
!
                JKNT = JKNT + 1
            END DO
!        
        END DO
!    
    END IF
!
    RETURN
!
    END SUBROUTINE LOC2GLB
