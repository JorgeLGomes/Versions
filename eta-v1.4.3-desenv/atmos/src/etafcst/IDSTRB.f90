    SUBROUTINE IDSTRB(ARRG, ARRL)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE IDSTRB
!> 
!> SUBPROGRAM: IDSTRB - DISTRIBUTE INTEGER GLOBAL ARRAY TO LOCAL ONES
!> PROGRAMMER: BLACK
!> ORG: W/NP2
!> DATE: 97-08-29
!>
!> ABSTRACT:
!> IDSTRB DISTRIBUTES THE ELEMENTS OF INTEGER GLOBAL ARRAY ARRG TO THE INTEGER LOCAL ARRAYS ARRL. 
!> LG IS THE VERTICAL DIMENSION OF THE GLOBAL ARRAY. 
!> LL IS THE VERTICAL DIMENSION OF THE LOCAL ARRAY. 
!> L1 IS THE SPECIFIC LEVEL OF ARRL THAT IS BEING FILLED DURING THIS CALL 
!> (PERTINENT WHEN LG=1 AND LL>1).
!>
!> PROGRAM HISTORY LOG:
!> 97-08-29  BLACK      - ORIGINATOR
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY 
!>
!> INPUT ARGUMENT LIST:
!> ARRG - GLOBAL ARRAY
!>
!> OUTPUT ARGUMENT LIST:
!> ARRL - LOCAL ARRAY
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: GLB_TABLE
!>              F77KINDS
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : READ_NHB
!>              READ_RESTRT
!>
!> CALLS      : MPI_BARRIER
!>              MPI_RECV
!>              MPI_SEND
!>--------------------------------------------------------------------------------------------------
    USE GLB_TABLE
    USE F77KINDS
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!--------------------------------------- 
! DISTRIBUTE ARRAYS FROM GLOBAL TO LOCAL
!---------------------------------------
    INCLUDE "mpif.h"
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), DIMENSION(IM, JM)                               , INTENT(IN)          ::&
    & ARRG     
!
    INTEGER(KIND=I4KIND), DIMENSION(IM*JM)                                                      ::&
    & ARRX   
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2,JDIM1:JDIM2)              , INTENT(INOUT)       ::&
    & ARRL
!
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , JGLB    , LOCJ    , IGLB    , LOCI    , IPE     , KNT     , ISEND   ,   &
    & NUMVALS , IRECV   , IRTN 
!--------------------------------------------------------------------------------------- 
! PE0 FILLS ITS OWN LOCAL DOMAIN THEN PARCELS OUT ALL THE OTHER PIECES TO THE OTHER PEs.
!--------------------------------------------------------------------------------------- 
    IF (MYPE == 0) THEN
!    
        DO JGLB=JS_GLB_TABLE(0),JE_GLB_TABLE(0)
            LOCJ = G2LJ(JGLB)
            DO IGLB=IS_GLB_TABLE(0),IE_GLB_TABLE(0)
                LOCI = G2LI(IGLB)
                ARRL(LOCI,LOCJ) = ARRG(IGLB,JGLB)
            END DO
        END DO
!    
        DO IPE=1,NPES-1
!
            KNT = 0
!        
            DO JGLB=JS_GLB_TABLE(IPE),JE_GLB_TABLE(IPE)
                DO IGLB=IS_GLB_TABLE(IPE),IE_GLB_TABLE(IPE)
                    KNT = KNT + 1
                    ARRX(KNT) = ARRG(IGLB,JGLB)
                END DO
            END DO
!        
            CALL MPI_SEND(ARRX, KNT, MPI_INTEGER, IPE, IPE, MPI_COMM_COMP, ISEND)
        END DO
!----------------------------------------------------------------------------  
! ALL OTHER PEs RECEIVE THEIR PIECE FROM PE0 AND THEN FILL THEIR LOCAL ARRAY.
!---------------------------------------------------------------------------- 
    ELSE
        NUMVALS = (IE_GLB_TABLE(MYPE) - IS_GLB_TABLE(MYPE) + 1)                                   &
    &           * (JE_GLB_TABLE(MYPE) - JS_GLB_TABLE(MYPE) + 1)
! 
        CALL MPI_RECV(ARRX, NUMVALS, MPI_INTEGER, 0, MYPE, MPI_COMM_COMP, ISTAT, IRECV)
!    
        KNT = 0
!
        DO J=MY_JS_LOC,MY_JE_LOC
            DO I=MY_IS_LOC,MY_IE_LOC
                KNT = KNT + 1
                ARRL(I,J) = ARRX(KNT)
            END DO
        END DO
!    
    END IF
!
    CALL MPI_BARRIER(MPI_COMM_COMP, IRTN)
!
    RETURN
!
    END SUBROUTINE IDSTRB
