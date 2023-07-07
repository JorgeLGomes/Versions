    SUBROUTINE DSTRB(ARRG, ARRL, LG, LL, L1)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE DSTRB
!>
!> SUBPROGRAM: DSTRB - DISTRIBUTE GLOBAL ARRAY TO LOCAL ARRAYS
!> PROGRAMMER: BLACK
!> ORG: W/NP2
!> DATE: 97-08-29
!>
!> ABSTRACT:
!> DSTRB DISTRIBUTES THE ELEMENTS OF REAL GLOBAL ARRAY ARRG TO THE REAL LOCAL ARRAYS ARRL.  
!> LG IS THE VERTICAL DIMENSION OF THE  GLOBAL ARRAY.
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
!> LG   - VERTICAL DIMENSION OF GLOBAL ARRAY
!> LL   - VERTICAL DIMENSION OF LOCAL ARRAY
!> L1   - VERTICAL LEVEL OF ARRL BEING FILLED IN THIS CALL (USED ONLY WHEN LG=1 AND LL>1,
!>        I.E. WHEN THE GLOBAL ARRAY IS ACTUALLY JUST ONE LEVEL OF A MULTI_LEVEL ARRAY)
!>
!> OUTPUT ARGUMENT LIST:
!> ARRL - LOCAL ARRAY
!>
!> INPUT/OUTPUT ARGUMENT LIST:
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
!>              INIT
!>              INITS
!>              READ_NHB
!>
!> CALLS      : MPI_BARRIER
!>              MPI_RECV 
!>              MPI_SEND
!>--------------------------------------------------------------------------------------------------
!
!--------------------------------------- 
! DISTRIBUTE ARRAYS FROM GLOBAL TO LOCAL
!--------------------------------------- 
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
    REAL   (KIND=R4KIND), DIMENSION(IM, JM, LG)                           , INTENT(IN)          ::&
    & ARRG
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LL)         , INTENT(INOUT)       ::&
    & ARRL
!
    REAL   (KIND=R4KIND), DIMENSION(IM*JM*LG)                                                 ::&
    & ARRX
!
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LL      , LG      , L1 
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IGLB    , JGLB    , LOCI    , LOCJ    , I       , J       , L       , IPE     , KNT     ,   &
    & ISEND   , NUMVALS , IRECV   , IRTN    
!---------------------------------------------------------------------------------------
! PE0 FILLS ITS OWN LOCAL DOMAIN THEN PARCELS OUT ALL THE OTHER PIECES TO THE OTHER PEs.
!---------------------------------------------------------------------------------------
    IF (MYPE == 0) THEN
 !   
        IF (LG == 1 ) THEN
            DO JGLB=JS_GLB_TABLE(0),JE_GLB_TABLE(0)
                LOCJ = G2LJ(JGLB)
                DO IGLB=IS_GLB_TABLE(0),IE_GLB_TABLE(0)
                    LOCI = G2LI(IGLB)
                    ARRL(LOCI,LOCJ,L1) = ARRG(IGLB,JGLB,1)
                END DO
            END DO
 !       
        ELSE
 !       
            DO L=1,LG
                DO JGLB=JS_GLB_TABLE(0),JE_GLB_TABLE(0)
                    LOCJ = G2LJ(JGLB)
                    DO IGLB=IS_GLB_TABLE(0),IE_GLB_TABLE(0)
                        LOCI = G2LI(IGLB)
                        ARRL(LOCI,LOCJ,L) = ARRG(IGLB,JGLB,L)
                    END DO
                END DO
            END DO
        END IF
!    
        DO IPE=1,NPES-1
            KNT = 0
!        
            DO L=1,LG
                DO JGLB=JS_GLB_TABLE(IPE),JE_GLB_TABLE(IPE)
                    DO IGLB=IS_GLB_TABLE(IPE),IE_GLB_TABLE(IPE)
                        KNT = KNT + 1
                        ARRX(KNT) = ARRG(IGLB,JGLB,L)
                    END DO
                END DO
            END DO
 !       
            CALL MPI_SEND(ARRX, KNT, MPI_REAL, IPE, IPE, MPI_COMM_COMP, ISEND)
        END DO
!----------------------------------------------------------------------------
! ALL OTHER PEs RECEIVE THEIR PIECE FROM PE0 AND THEN FILL THEIR LOCAL ARRAY.
!----------------------------------------------------------------------------
    ELSE
        NUMVALS = (IE_GLB_TABLE(MYPE) - IS_GLB_TABLE(MYPE)+1)                                     &
    &           * (JE_GLB_TABLE(MYPE) - JS_GLB_TABLE(MYPE)+1)                                     &
    &           * LG
!
        CALL MPI_RECV(ARRX, NUMVALS, MPI_REAL, 0, MYPE, MPI_COMM_COMP, ISTAT, IRECV)
!    
        KNT = 0
        IF (LG == 1) THEN
            DO J=MY_JS_LOC,MY_JE_LOC
                DO I=MY_IS_LOC,MY_IE_LOC
                    KNT = KNT + 1
                    ARRL(I,J,L1) = ARRX(KNT)
                END DO
            END DO     
        ELSE
            DO L=1,LG
                DO J=MY_JS_LOC,MY_JE_LOC
                    DO I=MY_IS_LOC,MY_IE_LOC
                        KNT = KNT + 1
                        ARRL(I,J,L) = ARRX(KNT)
                    END DO
                END DO
            END DO
        END IF
 !   
    END IF
!
    CALL MPI_BARRIER(MPI_COMM_COMP, IRTN)
!
    RETURN
!
    END SUBROUTINE DSTRB
