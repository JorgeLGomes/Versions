    SUBROUTINE VWR
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE VWR
!>
!> SUBPROGRAM: VWR - ?????
!> PROGRAMMER: ?????   
!> ORG: ?????
!> DATE: ??-??-??
!> 
!> ABSTRACT: 
!> ?????
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????       - ORIGINATOR
!> 18-01-15  LUCCI       - MODERNIZATION OF THE CODE, INCLUDING:
!>                         * F77 TO F90/F95
!>                         * INDENTATION & UNIFORMIZATION CODE
!>                         * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                         * DOCUMENTATION WITH DOXYGEN
!>                         * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NONE
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
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
!>              VRBLS
!>
!> DRIVER     : -----
!>
!> CALLS      : MPI_BARRIER
!>              MPI_RECV
!>              MPI_SEND
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & JSTAT
!
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE, 4)                                         ::&
    & STATUS_ARRAY
!
    REAL   (KIND=R4KIND), DIMENSION(IM, JM)                                                     ::&
    & VWRITE
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , L       , IOUT    , IRTN    , IPE     , ISTAT   , IRECV   , ISEND   ,   &
    & IENDX
!
    IOUT = 80
!
    DO 500 L=1,LM
!    
        DO J=JS_LOC_TABLE(MYPE),JE_LOC_TABLE(MYPE)
            DO I=IS_LOC_TABLE(MYPE),IE_LOC_TABLE(MYPE)
                VWRITE(I+MY_IS_GLB-1,J+MY_JS_GLB-1) = V(I,J,L)
            END DO
        END DO
!    
        CALL MPI_BARRIER(MPI_COMM_COMP, IRTN)
!    
        IF (MYPE == 0) THEN
            DO IPE=1,NPES-1
                DO J=JS_GLB_TABLE(IPE),JE_GLB_TABLE(IPE)
                    CALL MPI_RECV(VWRITE(IS_GLB_TABLE(IPE),J), ICHUNKTAB(IPE),                    &
    &                             MPI_REAL, IPE, 99, MPI_COMM_COMP, ISTAT, IRECV)
                END DO
            END DO
!        
        ELSE
!
            DO J=MY_JS_GLB,MY_JE_GLB
                CALL MPI_SEND(VWRITE(MY_IS_GLB,J), ICHUNKTAB(MYPE),                               &
    &                                MPI_REAL, 0, 99, MPI_COMM_COMP, ISEND)
            END DO
        END IF
!    
        IF (MYPE == 0) THEN
            DO J=1,JM
                IENDX=IM
                IF (MOD(J,2) == 1) IENDX = IM-1
                WRITE(IOUT) (VWRITE(I,J),I=1,IENDX)
            END DO
        END IF
!
500 END DO
!
    STOP 555
!
    RETURN
    END SUBROUTINE VWR
