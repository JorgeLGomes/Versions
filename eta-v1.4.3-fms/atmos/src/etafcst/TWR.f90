    SUBROUTINE TWR
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE TWR
!>
!> SUBPROGRAM: TWR - ?????
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
!>              PVRBLS
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
    USE PVRBLS
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
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE,4)                                          ::&
    & STATUS_ARRAY
!
    REAL   (KIND=R4KIND), DIMENSION(IM, JM)                                                     ::&
    & TWRITE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    &  ARRAY_LOC
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IOUT    , I       , J       , L       , IPE     , MAXVALS , IRECV   , JOFFSET , JGLB    ,   &
    & JLOC    , IOFFSET , IGLB    , NUMVALS , ISEND   , IRTN    , IENDX
!
    IOUT = 80
!
    REWIND IOUT
!
    DO 500 L=1,LM
!    
        IF (MYPE == 0) THEN
            DO J=MY_JS_LOC,MY_JE_LOC
                DO I=MY_IS_LOC,MY_IE_LOC
                    TWRITE(I + MY_IS_GLB - 1, J + MY_JS_GLB - 1) = T(I,J,L)
                END DO
            END DO
!        
            DO IPE=1,NPES-1
                MAXVALS = (IDIM2 - IDIM1 + 1) * (JDIM2 - JDIM1 + 1)
!            
                CALL MPI_RECV(ARRAY_LOC, MAXVALS, MPI_REAL, IPE, IPE, MPI_COMM_COMP, JSTAT, IRECV)
!            
                JOFFSET = 0
!
                DO JGLB=JS_GLB_TABLE(IPE),JE_GLB_TABLE(IPE)
                    JLOC    = JS_LOC_TABLE(IPE) + JOFFSET
                    IOFFSET = 0
!
                    DO IGLB=IS_GLB_TABLE(IPE),IE_GLB_TABLE(IPE)
                        TWRITE(IGLB,JGLB) = ARRAY_LOC(IS_LOC_TABLE(IPE) + IOFFSET, JLOC)
                        IOFFSET           = IOFFSET                     + 1
                    END DO
!
                    JOFFSET = JOFFSET + 1
!
                END DO
            END DO
!        
        ELSE
!
            NUMVALS = (IDIM2 - IDIM1 + 1) * (JDIM2 - JDIM1 + 1)
!												   
            CALL MPI_SEND(T(IDIM1,JDIM1,L), NUMVALS, MPI_REAL, 0, MYPE,MPI_COMM_COMP, ISEND)
!
        END IF
!    
        CALL MPI_BARRIER(MPI_COMM_COMP, IRTN)
!    
        IF (MYPE == 0) THEN
            DO J=1,JM
                IENDX = IM
!
                IF (MOD(J,2) == 0) IENDX = IM - 1
!
                WRITE(IOUT) (TWRITE(I,J), I=1,IENDX)
            END DO
        END IF
!
500 END DO
!
    RETURN
!
    END SUBROUTINE TWR
