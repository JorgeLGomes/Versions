    SUBROUTINE DIST(A, B)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE DIST
!> 
!> SUBROUTINE: DIST - DISTRIBUTES DATA FROM TASK 0
!> PROGRAMMER: TUCCILLO
!> ORG: IBM
!> DATE: 00-01-20
!>
!> ABSTRACT: DISTRIBUTES DATA FROM TASK 0
!>
!> PROGRAM HISTORY LOG:
!> 00-01-20  TUCCILLO - ORIGINATOR
!> 18-01-15  LUCCI    - MODERNIZATION OF THE CODE, INCLUDING:
!>                      * F77 TO F90/F95
!>                      * INDENTATION & UNIFORMIZATION CODE
!>                      * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                      * DOCUMENTATION WITH DOXYGEN
!>                      * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> A - ARRAY TO BE DISTRIBUTED
!>
!> OUTPUT ARGUMENT LIST:
!> B - DISTRIBUTED DATA
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MPPCOM
!>              PARA
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>  
!> DRIVER     : -----
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARA
    USE PARMETA
    USE TEMPCOM
    USE TOPO
! 
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!
    REAL   (KIND=R4KIND), DIMENSION(IM, JM)                               , INTENT(IN)          ::&
    & A
!
    REAL   (KIND=R4KIND), DIMENSION(IM, MY_JSD:MY_JED)                    , INTENT(INOUT)       ::&
    & B
!
    REAL   (KIND=R4KIND), DIMENSION(IM*JM)                                                      ::&
    & BUF
!
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & STATUS
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , IERR
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & K       , II      , III    
!
    IF (ME == 0) THEN
 !   
        DO K=JSTA(ME),JEND(ME)
            DO J=MY_JS_GLB_A(K),MY_JE_GLB_A(K)
                DO I=MY_IS_GLB_A(K),MY_IE_GLB_A(K)
                    B(I,J) = A(I,J)
                END DO
            END DO 
        END DO
!----------------------
! SEND TO EVERYONE ELSE
!----------------------
        DO II=1,NUM_PROCS-1
            III = 0
            DO K=JSTA(II),JEND(II)
                DO J=MY_JS_GLB_A(K),MY_JE_GLB_A(K)
                    DO I=MY_IS_GLB_A(K),MY_IE_GLB_A(K)
                        III = III + 1
                        BUF(III) = A(I,J)
                    END DO
                END DO
            END DO
!
            CALL MPI_SEND(BUF, III  , MPI_REAL, II, II, MPI_COMM_COMP, IERR)
!
        END DO
!
    ELSE
!
        CALL MPI_RECV(BUF, IM*JM, MPI_REAL, 0 , ME, MPI_COMM_COMP, STATUS, IERR)
!
        III = 0
        DO K=JSTA(ME),JEND(ME)
            DO J=MY_JS_GLB_A(K),MY_JE_GLB_A(K)
                DO I=MY_IS_GLB_A(K),MY_IE_GLB_A(K)
                    III = III + 1
                    B(I,J) = BUF(III)
                END DO
            END DO
        END DO
!    
    END IF
!
    END SUBROUTINE DIST
