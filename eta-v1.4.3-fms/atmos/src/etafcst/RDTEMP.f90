    SUBROUTINE RDTEMP
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE RDTEMP
!>
!> SUBROUTINE: RDTEMP - RADIATIVE TEMPERATURE CHANGE
!> PROGRAMMER: BLACK
!> ORG: W/NP22
!> DATE: 93-12-29
!>
!> ABSTRACT:
!> RDTEMP APPLIES THE TEMPERATURE TENDENCIES DUE TO RADIATION AT ALL LAYERS AT EACH ADJUSTMENT TIME
!> STEP
!>
!> PROGRAM HISTORY LOG:
!> 87-09-??  BLACK      - ORIGINATOR
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 95-11-20  ABELES     - PARALLEL OPTIMIZATION
!> 98-10-30  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
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
!> NONE 
!>
!> USE MODULES: ACMRDL
!>              ACMRDS
!>              CTLBLK
!>              F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              PARMETA
!>              PARM_TBL
!>              PHYS
!>              PVRBLS
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>
!> DRIVER     : EBU
!>              NEWFLT
!>
!> CALLS      : ZENITH              
!>--------------------------------------------------------------------------------------------------
    USE ACMRDL
    USE ACMRDS
    USE CTLBLK
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE PARMETA
    USE PARM_TBL
    USE PHYS
    USE PVRBLS
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
#include "sp.h"
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & FACTR
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K  
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & TIMES   , DAYI    , HOUR    , TTNDKL 
!---------------------------------------- 
!  GET CURRENT VALUE OF COS(ZENITH ANGLE)
!----------------------------------------
    TIMES = (NTSD-1) * DT
!
    CALL ZENITH(TIMES, DAYI, HOUR)
!
!$omp parallel do 
!
    DO 50 J=MYJS,MYJE
        DO 50 I=MYIS,MYIE
            IF (CZMEAN(I,J) > 0.) THEN
                FACTR(I,J) = CZEN(I,J) / CZMEAN(I,J)
            ELSE
                FACTR(I,J) = 0.
            END IF
 50 END DO
!
!$omp parallel do private (TTNDKL)
!
    DO 100 K=1,LM
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                TTNDKL = RSWTT(I,J,K) * FACTR(I,J) + RLWTT(I,J,K)
                T(I,J,K) = T(I,J,K) + TTNDKL * DT * HTM(I,J,K) * HBM2(I,J)
                IF (MYPE == 4 .AND. I == 16 .AND. J == 11 .AND. K >= 30) THEN
                END IF
                633 FORMAT('K,TNEW,INCR,SW,LW: ',I2,1X,F8.3,1X,F8.5,1X,E12.6,1X,E12.6)
            END DO
        END DO
100 END DO
!
    RETURN
!
    END SUBROUTINE RDTEMP
