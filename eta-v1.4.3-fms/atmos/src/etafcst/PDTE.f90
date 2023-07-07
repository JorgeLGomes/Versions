    SUBROUTINE PDTE
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE PDTE
!>
!> SUBROUTINE: PDTE - SURFACE PRESSURE TENDENCY CALC
!> PROGRAMMER: JANJIC
!> ORG: W/NP22
!> DATE: 94-03-08
!>
!> ABSTRACT:
!> PDTE VERTICALLY INTEGRATES THE MASS FLUX DIVERGENCE TO OBTAIN THE SURFACE PRESSURE TENDENCY AND 
!> ETADOT ON THE LAYER INTERFACES.
!>
!> PROGRAM HISTORY LOG:
!> 87-06-??  JANJIC     - ORIGINATOR
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 95-11-20  ABELES     - PARALLEL OPTIMIZATION
!> 96-03-29  BLACK      - REMOVED SCRCH
!> 98-10-30  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
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
!> OUTPUT FILES:
!> NONE
!>
!> USE MODULES: CLDWTR
!>              CONTIN
!>              CTLBLK
!>              DYNAM
!>              F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              PARMETA
!>              PVRBLS
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>
!> DRIVER     : NEWFLT
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE CLDWTR
    USE CONST
    USE CONTIN
    USE CTLBLK
    USE DYNAM
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE PARMETA
    USE PVRBLS
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND), PARAMETER :: CP   = 1004.6
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM = IM * JM -  JM /  2
    INTEGER(KIND=I4KIND), PARAMETER :: JAM  =  6 +  2 * (JM - 10)
    INTEGER(KIND=I4KIND), PARAMETER :: LM1  = LM -  1
    INTEGER(KIND=I4KIND), PARAMETER :: LP1  = LM + 1
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & RUN     , FIRST   , RESTRT  , SIGMA
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & PRET    , RPSL
!--------------------------------------------------
! COMPUTATION OF PRESSURE TENDENCY AND PREPARATIONS
!--------------------------------------------------
    DO 100 L=2,LM
!
!$omp parallel do
!
        DO 100 J=MYJS,MYJE2
            DO 100 I=MYIS,MYIE
                DIV(I,J,L) = DIV(I,J,L-1) + DIV(I,J,L)
100 END DO
!
!$omp parallel do 
!
    DO 110 J=MYJS2,MYJE2
        DO 110 I=MYIS,MYIE
            PSDT(I,J) = -DIV(I,J,LM)
            PRET(I,J) = PSDT(I,J) * RES(I,J)
            RPSL(I,J) = 1. / PDSL(I,J)
110 END DO
!--------------------- 
! COMPUTATION OF ETADT 
!--------------------- 
!
!$omp parallel do
!
    DO 120 L=1,LM1
        DO 120 J=MYJS2,MYJE2
            DO 120 I=MYIS,MYIE
                ETADT(I,J,L) = -(PRET(I,J) * ETA(L+1) + DIV(I,J,L)) * HTM(I,J,L+1) * RPSL(I,J)
120 END DO
!---------------------------------------------- 
! KINETIC ENERGY GENERATION TERMS IN T EQUATION 
!----------------------------------------------
!
!$omp parallel do 
!
    DO 130 J=MYJS2,MYJE2
        DO 130 I=MYIS,MYIE
            OMGALF(I,J,1) = OMGALF(I,J,1) - DIV(I,J,1) * RTOP(I,J,1) * EF4T
                 T(I,J,1) =      T(I,J,1) - DIV(I,J,1) * RTOP(I,J,1) * EF4T
130 END DO
!
!$omp parallel do
!
    DO 140 J=MYJS2,MYJE2
        DO 140 I=MYIS,MYIE
            OMGALF(I,J,L) = OMGALF(I,J,L) - (DIV(I,J,L-1) + DIV(I,J,L)) * RTOP(I,J,L) * EF4T      &
    &                     *    HTM(I,J,L)
!
                 T(I,J,L) =      T(I,J,L) - (DIV(I,J,L-1) + DIV(I,J,L)) * RTOP(I,J,L) * EF4T
140 END DO
145 CONTINUE
!
!$omp parallel do
!
    DO 150 J=MYJS2,MYJE2
        DO 150 I=MYIS,MYIE
            OMGALF(I,J,LM) = OMGALF(I,J,LM) + (PRET(I,J) - DIV(I,J,LM1)) * RTOP(I,J,LM) * EF4T    &
    &                      *    HTM(I,J,LM)
!
                 T(I,J,LM) =      T(I,J,LM) + (PRET(I,J) - DIV(I,J,LM1)) * RTOP(I,J,LM) * EF4T
150 END DO
!
    DO 160 L=LM,2,-1
!
!$omp parallel do
!
        DO 160 J=MYJS,MYJE2
            DO 160 I=MYIS,MYIE
                DIV(I,J,L) = DIV(I,J,L) - DIV(I,J,L-1)
160 END DO
!
    RETURN
!
    END SUBROUTINE PDTE
