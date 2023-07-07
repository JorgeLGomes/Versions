    SUBROUTINE PDNEW
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE PDNEW
!>
!> SUBROUTINE: PDNEW - UPDATE SURFACE PRESSURE
!> PROGRAMMER: JANJIC
!> ORG: W/NP22
!> DATE: 94-03-08
!>
!> ABSTRACT:
!> PDNEW UPDATES THE SURFACE PRESSURE FROM THE TENDENCY COMPUTED IN PDTE.
!>
!> PROGRAM HISTORY LOG:
!> 87-06-??  JANJIC     - ORIGINATOR
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
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
!> OUTPUT FILES:
!> NONE
!>
!> USE MODULES: CONTIN
!>              CTLBLK
!>              F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              TEMPCOM
!>              TOPO     
!>              VRBLS     
!>
!> DRIVER     : DIGFLT
!>              EBU
!>              NEWFLT       
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE CONTIN
    USE CTLBLK
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
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J
!
#include "sp.h"
!----------------------------- 
! UPDATING PRESSURE DIFFERENCE 
!----------------------------- 
    DO J=MYJS2,MYJE2
        DO I=MYIS,MYIE
            PD(I,J) = PSDT(I,J) * DT + PD(I,J)
        END DO
    END DO
!
    RETURN
!
    END SUBROUTINE PDNEW
