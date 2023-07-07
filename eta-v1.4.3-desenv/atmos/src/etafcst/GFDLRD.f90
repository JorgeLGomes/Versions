    BLOCK DATA  GFDLRD
!>--------------------------------------------------------------------------------------------------
!> BLOCK DATA GFDLRD
!> 
!> BLOCK DATA: GFDLRD - ?????
!> PROGRAMMER: ?????
!> ORG: ?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> ?????
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????      - ORIGINATOR
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                       * F77 TO F90/F95
!>                       * INDENTATION & UNIFORMIZATION CODE
!>                       * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                       * DOCUMENTATION WITH DOXYGEN
!>                       * OPENMP FUNCTIONALITY
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
!> USE MODULES: CO2DTA
!>              DIUCON
!>              HCON
!>              INPUT 
!>              PARMETA
!>              PHYCON
!>              RDPARM
!>              RNDDTA
!>              SSALB
!>              SWRSAV
!>              TBLTMP
!>
!> DRIVER     : -----
!>
!> CALLS      : -----  
!--------------------------------------------------------------------------------------------------
    USE CO2DTA
    USE DIUCON
    USE HCON
    USE INPUT
    USE PARMETA
    USE PHYCON
    USE RDPARM
    USE RNDDTA
    USE SSALB
    USE SWRSAV
    USE TBLTMP
!
    IMPLICIT NONE
!
#include "sp.h"
!
    END BLOCK DATA GFDLRD
