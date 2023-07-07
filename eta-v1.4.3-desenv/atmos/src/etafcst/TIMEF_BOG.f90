    FUNCTION TIMEF()
!>--------------------------------------------------------------------------------------------------
!> FUNCTION TIMEF
!>
!> FUNCTION: TIMEF - DESIGNED TO DUPLICATE TIMEF
!> PROGRAMMER: ?????   
!> ORG: ?????
!> DATE: ??-??-??
!> 
!> ABSTRACT:
!> DESIGNED TO DUPLICATE TIMEF
!>
!> PROGRAM HISTORY LOG:
!> 99-12-??  M. PYLE     - ORIGINATOR
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
!>
!> DRIVER     : CHKOUT
!>              EBU
!>              INIT
!>              INITS
!>              PDTEDT
!>              QUILT
!>              TURBL
!>-------------------------------------------------------------------------------------------------- 
    USE F77KINDS
!
    REAL   (KIND=R4KIND), DIMENSION(2)                                                          ::&
    & ET
!
    REAL   (KIND=R8KIND)                                                                        ::&
    & TIMEF
!
    TIMEF = ETIME(ET)
    TIMEF = TIMEF * 1.E3
!
    END FUNCTION TIMEF

