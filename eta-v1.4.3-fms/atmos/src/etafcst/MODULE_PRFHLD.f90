    MODULE PRFHLD
!>--------------------------------------------------------------------------------------------------
!> MODULE PRFHLD
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!> 
!> DRIVER     : CHKOUT
!>              READ_RESTRT
!>              READ_RESTRT2
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : IDIM1, IDIM2, JDIM1, JDIM2
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & TLMIN   , TLMAX   , MAXWIND   , MAXWU   , MAXWV
!
    END MODULE PRFHLD
