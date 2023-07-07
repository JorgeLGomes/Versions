    MODULE CNVCLD
!>--------------------------------------------------------------------------------------------------
!> MODULE CNVCLD
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : CHKOUT
!>              CUCNVC
!>              GOSSIP
!>              INIT
!>              INITS
!>              RADTN
!>              READ_NHB
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
    & CUPPT   , CFRACL  , CFRACM  , CFRACH
!
    END MODULE CNVCLD

