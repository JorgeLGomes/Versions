    MODULE ACMCLD
!>--------------------------------------------------------------------------------------------------
!> MODULE ACMCLD
!>
!> USE MODULES: F77KINDS
!>              PARMETA       
!> 
!> DRIVER     : CHKOUT
!>              GOSSIP
!>              INIT
!>              INITS
!>              RADTN
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>-------------------------------------------------------------------------------------------------- 
    USE F77KINDS
    USE PARMETA
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NCLOD
!                                                                     
    REAL   (KIND=R4KIND)                                                                        ::&
    & TCLOD
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & NCFRCV  , NCFRST
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & ACFRCV  , ACFRST
!
    END MODULE ACMCLD
