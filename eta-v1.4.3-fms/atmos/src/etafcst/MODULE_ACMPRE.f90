    MODULE ACMPRE
!>--------------------------------------------------------------------------------------------------
!> MODULE ACMPRE
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : CHKOUT
!>              GOSSIP
!>              INIT
!>              INITS
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>              SURFCE
!>-------------------------------------------------------------------------------------------------- 
    USE F77KINDS
    USE PARMETA
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & TPREC
!   
    REAL   (KIND=R4KIND), DIMENSION (IDIM1:IDIM2, JDIM1:JDIM2)                                  ::&
    & ACSNOW  , ACSNOM  , SSROFF  , BGROFF
!
    END MODULE ACMPRE
