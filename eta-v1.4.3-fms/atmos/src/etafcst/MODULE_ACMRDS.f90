    MODULE ACMRDS
!>--------------------------------------------------------------------------------------------------
!> MODULE ACMRDL
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : CHKOUT
!>              GOSSIP
!>              INIT
!>              INITS
!>              RADTN
!>              RDTEMP
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
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NRDSW
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ARDSW   , TRDSW 
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & RSWIN   , RSWOUT  , RSWTOA  ,                                                               &
    & ASWIN   , ASWOUT  , ASWTOA
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & RSWTT
!
    END MODULE ACMRDS

