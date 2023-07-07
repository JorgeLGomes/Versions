    MODULE ACMRDL
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
    & NRDLW
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ARDLW   , TRDLW 
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & RLWIN   , RLWOUT  , RLWTOA  ,                                                               &
    & ALWIN   , ALWOUT  , ALWTOA
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & RLWTT
!
    END MODULE ACMRDL
