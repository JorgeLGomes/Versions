    MODULE PARMSOIL
!>--------------------------------------------------------------------------------------------------
!> MODULE PARMSOIL
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : CHKOUT
!>              EXIT
!>              GOSSIP
!>              INIT
!>              INITS
!>              QUILT
!>              RADTN
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>              SURFCE
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND), PARAMETER :: NSOIL = 8
    INTEGER(KIND=I4KIND), PARAMETER :: NROOT = 3
!
    END MODULE PARMSOIL
