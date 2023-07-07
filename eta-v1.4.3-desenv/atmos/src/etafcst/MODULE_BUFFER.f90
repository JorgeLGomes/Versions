    MODULE BUFFER
!>--------------------------------------------------------------------------------------------------
!> MODULE BUFFER
!>
!> USE MODULES: F77KINDS
!>              PARMBUF 
!>
!> DRIVER     : CHKOUT
!>              QUILT    
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMBUF
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IP
!
    REAL   (KIND=R4KIND), DIMENSION(IBUFMAX)                                                    ::&
    & BUF
!
   END MODULE BUFFER
