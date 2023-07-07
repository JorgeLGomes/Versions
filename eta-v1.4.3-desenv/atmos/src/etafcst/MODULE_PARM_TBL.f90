    MODULE PARM_TBL
!>--------------------------------------------------------------------------------------------------
!> MODULE PARM_TBL
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : CHKOUT
!>              CUCNVC
!>              DIGFLT
!>              EXIT
!>              GOSSIP
!>              GSCOND
!>              GSMDRIVE     
!>              HDIFF
!>              HDIFFS
!>              INIT
!>              INITS
!>              NEWFLT
!>              PRECPD
!>              RADTN
!>              RDTEMP
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>              SURFCE
!>              WRTRST
!>              ZENITH
!>--------------------------------------------------------------------------------------------------  
    USE F77KINDS
! 
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND), PARAMETER :: ITB  = 76 
    INTEGER(KIND=I4KIND), PARAMETER :: JTB  = 134
    INTEGER(KIND=I4KIND), PARAMETER :: ITBQ = 152
    INTEGER(KIND=I4KIND), PARAMETER :: JTBQ = 440
!
    END MODULE PARM_TBL
