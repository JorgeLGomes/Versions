    MODULE LOOPS
!>--------------------------------------------------------------------------------------------------
!> MODULE LOOPS
!>
!> USE MODULES: F77KINDS
!>              PARAMETER
!>
!> DRIVER     : CHKOUT
!>              CUCNVC
!>              DIVHOA
!>              DIVHOAST
!>              EPS
!>              EXIT
!>              GOSSIP
!>              GSCOND
!>              GSMDRIVE
!>              HADZ
!>              HZADV
!>              HZADVS
!>              HZADV2
!>              HZADV_LM1
!>              INIT
!>              INITS
!>              PDTEDT
!>              PGCOR
!>              PRECPD
!>              RADTN
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>              SLADVT
!>              SURFCE
!>              TURBL
!>              VADZ
!>              WRTRST
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : JM, IDIM1, IDIM2, JDIM1, JDIM2
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND), PARAMETER :: JAM = 6 + 2 * (JM - 10)

    INTEGER(KIND=I4KIND), DIMENSION(JAM)                                                        ::&
    & IHLA    , IHHA    ,                                                                         &
    & IVLA    , IVHA    ,                                                                         &
    & JRA
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & LMH     , LMV
!
    END MODULE LOOPS

