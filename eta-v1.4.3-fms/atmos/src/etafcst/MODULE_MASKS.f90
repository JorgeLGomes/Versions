    MODULE MASKS
!>--------------------------------------------------------------------------------------------------
!> MODULE MASKS
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!> 
!> DRIVER     : BOCOH
!>              BOCOHF
!>              BOCOV
!>              CHKOUT
!>              CLTEND
!>              CUCNVC
!>              DDAMP
!>              DIGFLT
!>              DIVHOA
!>              DIVHOAST
!>              EPS
!>              EXIT
!>              GOSSIP
!>              GSCOND
!>              GSMDRIVE
!>              HADZ
!>              HDIFF
!>              HDIFFS
!>              HZADV
!>              HZADV2
!>              HZADV_LM1
!>              HZADVS
!>              INIT
!>              INITS
!>              NEWFLT
!>              PDTE
!>              PDTEDT
!>              PGCOR
!>              PRECPD
!>              RADTN
!>              RDTEMP
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>              SLADVT
!>              SURFCE
!>              TURBL
!>              VADZ
!>              VTADV
!>              VTADVF
!>              WRTRST
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : IDIM1, IDIM2, JDIM1, JDIM2, LM
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & VBM2    , VBM3    ,                                                                         &
    & SM      , SICE    ,                                                                         &
    & HBM2    , HBM3
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & HTM     , VTM
!
    END MODULE MASKS
