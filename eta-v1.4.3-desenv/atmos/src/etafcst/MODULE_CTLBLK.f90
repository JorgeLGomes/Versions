    MODULE CTLBLK
!>--------------------------------------------------------------------------------------------------
!> MODULE CTLBLK
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>              RDPARM
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
!>              EBU
!>              EPS
!>              GSCOND
!>              GSMDRIVE
!>              HADZ
!>              HDIFF
!>              HDIFFS
!>              HZADV
!>              HZADVS
!>              HZADV2
!>              HZADV_LM1
!>              INIT
!>              INITS
!>              NEWFLT
!>              PDNEW
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
!>              SOLARD
!>              SURFCE
!>              TURBL
!>              VADZ
!>              VTADV
!>              VTADVF
!>              ZENITH
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA
    USE RDPARM
!
    IMPLICIT NONE
!
    SAVE
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & RUN     , FIRST   , RESTRT  , SIGMA   , NEST    , SINGLRST, SUBPOST
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IHRST   , NFCST   , NBC     , LIST    , IOUT    , NTSD    , NTSTM   , NSTART  , NTDDMP  ,   &
    & NPREC   , IDTAD   , NBOCO   , NSHDE   , NCP     , NPHS    , NCNVC   , NRADS   , NRADL 
!
    INTEGER(KIND=I4KIND), DIMENSION(3)                                                          ::&
    & IDAT
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & DT
! 
    END MODULE CTLBLK

