    MODULE DYNAM
!>--------------------------------------------------------------------------------------------------
!> MODULE DYNAM
!>
!> USE MODULES: F77KINDS
!>              LOOPS
!>              PARMETA
!>              RDPARM
!>
!> DRIVER     : CHKOUT
!>              CUCNVC
!>              DDAMP
!>              DIGFLT
!>              DIVHOA
!>              DIVHOAST
!>              EPS
!>              GOSSIP
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
!>              VDIFH
!>              VTADV
!>              VTADVF
!>              ZENITH
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE LOOPS   , ONLY : JAM
    USE PARMETA , ONLY : LM , IDIM1, IDIM2, JDIM1, JDIM2
    USE RDPARM  , ONLY : LP1
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & PT      , R       ,                                                                         &
    & DY      , CPGFV   , EN      , ENT     ,                                                     &
    & F4D     , F4Q     ,                                                                         &
    & EF4T    
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                                         ::&
    & AETA    , DETA    ,                                                                         &                                                                             
    & RDETA   , F4Q2       
!
    REAL   (KIND=R4KIND), DIMENSION(LP1)                                                        ::&
    & ETA     ,                                                                                   &
    & DFL
!
    REAL   (KIND=R4KIND), DIMENSION(JAM)                                                        ::&
    & EM      , EMT
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & DX      , WPDAR   ,                                                                         & 
    & CPGFU   , CURV    ,                                                                         &
    & FCP     , FDIV    ,                                                                         &
    & F       ,                                                                                   &
    & DDMPU   , DDMPV   ,                                                                         &
    & FAD
!
    END MODULE DYNAM
