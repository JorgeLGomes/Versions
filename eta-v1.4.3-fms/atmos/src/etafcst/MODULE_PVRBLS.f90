    MODULE PVRBLS
!>--------------------------------------------------------------------------------------------------
!> MODULE PVRBLS
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : BOCOH
!>              BOCOHF
!>              CHKOUT
!>              CUCNVC
!>              DDAMP
!>              DIGFLT
!>              EBU
!>              EPS
!>              EXIT
!>              GOSSIP
!>              GSCOND
!>              GSMDRIVE
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
!>              PRECPD
!>              RADTN
!>              RDTEMP
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>              SLADVT
!>              SURFCE
!>              TURBL
!>              TWR
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
    & Z0      , USTAR   ,                                                                         &
    & UZ0     , VZ0     ,                                                                         &
    & THZ0    , QZ0     ,                                                                         &
    & THS     , QS      ,                                                                         &
    & AKMS    , AKHS    ,                                                                         &
    & RF      ,                                                                                   &
    & TWBS    , QWBS    ,                                                                         &
    & SNO     , SI      ,                                                                         &
    & CLDEFI  ,                                                                                   &
    & PREC    , ACPREC  ,                                                                         &
    & ACCLIQ  , CUPREC  ,                                                                         &
    & TH100   , Q100    ,                                                                         &
    & U100    , V100    ,                                                                         &
    & TH10    , Q10     ,                                                                         &
    & U10     , V10     ,                                                                         &
    & TSHLTR  , QSHLTR  ,                                                                         &
    & TSHLTR1 ,                                                                                   &
    & PSHLTR  ,                                                                                   &
    & XMOMFLUX,                                                                                   &
    & YMOMFLUX,                                                                                   &
    & AKM10   , AKM10V                                                                          
!    & PLM                                                                                         
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & Q2
!
    END MODULE PVRBLS
