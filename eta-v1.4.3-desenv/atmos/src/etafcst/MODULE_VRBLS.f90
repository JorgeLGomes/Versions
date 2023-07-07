    MODULE VRBLS
!>-------------------------------------------------------------------------------------------------- 
!> MODULE VRBLS
!>
!> USE MODULES: PARMETA
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
!>              EXIT
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
!>              PDNEW
!>              PDTE
!>              PDTEDT
!>              PGCOR
!>              PRECPD
!>              RADTN
!>              RDTEMP
!>              READ_NHB
!>              READ_RESTRT
!>              RES_RESTRT2
!>              SLADVT
!>              SURFCE
!>              TURBL
!v              TWR
!>              VADZ
!>              VTADV
!>              VTADVF
!>              VWR
!>              WRTRST
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : LM, IDIM1, IDIM2, JDIM1, JDIM2
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & PD      , FIS     , RES 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & T       , U       , V       , Q
!
    END MODULE VRBLS

