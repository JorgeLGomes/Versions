    MODULE CONTIN
!>--------------------------------------------------------------------------------------------------
!> MODULE CONTIN
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : CHKOUT
!>              DDAMP
!>              DIGFLT
!>              DIVHOA
!>              DIVHOAST
!>              EBU
!>              EPS
!>              GOSSIP
!>              HADZ
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
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>              SLADVT
!>              VADZ
!>              VTADV
!>              VTADVF
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : IDIM1, IDIM2, JDIM1, JDIM2, LM
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & PDSL    , PDSLO   , PSDT
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM  )                             ::&
    & RTOP    , OMGALF  , DIV
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM-1)                             ::&
    & ETADT
! 
    END MODULE CONTIN
