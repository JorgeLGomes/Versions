    MODULE CLDWTR
!>--------------------------------------------------------------------------------------------------
!> MODULE CLDWTR
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : BOCOH
!>              BOCOHF
!>              CHKOUT
!>              DDAMP
!>              DIVHOA
!>              EBU
!>              EXIT
!>              GOSSIP
!>              GSCOND
!>              GSMDRIVE
!>              HDIFF
!>              HDIFFS
!>              HZADV
!>              HZADV2
!>              HZADVS
!>              HZADV_LM1
!>              INIT
!>              INITS
!>              NEWFLT
!>              PDTE
!>              PRECPD
!>              RADTN
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>              VADZ
!>              VTADV
!>              WRTRST
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : IDIM1, IDIM2, JDIM1, JDIM2, LM
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & CWM
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & LC
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & U00     , SR
!
    REAL   (KIND=R4KIND), DIMENSION(2*LM)                                                       ::&
    & UL
!
    END MODULE CLDWTR
