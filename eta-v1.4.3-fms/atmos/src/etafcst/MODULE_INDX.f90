    MODULE INDX
!>--------------------------------------------------------------------------------------------------
!> MODULE INDX
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : CHKOUT
!>              CUCNVC
!>              DDAMP
!>              DIVHOA
!>              DIVHOAST
!>              EPS
!>              EXIT
!>              FILT25
!>              GOSSIP
!>              HADZ
!>              HDIFF
!>              HDIIFS
!>              HZADV
!>              HZADVS
!>              HZADV2
!>              HZADV_LM1
!>              INIT
!>              INITS
!>              PDTEDT
!>              PGCOR
!>              RADTN
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>              SLADVT
!>              TURBL
!>              VADZ
!>              VTADV
!>              VTADVF     
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : IM, JM, IDIM1, IDIM2, JDIM1, JDIM2
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND), DIMENSION(JDIM1:JDIM2)                                                ::&
    & IHE     , IHW     ,                                                                         &
    & IVE     , IVW
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & IRAD
!
    INTEGER(KIND=I4KIND), DIMENSION(JM)                                                         ::&
    & IHEG    , IHWG    ,                                                                         &
    & IVEG    , IVWG
!
    INTEGER(KIND=I4KIND), DIMENSION(2*IM-1)                                                     ::&
    & IRADG
!
    END MODULE INDX
