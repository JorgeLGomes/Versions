    MODULE PHYS
!>--------------------------------------------------------------------------------------------------
!> MODULE PHYS
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>              PARM_TBL
!>              RDPARM
!>
!> DRIVER     : CHKOUT
!>              GOSSIP
!>              INIT
!>              INITS
!>              READ_NHB
!>              TURBL 
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : LM ,IDIM1, IDIM2, JDIM1, JDIM2
    USE PARM_TBL, ONLY : ITB, JTB , ITBQ , JTBQ
    USE RDPARM  , ONLY : LP1
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & KTM
!
    REAL   (KIND=R4KIND)                                                                        ::&    
    & DTQ2    , TDTQ2   , DTD     , TDTD    ,                                                     &
    & ROS     , CS      , DS      , ROI     , CI      , DI      ,                                 &
    & PL      , THL     , RDQ     , RDTH    , RDP     , RDTHE   ,                                 &
    & PLQ     , RDPQ    , RDTHEQ
!
    REAL   (KIND=R4KIND), DIMENSION(LP1)                                                        ::&    
    & DFRLG
!
    REAL   (KIND=R4KIND), DIMENSION(JTB)                                                        ::&
    & QS0     , SQS
!
    REAL   (KIND=R4KIND), DIMENSION(ITB)                                                        ::&
    & THE0    , STHE
!
    REAL   (KIND=R4KIND), DIMENSION(ITBQ)                                                       ::&
    & THE0Q   , STHEQ
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & MXSNAL  , EPSR    ,                                                                         &
    & RADIN   , RADOT   ,                                                                         &
    & GLAT    , GLON    ,                                                                         &
    & CZEN    ,                                                                                   &
    & HTOP    , HBOT    ,                                                                         &
    & CNVTOP  , CNVBOT  ,                                                                         &
    & TG      , GFFC    ,                                                                         &
    & SST     , ALBASE  ,                                                                         &
    & ALBEDO  ,                                                                                   &
    & HDAC    , HDACV   ,                                                                         &
    & CZMEAN  , SIGT4
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, 565)                              ::&
    & SSTM
!
   REAL   (KIND=R4KIND), DIMENSION(ITB, JTB)                                                   ::&
    & PTBL
!
    REAL   (KIND=R4KIND), DIMENSION(JTB, ITB)                                                   ::&
    & TTBL
!
    REAL   (KIND=R4KIND), DIMENSION(JTBQ, ITBQ)                                                 ::&
    & TTBLQ
!
    END MODULE PHYS
