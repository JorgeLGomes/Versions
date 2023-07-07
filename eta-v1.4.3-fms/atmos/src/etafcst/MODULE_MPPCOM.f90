    MODULE MPPCOM
!>--------------------------------------------------------------------------------------------------
!> MODULE MPPCOM
!>
!> ABSTRACT:
!> IT WAS CREATED FROM MPP.h
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : BOCOH
!>              BOCOHF
!>              BOCOV
!>              CHKOUT
!>              CLO89
!>              CLTEND
!>              COLLECT
!>              CONRAD
!>              CUCNVC
!>              DDAMP
!>              DIFCOF
!>              DIGFLT
!>              DIST
!>              DIVHOA
!>              DIVHOAST
!>              DSTRB
!>              E1E290
!>              E290
!>              E2SPEC
!>              E3V88
!>              EBU
!>              EPS
!>              FILT25
!>              FST88
!>              GOSSIP
!>              GSCOND
!>              GSMCOLUMN
!>              GSMDRIVE
!>              HADZ
!>              HDIFF
!>              HDIFFS
!>              HZADV
!>              HZADVS
!>              HZADV2
!>              HZADV_LM1
!>              IDSTRB
!>              INIT
!>              INITS
!>              LOC2GLB
!>              LWR88
!>              MIXLEN
!>              MODULE_EXCHM
!>              MPI_FIRST
!>              MPPINIT
!>              NEWFLT
!>              OZON2D
!>              PDNEW
!>              PDTE
!>              PDTEDT
!>              PGCOR
!>              PRECPD
!>              QUILT
!>              RADFS
!>              RADTN
!>              RDTEMP
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>              SFCDIF
!>              SFCDIF_BHS
!>              SLP
!>              SLPSIG
!>              SLPSIGSPLINE
!>              SPA88
!>              SURFCE
!>              SWR93
!>              TABLE
!>              TTBLEX
!>              TURBL
!>              TWR
!>              UPDATE
!>              VADZ
!>              VDIFH
!>              VDIFV
!>              VTADV
!>              VTADVF
!>              VWR
!>              ZENITH
!>              ZERO2
!>              ZERO3
!>              ZERO3_T 
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : INPES, JNPES
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & MYPE              , NPES              ,                                                     &
    & MY_IS_GLB         , MY_IE_GLB         ,                                                     &
    & MY_JS_GLB         , MY_JE_GLB         ,                                                     &
    & MY_IS_LOC         , MY_IE_LOC         ,                                                     &
    & MY_JS_LOC         , MY_JE_LOC         ,                                                     &
    & MYIS              , MYIS1             ,                                                     &
    & MYIS2             , MYIE              ,                                                     &
    & MYIE1             , MYIE2             ,                                                     &
    & MYIS_P1           , MYIS_P2           ,                                                     &
    & MYIS_P3           , MYIS_P4           ,                                                     &
    & MYIS_P5           ,                                                                         &
    & MYIS1_P1          , MYIS1_P2          ,                                                     &
    & MYIS1_P3          , MYIS1_P4          ,                                                     &
    & MYIE_P1           , MYIE_P2           ,                                                     &
    & MYIE_P3           , MYIE_P4           ,                                                     &
    & MYIE_P5           ,                                                                         &
    & MYIE1_P1          , MYIE1_P2          ,                                                     &
    & MYIE1_P3          , MYIE1_P4          ,                                                     &
    & MYIE2_P1          ,                                                                         &
    & MYJS              , MYJS1             ,                                                     &
    & MYJS2             , MYJS3             ,                                                     &
    & MYJS4             , MYJS5             ,                                                     &
    & MYJS_P1           , MYJS_P2           ,                                                     &
    & MYJS_P3           , MYJS_P4           ,                                                     &
    & MYJS_P5           ,                                                                         &
    & MYJS1_P1          , MYJS1_P2          ,                                                     &
    & MYJS1_P3          , MYJS1_P4          ,                                                     &
    & MYJS2_P1          , MYJS2_P2          ,                                                     &
    & MYJS2_P3          , MYJS2_P4          ,                                                     &
    & MYJS3_P4          ,                                                                         &
    & MYJS4_P1          , MYJS4_P4          ,                                                     &
    & MYJS5_P1          , MYJS5_P2          ,                                                     &
    & MYJE              , MYJE1             ,                                                     &
    & MYJE2             , MYJE3             ,                                                     &
    & MYJE4             , MYJE5             ,                                                     &
    & MYJE_P1           , MYJE_P2           ,                                                     &
    & MYJE_P3           , MYJE_P4           ,                                                     &
    & MYJE_P5           ,                                                                         &
    & MYJE1_P1          , MYJE1_P2          ,                                                     &
    & MYJE1_P3          , MYJE1_P4          ,                                                     &
    & MYJE2_P1          , MYJE2_P2          ,                                                     &
    & MYJE2_P3          , MYJE2_P4          ,                                                     &
    & MYJE3_P4          ,                                                                         &
    & MYJE4_P1          , MYJE4_P4          ,                                                     &
    & MYJE5_P1          , MYJE5_P2          ,                                                     &
    & MY_N              , MY_E              ,                                                     &
    & MY_S              , MY_W              ,                                                     &
    & MY_NE             , MY_SE             ,                                                     &
    & MY_SW             , MY_NW
!
    INTEGER(KIND=I4KIND), DIMENSION(8)                                                          ::&
    & MY_NEB
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ILCOL   , IRCOL   ,                                                                         &
    & IBROW   , ITROW   ,                                                                         &
    & ILPAD1  , ILPAD2  ,                                                                         &
    & ILPAD3  , ILPAD4  ,                                                                         &
    & ILPAD5  ,                                                                                   &
    & IRPAD1  , IRPAD2  ,                                                                         &
    & IRPAD3  , IRPAD4  ,                                                                         &
    & IRPAD5  ,                                                                                   &
    & JBPAD1  , JBPAD2  ,                                                                         &
    & JBPAD3  , JBPAD4  ,                                                                         &
    & JBPAD5  ,                                                                                   &
    & JTPAD1  , JTPAD2  ,                                                                         &
    & JTPAD3  , JTPAD4  ,                                                                         &
    & JTPAD5
!
    INTEGER(KIND=I4KIND), DIMENSION(0:INPES*JNPES)                                              ::&
    & IS_LOC_TABLE      , JS_LOC_TABLE      ,                                                     &
    & IE_LOC_TABLE      , JE_LOC_TABLE      ,                                                     &
    & ICHUNKTAB
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & MPI_COMM_COMP     , MPI_COMM_INTER    ,                                                     &
    & IQUILT_GROUP 
!
    INTEGER(KIND=I4KIND), DIMENSION(100)                                                        ::&
    & MPI_COMM_INTER_ARRAY        , INUMQ
!
    END MODULE MPPCOM





