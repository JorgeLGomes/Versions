    MODULE GLB_TABLE
!>--------------------------------------------------------------------------------------------------
!> MODULE GLB_TABLE 
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
    INTEGER(KIND=I4KIND), DIMENSION(0:INPES*JNPES)                                              ::&
    & IS_GLB_TABLE      , IE_GLB_TABLE      ,                                                     &
    & JS_GLB_TABLE      , JE_GLB_TABLE
!
    END MODULE GLB_TABLE 
