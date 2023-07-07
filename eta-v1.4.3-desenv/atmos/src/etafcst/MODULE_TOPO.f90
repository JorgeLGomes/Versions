    MODULE TOPO
!>--------------------------------------------------------------------------------------------------
!> MODULE TOPO
!>
!> ABSTRACT:
!> IT WAS CREATED FROM MPP.h
!>
!> USE MODULES: PARMETA
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
    USE PARMETA , ONLY : IM, JM, LM
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(0:IM+1, 0:JM+1)                                             ::&
    & TEMP2X  ,                                                                                   &
    & TTVG 
!
    REAL   (KIND=R4KIND), DIMENSION(0:IM+1, 0:JM+1, LM)                                         ::&
    & HTMG
!
    END MODULE TOPO
