    MODULE ACMCLH
!>--------------------------------------------------------------------------------------------------
!> MODULE ACMCLH
!>
!> USE MODULES: F77KINDS
!>              PARMETA        
!>
!> DRIVER     : CHKOUT
!>              CUCNVC
!>              GOSSIP
!>              GSCOND
!>              GSMDRIVE
!>              INIT
!>              INITS
!>              PRECPD
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NHEAT
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & THEAT   , AVRAIN  , AVCNVC  , ARATIM  , ACUTIM 
! 
    REAL   (KIND=R4KIND), DIMENSION (IDIM1:IDIM2, JDIM1:JDIM2, LM)                              ::&
    & TRAIN   , TCUCN
!
    END MODULE ACMCLH
