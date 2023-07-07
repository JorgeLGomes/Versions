    MODULE SOIL
!>--------------------------------------------------------------------------------------------------
!> MODULE SOIL
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>              PARMSOIL
!>
!> DRIVER     : CHKOUT
!>              EXIT
!>              GOSSIP
!>              INIT
!>              INITS
!>              RADTN
!>              READ_NHB
!>              READ_RESTRT
!>              READ_RESTRT2
!>              SURFCE
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : IDIM1, IDIM2, JDIM1, JDIM2
    USE PARMSOIL, ONLY : NSOIL 
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & IVGTYP  , ISLTYP  ,                                                                         &
    & ISLOPE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & SOILTB  , SFCEXC  ,                                                                         &
    & SMSTAV  , SMSTOT  ,                                                                         &
    & GRNFLX  , PCTSNO  ,                                                                         &
    & VEGFRC  , CMC
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, NSOIL)                            ::&
    & SMC     ,                                                                                   &
    & STC     ,                                                                                   &
    & SH2O
!
    REAL   (KIND=R4KIND), DIMENSION(NSOIL)                                                      ::&
    & SLDPTH  , RTDPTH
!
    END MODULE SOIL


