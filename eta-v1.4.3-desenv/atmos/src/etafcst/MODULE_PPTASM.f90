    MODULE PPTASM
!>--------------------------------------------------------------------------------------------------
!> MODULE PPTASM
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!> 
!> DRIVER     : CUCNVC
!>              GSMDRIVE
!>              INIT
!>              INITS
!>              PRECPD
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : IDIM1, IDIM2, JDIM1, JDIM2, LM
!
    IMPLICIT NONE
!
    SAVE
!--------------------------------------------------------------------------------
! ITSTLOC, JTSTLOC, MTSTPE: LOCAL(ITEST, JTEST) POINT, AND THE NODE IT BELONGS TO
!--------------------------------------------------------------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & MTSTPE  , ITSTLOC , JTSTLOC 
!-----------------------------------------------------------------------------------------------
! APREC - GRID-SCALE PRECIP, CALCULATED IN PRECPD. DID NOT SEEM TO SERVE ANY PARTICULAR PURPOSE.
!-----------------------------------------------------------------------------------------------
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & PHOUR   ,                                                                                   &
    & APREC
!-----------------------------------------------------------------------------
! TLATCU: CUCNVC LATENT HEAT
! TLATGS: GRID-SCALE LATENT HEAT, HEAT WHICH IS KEPT TRACK OF BY BRAD IN TMOD.
!-----------------------------------------------------------------------------
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & TLATCU  ,                                                                                   &
    & TLATGS
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, 3)                                ::&
    & PPTDAT
!
    END MODULE PPTASM
    

