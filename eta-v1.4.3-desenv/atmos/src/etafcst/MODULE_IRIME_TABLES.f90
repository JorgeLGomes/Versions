    MODULE IRIME_TABLES
!>--------------------------------------------------------------------------------------------------
!> MODULE IRIME_TABLES
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : GSMCOLUMN
!>              GSMCONST
!>--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
! VEL_RF - VELOCITY INCREASE OF RIMED PARTICLES AS FUNCTIONS OF CRUDE PARTICLE SIZE CATEGORIES 
! (AT 0.1 MM INTERVALS OF MEAN ICE PARTICLE SIZES) AND RIME FACTOR (DIFFERENT VALUES OF RIME FACTOR
! OF 1.1 ** N, WHERE N = 0 TO NRIME).
!--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND), PARAMETER :: NRIME = 40
!
    REAL   (KIND=R4KIND), DIMENSION(2:9, 0:NRIME)                                               ::&
    & VEL_RF
!
    END MODULE IRIME_TABLES
