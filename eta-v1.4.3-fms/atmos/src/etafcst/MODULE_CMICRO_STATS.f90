    MODULE CMICRO_STATS
!>--------------------------------------------------------------------------------------------------
!> MODULE CMICRO_STATS
!>
!> USE MODULES: F77KINDS
!> 
!> DRIVER     : GSMCOLUMN
!>              GSMDRIVE
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    SAVE
!--------------------------------------------------------- 
! THE FOLLOWING VARIABLES ARE FOR MICROPHYSICAL STATISTICS
!--------------------------------------------------------- 
    INTEGER(KIND=I4KIND), PARAMETER :: ITLO = -60
    INTEGER(KIND=I4KIND), PARAMETER :: ITHI =  40
!
    INTEGER(KIND=I4KIND), DIMENSION(ITLO:ITHI,  4)                                              ::&
    & NSTATS
!
    REAL   (KIND=R4KIND), DIMENSION(ITLO:ITHI,  5)                                              ::&
    & QMAX
!
    REAL   (KIND=R4KIND), DIMENSION(ITLO:ITHI, 22)                                              ::&
    & QTOT
!
    END MODULE CMICRO_STATS
