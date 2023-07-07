    MODULE SEASO3
!>--------------------------------------------------------------------------------------------------
!> MODULE SEASO3
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : O3CLIM
!>              OZON2D
!>--------------------------------------------------------------------------------------------------  
    USE F77KINDS
!
    IMPLICIT NONE
!
    SAVE
!----------------------- 
! XDUO3N(37,NL) - WINTER
! XDO3N2(37,NL) - SPRING
! XDO3N3(37,NL) - SUMMER
! XDO3N4(37,NL) - FALL
!----------------------- 
    INTEGER(KIND=I4KIND), PARAMETER :: NL = 81
    INTEGER(KIND=I4KIND), PARAMETER :: LNGTH = 37 * NL
!
    REAL   (KIND=R4KIND), DIMENSION(37, NL)                                                     ::&
    & XDUO3N  , XDO3N2  , XDO3N3  , XDO3N4
!
    REAL   (KIND=R4KIND), DIMENSION(NL)                                                         ::&
    & PRGFDL
!
    REAL   (KIND=R4KIND), DIMENSION(37,NL,4)                                                    ::&
    & O3O3
!
    REAL   (KIND=R4KIND), DIMENSION(LNGTH)                                                      ::&
    & XRAD1   , XRAD2   , XRAD3   , XRAD4 
!         
    EQUIVALENCE (XRAD1(1),XDUO3N(1,1),O3O3(1,1,1)),                                               &
    &           (XRAD2(1),XDO3N2(1,1),O3O3(1,1,2)),                                               &
    &           (XRAD3(1),XDO3N3(1,1),O3O3(1,1,3)),                                               &
    &           (XRAD4(1),XDO3N4(1,1),O3O3(1,1,4))
!
    END MODULE SEASO3
