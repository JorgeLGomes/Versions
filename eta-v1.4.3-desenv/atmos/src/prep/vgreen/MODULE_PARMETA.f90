    MODULE PARMETA
!>--------------------------------------------------------------------------------------------------
!> MODULE PARMETA
!>
!> USE MODULES: F77KINDS
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!------------------------------------------------------- 
! SET PRIMARY GRID DIMENSIONS AND PRESSURE OUTPUT LEVELS
!------------------------------------------------------- 
    INTEGER(KIND=I4KIND), PARAMETER :: IM    = 285
    INTEGER(KIND=I4KIND), PARAMETER :: JM    = 601
    REAL   (KIND=R4KIND), PARAMETER :: DLMD  = .14423076925000000000 
    REAL   (KIND=R4KIND), PARAMETER :: DPHD  = .13461538475000000000 
    REAL   (KIND=R4KIND), PARAMETER :: TLM0D = -65.0
    REAL   (KIND=R4KIND), PARAMETER :: TPH0D = -12.0
!
    END MODULE PARMETA


