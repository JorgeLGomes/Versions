    MODULE C_FRACN
!>--------------------------------------------------------------------------------------------------
!> MODULE C_FRACN
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : GSMDRIVE
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(LM, IDIM1:IDIM2, JDIM1:JDIM2)                               ::&
    & F_ICE   ,                                                                                   &
    & F_RAIN  ,                                                                                   &
    & F_RIMEF
! 
    END MODULE C_FRACN
