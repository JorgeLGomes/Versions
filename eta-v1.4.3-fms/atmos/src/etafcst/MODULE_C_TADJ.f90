    MODULE C_TADJ
!>--------------------------------------------------------------------------------------------------
!> MODULE C_TADJ
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : CLTEND
!>              INIT
!>              INITS
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : IDIM1, IDIM2, JDIM1, JDIM2, LM
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & T_ADJ , T_OLD
!
    END MODULE C_TADJ
