    MODULE SLOPES
!>--------------------------------------------------------------------------------------------------
!> MODULE SLOPES
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!> 
!> DRIVER     : DIVHOAST
!>              HDIFFS
!>              HZADVS
!>              INITS
!>              READ_NHB
!>              SLADVT
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : IDIM1, IDIM2, JDIM1, JDIM2, LM
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & ISLD
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & VTMS 
!
    END MODULE SLOPES
