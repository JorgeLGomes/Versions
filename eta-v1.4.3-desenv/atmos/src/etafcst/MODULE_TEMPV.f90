    MODULE TEMPV
!>--------------------------------------------------------------------------------------------------
!> MODULE TEMPV
!>
!> CALLS : F77KINDS
!>         PARMETA
!>
!> DRIVER: CHKOUT
!>         GSCOND
!>         INIT
!>         INITS
!>-------------------------------------------------------------------------------------------------- 
    USE F77KINDS
    USE PARMETA , ONLY : IDIM1, IDIM2, JDIM1, JDIM2, LM
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & P0
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & T0      , Q0
!
    END MODULE TEMPV
