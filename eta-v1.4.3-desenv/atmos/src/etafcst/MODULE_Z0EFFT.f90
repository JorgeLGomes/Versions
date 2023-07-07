    MODULE Z0EFFT
!>-------------------------------------------------------------------------------------------------- 
!> MODULE Z0EFFT
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : GOSSIP
!>              INIT
!>              INITS
!>              TURBL
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, 4)                                ::&
    & ZEFFIJ
!
    END MODULE Z0EFFT

