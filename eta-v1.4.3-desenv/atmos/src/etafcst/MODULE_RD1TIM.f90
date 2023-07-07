    MODULE RD1TIM
!>--------------------------------------------------------------------------------------------------
!> MODULE RD1TIM
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : INIT
!>              INITS
!>              RADTN
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : IDIM1, IDIM2, JDIM1, JDIM2
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & K400    , RAD1
!
    REAL   (KIND=R4KIND), DIMENSION(3)                                                          ::&
    & CTHK    , TAUCV
!
    REAL   (KIND=R4KIND), DIMENSION(4)                                                          ::&
    & PTOPC
!
    INTEGER(KIND=R4KIND), DIMENSION(3)                                                          ::&
    & LTOP
!
    INTEGER(KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & LVL
!
    DATA                                                                                          &
    & CTHK  /20000.0, 20000.0, 20000.0/
!
    DATA                                                                                          &
    & TAUCV /0.16   , 0.14   , 0.12   /
!
    DATA                                                                                          &
    & LTOP  /0      , 0      , 0      /
!
    END MODULE RD1TIM
