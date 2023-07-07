    MODULE TIMCHK
!>--------------------------------------------------------------------------------------------------
!> MODULE TIMCHK
!>
!> USE MODULES: F77KINDS
!> 
!> DRIVER     : CHKOUT
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
! 
    IMPLICIT NONE
!
    SAVE
! 
    REAL   (KIND=R4KIND)                                                                        ::&
    & SLP_TIM , GATH_TIM, WRT_TIM , PROF_TIM, BCEX_TIM, STAT_TIM
!
    END MODULE TIMCHK
