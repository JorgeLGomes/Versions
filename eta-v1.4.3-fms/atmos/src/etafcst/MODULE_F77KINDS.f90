    MODULE F77KINDS
!>--------------------------------------------------------------------------------------------------
!> MODULE F77KINDS
!>
!> USE MODULES: ----- 
!>        
!> DRIVER: ALL SUBROUTINES AND MODULES
!>--------------------------------------------------------------------------------------------------
    IMPLICIT NONE
!
    SAVE
!
    INTEGER, PARAMETER :: I1KIND  =  1  ! INTEGER*1	  
    INTEGER, PARAMETER :: I2KIND  =  2  ! INTEGER*2
    INTEGER, PARAMETER :: I4KIND  =  4  ! INTEGER*4
    INTEGER, PARAMETER :: I8KIND  =  8  ! INTEGER*8
!
    INTEGER, PARAMETER :: L1KIND  =  1  ! LOGICAl*1
    INTEGER, PARAMETER :: L4KIND  =  4  ! LOGICAL*4
    INTEGER, PARAMETER :: L8KIND  =  8  ! LOGICAL*8
!
    INTEGER, PARAMETER :: R4KIND  =  4  ! REAL*4
    INTEGER, PARAMETER :: R8KIND  =  8  ! REAL*8
!
    INTEGER, PARAMETER :: CX8KIND =  8  ! COMPLEX*8
    INTEGER, PARAMETER :: CX16KIND=  16 ! COMPLEX*16
!
  END MODULE F77KINDS
