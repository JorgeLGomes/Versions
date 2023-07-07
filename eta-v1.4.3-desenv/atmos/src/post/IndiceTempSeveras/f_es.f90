!
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION F_ES                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
FUNCTION f_es( p, t )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the saturation specific humidity using enhanced Teten's
!  formula.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  01/08/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    p        Pressure (Pascal)
!    t        Temperature (K)
!
!  OUTPUT:
!
!    f_es     Saturation water vapor pressure (Pa)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE

 REAL :: p         ! Pressure (Pascal)
 REAL :: t         ! Temperature (K)
 REAL :: f_es      ! Saturation water vapor pressure (Pa)
!
!-----------------------------------------------------------------------
!
!  Function f_es and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
 REAL :: f_esl, f_esi

!fpp$ expand (f_esl)
!fpp$ expand (f_esi)
!dir$ inline always f_esl, f_esi
!*$*  inline routine (f_esl, f_esi)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
 IF ( t >= 273.15 ) THEN      ! for water
   f_es = f_esl( p,t )
 ELSE                            ! for ice
   f_es = f_esi( p,t )
 END IF

 RETURN
END FUNCTION f_es
