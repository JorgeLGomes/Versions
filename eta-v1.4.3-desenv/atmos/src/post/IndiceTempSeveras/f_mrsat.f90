!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION F_MRSAT                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
FUNCTION f_mrsat( p, t )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the saturation water vapor mixing ratio using enhanced
!  Teten's formula.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  01/08/1998
!
!  MODIFICATION HISTORY:
!
!  16 FEB 2004
!  Ernani L. Nascimento: modified to include constant rddrv without
!  INCLUDE command.
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
!    f_mrsat  Saturation water vapor mixing ratio (kg/kg).
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
 REAL :: f_mrsat   ! Saturation water vapor mixing ratio (kg/kg)

 REAL :: rd        ! Gas constant for dry air  (m**2/(s**2*K))
 PARAMETER( rd     = 287.0 )

 REAL :: rv        ! Gas constant for water vapor  (m**2/(s**2*K)).
 PARAMETER( rv     = 461.0 )

 REAL :: rddrv
 PARAMETER( rddrv  = rd/rv )
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
 REAL :: fes
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
! Removed line below. ERNANI L. NASCIMENTO
!  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Function f_es and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
 REAL :: f_es
!fpp$ expand (f_es)
!dir$ inline always f_es
!*$*  inline routine (f_es)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
 fes = f_es( p,t )
 f_mrsat = rddrv * fes / (p-fes)

 RETURN
END FUNCTION f_mrsat
