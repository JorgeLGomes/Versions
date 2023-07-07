!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
FUNCTION f_esi( p, t )
!-----------------------------------------------------------------------
!  Calculate the saturation water vapor over ice using enhanced
!  Teten's formula.
!
!  Feb 16 2004:
!  Modified by Ernani L. Nascimento to include constants satfia,satfib,
!  sateia,sateib,sateic without INCLUDE command.
!-----------------------------------------------------------------------
 IMPLICIT NONE

 REAL :: p         ! Pressure (Pascal)
 REAL :: t         ! Temperature (K)
 REAL :: f_esi     ! Saturation water vapor pressure over ice (Pa)

 REAL :: f

 REAL :: satfia, satfib
 PARAMETER ( satfia = 1.0003 )
 PARAMETER ( satfib = 4.18E-8 )  ! for p in Pa

 REAL :: sateia, sateib, sateic
 PARAMETER ( sateia = 611.15 )   ! es in Pa
 PARAMETER ( sateib = 22.452 )
 PARAMETER ( sateic = 0.6 )

!  Removed line below. ERNANI L. NASCIMENTO
!  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
 f = satfia + satfib * p
 f_esi = f * sateia * EXP( sateib*(t-273.15)/(t-sateic) )

 RETURN
END FUNCTION f_esi
