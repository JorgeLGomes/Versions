!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
FUNCTION f_esl( p, t )
!-----------------------------------------------------------------------
!  Calculate the saturation water vapor over liquid water using
!  enhanced Teten's formula.
!
!  Feb 16 2004:
!  Modified by Ernani L. Nascimento to include constants satfwa,satfwb,
!  satewa,satewb,satewc without INCLUDE command.
!-----------------------------------------------------------------------
 IMPLICIT NONE

 REAL :: p         ! Pressure (Pascal)
 REAL :: t         ! Temperature (K)
 REAL :: f_esl     ! Saturation water vapor pressure over liquid water

 REAL :: f

 REAL :: satfwa, satfwb
 PARAMETER ( satfwa = 1.0007 )
 PARAMETER ( satfwb = 3.46E-8 )  ! for p in Pa

 REAL :: satewa, satewb, satewc
 PARAMETER ( satewa = 611.21 )   ! es in Pa
 PARAMETER ( satewb = 17.502 )
 PARAMETER ( satewc = 32.18 )

!  Removed line below. ERNANI L. NASCIMENTO
!  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
 f = satfwa + satfwb * p
 f_esl = f * satewa * EXP( satewb*(t-273.15)/(t-satewc) )

 RETURN
END FUNCTION f_esl
