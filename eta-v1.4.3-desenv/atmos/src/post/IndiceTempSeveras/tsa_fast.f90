!
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                FUNCTION TSA_FAST                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
FUNCTION tsa_fast(os,p)
!
!   THIS FUNCTION RETURNS THE TEMPERATURE TSA (CELSIUS) ON A SATURATION
!   ADIABAT AT PRESSURE P (MILLIBARS). OS IS THE EQUIVALENT POTENTIAL
!   TEMPERATURE OF THE PARCEL (CELSIUS). SIGN(A,B) REPLACES THE
!   ALGEBRAIC SIGN OF A WITH THAT OF B.
!
!    BAKER,SCHLATTER 17-MAY-1982     Original version
!    Modification for better convergence, Keith Brewster, Feb 1994.
!
!   B IS AN EMPIRICAL CONSTANT APPROXIMATELY EQUAL TO THE LATENT HEAT
!   OF VAPORIZATION FOR WATER DIVIDED BY THE SPECIFIC HEAT AT CONSTANT
!   PRESSURE FOR DRY AIR.
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate Convective Available Potential Energy (CAPE)
!  in each column of ARPS grid.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  Cleaned-up, removed OLAPS artifacts.
!
!-----------------------------------------------------------------------
!
 REAL :: b
!  PARAMETER (B=2.6518986)
 PARAMETER (b=2651.8986)
 a= os+273.15
!
!   Above 200 mb figure all the moisture is wrung-out, so
!   the temperature is that which has potential temp of theta-e.
!   Otherwise iterate to find combo of moisture and temp corresponding
!   to thetae.
!
 IF( p < 200.) THEN
   tq=a*((p/1000.)**.286)
 ELSE
!   D IS AN INITIAL VALUE USED IN THE ITERATION BELOW.
   d= 120.
!   TQ IS THE FIRST GUESS FOR TSA.
   tq= 253.15
   x = 0.
!
!   ITERATE TO OBTAIN SUFFICIENT ACCURACY....SEE TABLE 1, P.8
!   OF STIPANUK (1973) FOR EQUATION USED IN ITERATION.
   DO i= 1,25
     d= 0.5*d

     x_last = x

     x= a*EXP(-b*f_mrsat(p*100.,tq)/tq)-tq*((1000./p)**.286)

     IF (ABS(x) < 1E-3) GO TO 2
!
     IF (x_last * x < 0.) THEN
       slope = (x-x_last) / (tq - tq_last)
       delta = - x / slope
       ad = AMIN1(ABS(delta),d)
       tq_last = tq
       tq = tq + SIGN(ad,delta)
     ELSE
       tq_last = tq
       tq= tq+SIGN(d,x)
     END IF

   END DO
 END IF
2 tsa_fast = tq-273.15
 RETURN
END FUNCTION tsa_fast
