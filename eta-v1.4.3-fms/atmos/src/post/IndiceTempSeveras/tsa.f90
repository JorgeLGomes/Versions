!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION tsa(os,p)
!
!    g.s. stipanuk     1973           original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982

!   this function returns the temperature tsa (celsius) on a saturation
!   adiabat at pressure p (millibars). os is the equivalent potential
!   temperature of the parcel (celsius). sign(a,b) replaces the
!   algebraic sign of a with that of b.
!   b is an empirical constant approximately equal to 0.001 of the latent
!   heat of vaporization for water divided by the specific heat at constant
!   pressure for dry air.

 DATA b/2.6518986/
 a= os+273.15

!   tq is the first guess for tsa.

 tq= 253.15

!   d is an initial value used in the iteration below.

 d= 120.

!   iterate to obtain sufficient accuracy....see table 1, p.8
!   of stipanuk (1973) for equation used in iteration.

 DO i= 1,12
   tqk= tq-273.15
   d= d/2.
   x= a*EXP(-b*w(tqk,p)/tq)-tq*((1000./p)**.286)
   IF (ABS(x) < 1E-7) GOTO 2
   tq= tq+SIGN(d,x)
 END DO
2 tsa= tq-273.15
 RETURN
 END FUNCTION tsa
