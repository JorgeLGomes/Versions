!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION tw(t,td,p)

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

!   this function returns the wet-bulb temperature tw (celsius)
!   given the temperature t (celsius), dew point td (celsius)
!   and pressure p (mb).  see p.13 in stipanuk (1973), referenced
!   above, for a description of the technique.
!
!
!   determine the mixing ratio line thru td and p.

 aw = w(td,p)
!
!   determine the dry adiabat thru t and p.

 ao = o(t,p)
 pi = p

!   iterate to locate pressure pi at the intersection of the two
!   curves .  pi has been set to p for the initial guess.

 DO i= 1,10
   x= .02*(tmr(aw,pi)-tda(ao,pi))
   IF (ABS(x) < 0.01) EXIT
   pi= pi*(2.**(x))
 END DO

!   find the temperature on the dry adiabat ao at pressure pi.

 ti= tda(ao,pi)

!   the intersection has been located...now, find a saturation
!   adiabat thru this point. function os returns the equivalent
!   potential temperature (c) of a parcel saturated at temperature
!   ti and pressure pi.

 aos= os(ti,pi)

!   function tsa returns the wet-bulb temperature (c) of a parcel at
!   pressure p whose equivalent potential temperature is aos.

 tw = tsa(aos,p)
 RETURN
 END FUNCTION tw
