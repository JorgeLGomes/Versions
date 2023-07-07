!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION oe(t,td,p)
!
!    g.s. stipanuk     1973          original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982

!   this function returns equivalent potential temperature oe (celsius)
!   of a parcel of air given its temperature t (celsius), dew point
!   td (celsius) and pressure p (millibars).
!   find the wet bulb temperature of the parcel.

 atw = tw(t,td,p)

!   find the equivalent potential temperature.

 oe = os(atw,p)
 RETURN
 END FUNCTION oe
