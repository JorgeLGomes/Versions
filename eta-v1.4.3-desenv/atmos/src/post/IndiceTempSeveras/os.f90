!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION os(t,p)
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

!   this function returns the equivalent potential temperature os
!   (celsius) for a parcel of air saturated at temperature t (celsius)
!   and pressure p (millibars).
 DATA b/2.6518986/
!   b is an empirical constant approximately equal to the latent heat
!   of vaporization for water divided by the specific heat at constant
!   pressure for dry air.

 tk = t+273.15
 osk= tk*((1000./p)**.286)*(EXP(b*w(t,p)/tk))
 os= osk-273.15
 RETURN
 END FUNCTION os
