!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION w(t,p)
!
!    g.s. stipanuk     1973              original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982
!
!  this function returns the mixing ratio (grams of water vapor per
!  kilogram of dry air) given the dew point (celsius) and pressure
!  (millibars). if the temperture  is input instead of the
!  dew point, then saturation mixing ratio (same units) is returned.
!  the formula is found in most meteorological texts.

 x= esat(t)
 w= 622.*x/(p-x)
 RETURN
 END FUNCTION w
