!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION tmr(w,p)
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

!   this function returns the temperature (celsius) on a mixing
!   ratio line w (g/kg) at pressure p (mb). the formula is given in
!   table 1 on page 7 of stipanuk (1973).
!
!   initialize constants

 DATA c1/.0498646455/,c2/2.4082965/,c3/7.07475/
 DATA c4/38.9114/,c5/.0915/,c6/1.2035/

 x= ALOG10(w*p/(622.+w))
 tmrk= 10.**(c1*x+c2)-c3+c4*((10.**(c5*x)-c6)**2.)
 tmr= tmrk-273.15
 RETURN
 END FUNCTION tmr
