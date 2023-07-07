!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION tda(o,p)
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

!   this function returns the temperature tda (celsius) on a dry adiabat
!   at pressure p (millibars). the dry adiabat is given by
!   potential temperature o (celsius). the computation is based on
!   poisson's equation.

 ok= o+273.15
 tdak= ok*((p*.001)**.286)
 tda= tdak-273.15
 RETURN
 END FUNCTION tda
