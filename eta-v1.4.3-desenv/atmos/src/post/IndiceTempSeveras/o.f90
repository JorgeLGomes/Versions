!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION o(t,p)
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

!   this function returns potential temperature (celsius) given
!   temperature t (celsius) and pressure p (mb) by solving the poisson
!   equation.

 tk= t+273.15
 ok= tk*((1000./p)**.286)
 o= ok-273.15
 RETURN
 END FUNCTION o
