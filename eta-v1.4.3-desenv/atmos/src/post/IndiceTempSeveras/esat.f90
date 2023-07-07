!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION esat(t)
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
!
!   this function returns the saturation vapor pressure over
!   water (mb) given the temperature (celsius).
!   the algorithm is due to nordquist, w.s.,1973: "numerical approxima-
!   tions of selected meteorlolgical parameters for cloud physics prob-
!   lems," ecom-5475, atmospheric sciences laboratory, u.s. army
!   electronics command, white sands missile range, new mexico 88002.

 tk = t+273.15
 p1 = 11.344-0.0303998*tk
 p2 = 3.49149-1302.8844/tk
 c1 = 23.832241-5.02808*ALOG10(tk)
 esat = 10.**(c1-1.3816E-7*10.**p1+8.1328E-3*10.**p2-2949.076/tk)
 RETURN
 END FUNCTION esat
