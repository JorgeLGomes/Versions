!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE ptlcl(p,t,td,pc,tc)
!
!   this subroutine estimates the pressure pc (mb) and the temperature
!   tc (celsius) at the lifted condensation level (lcl), given the
!   initial pressure p (mb), temperature t (celsius) and dew point
!   (celsius) of the parcel.  the approximation is that lines of
!   constant potential temperature and constant mixing ratio are
!   straight on the skew t/log p chart.
!
!    baker,schlatter   17-may-1982   original version
!
!   teten's formula for saturation vapor pressure as a function of
!   pressure was used in the derivation of the formula below.  for
!   additional details, see math notes by t. schlatter dated 8 sep 81.
!   t. schlatter, noaa/erl/profs program office, boulder, colorado,
!   wrote this subroutine.
!
!   akap = (gas constant for dry air) / (specific heat at constant
!       pressure for dry air)
!   cta = difference between kelvin and celsius temperatures
!
 DATA akap,cta/0.28541,273.16/
 c1 = 4098.026/(td+237.3)**2
 c2 = 1./(akap*(t+cta))
 pc = p*EXP(c1*c2*(t-td)/(c2-c1))
 tc = t+c1*(t-td)/(c2-c1)
 RETURN
END SUBROUTINE ptlcl
