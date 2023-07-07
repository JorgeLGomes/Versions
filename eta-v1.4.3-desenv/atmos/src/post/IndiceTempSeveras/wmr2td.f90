!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION WMR2TD                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
 FUNCTION wmr2td(pres,wmr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate Convective Available Potential Energy (CAPE)
!  in each column of ARPS grid.  => WRONG!!
!
!  CONVERTS WATER VAPOR MIXING RATIO INTO DEW POINT TEMPERATURE
!  (Ernani L. Nascimento)
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on GEMPAK routine of same name.
!
!  MODIFICATION HISTORY:
!
!
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
 REAL :: wmr2td
 REAL :: pres
 REAL :: wmr
 REAL :: wkgkg,e,evap
!
 wkgkg = 0.001 * wmr
 wkgkg = AMAX1(wmr,0.00005)
 e= (pres*wkgkg) / (0.62197 + wkgkg)
 evap = e /(1.001 + (( pres - 100.) /900.) * 0.0034)
 wmr2td = ALOG(evap/6.112) * 243.5 /( 17.67 - ALOG (evap/6.112))

 RETURN
 END FUNCTION wmr2td
