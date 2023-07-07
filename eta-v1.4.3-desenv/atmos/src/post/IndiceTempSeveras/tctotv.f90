!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION TCTOTV                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
 FUNCTION tctotv(tt,ww)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Virtual Temperature
!
!  Given T in Celcius and mixing ratio in g/kg
!  find the virtual temperature.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
!
 REAL :: tctotv,tt,ww
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
 tctotv=(tt+273.15)*(1.+0.0006*ww)
 RETURN
 END FUNCTION tctotv
