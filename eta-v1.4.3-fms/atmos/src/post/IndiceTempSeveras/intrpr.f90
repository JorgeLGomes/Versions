!
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INTRPR                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######   and modified by Ernani L. Nascimento at Simepar    ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE intrpr(nlev,p,var,plvl,varatp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate variable "var" linearly in log-pressure (p).
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  Cleaned-up, removed OLAPS artifacts.
!
!  Feb 13 2004  Ernani L. Nascimento
!  Removed variable maxlev
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
 INTEGER :: nlev
 REAL :: p(nlev+2),var(nlev+2)
 REAL :: plvl
 REAL :: varatp
!
 INTEGER :: k
 REAL :: w1
!
 DO k=2,nlev
   IF(p(k) < plvl) EXIT
 END DO
!
 w1=ALOG(p(k)/plvl)/ALOG(p(k)/p(k-1))
 varatp = w1*var(k-1) + (1.-w1)*var(k)
!
 RETURN
END SUBROUTINE intrpr
