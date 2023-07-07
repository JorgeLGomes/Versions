!
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CALCLI                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######  and modified by Ernani L. Nascimento at Simepar     ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE calcli(nlevel,thepcl,p,t,li)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate Convective Available Potential Energy (CAPE)
!  in each column of ARPS grid.
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
!  Fev/14/2004  Ernani L. Nascimento: removed variable maxlev
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
!
!  Input variables
!
 INTEGER :: nlevel
 REAL :: thepcl
 REAL :: p(nlevel+2),t(nlevel+2)
!
!  Output variable
!
 REAL :: li
!
!  Functions
!
 REAL :: tsa_fast
!
!  Misc internal variables
!
 INTEGER :: n
 REAL :: dp,wlow,t500,par500
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
 DO n=2,nlevel-1
   IF( p(n) <= 500.) EXIT
 END DO
!  101 CONTINUE
 dp=ALOG(p(n-1)/p(n))
 wlow=ALOG(500./p(n))/dp
 t500=t(n)*(1.-wlow) + t(n-1)*wlow
 par500=tsa_fast(thepcl,500.)
!
 li=t500-par500
!
 RETURN
END SUBROUTINE calcli
