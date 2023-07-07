!
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE THWXN                     ######
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
SUBROUTINE thwxn(nlevel,p,ht,t,td,elev,twdf)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
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
!  Feb/14/2004 Ernani L. Nascimento: removed variable maxlev
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
!
!  Input variables
!
 INTEGER :: nlevel
 REAL :: p(nlevel+2),ht(nlevel+2),t(nlevel+2),td(nlevel+2)
 REAL :: elev
!
!  Output variables
!
 REAL :: twdf
!
!  Functions
!
 REAL :: oe,tsa_fast
!
!  Misc internal variables
!
 INTEGER :: n
 REAL :: h3km,thaec,thw,twx,twn
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
 twx=-999.
 twn=999.
 h3km=elev+3000.
 DO n=1,nlevel
   IF(ht(n) >= elev) THEN
     IF(ht(n) > h3km) EXIT
     thaec=oe(t(n),td(n),p(n))
     thw=tsa_fast(thaec,1000.)
     twx=AMAX1(twx,thw)
     twn=AMIN1(twn,thw)
   END IF
 END DO
!  101 CONTINUE
!
!  Find difference between max and min
!
 twdf=twx-twn
 RETURN
END SUBROUTINE thwxn
