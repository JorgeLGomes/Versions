!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SINDEX                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######    and modified by Ernank L. Nascimento at Simepar   ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE sindex(nlevel,p,ht,t,tv,td,w,                                &
          partem,buoy,wload,mbuoy,pbesnd,mbesnd,                       &
          lcl_pbe,lfc,el,twdf,li,cape,mcape,cin,tcap)
!
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
!  5/13/1996  Added cap strength.
!  Feb 13 2004 ERNANI L. NASCIMENTO: removed variable maxlev
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
!
 INTEGER :: nlevel
!
 REAL :: p(nlevel+2),ht(nlevel+2),t(nlevel+2),tv(nlevel+2),td(nlevel+2),w(nlevel+2)
 REAL :: partem(nlevel+2),buoy(nlevel+2),wload(nlevel+2),mbuoy(nlevel+2)
 REAL :: pbesnd(nlevel+2),mbesnd(nlevel+2)
!
!  Returned from sindex
!
 REAL :: lfc,el,twdf,li,cape,mcape,cin,mcin,tcap
!
!  Potbe variables
!
 REAL :: plcl_pbe,tlcl_pbe,lcl_pbe,thepcl
 REAL :: velneg,mvelneg
!
!  Functions
!
 REAL :: oe
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
 thepcl=oe(t(1),td(1),p(1))
!
!  print *, ' theta-e of parcel: ',thepcl
!
 CALL ptlcl(p(1),t(1),td(1),plcl_pbe,tlcl_pbe)
!  print *, ' press and temp at LCL: ',plcl_pbe,tlcl_pbe
!
!  Find height of LCL
!

 CALL intrpr(nlevel,p,ht,plcl_pbe,lcl_pbe)
!  print *, ' NCL: ', lcl_pbe
!
!  Calculate the CAPE and such
!
 CALL potbe(nlevel,p(1),t(1),w(1),                                     &
            thepcl,plcl_pbe,tlcl_pbe,lcl_pbe,                          &
            p,ht,t,tv,td,w,                                            &
            partem,buoy,wload,mbuoy,pbesnd,mbesnd,                     &
            cin,velneg,cape,mcin,mvelneg,mcape,lfc,el,tcap)
!
!  Calculate Lifted Index
!
 CALL calcli(nlevel,thepcl,p,t,li)
!
!  Calculate max and min wet bulb potential temperature
!
 CALL thwxn(nlevel,p,ht,t,td,ht(1),twdf)
!
 RETURN
END SUBROUTINE sindex
