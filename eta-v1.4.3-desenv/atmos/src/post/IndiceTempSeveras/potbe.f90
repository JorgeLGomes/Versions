!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE POTBE                      ######
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
SUBROUTINE potbe(nlevel,pmean,tmean,wmean,                              &
          blthte,plcl,tlcl,lcl,                                        &
          p,ht,t,tv,td,w,                                              &
          partem,buoy,wload,mbuoy,pbesnd,mbesnd,                       &
          pbeneg,velneg,pos_max,                                       &
          mbeneg,mvelneg,mpos_max,lfc,el,tcap)
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
!  February, 1994  Based on OLAPS, hence LAPS, version of same.
!                  from FSL, by Steve Albers 1991
!
!  MODIFICATION HISTORY:
!  3/11/1996  Cleaned-up, removed OLAPS artifacts.
!
!  Feb/14/2004  Ernani L. Nascimento: removed variable maxlev
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
!
 INTEGER :: nlevel
 REAL :: pmean,tmean,wmean,blthte,plcl,tlcl,lcl
 REAL :: p(nlevel+2),ht(nlevel+2)
 REAL :: t(nlevel+2),tv(nlevel+2),td(nlevel+2),w(nlevel+2)
 REAL :: partem(nlevel+2),buoy(nlevel+2),wload(nlevel+2),mbuoy(nlevel+2)
 REAL :: pbesnd(nlevel+2),mbesnd(nlevel+2)
 REAL :: pbeneg,velneg,pos_max
 REAL :: mbeneg,mvelneg,mpos_max,lfc,el,tcap
!
!  Parameters
!
 REAL :: g,gamma
 PARAMETER (g=9.80665,                                                 &
            gamma = .009760)   ! Dry Adiabatic Lapse Rate Deg/m
!
!  Functions
!
 REAL :: tsa_fast,tctotv
!
!  Misc internal variables
!
 INTEGER :: n,nel
 REAL :: deltah,delta_ht_dry,delta_ht_wet
 REAL :: sntlcl,buoy_lcl,wsat,partv
 REAL :: nbe_min,pbe_wet,pbe_dry,pos_area
 REAL :: wlow,htzero,adjeng
!
!-----------------------------------------------------------------------
!
!  Function f_mrsat and inline directive for Cray PVP
!
!  Note from Ernani L. Nascimento: this directive is ignored whenever
!  this code is compiled in a machine different from Cray. Thus, no need to
!  delete the piece of code !fpp$, !dir$, !*$* below.
!
!-----------------------------------------------------------------------
!
 REAL :: f_mrsat

!fpp$ expand (f_mrsat)
!dir$ inline always f_mrsat
!*$*  inline routine (f_mrsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!  Reset output variables.
!
!  These should be the same as what is assigned
!  no positive area found.  They are limited in
!  range to allow for contouring when positive
!  areas exist in some columns of a domain and
!  not in others.
!
 pbeneg=-400.
 velneg=20.
 pos_max=0.
 mbeneg=-400.
 mvelneg=20.
 mpos_max=0.
 lfc=10000.
 el=0.
 tcap=0.
!
!  Initialize parcel path arrays
!
 partem(1) = t(1)
 buoy(1) = 0.
 wload(1) = 0.
 mbuoy(1) = 0.
 pbesnd(1) = 0.
 mbesnd(1) = 0.

!  WRITE(6,810)pmean,tmean,wmean,plcl,tlcl,lcl
! 810 format(' pmean,tmean,wmean,plcl,tlcl,lcl',2F10.2,F10.5,2F10.2
!    +   ,F5.1)

 DO n=2,nlevel
   deltah = ht(n) - ht(n-1)
   IF(plcl < p(n-1))THEN ! lower level is below LCL
     IF(plcl < p(n))THEN ! upper level is below LCL
!        WRITE(6,*)' DRY CASE'
       partem(n)=partem(n-1)-gamma*deltah
       partv=tctotv(partem(n),w(1))
       buoy(n)=(partv-tv(n))/tv(n)
       pbesnd(n)=pbesnd(n-1)+g*0.5*(buoy(n)+buoy(n-1))*deltah
       wload(n)=0.
       mbuoy(n)=buoy(n)
       mbesnd(n)=pbesnd(n)
       IF((p(1)-p(n)) < 300.) tcap=AMAX1(tcap,(tv(n)-partv))

     ELSE ! Upper level is above LCL
!
!  BRACKETING CASE AROUND lcl - DRY ADIABATIC PART
!
!        WRITE(6,*)' DRY ADIABATIC PART'
       delta_ht_dry = lcl - ht(n-1)
!        WRITE(6,307)tlcl
!307        format(' PARCEL TEMP AT lcl= ',F10.3)
       CALL intrpr(nlevel,p,tv,plcl,sntlcl)
       partv=tctotv(tlcl,w(1))
       buoy_lcl=(partv-sntlcl)/sntlcl
       pbe_dry=g*0.5*(buoy_lcl+buoy(n-1))*delta_ht_dry
       IF((p(1)-plcl) < 300.) tcap=AMAX1(tcap,(sntlcl-partv))
!        WRITE(6,777)N,P(N),tlcl,sntlcl,buoy_lcl
!#          ,buoy(N-1),delta_ht_dry,HT(N),pbesnd(N-1)+pbe_dry
!
!        MOIST ADIABATIC PART
!
!        WRITE(6,*)' MOIST ADIABATIC PART'
       delta_ht_wet=deltah-delta_ht_dry

       partem(n) = tsa_fast(blthte,p(n))
       wsat=1000.*f_mrsat( p(n)*100., partem(n)+273.15 )
       partv=tctotv(partem(n),wsat)
       buoy(n)=(partv-tv(n))/tv(n)
       pbe_wet = g*0.5*(buoy(n)+buoy_lcl)*delta_ht_wet
       pbesnd(n)=pbesnd(n-1) + pbe_dry + pbe_wet
!
       wload(n)=0.001*(w(1)-wsat)
       mbuoy(n)=buoy(n) - wload(n)
       pbe_wet = g*0.5*(mbuoy(n)+buoy_lcl)*delta_ht_wet
       mbesnd(n)=mbesnd(n-1) + pbe_dry + pbe_wet
       IF((p(1)-plcl) < 300.) tcap=AMAX1(tcap,(tv(n)-partv))

     END IF ! Upper level below LCL (Dry or bracket)
   ELSE ! Lower Level is above LCL
!      WRITE(6,*)' GETTING PARCEL TEMPERATURE FOR MOIST CASE'
     partem(n) = tsa_fast(blthte,p(n))

     wsat=1000.*f_mrsat( p(n)*100., partem(n)+273.15 )
     partv=tctotv(partem(n),wsat)
     buoy(n)=(partv-tv(n))/tv(n)
     pbesnd(n)=pbesnd(n-1)+g*0.5*(buoy(n)+buoy(n-1))*deltah
!
     wload(n)=0.001*(w(1)-wsat)
     mbuoy(n)=buoy(n) - wload(n)
     mbesnd(n)=mbesnd(n-1)+g*0.5*(mbuoy(n)+mbuoy(n-1))*deltah
     IF((p(1)-p(n)) < 300.) tcap=AMAX1(tcap,(tv(n)-partv))

   END IF

!    WRITE(6,777)N,P(N),partem(N),T(N),(buoy(n)*1000.),pbesnd(n)
!777    format(' PBE: P,partem,t,b,pbe=',I3,F6.1,4F8.2)
 END DO
!
!  DETERMINE ENERGY EXTREMA
!  Find heights with nuetral buoyancy
!
 pos_area=0.
 nbe_min=0.
 DO n=2,nlevel
!    WRITE(6,940)N
!940    format(
!    :' LOOKING FOR NEUTRAL BUOYANCY - ENERGY EXTREMUM, LEVEL',I3)

   IF((buoy(n)*buoy(n-1)) < 0.)THEN
     wlow=buoy(n)/(buoy(n)-buoy(n-1))
     htzero=ht(n)*(1.-wlow) + wlow*ht(n-1)
     deltah=htzero-ht(n-1)
     adjeng=pbesnd(n-1)+g*0.5*buoy(n-1)*deltah
!
     IF (p(n) >= 500.)  THEN
       nbe_min=AMIN1(adjeng,nbe_min)
     END IF
!
     pos_area=adjeng-nbe_min
     pos_max=AMAX1(pos_area,pos_max)
   END IF
 END DO

!  WRITE(6,464)ICP,ICT,N1,NLEVEL
!464  format(' ICP,ICT,N1,NLEVEL',4I5)
!
!  Case when equlibrium level is above top of domain
!
 pos_area=pbesnd(nlevel)-nbe_min
 pos_max=AMAX1(pos_area,pos_max)
!
!  At least one region of positive area in sounding
!  Make sure there is at least 1 J/kg to avoid some
!  round-off errors esp near LCL.
!
 IF(pos_max > 1.0)THEN
   pbeneg=AMAX1(nbe_min,-400.)
   velneg=SQRT(2.0*ABS(pbeneg))
   velneg=AMIN1(velneg,20.)
 ELSE ! Case when no positive area exists anywhere in sounding
   pos_max=0.0
   pbeneg =-400.
   velneg = 20.
 END IF
!  WRITE(6,485)pos_max,PBENEG,VELNEG
!485  format(' pos_max',F10.1,' PBENEG',F10.1,' VELNEG',F10.1)

!
!  DETERMINE ENERGY EXTREMA FOR MOIST BUOYANCY
!  Find heights with nuetral buoyancy
!
 pos_area=0.
 nbe_min=0.
 DO n=2,nlevel
!    WRITE(6,940)N

   IF((mbuoy(n)*mbuoy(n-1)) < 0.)THEN
     wlow=mbuoy(n)/(mbuoy(n)-mbuoy(n-1))
     htzero=ht(n)*(1.-wlow) + wlow*ht(n-1)
     deltah=htzero-ht(n-1)
     adjeng=mbesnd(n-1)+g*0.5*mbuoy(n-1)*deltah
!
     IF (p(n) >= 500.)  THEN
       nbe_min=AMIN1(adjeng,nbe_min)
     END IF
!
     pos_area=adjeng-nbe_min
     mpos_max=AMAX1(pos_area,mpos_max)
   END IF
 END DO

!  WRITE(6,464)ICP,ICT,N1,NLEVEL
!
!  Case when equlibrium level is above top of domain
!
 pos_area=mbesnd(nlevel)-nbe_min
 mpos_max=AMAX1(pos_area,mpos_max)
!
!  At least one region of positive area in sounding
!  Make sure there is at least 1 J/kg to
!  spurious pos energy due to round off.
!
 IF(mpos_max > 1.0)THEN
   mbeneg=AMAX1(nbe_min,-400.)
   mvelneg=SQRT(2.0*ABS(pbeneg))
   mvelneg=AMIN1(mvelneg,20.)
 ELSE ! Case when no positive area exists anywhere in sounding
   mpos_max=0.0
   mbeneg =-400.
   mvelneg = 20.
 END IF
!  WRITE(6,486)mpos_max,PBENEG,VELNEG
!486  format(' Mpos_max',F10.1,' MBENEG',F10.1,' mVELNEG',F10.1)
!
!    Case when equlibrium level is above top of domain
!
 mpos_max = MAX(mpos_max,(mbesnd(nlevel) - nbe_min))
!
!  Find EL and LFC
!  Unxts are set to km ASL
!
 IF(pos_max > 1.0) THEN
   IF(buoy(nlevel) > 0.) THEN
     nel=nlevel
     el=0.001*ht(nlevel)
   ELSE
     DO  n=nlevel-1,2,-1
       IF(buoy(n) > 0.) EXIT
     END DO
!      1201     CONTINUE
     nel=n
     wlow=buoy(n+1)/(buoy(n+1)-buoy(n))
     el=0.001 * (ht(n+1)*(1.-wlow) + ht(n)*wlow)
   END IF
!
   DO n=nel,1,-1
     IF(buoy(n) < 0.) EXIT
   END DO
!    1301   CONTINUE
   IF(n > 0) THEN
     wlow=buoy(n+1)/(buoy(n+1)-buoy(n))
     lfc=ht(n+1)*(1.-wlow) + ht(n)*wlow
   ELSE
     lfc=ht(1)
   END IF
 ELSE
   el=0.
   lfc=10000.
 END IF
 lfc=AMIN1(lfc,10000.)
 RETURN
END SUBROUTINE potbe
