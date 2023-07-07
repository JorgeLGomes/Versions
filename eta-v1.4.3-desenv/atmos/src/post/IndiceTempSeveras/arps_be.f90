!
!
! BELOW: FUNCTIONS AND SUBROUTINES USED BY PROGRAM ÍNDICES SEVEROS
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ARPS_BE                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######    and modified by Ernani L. Nascimento at Simepar   ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE arps_be(nrad,lm,pres,hgt,tc,tdc,wmr,lcl,lfc,el,twdf,li,cape,mcape,   &
          cin,tcap,p1d,ht1d,t1d,tv1d,td1d,wmr1d,partem,buoy,wload,      &
          mbuoy,pbesnd,mbesnd,nlev,estacao,hora,dia,mes,ano)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the lifting condensation level (lcl), level of free
!  convection (lfc), equilibrium level (el), max wet-bulb potential
!  temperature difference (twdf), lifted index (LI), Convective
!  Available Potential Energy (CAPE), Moist Convective Potential
!  Energy (MCAPE, includes water loading), convective inhibition
!  (CIN) and lid strength (tcap)  over the ARPS grid.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  (Keith Brewster)
!  Cleaned-up, removed OLAPS artifacts.
!
!  FEB 13-14 2004 (ERNANI L. NASCIMENTO)
!  Modified from ARPS to compute indices from rawinsonde data.
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
 INTEGER :: nrad,lm
!
!-----------------------------------------------------------------------
!
!  "1-D" input variables (Ernani L. Nascimento)
!
!-----------------------------------------------------------------------
!
 INTEGER :: nlev(nrad),estacao(nrad),hora(nrad),dia(nrad),mes(nrad),   &
            ano(nrad)
!
!-----------------------------------------------------------------------
!
!  "2-D" input variables (Ernani L. Nascimento)
!
!-----------------------------------------------------------------------
!
 REAL :: pres(lm,nrad)
! hgt is integer! Different from the original code: Ernani L. Nascimento
 INTEGER :: hgt(lm,nrad)
 REAL :: tc(lm,nrad)
 REAL :: wmr(lm,nrad)
 REAL :: tdc(lm,nrad)
!
!-----------------------------------------------------------------------
!
!  Output variables ("1-D") (Ernani L. Nascimento)
!
!-----------------------------------------------------------------------
!
 REAL :: lcl(nrad)
 REAL :: lfc(nrad)
 REAL :: el(nrad)
 REAL :: twdf(nrad)
 REAL :: li(nrad)
 REAL :: cape(nrad)
 REAL :: mcape(nrad)
 REAL :: cin(nrad)
 REAL :: tcap(nrad)
!
!-----------------------------------------------------------------------
!
!  Scratch space for calculations
!
!-----------------------------------------------------------------------
!
 REAL :: p1d(lm),ht1d(lm),tv1d(lm),td1d(lm)
 REAL :: partem(lm),buoy(lm),wload(lm)
 REAL :: mbuoy(lm),pbesnd(lm),mbesnd(lm)
 REAL :: wmr1d(lm),t1d(lm)
!
!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------
!
 REAL :: wmr2td,tctotv
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
 INTEGER :: i,j,k,n,nlevel,bla
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!  Loop over all soundings (Ernani L. Nascimento)
!-----------------------------------------------------------------------
!  print *, 'hgt(m) depois: ',hgt(1,1)
 DO k=1,nrad
  if (nlev(k)==0) then
   lcl(k)=-999.9
   lfc(k)=-999.9
   el(k)=-999.9
   twdf(k)=-999.9
   li(k)=-999.9
   cape(k)=-999.9
   mcape(k)=-999.9
   cin(k)=-999.9
   tcap(k)=-999.9
  else
   j=nlev(k)
   if (pres(j,k)>300.0) then
    write(*,*) 'Sondagem abaixo não chegou aos 300hPa (estacao,hora,dia,mes,ano):'
    write(*,*) estacao(k),hora(k),dia(k),mes(k),ano(k)
    write(*,*) pres(j,k),nlev(k)
    write(*,*) 'Buscando próxima sondagem.'
    cycle
   endif
!-----------------------------------------------------------------------
! Loop over all levels of the sounding
!-----------------------------------------------------------------------
   DO i=1,nlev(k)
     p1d(i) = pres(i,k)
     ht1d(i) = float(hgt(i,k))
     wmr1d(i) = wmr(i,k)
     t1d(i) = tc(i,k)
!      td1d(i) = wmr2td(p1d(i),wmr1d(i))
     td1d(i) = tdc(i,k)
     tv1d(i) = tctotv(t1d(i),wmr1d(i))
   END DO
!
!    print *, 'ht1d(1): ', ht1d(1)
   nlevel = nlev(k)
!
   CALL sindex(nlevel,p1d,ht1d,t1d,tv1d,td1d,wmr1d,partem,buoy,wload,    &
               mbuoy,pbesnd,mbesnd,lcl(k),lfc(k),el(k),twdf(k),li(k),    &
               cape(k),mcape(k),cin(k),tcap(k))
  endif
 END DO
RETURN
END SUBROUTINE arps_be
