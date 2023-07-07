      SUBROUTINE SIG2PSPLINE(IMOUT,JMOUT)
C     ******************************************************************
C     *                                                                *
C     *  THIS IS THE PROGRAM FOR CONVERSION FROM SIGMA TO P SYSTEM     *
C     *  USING SPLINE FITTING                                          *
C     *  02-1981 TO 03- 1998 Z. JANJIC, D. JOVIC, S. NICKOVIC          *
C     *  09-2000  H CHUANG MODIFIED THE PROGRAM TO BE INCLUDED IN      *
C     *              OPERATIONAL POST                                  *
C     *  NOTE: CURRENT POST ONLY PROCESS UP TO 1000 MB PRESSURE LEVEL  *
C     *        2 ADDITIONAL LEVELS (1025 1050) ARE ADDED DURING SPLINE *
C     *        FITTING COMPUTATION BUT THE FIELDS ON THESE TWO LEVELS  *
C     *        ARE NOT OUTPUT TO GRID FILES                            *   
C     ******************************************************************
!-----------------------------------------------------------------------
      INCLUDE "parmeta"
      INCLUDE "parmout"
      INCLUDE "params"    
!-----------------------------------------------------------------------
Cmptest      parameter(nsmud=100,lp2=lm+2,LSMP2=LSM+2)
      parameter(nsmud=25,lp2=lm+2,LSMP2=LSM+2)
      parameter(zero=1.e-5)
      parameter(iii=2,jjj=2)
c-----------------------------------------------------------------------
      LOGICAL IOOMG,IOALL
                             d i m e n s i o n
     &ihw(jm),ihe(jm),ivw(jm),ive(jm)
c
     &,ztt   (lm+1)
     &,tcol  (lm+1),ocol  (lm+1),qcol  (lm+1),cwcl  (lm+1)
     &,zth   (lm+2)
     &,hcol  (lm+2),q2cl  (lm+2)
     &,ztu   (lm+1),ztv   (lm+1)
     &,ucol  (lm+1),vcol  (lm+1)
c
     &,y2    (lm+2),phld  (lm+2),qhld  (lm+2)
c
     &,ztsl  (lsmp2),ovrlx (lsmp2)
     &,tcolsl(lsmp2),ocolsl(lsmp2),qcolsl(lsmp2),cwclsl(lsmp2)
     &,hcolsl(lsmp2),q2clsl(lsmp2)
     &,ucolsl(lsmp2),vcolsl(lsmp2)
C
     &,hcol3(3),pcol3(3),y23(3)
c     
     &,tmask (im,jm),hs    (im,jm)
c
     &,tsll  (im,jm),osll  (im,jm),qsll  (im,jm),qcll  (im,jm)
     &,fsll  (im,jm),q2ll  (im,jm)
c
     &,zet   (im,jm,lm+1)
c
     &,tprs(im,jm,lsmp2) ,oprs(im,jm,lsmp2),qprs(im,jm,lsmp2)
     &,fprs(im,jm,lsmp2) ,qcprs(im,jm,lsmp2),q2prs(im,jm,lsmp2)
     &,uprs(im,jm,lsmp2) ,vprs(im,jm,lsmp2)
c
      real iw(im,jm,lm),IWU,IWL,icecl(lm+1),icelsl(lsmp2),icell(im,jm)
     &,iceprs(im,jm,lsmp2) 
      REAL EGRID1(IM,JM),EGRID2(IM,JM)
      REAL GRID1(IMOUT,JMOUT),GRID2(IMOUT,JMOUT)
cc      REAL,ALLOCATABLE,DIMENSION(:,:,:)::tprs,oprs,qprs,fprs
cc     +,qcprs,q2prs,uprs,vprs 
!-----------------------------------------------------------------------
c      equivalence
c     & (t(1),tsl(1))
c     &,(q(1),qsl(1))
c     &,(cwm(1),qcsl(1))
c     &,(q2(1),q2sl(1))
c     &,(omgalf(1),osl(1))
c     &,(div(1),dsl(1))
c     &,(u(1),usl(1))
c     &,(v(1),vsl(1))
c     &,(dwdt(1),dwsl(1))
c     &,(w(1),wsl(1))
!-----------------------------------------------------------------------
C     INCLUDE COMMON BLOCKS.
      INCLUDE "CTLBLK.comm"
      INCLUDE "OMGAOT.comm"
      INCLUDE "LOOPS.comm"
      INCLUDE "MASKS.comm"
      INCLUDE "MAPOT.comm"
      INCLUDE "VRBLS.comm"
      INCLUDE "PVRBLS.comm"
      INCLUDE "RQSTFLD.comm"
      INCLUDE "EXTRA.comm"
      INCLUDE "CLDWTR.comm"
      INCLUDE "E2PFLG.comm"
      INCLUDE "DYNAMD.comm"
!-----------------------------------------------------------------------
                              d a t a
     & y2/lp2*0./
!-----------------------------------------------------------------------
                              d a t a
     & ovrlx/lsmp2*0.175/
!-----------------------------------------------------------------------
C     VERTICAL INTERPOLATION.  EXECUTE ONLY
C     IF THERE'S SOMETHING WE WANT.
C
      IF((IGET(012).GT.0).OR.(IGET(013).GT.0).OR.
     X   (IGET(014).GT.0).OR.(IGET(015).GT.0).OR.
     X   (IGET(016).GT.0).OR.(IGET(017).GT.0).OR.
     X   (IGET(018).GT.0).OR.(IGET(019).GT.0).OR.
     X   (IGET(020).GT.0).OR.(IGET(030).GT.0).OR.
     X   (IGET(021).GT.0).OR.(IGET(022).GT.0).OR.
     X   (IGET(153).GT.0).OR.(IGET(166).GT.0).OR.
     X   (IGET(23).GT.0)) THEN
C
C  SET UP UTIM FOR THIS TIME STEP
C
       UTIM=1.
       CLIMIT=1.E-20
c
       DO 75 L=1,LM
        IF(L.EQ.1)THEN
         DO J=JSTA,JEND
          DO I=1,IM
           IW(I,J,L)=0.
          ENDDO
         ENDDO
         GO TO 75
        ENDIF
        doout70: DO J=JSTA,JEND
         doin70: DO I=1,IM
          LML=LM-LMH(I,J)
          HH=HTM(I,J,L)*HBM2(I,J)
          TKL=T(I,J,L)
          QKL=Q(I,J,L)
          CWMKL=CWM(I,J,L)
          TMT0=(TKL-273.16)*HH
          TMT15=AMIN1(TMT0,-15.)*HH    
          PP=PDSL(I,J)*AETA(L)+PT
          QW=HH*PQ0/PP*EXP(HH*A2*(TKL-A3)/(TKL-A4))
          QI=QW*(1.+0.01*AMIN1(TMT0,0.))     
          U00KL=U00(I,J)+UL(L+LML)*(0.95-U00(I,J))*UTIM
C
          IF(TMT0.LT.-15.0)THEN
           FIQ=QKL-U00KL*QI
           IF(FIQ.GT.0..OR.CWMKL.GT.CLIMIT) THEN
            IW(I,J,L)=1.
           ELSE
            IW(I,J,L)=0.
           ENDIF
          ENDIF
C
          IF(TMT0.GE.0.0)IW(I,J,L)=0.
          IF(TMT0.LT.0.0.AND.TMT0.GE.-15.0)THEN
           IW(I,J,L)=0.
           IF(IW(I,J,L-1).EQ.1..AND.CWMKL.GT.1.E-20)IW(I,J,L)=1.
          ENDIF
C
        END DO doin70
        END DO doout70
   75  CONTINUE
c       	  
       do 100 l=1,lsm
        ztsl(l)=alog(spl(l))**2
 100   continue
       ztsl(lsm+1)=alog(102500.)**2  
       ztsl(lsmp2)=alog(105000.)**2	
c
       ztbot=ztsl(lsmp2)
c
       do 110 j=jsta,jend
        ihw(j)=  -mod(j,2)
        ihe(j)= 1-mod(j,2)
        ivw(j)=-1+mod(j,2)
        ive(j)=   mod(j,2)
 110   continue
c
cc       allocate(tprs(im,jm,lsmp2)) 
cc       allocate(oprs(im,jm,lsmp2))
cc       allocate(qprs(im,jm,lsmp2))
cc       allocate(fprs(im,jm,lsmp2))
cc       allocate(qcprs(im,jm,lsmp2))
cc       allocate(q2prs(im,jm,lsmp2))
cc       allocate(uprs(im,jm,lsmp2))
cc       allocate(vprs(im,jm,lsmp2))
c--------------mass point variables at pressure levels------------------
       do 200 j=jsta,jend
        ihl=1
        ihh=im-1+mod(j,2)
        do 201 i=ihl,ihh
	 if(abs(pdsl(i,j)-pd(i,j)).gt. 1)print*,'inconsistent pd',i,j
     +,  pd(i,j),pdsl(i,j) 	 
         pdp=pdsl(i,j)
c
         hsp=fis(i,j)/g
         hcol(lm+1)=hsp
         hs(i,j)=hsp
c
c         alp1l=alog(pint(i,j,lm+1))   ! for nonhydrostatic version
         alp1l=alog(pt+pdp)
         alp2l=(alog(pdp+pt))**2
         zth(lm+1)=alp2l
c
c         pcol(1)=pt
         q2cl(1)=0.
c--------------loading values at mid-layers and interfaces--------------
         do 220 ivi=1,lm
          l=lm+1-ivi
c
c          alp1u=alog(pint(i,j,l))   ! for nonhydrostatic version
          alp1u=alog(eta(l)*pdp+pt)
          alp2u=(alog(eta(l)*pdp+pt))**2
c
          dpd=(eta(l+1)-eta(l))*pdp
c
c         dh=(q(i,j,l)*0.608+1.)*(alp1l-alp1u)*t(i,j,l)*RD/g
c     2   *dpd/(pint(i,j,l+1)-pint(i,j,l))
          dh=(q(i,j,l)*0.608+1.)*(alp1l-alp1u)*t(i,j,l)*RD/g
c
          alpcp=dh*g/(dpd*cp)
c
          zth(l)=alp2u
czj!!      ztt(l)=(alp2l+alp2u)*0.5
          ztt(l)=alog((eta(l)+eta(l+1))*0.5*pdp+pt)**2
c	  if(i.eq.120.and.j.eq.279)print*,'l,eta(l),eta(l+1),ztt='
c     +,   l,eta(l),eta(l+1),ztt(l),pdp,pt	  
c
          tcol(l)=t(i,j,l)
          ocol(l)=omga(i,j,l)/alpcp
          qcol(l)=q(i,j,l)
          cwcl(l)=cwm(i,j,l)*(1-iw(i,j,l))
	  icecl(l)=cwm(i,j,l)*iw(i,j,l)
c
c          pcol(l+1)=pint(i,j,l+1)
          q2cl(l+1)=q2(i,j,l+1)     
          hcol(l  )=hcol(l+1)+dh
c
          alp1l=alp1u
          alp2l=alp2u
 220     continue
c-------------extrapolation underground for mid-layer variables---------
         if((ztt(lm)-ztbot) .ge. 0.0)then
          lmd=lm
         else
          lmd=lm+1
c
          ztt(lm+1)=ztbot
c
          dum=(tcol(lm)-tcol(lm-4))/(ztt(lm)-ztt(lm-4))
          x=ztbot-ztt(lm)
c
          tcol(lm+1)=dum*x+tcol(lm)
          ocol(lm+1)=0.
          qcol(lm+1)=0.
          cwcl(lm+1)=0.
	  icecl(lm+1)=0.
         endif
         if(i.eq.iii.and.j.eq.jjj)then
	  print*,'on sigma at ',iii,jjj
	  print*,'terrain= ',fis(i,j)/g
	  do l=1,lmd
	   print*,'l,P,T= ',l,exp(ztt(l)**0.5),tcol(l) 	
          end do
	 end if
c         if(i.eq.120.and.j.eq.279)then
c	  print*,'lmd, ztbot= ',lmd,ztbot
c	  print*,'ztt=',ztt
c	 end if 
c-------------interpolation of mid-layer variables to pressure levels---
         y2(lmd)=0.
c	
         call splinef(lmd,ztt,tcol,y2,lsmp2,ztsl,tcolsl,phld,qhld)
         call splinef(lmd,ztt,ocol,y2,lsmp2,ztsl,ocolsl,phld,qhld)
         call splinef(lmd,ztt,qcol,y2,lsmp2,ztsl,qcolsl,phld,qhld)
         call splinef(lmd,ztt,cwcl,y2,lsmp2,ztsl,cwclsl,phld,qhld)
	 call splinef(lmd,ztt,icecl,y2,lsmp2,ztsl,icelsl,phld,qhld)
c
         do 230 l=1,lsmp2
          tprs (i,j,l)=tcolsl(l)
	  oprs (i,j,l)=ocolsl(l)
          qprs (i,j,l)=qcolsl(l)
          qcprs(i,j,l)=cwclsl(l)
	  iceprs(i,j,l)=icelsl(l)
 230     continue
c-------------extrapolation underground for interface variables---------
         if(zth(lm+1).ge.ztbot)    then
          lmd=lm+1      
         else
          lmd=lm+2
c
          zth(lm+2)=ztbot
          x=ztbot-zth(lm+1)
c          if(abs(x/ztbot).lt.zero)x=0.
c
          dum=(hcol(lm+1)-hcol(lm-5))/(zth(lm+1)-zth(lm-5))
          d2=0.
c
c          if(x.lt.0.0000001)print*,'i,j,zth(lm+1),ztbot,x= ',i,j
c     +	  ,zth(lm+1),ztbot,x
          hcol(lm+2)=d2*x*x+dum*x+hcol(lm+1)
c          pcol(lm+2)=spl(lsmp2)
          q2cl(lm+2)=0.
         endif
c-------------interpolation of interface variables to pressure levels---
         if(i.eq.iii.and.j.eq.jjj)then
	  print*,'on sigma at ',iii,jjj
	  do l=1,lmd
	   print*,'l,H,P= ',l,hcol(l),exp(zth(l)**0.5) 	
          end do
	 end if
c
         y2(lmd)=0.
c
         call splinef(lmd,zth,hcol,y2,lsmp2,ztsl,hcolsl,phld,qhld)
c         call splinef(lmd,zth,pcol,y2,lsmp2,ztsl,pcolsl,phld,qhld)
         call splinef(lmd,zth,q2cl,y2,lsmp2,ztsl,q2clsl,phld,qhld)
c
         if(i.eq.iii.and.j.eq.jjj)then
	  print*,'on pressure before relaxation at ',iii,jjj
	  do l=1,lsmp2
	   print*,'l,H,P,T= ',l,hcolsl(l),exp(ztsl(l)**0.5)
     +	   ,tcolsl(l) 	
          end do
	 end if
c	 
         do 240 l=1,lsmp2
          fprs (i,j,l)=hcolsl(l)
c          psl (i,j,l)=pcolsl(l)
          q2prs(i,j,l)=max(q2clsl(l),0.)
 240     continue
c
         do 250 l=1,lm+1
          zet(i,j,l)=zth(l)
 250     continue
 201    continue  
 200   continue
c--------------filling boundaries--------------------------------------
       do 260 l=1,lm+1
        do 261 j=jsta,jend
         zet(im,j,l)=zet(im-1,j,l)
 261    continue
 260   continue
c--------------end of mass point interpolations------------------------
c
c
       print*,'in sig2pspline before relaxation of mass points'
       do 300 l=1,lsmp2
c--------------filling remaining undefined mass point boundaries-------
        if(mod(jsta,2) .lt. 1)then   !even jsta
         do 351 j=jsta,jend,2
          tprs(im,j,l)=tprs(im-1,j,l)
          oprs(im,j,l)=oprs(im-1,j,l)
          qprs(im,j,l)=qprs(im-1,j,l)
          qcprs(im,j,l)=qcprs(im-1,j,l)
	  iceprs(im,j,l)=iceprs(im-1,j,l)
          fprs(im,j,l)=fprs(im-1,j,l)
          q2prs(im,j,l)=q2prs(im-1,j,l)
 351     continue
        else
	 do 352 j=jsta+1,jend,2
          tprs(im,j,l)=tprs(im-1,j,l)
          oprs(im,j,l)=oprs(im-1,j,l)
          qprs(im,j,l)=qprs(im-1,j,l)
          qcprs(im,j,l)=qcprs(im-1,j,l)
	  iceprs(im,j,l)=iceprs(im-1,j,l)
          fprs(im,j,l)=fprs(im-1,j,l)
          q2prs(im,j,l)=q2prs(im-1,j,l)
 352     continue
        end if    
c	if(l.eq.lsm)then
c	 print*,'1000 mb enormal height'
c         do i=1,im
c	  do j=jsta,jend
c	   if(fprs(i,j,l).gt.200.)print*,i,j,fprs(i,j,l)
c          end do
c         end do
c	end if 
c--------------smoothing pressure levels ruptured by topography--------
        rlx=ovrlx(l)
        yes=0.
c
        avt=0.
        avh=0.
        kt=0
        kh=0 
        do 310 j=jsta,jend
         ihl=1
         ihh=im-1+mod(j,2)
         do 311 i=ihl,ihh
          href=hs(i,j)
          if(href.gt.300.) then
           href=href+1500.
          endif
          if(href.ge.fprs(i,j,l))then
           yes=1.
           tmask(i,j)=rlx
          else
           tmask(i,j)=0.
c
c      avt=tprs(i,j,l)+avt
c      avh=fprs(i,j,l)+avh
c      kt=kt+1
c      kh=kh+1
c
          endif
c          if(href.ge.fprs(i,j,l).and.l.lt.lsm)then
c	   tprs(i,j,l)=(((fprs(i,j,l-1)-fprs(i,j,l)))
c     &           /((alog(spl(l))-alog(spl(l-1))))
c     &           +((fprs(i,j,l)-fprs(i,j,l+1)))
c     &           /((alog(spl(l+1))-alog(spl(l)))))*0.5*g/rd
c          end if
          if(href.ge.fprs(i,j,l))then
           if(l.lt.lsmp2)then
            tprs(i,j,l)=(((fprs(i,j,l-1)-fprs(i,j,l)))
     &	         /(ztsl(l)**0.5-ztsl(l-1)**0.5)
     &           +((fprs(i,j,l)-fprs(i,j,l+1)))
     &           /(ztsl(l+1)**0.5-ztsl(l)**0.5))*0.5*g/rd
           else
            tprs(i,j,l)=tprs(i,j,l-1)+(tprs(i,j,l-1)-tprs(i,j,l-2))
     &       *(ztsl(l)-ztsl(l-1))/(ztsl(l-1)-ztsl(l-2))  
           end if
          end if 
 311     continue
 310    continue
c
c??        if(yes.gt.0.)then
c
c              if(kt*kh.ne.0)    then
c          avt=avt/kt
c          avh=avh/kh
c              endif
c
         do 312 j=jsta,jend
          ihl=1
          ihh=im-1+mod(j,2)
          do 313 i=ihl,ihh
           if(hs(i,j).ge.fprs(i,j,l))then
            oprs(i,j,l)=0.
            qprs(i,j,l)=0.
            qcprs(i,j,l)=0.
	    iceprs(i,j,l)=0.
            q2prs(i,j,l)=0.
           endif
 313      continue
 312     continue 
         do 320 n=1,nsmud
	  call exch(tprs(1,1,l))
	  call exch(oprs(1,1,l))
	  call exch(qprs(1,1,l))
	  call exch(qcprs(1,1,l))
	  call exch(iceprs(1,1,l))
	  call exch(fprs(1,1,l))
	  call exch(q2prs(1,1,l))
          do 330 j=jsta_m,jend_m
           ihl=1+mod(j,2)
           ihh=im-1
           do 331 i=ihl,ihh
            if(tmask(i,j).gt.0.05) then
             tsll(i,j)=tprs(i+ihw(j),j-1,l)+tprs(i+ihe(j),j-1,l)
     2         +tprs(i+ihw(j),j+1,l)+tprs(i+ihe(j),j+1,l)-tprs(i,j,l)*4.
             osll(i,j)=oprs(i+ihw(j),j-1,l)+oprs(i+ihe(j),j-1,l)
     2         +oprs(i+ihw(j),j+1,l)+oprs(i+ihe(j),j+1,l)-oprs(i,j,l)*4.
             qsll(i,j)=qprs(i+ihw(j),j-1,l)+qprs(i+ihe(j),j-1,l)
     2         +qprs(i+ihw(j),j+1,l)+qprs(i+ihe(j),j+1,l)-qprs(i,j,l)*4.
             qcll(i,j)=qcprs(i+ihw(j),j-1,l)+qcprs(i+ihe(j),j-1,l)
     2       +qcprs(i+ihw(j),j+1,l)+qcprs(i+ihe(j),j+1,l)-qcprs(i,j,l)*4.
             icell(i,j)=iceprs(i+ihw(j),j-1,l)+iceprs(i+ihe(j),j-1,l)
     2       +iceprs(i+ihw(j),j+1,l)+iceprs(i+ihe(j),j+1,l)
     3       -iceprs(i,j,l)*4.
             fsll(i,j)=fprs(i+ihw(j),j-1,l)+fprs(i+ihe(j),j-1,l)
     2         +fprs(i+ihw(j),j+1,l)+fprs(i+ihe(j),j+1,l)-fprs(i,j,l)*4.
             q2ll(i,j)=q2prs(i+ihw(j),j-1,l)+q2prs(i+ihe(j),j-1,l)
     2       +q2prs(i+ihw(j),j+1,l)+q2prs(i+ihe(j),j+1,l)-q2prs(i,j,l)*4.
            endif
 331       continue
 330      continue
          do 340 j=jsta_m,jend_m
           ihl=1+mod(j,2)
           ihh=im-1
           do 341 i=ihl,ihh
            if(tmask(i,j).gt.0.05) then
             tprs (i,j,l)=tsll(i,j)*tmask(i,j)+tprs (i,j,l)
             oprs (i,j,l)=osll(i,j)*tmask(i,j)+oprs (i,j,l)
             qprs (i,j,l)=qsll(i,j)*tmask(i,j)+qprs (i,j,l)
             qcprs(i,j,l)=qcll(i,j)*tmask(i,j)+qcprs(i,j,l)
	     iceprs(i,j,l)=icell(i,j)*tmask(i,j)+iceprs(i,j,l)
             fprs (i,j,l)=fsll(i,j)*tmask(i,j)+fprs (i,j,l)
             q2prs(i,j,l)=max(q2ll(i,j)*tmask(i,j)+q2prs(i,j,l),0.)
            endif
 341       continue
 340      continue
 320     continue
c??        endif
 300   continue
       print*,'in sig2pspline after relaxation of mass points'
ccc Recalculate QC so that RHs underground are the same as that at lowest
ccc level above ground
       doout346: do j=jsta,jend
        ihl=1
        ihh=im-1+mod(j,2)
        doin346: do i=ihl,ihh 
         do347: do l=1,lsmp2
	  gar1=(eta(lm)+eta(lm+1))*0.5*pdsl(i,j)+pt
	  gar2=exp(ztsl(l)**0.5)
	  if(gar2.gt.gar1)then
	   qw=PQ0/exp(ztsl(l-1)**0.5)*EXP(A2*(tprs(i,j,l-1)-A3)
     +	   /(tprs(i,j,l-1)-A4))
	   rhl=qprs(i,j,l-1)/qw  
	   if(rhl.gt.1)then
cc	    print*,'enormal rh',i,j,rhl
	    rhl=1.
	   end if 
	   if(i.eq.iii.and.j.eq.jjj)print*,'sample rh',rhl
	   do ll=l,lsmp2
	    qw=PQ0/exp(ztsl(ll-1)**0.5)*EXP(A2*(tprs(i,j,ll-1)-A3)
     +  	    /(tprs(i,j,ll-1)-A4))
	    qprs(i,j,ll)=rhl*qw
	   end do
	   exit do347  
          end if
         end do do347
         end do doin346
         end do doout346
c
c--------------velocity point variables at pressure levels--------------
       do l=1,lm+1
        call exch(zet(1,1,l))
       end do
       print*,'in sig2pspline after relaxation of zet'
       do 400 j=jsta_m,jend_m
        ivl=2-mod(j,2)
        ivh=im-1
        do 401 i=ivl,ivh
         zup=(zet(i,j-1,1)+zet(i+ivw(j),j,1)
     2    +zet(i+ive(j),j,1)+zet(i,j+1,1))
         zupu=(zet(i+ivw(j),j,1)+zet(i+ive(j),j,1)+zup)/12.
         zupv=(zet(i,j-1,1)+zet(i,j+1,1)+zup)/12.
c--------------loading variables into vertical columns------------------
         do 410 l=1,lm
          zlo=(zet(i,j-1,l+1)+zet(i+ivw(j),j,l+1)
     2    +zet(i+ive(j),j,l+1)+zet(i,j+1,l+1))
          zlou=(zet(i+ivw(j),j,l+1)+zet(i+ive(j),j,l+1)+zup)/12.
          zlov=(zet(i,j-1,l+1)+zet(i,j+1,l+1)+zup)/12.
c
          ztu(l)=zlou+zupu
          ztv(l)=zlov+zupv
c
          ucol(l)=u(i,j,l)
          vcol(l)=v(i,j,l)
cc	  if(ucol(l).gt.100.)print*,'large I u ',i,j,l,ucol(l)
cc	  if(vcol(l).gt.100.)print*,'large I v ',i,j,l,vcol(l)
c
          zup=zlo
          zupu=zlou
          zupv=zlov
 410     continue 
c
cc         if(amin1(ztu(lm),ztv(lm)).ge.ztbot)then
cc          lmd=lm
cc         else
cc          lmd=lm+1
cc          ztu(lm+1)=ztbot
cc          ztv(lm+1)=ztbot
cc          ucol(lm+1)=ucol(lm)
cc          vcol(lm+1)=vcol(lm)
cc         endif
c
         if(ztu(lm).ge.ztbot)then
          lmd=lm
         else
          lmd=lm+1
          ztu(lm+1)=ztbot
          ucol(lm+1)=ucol(lm)
         endif
c	 
         y2(lmd)=0.
c	 
         if(i.eq.iii.and.j.eq.jjj)then
	  print*,'large input wind input at',i,j
	  do l=1,lmd
	   print*,l,ztu(l),ztv(l),ucol(l),vcol(l)
	  end do
	 end if
c	 
         call splinef(lmd,ztu,ucol,y2,lsmp2,ztsl,ucolsl,phld,qhld)
c	 
         if(ztv(lm).ge.ztbot)then
          lmd=lm
         else
          lmd=lm+1
          ztv(lm+1)=ztbot
          vcol(lm+1)=vcol(lm)
         endif
	 y2(lmd)=0.
         call splinef(lmd,ztv,vcol,y2,lsmp2,ztsl,vcolsl,phld,qhld) 
c
         do 420 l=1,lsmp2
cc	  if(ucolsl(l).gt.100.)print*,'large O u ',i,j,l,ucolSL(l)
cc	  if(vcolSL(l).gt.100.)print*,'large O v ',i,j,l,vcolSL(l)
          uprs(i,j,l)=ucolsl(l)
          vprs(i,j,l)=vcolsl(l)
 420     continue
         if(i.eq.iii.and.j.eq.jjj)then
	  print*,'output wind before geostrophic adjus',i,j
	  do l=1,lsmp2
	   print*,l,ztsl(l),uprs(i,j,l),vprs(i,j,l)
	  end do
	 end if 
cc         if(hs(i,j).gt.300.) then
cc          href=hs(i,j)+1500.
cc         endif
	 do l=1,lsmp2
          if(ztsl(l).ge.amin1(ztu(lm),ztv(lm)))then
	   do ll=l,lsmp2
            uprs(i,j,ll)=uprs(i,j,l-1)
	    vprs(i,j,ll)=vprs(i,j,l-1) 
cc	    ihve=i+mod(j,2)
cc	    ihvw=ihve-1
cc	    TPHI=(J-JMT)*DPHD*DTR
cc	    favg=0.25*(F(I+IHVE,J)+F(I+IHVW,J)
cc     1           +F(I,J+1)+F(I,J-1))*2./dt
cc     2           +uprs(i,j,l-1)*tan(tphi)/erad 
cc            if(i.eq.iii.and.j.eq.jjj)print*,'sample F',favg
cc	    uprs(i,j,ll)=-g/favg*(fprs(i,j+1,ll)-fprs(i,j-1,ll))
cc     1          	    /dy/2.
cc	    vprs(i,j,ll)=g/favg*(fprs(i+ihve,j,ll)-fprs(i+ihvw,j,ll))
cc     1          	    /dx(i,j)/2.
	   end do
	   go to 403
	  end if  
	 end do 
 403     continue
c 
         do l=1,lsmp2 
         if(uprs(i,j,l).gt.100.)print*,'large O u ',i,j,l,uprs(i,j,l)
	 if(vprs(i,j,l).gt.100.)print*,'large O v ',i,j,l,vprs(i,j,l) 
         end do 
         if(i.eq.iii.and.j.eq.jjj)then
	  print*,'output wind after geostrophic adjus',i,j
	  do l=1,lsmp2
	   print*,l,ztsl(l),uprs(i,j,l),vprs(i,j,l)
	  end do
	 end if  
 401    continue
 400   continue
c--------------end of wind point interpolations------------------------
       do l=1,lsmp2
c--------------filling remaining undefined wind point boundaries-------
        do j=jsta_m,jend_m
         uprs(1,j,l)=uprs(2,j,l)
         uprs(im,j,l)=uprs(im-1,j,l)
         vprs(1,j,l)=vprs(2,j,l)
         vprs(im,j,l)=vprs(im-1,j,l)
        enddo
	if(me.eq.0)then
	 do i=1,im	 
          uprs(i,1,l)=uprs(i,3,l)
          vprs(i,1,l)=vprs(i,3,l)
         enddo
	end if
	if(me.eq.(NUM_PROCS-1))then
         do i=1,im	 
          uprs(i,jm,l)=uprs(i,jm-2,l)
          vprs(i,jm,l)=vprs(i,jm-2,l)
         enddo
	end if
c-----------------------------------------------------------------------
       enddo
c--------------sea level pressure---------------------------------------
CCC COMMENT OUT SEA LEVEL PRESSURE REDUCTION BECAUSE IT'S DONE IN QUILT NOW
cc       do 600 j=jsta,jend
cc        ihl=1
cc        ihh=im-1+mod(j,2)
cc        do 601 i=ihl,ihh
cc         if(fis(i,j).gt.-1..and.fis(i,j).lt.1.) then
c          pslp(i,j)=pint(i,j,lm+1)
cc          pslp(i,j)=pt+pd(i,j) 
cc         else
cc	  if(fprs(i,j,lsmp2).gt.0)then
cc	   gar=exp(ztsl(lsmp2)**0.5)+5000.
cc	   hcol3(1)=fprs(i,j,lsmp2)+(fprs(i,j,lsmp2)-fprs(i,j,lsmp2-5))
cc     +     *(alog(gar)**2.-ztsl(lsmp2))/(ztsl(lsmp2)-ztsl(lsmp2-5))
cc           hcol3(2)=fprs(i,j,lsmp2)
cc	   hcol3(3)=fprs(i,j,lsmp2-1) 
cc	   pcol3(1)=gar      
cc           pcol3(2)=exp(ztsl(lsmp2)**0.5)  
cc           pcol3(3)=exp(ztsl(lsmp2-1)**0.5)	 
cc	  else   
cc	   do l=lsmp2,2,-1 
cc	    if(fprs(i,j,l).le.0. .and. fprs(i,j,l-1).gt.0.)then 
cc	     lll=l
cc	     go to 603
cc	    end if 
cc	   end do  
cc 603       hcol3(1)=fprs(i,j,lll)      
cc           hcol3(2)=fprs(i,j,lll-1)    
cc           hcol3(3)=fprs(i,j,lll-2)    
cc           pcol3(1)=exp(ztsl(lll)**0.5)    
cc           pcol3(2)=exp(ztsl(lll-1)**0.5)     
cc           pcol3(3)=exp(ztsl(lll-2)**0.5)     
cc          end if
cc          do l=1,3
cc           y23(l)=0.0
cc          end do
cc          call splinef(3,hcol3,pcol3,y23,1,0.0,slp1,phld,qhld)
cc          pslp(i,j)=slp1
cc         endif
cc         if(pslp(i,j).lt.65000.OR.PSLP(I,J).GT.105000.) print*,
cc     +    'alert, pslp lt 650 or gt 1050 at ',i,j,pslp(i,j)
c     +,  hcol3(2),hcol3(3),pcol3(1),pcol3(2),pcol3(3) 
cc         if(i.eq.iii.and.j.eq.jjj)then
cc	  print*,'on P after relaxation at ',iii,jjj
cc	  do l=1,lsmp2
cc	   print*,'l,H,P,T,U,V= ',l,fprs(i,j,l),exp(ztsl(l)**0.5)
cc     +,    tprs(i,j,l),uprs(i,j,l),vprs(i,j,l) 	
cc          end do
cc	  print*,'pslp = ',pslp(i,j)
cc	 end if    
cc 601    continue
cc 600   continue
C
C*********  END OF VERTICAL INTERPOLATION       
C      
C        INTERPOLATE HORIZONTALLY/OUTPUT SELECTED FIELDS.
C
C---------------------------------------------------------------------
C     
C***  SPECIFIC HUMIDITY.
C
       DO 650 LP=1,LSM
c        print*,'horizontal interp at L=',lp       
        IF(IGET(016).GT.0)THEN
          IF(LVLS(LP,IGET(016)).GT.0)THEN
            CALL E2OUT(016,000,QPRS(1,1,LP),EGRID2,GRID1,GRID2
     +	            ,IMOUT,JMOUT)
            CALL BOUND(GRID1,H1M12,H99999,IMOUT,JMOUT)
            ID(1:25)=0
            print*,'calling output humidity from sig2pspline'
            CALL OUTPUT(IOUTYP,IGET(016),LP,GRID1,IMOUT,JMOUT)
            if(lp.eq.1)print*,'sample SH IOUTYP',IOUTYP,IGET(016)
          ENDIF
        ENDIF
C     
C***  OMEGA
C
        IF(IGET(020).GT.0)THEN
          IF(LVLS(LP,IGET(020)).GT.0)THEN
            CALL E2OUT(020,000,OPRS(1,1,LP),EGRID2,GRID1,GRID2
     +        	    ,IMOUT,JMOUT)
            ID(1:25)=0
            print*,'calling output omeg from sig2pspline'
            CALL OUTPUT(IOUTYP,IGET(020),LP,GRID1,IMOUT,JMOUT)
          ENDIF
        ENDIF
C     
C***  MOISTURE CONVERGENCE
C
        IF(IGET(085).GT.0)THEN
          IF(LVLS(LP,IGET(085)).GT.0)THEN
            CALL CALMCVG(QPRS(1,1,LP),UPRS(1,1,LP),VPRS(1,1,LP)
     1   	    ,-1,EGRID1)
            CALL E2OUT(085,000,EGRID1,EGRID2,
     1                 GRID1,GRID2,IMOUT,JMOUT)
C
C     CONVERT TO DIVERGENCE FOR GRIB UNITS
C
            CALL SCLFLD(GRID1,-1.0,IMOUT,JMOUT)
            ID(1:25)=0
            print*,'calling output moisture from sig2p'
            CALL OUTPUT(IOUTYP,IGET(085),LP,GRID1,IMOUT,JMOUT)
          ENDIF
        ENDIF
C     
C***  U AND/OR V WIND
C
        IF(IGET(018).GT.0.OR.IGET(019).GT.0)THEN
          IF(LVLS(LP,IGET(018)).GT.0.OR.LVLS(LP,IGET(019)).GT.0)THEN
            CALL E2OUT(018,019,UPRS(1,1,LP),VPRS(1,1,LP),GRID1,GRID2
     1	            ,IMOUT,JMOUT)
            ID(1:25)=0
            IF(IGET(018).GT.0)THEN
              CALL OUTPUT(IOUTYP,IGET(018),LP,GRID1,IMOUT,JMOUT)
            ENDIF
            ID(1:25)=0
            IF(IGET(019).GT.0) 
     1       CALL OUTPUT(IOUTYP,IGET(019),LP,GRID2,IMOUT,JMOUT)
          ENDIF
        ENDIF
C     
C***  ABSOLUTE VORTICITY
C
         IF (IGET(021).GT.0) THEN
          IF (LVLS(LP,IGET(021)).GT.0) THEN
            CALL CALVOR(UPRS(1,1,LP),VPRS(1,1,LP),EGRID1)
            CALL E2OUT(021,000,EGRID1,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(021),LP,GRID1,IMOUT,JMOUT)
          ENDIF
         ENDIF
C     
C        GEOSTROPHIC STREAMFUNCTION.
         IF (IGET(086).GT.0) THEN
          IF (LVLS(LP,IGET(086)).GT.0) THEN
!$omp  parallel do
            DO J=JSTA,JEND
            DO I=1,IM
              EGRID2(I,J)=FPRS(I,J,LP)*GI
            ENDDO
            ENDDO
            CALL CALSTRM(EGRID2,EGRID1)
            CALL E2OUT(086,000,EGRID1,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(086),LP,GRID1,IMOUT,JMOUT)
          ENDIF
         ENDIF
C     
C***  TURBULENT KINETIC ENERGY
C
         IF (IGET(022).GT.0) THEN
          IF (LVLS(LP,IGET(022)).GT.0) THEN
            CALL E2OUT(022,000,Q2PRS(1,1,LP),EGRID2,GRID1,GRID2
     +  	    ,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(022),LP,GRID1,IMOUT,JMOUT)
          ENDIF
         ENDIF
C     
C***  TOTAL CLOUD WATER
C
         IF (IGET(153).GT.0) THEN
          IF (LVLS(LP,IGET(153)).GT.0) THEN
            CALL E2OUT(153,000,QCPRS(1,1,LP),EGRID2,GRID1,GRID2
     +  	    ,IMOUT,JMOUT)
            CALL BOUND(GRID1,H1M12,H99999,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(153),LP,GRID1,IMOUT,JMOUT)
          ENDIF
         ENDIF
C
C***  TOTAL CLOUD ICE 
C
         IF (IGET(166).GT.0) THEN
          IF (LVLS(LP,IGET(166)).GT.0) THEN
            CALL E2OUT(166,000,ICEPRS(1,1,LP),EGRID2,GRID1,GRID2
     +  	    ,IMOUT,JMOUT)
            CALL BOUND(GRID1,H1M12,H99999,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(166),LP,GRID1,IMOUT,JMOUT)
          ENDIF
         ENDIF
C
 650   continue
C***  OUTPUT SEA LEVEL PRESSURE IF REQUESTED.
cc        IF(IGET(023).GT.0)THEN
cc          CALL E2OUT(023,000,PSLP,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
cc          ID(1:25)=0
cc          CALL OUTPUT(IOUTYP,IGET(023),LVLS(1,IGET(023)),
cc     1                GRID1,IMOUT,JMOUT)
cc        ENDIF
C
C***  SECOND, STANDARD NGM SEA LEVEL PRESSURE.
C
        IF(IGET(105).GT.0)THEN
          CALL NGMSLP2
          CALL E2OUT(105,000,SLP,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
          ID(1:25)=0
          CALL OUTPUT(IOUTYP,IGET(105),LVLS(1,IGET(105)),
     1        GRID1,IMOUT,JMOUT)
        ENDIF
C
C---------------------------------------------------------------------
C***  OUTPUT THE TEMPERATURES AND QUANTITIES DERIVED FROM IT
C---------------------------------------------------------------------
C
C     
       DO 680 LP=1,LSM
C***  TEMPERATURE
C
        IF(IGET(013).GT.0) THEN
          IF(LVLS(LP,IGET(013)).GT.0)THEN
             CALL E2OUT(013,000,TPRS(1,1,LP),EGRID2,GRID1,GRID2
     1,                 IMOUT,JMOUT)
             ID(1:25)=0
             CALL OUTPUT(IOUTYP,IGET(013),LP,GRID1,IMOUT,JMOUT)
             if(lp.eq.1)print*,'sample T IOUTYP',IOUTYP,IGET(013)
          ENDIF
        ENDIF
C     
C***  POTENTIAL TEMPERATURE.
C
        IF(IGET(014).GT.0)THEN
          IF(LVLS(LP,IGET(014)).GT.0)THEN
!$omp  parallel do
            DO J=JSTA,JEND
            DO I=1,IM
              EGRID2(I,J)=SPL(LP)
            ENDDO
            ENDDO
C
            CALL CALPOT2(EGRID2,TPRS(1,1,LP),EGRID1,IM,JM)
            CALL E2OUT(014,000,EGRID1,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(014),LP,GRID1,IMOUT,JMOUT)
          ENDIF
        ENDIF
C     
C***  RELATIVE HUMIDITY.
C
        IF(IGET(017).GT.0)THEN
          IF(LVLS(LP,IGET(017)).GT.0)THEN
!$omp  parallel do
            DO J=JSTA,JEND
            DO I=1,IM
              EGRID2(I,J)=SPL(LP)
            ENDDO
            ENDDO
C
            CALL CALRH2(EGRID2,TPRS(1,1,LP),QPRS(1,1,LP),ICE,EGRID1
     1,                 IM,JM)
            CALL E2OUT(017,000,EGRID1,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            CALL SCLFLD(GRID1,H100,IMOUT,JMOUT)
            CALL BOUND(GRID1,H1,H100,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(017),LP,GRID1,IMOUT,JMOUT)
          ENDIF
        ENDIF
C     
C***  DEWPOINT TEMPERATURE.
C
        IF(IGET(015).GT.0)THEN
          IF(LVLS(LP,IGET(015)).GT.0)THEN
!$omp  parallel do
            DO J=JSTA,JEND
            DO I=1,IM
              EGRID2(I,J)=SPL(LP)
            ENDDO
            ENDDO
C
            CALL CALDWP2(EGRID2,QPRS(1,1,LP),EGRID1,TPRS(1,1,LP))
            CALL E2OUT(015,000,EGRID1,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(015),LP,GRID1,IMOUT,JMOUT)
          ENDIF
        ENDIF
C     
C---------------------------------------------------------------------
C***  CALCULATE 1000MB GEOPOTENTIALS CONSISTENT WITH SLP OBTAINED 
C***  FROM THE MESINGER OR NWS SHUELL SLP REDUCTION.
C---------------------------------------------------------------------
C     
C***  FROM MESINGER SLP
C
cc        IF(IGET(023).GT.0.AND.ABS(SPL(LP)-1.E5).LE.0.1)THEN
cc          ALPTH=ALOG(1.E5)
cc!$omp  parallel do private(i,j)
cc          DO J=JSTA,JEND
cc         DO I=1,IM
cc            PSLPIJ=PSLP(I,J)
cc            ALPSL=ALOG(PSLPIJ)
cc            PSFC=PD(I,J)+PT
cc            IEND=IM-MOD(J+1,2)
cc            IF(ABS(PSLPIJ-PSFC).LT.5.E2.OR.
cc     1        ((I.EQ.1.OR.I.LE.IEND.OR.J.LE.2.OR.J.GE.JM-1)
cc     2         .AND.SM(I,J).GT.0.5))THEN
cc              FPRS(I,J,LP)=R*TPRS(I,J,LSL)*(ALPSL-ALPTH)
cc            ELSE
cc              FPRS(I,J,LP)=FIS(I,J)/(ALPSL-ALOG(PD(I,J)+PT))*
cc     1                              (ALPSL-ALPTH)
cc            ENDIF
cc            Z1000(I,J)=FPRS(I,J,LP)*GI
cc          ENDDO
cc          ENDDO
C     
C***  FROM NWS SHUELL SLP. NGMSLP2 COMPUTES 1000MB GEOPOTENTIAL.
C
cc        ELSEIF(IGET(023).LE.0.AND.ABS(SPL(LP)-1.E5).LE.0.1)THEN
cc!$omp  parallel do private(i,j)
cc          DO J=JSTA,JEND
cc          DO I=1,IM
cc            FPRS(I,J,LP)=Z1000(I,J)*G
cc          ENDDO
cc          ENDDO
cc        ENDIF
C
C***  OUTPUT GEOPOTENTIAL (SCALE BY GI)
C
        IF(IGET(012).GT.0)THEN
          IF(LVLS(LP,IGET(012)).GT.0)THEN
!$omp  parallel do
            DO J=JSTA,JEND
            DO I=1,IM
              EGRID1(I,J)=FPRS(I,J,LP)
            ENDDO
            ENDDO
C
            if(me.eq.0)print*,'H at (3,3,lp)= ',lp,fprs(3,3,lp)
            CALL E2OUT(012,000,EGRID1,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(012),LP,GRID1,IMOUT,JMOUT)
          ENDIF
        ENDIF
C
  680  CONTINUE
C
        IOALL=.TRUE.
C
C***  ENDIF FOR IF TEST SEEING IF WE WANT ANY OTHER VARIABLES
C
      ENDIF
C     
C     END OF ROUTINE.
C
cc      deallocate(tprs) 
cc      deallocate(oprs)
cc      deallocate(qprs)
cc      deallocate(fprs)
cc      deallocate(qcprs)
cc      deallocate(q2prs)
cc      deallocate(uprs)
cc      deallocate(vprs)
      RETURN
      END
C
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SPLINEF(NOLD,XOLD,YOLD,Y2,NNEW,XNEW,YNEW,P,Q)           
C                                                                       
C     ******************************************************************
C     *                                                                *
C     *  THIS IS A ONE-DIMENSIONAL CUBIC SPLINE FITTING ROUTINE        *
C     *  PROGRAMED FOR A SMALL SCALAR MACHINE.                         *
C     *                                                                *
C     *  PROGRAMER[ Z. JANJIC, YUGOSLAV FED. HYDROMET. INST., BEOGRAD  *
C     *                                                                *
C     *                                                                *
C     *                                                                *
C     *  NOLD - NUMBER OF GIVEN VALUES OF THE FUNCTION.  MUST BE GE 3. *
C     *  XOLD - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE     *
C     *         FUNCTION ARE GIVEN.  MUST BE IN ASCENDING ORDER.       *
C     *  YOLD - THE GIVEN VALUES OF THE FUNCTION AT THE POINTS XOLD.   *
C     *  Y2   - THE SECOND DERIVATIVES AT THE POINTS XOLD.  IF NATURAL *
C     *         SPLINE IS FITTED Y2(1)=0. AND Y2(NOLD)=0. MUST BE      *
C     *         SPECIFIED.                                             *
C     *  NNEW - NUMBER OF VALUES OF THE FUNCTION TO BE CALCULATED.     *
C     *  XNEW - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE     *
C     *         FUNCTION ARE CALCULATED.  XNEW(K) MUST BE GE XOLD(1)   *
C     *         AND LE XOLD(NOLD).                                     *
C     *  YNEW - THE VALUES OF THE FUNCTION TO BE CALCULATED.           *
C     *  P, Q - AUXILIARY VECTORS OF THE LENGTH NOLD-2.                *
C     *                                                                *
C     ******************************************************************
C                                                                       
                             D I M E N S I O N                          
     2 XOLD(NOLD),YOLD(NOLD),Y2(NOLD),P(NOLD),Q(NOLD)
     3,XNEW(NNEW),YNEW(NNEW)
C-----------------------------------------------------------------------
      NOLDM1=NOLD-1                                                     
C                                                                       
      DXL=XOLD(2)-XOLD(1)                                               
      DXR=XOLD(3)-XOLD(2)                                               
      DYDXL=(YOLD(2)-YOLD(1))/DXL                                       
      DYDXR=(YOLD(3)-YOLD(2))/DXR                                       
      RTDXC=.5/(DXL+DXR)                                                
C                                                                       
      P(1)= RTDXC*(6.*(DYDXR-DYDXL)-DXL*Y2(1))                          
      Q(1)=-RTDXC*DXR                                                   
C                                                                       
      IF(NOLD.EQ.3) GO TO 700                                           
C-----------------------------------------------------------------------
      K=3                                                               
C                                                                       
 100  DXL=DXR                                                           
      DYDXL=DYDXR                                                       
      DXR=XOLD(K+1)-XOLD(K)  
c      if(i.eq.120.and.j.eq.279)then
c       print*,'in spline',k,XOLD(K),XOLD(K+1),dxr 
c      end if                                           
      DYDXR=(YOLD(K+1)-YOLD(K))/DXR                                     
      DXC=DXL+DXR                                                       
      DEN=1./(DXL*Q(K-2)+DXC+DXC)                                       
C                                                                       
      P(K-1)= DEN*(6.*(DYDXR-DYDXL)-DXL*P(K-2))                         
      Q(K-1)=-DEN*DXR                                                   
C                                                                       
      K=K+1                                                             
      IF(K.LT.NOLD) GO TO 100                                           
C-----------------------------------------------------------------------
 700  K=NOLDM1                                                          
C                                                                       
 200  Y2(K)=P(K-1)+Q(K-1)*Y2(K+1)                                       
C                                                                       
      K=K-1                                                             
      IF(K.GT.1) GO TO 200                                              
C-----------------------------------------------------------------------
      K1=1                                                              
C                                                                       
 300  XK=XNEW(K1)                                                       
C                                                                       
      DO 400 K2=2,NOLD                                                  
      IF(XOLD(K2).LE.XK) GO TO 400                                      
      KOLD=K2-1                                                         
      GO TO 450                                                         
 400  CONTINUE                                                          
      YNEW(K1)=YOLD(NOLD)                                               
      GO TO 600                                                         
C                                                                       
 450  IF(K1.EQ.1)   GO TO 500                                           
      IF(K.EQ.KOLD) GO TO 550                                           
C                                                                       
 500  K=KOLD                                                            
C                                                                       
      Y2K=Y2(K)                                                         
      Y2KP1=Y2(K+1)                                                     
      DX=XOLD(K+1)-XOLD(K)                                              
      RDX=1./DX                                                         
C                                                                       
      AK=.1666667*RDX*(Y2KP1-Y2K)                                       
      BK=.5*Y2K                                                         
      CK=RDX*(YOLD(K+1)-YOLD(K))-.1666667*DX*(Y2KP1+Y2K+Y2K)            
C                                                                       
 550  X=XK-XOLD(K)                                                      
      XSQ=X*X                                                           
C                                                                       
      YNEW(K1)=AK*XSQ*X+BK*XSQ+CK*X+YOLD(K)                             
C                                                                       
 600  K1=K1+1                                                           
      IF(K1.LE.NNEW) GO TO 300       
C-----------------------------------------------------------------------
      RETURN                                                            
      END
