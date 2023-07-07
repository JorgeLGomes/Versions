c===============================================================================
c
      subroutine table(ptbl,ttbl,pt
     .                ,rdq,rdth,rdp,rdthe,pl,thl,qs0,sqs,sthe,the0)
c
c     implicit none
c
c *** Generate values for look-up tables used in convection.
c
      integer*4 itb,jtb
      parameter (itb=76,jtb=134)
      real*4 thh,ph,pq0,a1,a2,a3,a4,r,cp,eliwv,eps
Cmp      parameter (thh=350.,ph=105000.
C	changed 3/8/99 to make more like operation grdeta stuff
C	appeared to have no change on underflow problem
      parameter (thh=365.,ph=105000.
     .          ,pq0=379.90516
     .          ,a1=610.78,a2=17.2693882,a3=273.16,a4=35.86
     .          ,r=287.04,cp=1004.6,eliwv=2.683e6,eps=1.e-10)
c
      real*4 ptbl(itb,jtb),ttbl(jtb,itb),qsold (jtb),pold(jtb)
     .      ,qs0 (jtb),sqs   (jtb),qsnew(jtb)
     .      ,y2p (jtb),app   (jtb),aqp  (jtb),pnew(jtb)
     .      ,told(jtb),theold(jtb),the0 (itb),sthe(itb)
     .      ,y2t (jtb),thenew(jtb),apt  (jtb),aqt (jtb),tnew(jtb)
c_______________________________________________________________________________
c
c *** Coarse look-up table for saturation point.
c
      kthm=jtb
      kpm=itb
      kthm1=kthm-1
      kpm1=kpm-1
c
      pl=pt
c
      dth=(thh-thl)/real(kthm-1)
      dp =(ph -pl )/real(kpm -1)
c
      rdth=1./dth
      rdp=1./dp
      rdq=kpm-1
c
      th=thl-dth
c
      do kth=1,kthm
         th=th+dth
         p=pl-dp
         do kp=1,kpm
            p=p+dp
            ape=(100000./p)**(r/cp)
            qsold(kp)=pq0/p*exp(a2*(th-a3*ape)/(th-a4*ape))
            pold(kp)=p
         enddo
c
         qs0k=qsold(1)
         sqsk=qsold(kpm)-qsold(1)
         qsold(1  )=0.
         qsold(kpm)=1.
c
         do kp=2,kpm1
            qsold(kp)=(qsold(kp)-qs0k)/sqsk
c
c ********* Fix due to cyber half prec. limitation.
c
            if (qsold(kp)-qsold(kp-1) .lt. eps) 
     .         qsold(kp)=qsold(kp-1)+eps
c
c ********* End fix.
c
         enddo
c
         qs0(kth)=qs0k
         sqs(kth)=sqsk
         qsnew(1  )=0.
         qsnew(kpm)=1.
         dqs=1./real(kpm-1)
c
         do kp=2,kpm1
            qsnew(kp)=qsnew(kp-1)+dqs
         enddo
c
         y2p(1   )=0.
         y2p(kpm )=0.
c
         call spline(jtb,kpm,qsold,pold,y2p,kpm,qsnew,pnew,app,aqp)
c
         do kp=1,kpm
            ptbl(kp,kth)=pnew(kp)
         enddo
c
      enddo
c
c *** Coarse look-up table for t(p) from constant the.
c
      p=pl-dp
      do kp=1,kpm
         p=p+dp
         th=thl-dth
         do kth=1,kthm
            th=th+dth
            ape=(100000./p)**(r/cp)
            qs=pq0/p*exp(a2*(th-a3*ape)/(th-a4*ape))
            told(kth)=th/ape
            theold(kth)=th*exp(eliwv*qs/(cp*told(kth)))
         enddo
c
         the0k=theold(1)
         sthek=theold(kthm)-theold(1)
         theold(1   )=0.
         theold(kthm)=1.
c
         do kth=2,kthm1
            theold(kth)=(theold(kth)-the0k)/sthek
c
            if (theold(kth)-theold(kth-1) .lt. eps)
     .         theold(kth)=theold(kth-1)+eps
c
         enddo
c
         the0(kp)=the0k
         sthe(kp)=sthek
         thenew(1  )=0.
         thenew(kthm)=1.
         dthe=1./real(kthm-1)
         rdthe=1./dthe
c
         do kth=2,kthm1
            thenew(kth)=thenew(kth-1)+dthe
         enddo
c
         y2t(1   )=0.
         y2t(kthm)=0.
c
         call spline(jtb,kthm,theold,told,y2t,kthm,thenew,tnew,apt,aqt)
c
         do kth=1,kthm
            ttbl(kth,kp)=tnew(kth)
         enddo
c
      enddo
c
      return
      end
