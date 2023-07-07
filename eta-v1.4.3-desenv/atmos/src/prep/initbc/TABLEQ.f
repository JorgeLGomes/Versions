c===============================================================================
c
      subroutine tableq(ttblq,rdp,rdthe,pl,thl,sthe,the0)
c
c *** Generate values for finer look-up tables used in convection.
c
      parameter (itb=152,jtb=440)
      parameter (thh=325.,ph=105000.
     .          ,pq0=379.90516
     .          ,a1=610.78,a2=17.2693882,a3=273.16,a4=35.86
     .          ,r=287.04,cp=1004.6,eliwv=2.683e6,eps=1.e-9)
Cmp	change eps from 1.e-10 to 1.e-9 (make more like operational)
c
      real*4 ttblq(jtb,itb)
     .      ,told (jtb),theold(jtb),the0(itb),sthe(itb)
     .      ,y2t  (jtb),thenew(jtb),apt (jtb),aqt (jtb),tnew(jtb)
c_______________________________________________________________________________
c
c *** Coarse look-up table for saturation point.
c
      kthm=jtb
      kpm=itb
      kthm1=kthm-1
      kpm1=kpm-1
c
      dth=(thh-thl)/real(kthm-1)
      dp =(ph -pl )/real(kpm -1)
c
      rdp=1./dp
      th=thl-dth
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
            if (theold(kth)-theold(kth-1) .lt. eps)
     .         theold(kth)=theold(kth-1)+eps
c
         enddo
c
         the0(kp)=the0k
         sthe(kp)=sthek
c
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
            ttblq(kth,kp)=tnew(kth)
         enddo
c
      enddo
c
      return
      end
c

