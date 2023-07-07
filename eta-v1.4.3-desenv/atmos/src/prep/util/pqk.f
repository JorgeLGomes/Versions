      subroutine pqk(almd,aphd,inside,p,q,k)
c
c *** Routine to transform alm, aph (longitude, latitude) into
c        transformed (rotated) longitude, latitude, check whether the 
c        values obtained are within the model region, and, if they 
c        are, calculate values of p,q (coordinates within a square
c        formed by connecting four neighboring height points) and k.
c
c *** Grid constants:
c
c        wbd,sbd - Western and southern boundaries on ll or tll grid
c                  in degrees                                          
c        tlm0d,tph0d - Angles of rotation of the ll coordinate system 
c                  in the direction of lambda and phi respectively   
c                  in order to obtain the tll coordinate system     
c        dlmd,dphd - Mesh sides in degrees                        
c
      implicit none
c
      include 'ecommons.h'
c
      real*4 almd,aphd,p,q
     .      ,rim1,rdlmd,rdphd
     .      ,ebd,nbd
     .      ,rlmd,srlm,crlm,saph,caph,cc
     .      ,anum,denom
     .      ,x1,x2,y1,y2,d2r


	parameter(d2r=.017453293)
c
      integer*4 k,i2,j2,i3,j3,jr
c
      logical inside
c_______________________________________________________________________________
c
c *** Define constants.
c
      rim1=im1
      rdlmd=1./dlmd
      rdphd=1./dphd
      ebd=wbd+2*(im1)*dlmd
      nbd=sbd+  (jm1)*dphd
c
c *** almd, aphd is lon, lat in degrees.
c *** Transform into rotated longitude, latitude.
c
      rlmd=almd-tlm0d

C      srlm=sind(rlmd)
C      crlm=cosd(rlmd)
C      saph=sind(aphd)
C      caph=cosd(aphd)

      srlm=sin(rlmd*d2r)
      crlm=cos(rlmd*d2r)
      saph=sin(aphd*d2r)
      caph=cos(aphd*d2r)

      cc=crlm*caph
c
      anum=srlm*caph
      denom=ctph0*cc+stph0*saph
c
C      almd=atan2d(anum,denom)
C      aphd=asind(ctph0*saph-stph0*cc)

      almd=atan2(anum,denom)
	almd=almd*(1./d2r)
      aphd=asin(ctph0*saph-stph0*cc)
	aphd=aphd*(1./d2r)
c
c *** Is almd, aphd inside the (one line reduced) model domain?
c
      if (almd .le. wbd .or. almd .ge. ebd .or.
     .    aphd .le. sbd .or. aphd .ge. nbd) then
             inside=.false.
             return
      endif
      inside=.true.
c
c *** x1, y1 is a coordinate system with dlmd, dphd as length units.
c
      x1=(almd-wbd)*rdlmd
      y1=(aphd-sbd)*rdphd
c
c *** x2, y2 rotated for +45 deg. & translated for im1.
c
      x2=0.5*( x1+y1)
      y2=0.5*(-x1+y1)+rim1
c
c *** i2, j2 are coordinates of southern points of grid boxes.
c
      i2=int(x2)
      j2=int(y2)
c
c *** Remaining parameters need to know the position
c        of the transformed point within the model region.
c
      p=x2-i2
      q=y2-j2
c
c *** Index k corresponds to i2, j2.
c
      jr=j2-im1
      i3=i2-jr
      j3=i2+jr
      k=min(imjm,max(1,j3*im-j3/2+(i3+2)/2))
c
      return
      end
