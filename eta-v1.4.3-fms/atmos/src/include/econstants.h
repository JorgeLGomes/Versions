c
c *** Universal constants.
c
      real*4 g,r,cp,kappa,gamma,prf0,t0
      parameter (g=9.80
     .          ,r=287.04
     .          ,cp=1004.6
     .          ,kappa=r/cp
     .          ,gamma=0.0065
     .          ,prf0=101325.
     .          ,t0=288.)
c
c *** Other constants.
c
      real*4 pi,dtr
      parameter (pi=3.141592654
     .          ,dtr=pi/180.)
c
c *** Stereographic map constants.
c
c     Unit of distance of the image surface x1,y1 and x,y systems is
c        the shorter side of the printing cell (2.54/10. cm)
c     x1p,y1p are the north pole coordinates.
c     dp30 is the distance, on the overlay map, between the
c        north pole and 30 deg N latitude circle, in cm.
c        dp30 = 18.02 * 2.54
c     cmld is the central map meridian longitude, in deg 
c        (which = tlm0d, and is set indirectly through the namelist).
c
      real*4 dp30,x1p,y1p
      integer*4 ixm,iym
      parameter (dp30=45.78
     .          ,x1p=57.,y1p=172.
     .          ,ixm=65,iym=82)
