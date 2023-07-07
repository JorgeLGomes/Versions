      subroutine rtll(tlm,tph,tlm0d,dtr,ctph0,stph0,almd,aphd)
!     ****************************************************************
!     *                                                              *
!     *  programer: z. janjic, shmz, feb. 1981                       *
!     *  ammended:  z. janjic, ncep, jan. 1996                       *
!     *                                                              *
!     *  transformation from rotated lat-lon to lat-lon coordinates  *
!     ****************************************************************
!     ****************************************************************
!     *  tlm   - transformed longitude, rad.                         *
!     *  tph   - transformed latitude, rad.                          *
!     *  tlm0d - the angle of rotation of the transformed lat-lon    *
!     *          system in the longitudinal direction, degs          *
!     *  ctph0 - cos(tph0), tph0 is the angle of rotation of the     *
!     *          transformed lat-lon systemn in the latitudinal      *
!     *          direction, precomputed                              *
!     *  stph0 - sin(tph0), tph0 is the angle of rotation of the     *
!     *          transformed lat-lon systemn in the latitudinal      *
!     *          direction, precomputed                              *
!     *  almd  - geographical longitude, degs, range -180.,180       *
!     *  aphd  - geographical latitude,  degs, range - 90., 90.,     *
!     *          poles are singular                                  *
!     ****************************************************************
!
      stlm=sin(tlm)
      ctlm=cos(tlm)
      stph=sin(tph)
      ctph=cos(tph)
!
      sph=ctph0*stph+stph0*ctph*ctlm
      aph=asin(sph)
      aphd=aph/dtr
      anum=ctph*stlm
      denom=(ctlm*ctph-stph0*sph)/ctph0
      relm=atan2(anum,denom)
      almd=relm/dtr+tlm0d
!
      if(almd.gt. 180.)    almd=almd-360.
      if(almd.lt.-180.)    almd=almd+360.
!
      return
      end
