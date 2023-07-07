!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CALCSHR                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE calcshr(nrad,lm,zp,sigma,                                       &
          p_pa,t3d,u3d,v3d,cape,                                       &
          shr37,ustrm,vstrm,srlfl,srmfl,helicity,brn,brnu,brnu2km,blcon, &
          tem2,tem3,nlev,estacao,hora,dia,mes,ano,strm_opt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate various wind shear parameters useful for gauging
!  the potential for severe storms.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  February, 1994   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  Keith Brewster
!  Added storm-relative flows, general clean-up to
!  meet ARPS coding standards.
!
!  3/22/1996  (Keith Brewster)
!  Fixed some bugs, added smoothing at the end.
!
!  06/20/2000 (Eric Kemp and Keith Brewster)
!  Changed BRN Shear to be the denominator of BRN, instead of wind
!  speed, now has units of speed squared.
!
!  04/19/2004 (Ernani L. Nascimento)
!  Modified the code to fit in índices_severos3.f90, and adapted
!  the computation of expected storm-motion for the Southern Hemisphere.
!
!-----------------------------------------------------------------------
!
!  Calculates some of the shear related variables from the
!  ARPS 3D wind field.
!
!  shr37         Magnitude of wind shear between 3 and 7 km AGL
!  ustrm,vstrm   Estimated storm motion (modified from Bob Johns)
!  srlfl         Low-level storm-relative wind
!  srmfl         Mid-level storm-relative wind
!  helicity      Helicity, storm relative
!  brn           Bulk Richardson Number (Weisman and Klemp)
!  brnu          Shear parameter of BRN, "U"
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
 INTEGER :: nrad,lm,strm_opt
!
!-----------------------------------------------------------------------
!
!  "1-D" input variables (Ernani L. Nascimento)
!
!-----------------------------------------------------------------------
!
 INTEGER :: nlev(nrad),estacao(nrad),hora(nrad),dia(nrad),mes(nrad),   &
            ano(nrad)
 REAL :: sigma(nrad)
 REAL :: cape(nrad)
!
!-----------------------------------------------------------------------
!
!  "2-D" input variables (Ernani L. Nascimento)
!
!-----------------------------------------------------------------------
!
 INTEGER :: zp(lm,nrad)
 REAL :: p_pa(lm,nrad)   ! Pressure in Pascals
 REAL :: t3d(lm,nrad)    ! Temperature in Kelvin
 REAL :: u3d(lm,nrad)
 REAL :: v3d(lm,nrad)

!
!-----------------------------------------------------------------------
!
!  Output variables
!
!-----------------------------------------------------------------------
!
 REAL :: shr37(nrad)      ! 7km - 3km wind shear
 REAL :: ustrm(nrad)
 REAL :: vstrm(nrad)  ! Estimated storm motion (Bob Johns)
 REAL :: srlfl(nrad)
 REAL :: srmfl(nrad)
 REAL :: helicity(nrad)   ! Helicity, storm relative
 REAL :: brn(nrad)        ! Bulk Richardson Number (Weisman and Klemp)
 REAL :: brnu(nrad)       ! Shear parameter of BRN, "U"
 REAL :: brnu2km(nrad)    ! 2km shear parameter of BRN
 REAL :: blcon(nrad)
!
!-----------------------------------------------------------------------
!
!  Temporary variables
!
!-----------------------------------------------------------------------
!

 REAL :: tem2(nrad)
 REAL :: tem3(nrad)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
!  INTEGER :: imid,jmid
 REAL :: elev(nrad),hgt3d(lm,nrad)
 INTEGER :: i,j,k,ksfc,k2km,k3km
 REAL :: h3km,u3km,v3km,h7km,u7km,v7km,u2,v2
 REAL :: p2km,t2km,h2km,u2km,v2km,p9km,t9km,h9km,u9km,v9km
 REAL :: p500m,t500m,h500m,u500m,v500m,p6km,t6km,h6km,u6km,v6km
 REAL :: sumu,sumv,sump,sumh,wlow,whigh,dx,dy,dz,dp,dx2,dy2
 REAL :: rhohi,rholo,rhoinv,arg,new_ustrm,new_vstrm
 REAL :: dirmean,spmean,ushr,vshr,ddir,perc,obs_ddstrm,obs_ffstrm
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
!
! write(*,*) 'Entrada da subrotina calcshr'
 DO k=1,nrad
  if (nlev(k)==0) then
   shr37(k)=-999.9
   ustrm(k)=-999.9
   vstrm(k)=-999.9
   srlfl(k)=-999.9
   srmfl(k)=-999.9
   helicity(k)=-999.9
   brn(k)=-999.9
   brnu(k)=-999.9
   brnu2km(k)=-999.9
   blcon(k)=-999.9
  else
!   print*, estacao(k),hora(k),dia(k),mes(k),ano(k)
   j=nlev(k)
   if ((p_pa(j,k)*0.01)>300.0) then
    cycle
   endif

   arg=0.0
   elev(k)=real(zp(1,k))
!    write(*,*) 'elev(k)= ',elev(k),' para k= ',k
   DO j=1,nlev(k)
     hgt3d(j,k)=real(zp(j,k))
   ENDDO
!   write(*,*) 'PASSEI AQUI 3, com k= ',k
!
!-----------------------------------------------------------------------
!
!  Find mass weighted mean wind in first 500m
!
!-----------------------------------------------------------------------
!
     sumu=0.
     sumv=0.
     sump=0.
     h500m=elev(k)+500.
     DO j=2,nlev(k)
       IF( hgt3d(j,k) < h500m ) THEN
         dp=p_pa(j-1,k)-p_pa(j,k)
         rhohi=p_pa(j,k)/t3d(j,k)
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u3d(j,k))
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v3d(j,k))
         sump=sump+dp
       ELSE
         dz=hgt3d(j,k)-hgt3d(j-1,k)
         wlow=(hgt3d(j,k)-h500m)/dz
         u500m=u3d(j,k)*(1.-wlow) + u3d(j-1,k)*wlow
         v500m=v3d(j,k)*(1.-wlow) + v3d(j-1,k)*wlow
         p500m=p_pa(j,k)*(1.-wlow) + p_pa(j-1,k)*wlow
         t500m=t3d(j,k)*(1.-wlow) + t3d(j-1,k)*wlow
         dp=p_pa(j-1,k)-p500m
         rhohi=p500m/t500m
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u500m)
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v500m)
         sump=sump+dp
!          print *, ' sumu,sumv,sump = ',sumu,sumv,sump,'para k= ',k
         EXIT
       END IF
     END DO
!      121    CONTINUE
     u500m=sumu/sump
     v500m=sumv/sump
!     write(*,*) 'PASSEI AQUI 4, com k= ',k
!      write(*,*) '(u500m,v500m)= ',u500m,v500m,' para k= ',k
!
!-----------------------------------------------------------------------
!
!  Find mass weighted mean wind sfc-2km AGL
!
!-----------------------------------------------------------------------
!
     sumu=0.
     sumv=0.
     sump=0.
     h2km=elev(k)+2000.
     DO j=2,nlev(k)
       IF( hgt3d(j,k) < h2km ) THEN
         dp=p_pa(j-1,k)-p_pa(j,k)
         rhohi=p_pa(j,k)/t3d(j,k)
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u3d(j,k))
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v3d(j,k))
         sump=sump+dp
       ELSE
         dz=hgt3d(j,k)-hgt3d(j-1,k)
         wlow=(hgt3d(j,k)-h2km)/dz
         u2km=u3d(j,k)*(1.-wlow) + u3d(j-1,k)*wlow
         v2km=v3d(j,k)*(1.-wlow) + v3d(j-1,k)*wlow
         p2km=p_pa(j,k)*(1.-wlow) + p_pa(j-1,k)*wlow
         t2km=t3d(j,k)*(1.-wlow) + t3d(j-1,k)*wlow
         dp=p_pa(j-1,k)-p2km
         rhohi=p2km/t2km
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u2km)
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v2km)
         sump=sump+dp
         EXIT
       END IF
     END DO
!      141    CONTINUE
!     write(*,*) 'PASSEI AQUI 4.5, com k= ',k
     u2km=sumu/sump
     v2km=sumv/sump
!     write(*,*) 'u2km,v2km: ',u2km,v2km
     tem2(k)=u2km
     tem3(k)=v2km
!     write(*,*) 'tem2(k),tem3(k): ',tem2(k),tem3(k)
!     write(*,*) 'PASSEI AQUI 5, com k= ',k
!
!-----------------------------------------------------------------------
!
!  Find mass weighted mean wind sfc-6km AGL
!
!-----------------------------------------------------------------------
!
     sumu=0.
     sumv=0.
     sump=0.
     h6km=elev(k)+6000.
     DO j=2,nlev(k)
       IF( hgt3d(j,k) < h6km ) THEN
         dp=p_pa(j-1,k)-p_pa(j,k)
         rhohi=p_pa(j,k)/t3d(j,k)
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u3d(j,k))
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v3d(j,k))
         sump=sump+dp
       ELSE
         dz=hgt3d(j,k)-hgt3d(j-1,k)
         wlow=(hgt3d(j,k)-h6km)/dz
         u6km=u3d(j,k)*(1.-wlow) + u3d(j-1,k)*wlow
         v6km=v3d(j,k)*(1.-wlow) + v3d(j-1,k)*wlow
         p6km=p_pa(j,k)*(1.-wlow) + p_pa(j-1,k)*wlow
         t6km=t3d(j,k)*(1.-wlow) + t3d(j-1,k)*wlow
         dp=p_pa(j-1,k)-p6km
         rhohi=p6km/t6km
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u6km)
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v6km)
         sump=sump+dp
         EXIT
       END IF
     END DO
     u6km=sumu/sump
     v6km=sumv/sump
!     write(*,*) 'PASSEI AQUI 6, com k= ',k,' e (u6km,v6km)= ',u6km,v6km
!
!-----------------------------------------------------------------------
!
!  Find mass weighted mean wind 2km-9km AGL
!
!-----------------------------------------------------------------------
!
     sumu=0.
     sumv=0.
     sump=0.
     h9km=elev(k)+9000.
     DO j=2,nlev(k)
       IF( hgt3d(j,k) > h2km ) EXIT
     END DO
!      181   CONTINUE
     k2km=j
     dz=hgt3d(j,k)-hgt3d(j-1,k)
     wlow=(hgt3d(j,k)-h2km)/dz
     u2=u3d(j,k)*(1.-wlow) + u3d(j-1,k)*wlow
     v2=v3d(j,k)*(1.-wlow) + v3d(j-1,k)*wlow
     p2km=p_pa(j,k)*(1.-wlow) + p_pa(j-1,k)*wlow
     t2km=t3d(j,k)*(1.-wlow) + t3d(j-1,k)*wlow
     dp=p2km-p_pa(j,k)
     rholo=p2km/t2km
     rhohi=p_pa(j,k)/t3d(j,k)
     rhoinv=1./(rhohi+rholo)
     sumu=sumu+dp*rhoinv*(rholo*u2+rhohi*u3d(j,k))
     sumv=sumv+dp*rhoinv*(rholo*v2+rhohi*v3d(j,k))
     sump=sump+dp
     DO j=k2km+1,nlev(k)
       IF( hgt3d(j,k) < h9km ) THEN
         dp=p_pa(j-1,k)-p_pa(j,k)
         rhohi=p_pa(j,k)/t3d(j,k)
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u3d(j,k))
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v3d(j,k))
         sump=sump+dp
       ELSE
         dz=hgt3d(j,k)-hgt3d(j-1,k)
         wlow=(hgt3d(j,k)-h9km)/dz
         u9km=u3d(j,k)*(1.-wlow) + u3d(j-1,k)*wlow
         v9km=v3d(j,k)*(1.-wlow) + v3d(j-1,k)*wlow
         p9km=p_pa(j,k)*(1.-wlow) + p_pa(j-1,k)*wlow
         t9km=t3d(j,k)*(1.-wlow) + t3d(j-1,k)*wlow
         dp=p_pa(j-1,k)-p9km
         rhohi=p9km/t9km
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u9km)
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v9km)
         sump=sump+dp
         EXIT
       END IF
     END DO
!      191    CONTINUE
     u9km=sumu/sump
     v9km=sumv/sump
!     write(*,*) 'PASSEI AQUI 7, com k= ',k
!
!-----------------------------------------------------------------------
!
!  Storm motion estimation
!  From Davies and Johns, 1993
!  "Some wind and instability parameters associated With
!  strong and violent tornadoes."
!  AGU Monograph 79, The Tornado...(Page 575)
!
!  Becuase of the discontinuity produced by that method
!  at the 15.5 m/s cutoff, their rules have been modified
!  to provide a gradual transition, and accomodate all the
!  data they mention in the article.
!
!  (04/19/2004) Modified by Ernani L. Nascimento. Adapting for the
!  Southern Hemisphere.
!
!-----------------------------------------------------------------------
!
     CALL uv2ddff(u6km,v6km,dirmean,spmean)
!      write(*,*) 'PASSAGEM 5.1, com k= ',k
!      write(*,*) 'PASSAGEM 5.1, com k= ',k,' e spmean= ',spmean
!      write(*,*) 'PASSAGEM 5.1, com k= ',k,' e dirmean= ',dirmean
     IF(spmean >= 20.0) THEN
!        write(*,*) 'Hei, passei aqui!'
       dirmean=dirmean-18.
       IF(dirmean <= 0.) dirmean=360.+dirmean
       spmean=spmean*0.89
     ELSE IF (spmean > 8.0) THEN
!        write(*,*) 'Não! Passei foi aqui!'
       whigh=(spmean - 8.0)/12.
       wlow =1.-whigh
       ddir=wlow*32.0 + whigh*18.0
       perc=wlow*0.75 + whigh*0.89
       dirmean=dirmean-ddir
       IF(dirmean <= 0.) dirmean=360.+dirmean
       spmean=spmean*perc
     ELSE
!        write(*,*) 'Na verdade, passei foi aqui!'
       dirmean=dirmean-32.
       IF(dirmean <= 0.) dirmean=360.+dirmean
       spmean=spmean*0.75
     END IF
!      write(*,*) 'PASSAGEM 5.2, com k= ',k
     arg = (dirmean * (3.141592654/180.))
     ustrm(k) = -spmean * sin(arg)
     vstrm(k) = -spmean * cos(arg)
! Utilizando movimento estimado da célula da esquerda em 9 de outubro de 2003
     IF (strm_opt == 1) then
      IF ((estacao(k)==83827).and.(hora(k)==00).and.(dia(k)==09).and. &
          (mes(k)==10).and.(ano(k)==2003)) THEN
       obs_ddstrm=270.0
       obs_ffstrm=10.7
       CALL ddff2uv(obs_ddstrm,obs_ffstrm,new_ustrm,new_vstrm)
       ustrm(k)=new_ustrm
       vstrm(k)=new_vstrm
      ENDIF
! Utilizando movimento estimado da célula da célula de 24 de maio de 2005
!       IF ((estacao(k)==99999).and.(hora(k)==18).and.(dia(k)==24).and. &
!           (mes(k)==05).and.(ano(k)==2005)) THEN
!        obs_ddstrm=297.0
!        obs_ffstrm=17.8
!        CALL ddff2uv(obs_ddstrm,obs_ffstrm,new_ustrm,new_vstrm)
!        ustrm(k)=new_ustrm
!        vstrm(k)=new_vstrm
!       ENDIF
     ENDIF
!      CALL ddff2uv(dirmean,spmean,ustrm(k),vstrm(k),5)
!      write(*,*) 'PASSAGEM 6, com k= ',k
!
!-----------------------------------------------------------------------
!
!  Storm-relative low-level flow
!
!-----------------------------------------------------------------------
!
     srlfl(k)=SQRT((ustrm(k)-u2km)*(ustrm(k)-u2km) +             &
                     (vstrm(k)-v2km)*(vstrm(k)-v2km))
!
!-----------------------------------------------------------------------
!
!  Storm relative mid-level flow
!
!-----------------------------------------------------------------------
!
     srmfl(k)=SQRT((ustrm(k)-u9km)*(ustrm(k)-u9km) +             &
                     (vstrm(k)-v9km)*(vstrm(k)-v9km))
!
!-----------------------------------------------------------------------
!
!  Shear parameter for Bulk Richardson number
!
!-----------------------------------------------------------------------
!
!    print *, ' density-weight mean 0-500 m ',u500m,v500m
!    print *, ' density-weight mean 0-6  km ',u6km,v6km
!
     brnu(k)=0.5*( (u6km-u500m)*(u6km-u500m) +                       &
                     (v6km-v500m)*(v6km-v500m) )
! Added the 2km BRNSHR (Ernani L. Nascimento, 03/feb/2005
     brnu2km(k)=0.5*( (u2km-u500m)*(u2km-u500m) +                    &
                     (v2km-v500m)*(v2km-v500m) )
!     write(*,*) 'PASSAGEM 7, com k= ',k
!
!-----------------------------------------------------------------------
!
!  Bulk Richardson number
!  A limit of 200 is imposed, since this could
!  go to inifinity.
!
!-----------------------------------------------------------------------
!
     IF(brnu(k) > 0.) THEN
       brn(k)=cape(k)/brnu(k)
       brn(k)=AMIN1(brn(k),200.)
     ELSE
       brn(k)=200.
     END IF
!    write(*,*) 'PASSEI AQUI 8, com k= ',k
!
!
!-----------------------------------------------------------------------
!
!  Calculate Helicity and 3km to 7km shear
!  since both involve the 3km wind.
!
!  For more efficient computation the Helicity is
!  computed for zero storm motion and the storm
!  motion is accounted for by adding a term at the end.
!  This is mathematically equivalent to accounting
!  for the storm motion at each level.
!
!-----------------------------------------------------------------------
!
     h3km=elev(k)+3000.
     h7km=elev(k)+7000.
!
!-----------------------------------------------------------------------
!
!  Find level just above 3 km AGL
!  Note, it is assumed here that there is at least
!  one level between the sfc and 3 km.
!
!-----------------------------------------------------------------------
!
     sumh=0.
     DO j=2,nlev(k)
       IF(hgt3d(j,k) > h3km) EXIT
       sumh=sumh +                                                     &
           ( u3d(j,k)*v3d(j-1,k) ) -                                   &
           ( v3d(j,k)*u3d(j-1,k) )
     END DO
!      240    CONTINUE
     k3km=j
     dz=hgt3d(j,k)-hgt3d(j-1,k)
     wlow=(hgt3d(j,k)-h3km)/dz
     u3km=u3d(j,k)*(1.-wlow) + u3d(j-1,k)*wlow
     v3km=v3d(j,k)*(1.-wlow) + v3d(j-1,k)*wlow
     sumh=sumh +                                                       &
           ( u3km*v3d(j-1,k) ) -                                       &
           ( v3km*u3d(j-1,k) )
     ushr=u3km-u3d(1,k)
     vshr=v3km-v3d(1,k)
     helicity(k)=sumh + vshr*ustrm(k) - ushr*vstrm(k)
!      write(*,*) 'PASSAGEM 9, com k= ',k
!
!-----------------------------------------------------------------------
!
!  Now Find 7km wind for 3-to-7km shear
!
!-----------------------------------------------------------------------
!
     DO j=k3km,nlev(k)
       IF(hgt3d(j,k) > h7km) EXIT
     END DO
!      260   CONTINUE
     dz=hgt3d(j,k)-hgt3d(j-1,k)
     wlow=(hgt3d(j,k)-h7km)/dz
     u7km=u3d(j,k)*(1.-wlow) + u3d(j-1,k)*wlow
     v7km=v3d(j,k)*(1.-wlow) + v3d(j-1,k)*wlow
     shr37(k)=(SQRT( (u7km-u3km)*(u7km-u3km) +                       &
                       (v7km-v3km)*(v7km-v3km) ))/4000.
  endif
 END DO
!  write(*,*) 'Saída da subrotina calcshr.'

!-----------------------------------------------------------------------
!
!  Calculate the low-level convergence
!
!-----------------------------------------------------------------------
!
!  dx=(x(2)-x(1))
!  dy=(y(2)-y(1))
!  dx2=2.*dx
!  dy2=2.*dy
!  DO j=2,ny-2
!    DO i=2,nx-2
!      blcon(i,j)=-1000.*(                                               &
!                 (tem2(i+1,j)-tem2(i-1,j))/(sigma(i,j)*dx2)+            &
!                 (tem3(i,j+1)-tem3(i,j-1))/(sigma(i,j)*dy2) )
!    END DO
!  END DO
!
!  DO j=2,ny-2
!    blcon(1,j)=-1000.*(                                                 &
!                (tem2(2,j)-tem2(1,j))/(sigma(1,j)*dx)+                  &
!                (tem3(1,j+1)-tem3(1,j-1))/(sigma(1,j)*dy2) )
!    blcon(nx-1,j)=-1000.*(                                              &
!           (tem2(nx-1,j)-tem2(nx-2,j))/(sigma(nx-1,j)*dx)+              &
!           (tem3(nx-1,j+1)-tem3(nx-1,j-1))/(sigma(nx-1,j)*dy2) )
!  END DO
!  DO i=2,nx-2
!    blcon(i,1)=-1000.*(                                                 &
!                (tem2(i+1,1)-tem2(i-1,1))/(sigma(i,1)*dx2)+             &
!                (tem3(i,2)-tem3(i,1))/(sigma(i,1)*dy) )
!    blcon(i,ny-1)=-1000.*(                                              &
!            (tem2(i+1,ny-1)-tem2(i-1,ny-1))/(sigma(i,ny-1)*dx2)+        &
!            (tem3(i,ny-1)-tem3(i,ny-2))/(sigma(i,ny-1)*dy) )
!  END DO
!  blcon(1,1)=blcon(2,2)
!  blcon(nx-1,ny-1)=blcon(nx-2,ny-2)
!
!-----------------------------------------------------------------------
!
!  Smooth the output arrays
!
!-----------------------------------------------------------------------
!
!  CALL smooth9p(   shr37,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(   ustrm,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(   vstrm,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(   srlfl,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(   srmfl,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(helicity,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(     brn,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(    brnu,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(   blcon,nx,ny,1,nx-1,1,ny-1,tem1)
!
 RETURN
END SUBROUTINE calcshr
