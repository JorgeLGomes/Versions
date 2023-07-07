      subroutine eta_commons
c
      implicit none

      real factor
      parameter(factor=57.2957795)
c
      include 'ecommons.h'
c
      namelist/model_grids/tlm0d,tph0d,im,jm,lm,ptinp,dlmd,dphd
     .                    ,dt,idtad,imonth,idate,iyear,istrtim
     .                    ,nsoil
     .                    ,ninit,tboco,init_in,init_gdsdir
     .                    ,init_out
c
      namelist/sfc_init/topo_in,topo_out,SEARES,FRACLK,GRIBSOIL,
     .         SIGMA,REAN,REAN_SFC

c
      namelist/coac_nml/coac,lcoac
c
c_______________________________________________________________________________
c
c *** Read ETA namelist.
c
      open(1,file='ETAIN',form='formatted',status='old',err=900)
      read(1,model_grids,end=901)
      read(1,sfc_init,end=902)
      close(1)

Cnew
       
        NHOUR=(NINIT-1)*TBOCO
        write(6,*) 'initializing for a model run of length: ', 
     +  NHOUR, 'hours'
Cnew
c
c *** Fill model grid constants.
c
Cmp      im=-wbd/dlmd+1.5
Cmp      jm=-2*sbd/dphd+1.5

      WBD=-(IM-1)*DLMD
      SBD=-((JM-1)/2.)*DPHD
Cmp
      imt=2*im-1
      jmt=jm/2+1
      imjm=im*jm-jm/2
      imjm1=imjm-1
      im1=im-1
      im2=im-2
      jm1=jm-1
      jm2=jm-2
      jm3=jm-3
      jm4=jm-4
      jm5=jm-5
      imm1=imt-1
      imm3=imt-3
      imm5=imt-5
      jmm1=jmt-1
      jmm2=jmt-2
      kb=im*2+jm-3
      lb=jmt*2+imt-3
      lmp1=lm+1
      jam=6+2*(jm-10)
c
C      ctph0=cosd(tph0d)
C      stph0=sind(tph0d)
      ctph0=cos(tph0d/factor)
      stph0=sin(tph0d/factor)
c
      print *,' '
      print *,'ETA grid parameters:'
      print *,'   im   =',im   ,' jm  =',jm
      print *,'   imt  =',imt  ,' jmt =',jmt
      print *,'   imjm =',imjm
      print *,' '
        print*, 'soil levels  in ecommons ', nsoil
Cmp	
        print*,' centered at ',TPH0D,TLM0D
C	print*, TPH0D, TLM0D
C	print*, 'cos/sin vals...', ctph0,stph0
Cmp
c

       write(6,*) 'seares in eta_commons: ', seares
Cnew
c_______________________________________________________________________________
c
c *** Read COAC namelist.
c
      open(3,file='COAC.nml',form='formatted',status='old',err=905)
      read(3,coac_nml,end=903)
      close(3)

      return
c
c *** Error trapping.
c
900   print *,'ETAIN namelist not found.'
      stop
c
901   print *,'Error reading namelist - model_grids'
      stop
c
902   print *,'Error reading namelist - surface'
      stop
c
903   print *,'Error reading namelist - plot_diag'
      stop
c
904   print *,'Error reading namelist - init_diag'
      stop
c
905   print *,'Error reading namelist - coac'
      stop
c
      end
