c
c *** Maximum allowable init files.
c
      integer*4 max_init_files
      parameter (max_init_files=50)
c
c *** Grid constants.
c
      real*4 ctph0,stph0
c
      integer*4  imjm1
     .         ,im2
     .         ,jm1,jm3,jm4,jm5
     .         ,imm1,imm3,imm5
     .         ,jmm1,jmm2
     .         ,kb,lb,lmp1
     .         ,jam

c
      common /grid_const/ctph0,stph0
C     .                  ,im,jm
C     .                  
     .                  ,im2
     .                  ,jm1,jm3,jm4,jm5
     .                  ,imm1,imm3,imm5
     .                  ,jmm1,jmm2
     .                  ,kb,lb,lmp1
     .                  ,jam
c
c *** Namelist variables - model_grids.
c
      real*4 tlm0d,tph0d,wbd,sbd,dlmd,dphd
     .      ,dt,w
     .      ,tboco
c
      integer*4
     .         idtad
     .         ,imonth,idate,iyear,istrtim
     .         ,nsoil
     .         ,ninit
     .         ,nhour
c
      character*256 init_in(max_init_files),init_out,fcst_out
c
      common /model_grids/tlm0d,tph0d,wbd,sbd,dlmd,dphd
C     .                   ,lm
     .                   ,dt,w,idtad
     .                   ,imonth,idate,iyear,istrtim
     .                   ,nsoil
     .                   ,ninit,init_in,init_out
     .                   ,tboco,nhour
     .                   ,fcst_out
c
c *** Namelist variables - surface.
c
      character*256 topo_in,topo_out
     .             ,soil_in,soil_out
     .             ,veg_in,veg_out
c
      logical soili,vegi,surfi,climsst,surft,surfw,gribsoil
	real seares
c
      common /surface/topo_in,topo_out
     .               ,soili,soil_in,soil_out
     .               ,vegi,veg_in,veg_out
     .               ,surfi,climsst,surft,surfw,gribsoil,seares
c
c *** Namelist variables - plot_diag.
c
      logical plot
c
      character*256 plot_out
c
      common /plot_diag/plot,plot_out
c
c *** Namelist variables - init_diag.
c
      logical idtbls,hmonly,print,fpmnts
     .       ,siluet,sigma,mrmsxl,nomnts
c
      common /init_diag/idtbls,hmonly,print,fpmnts
     .                 ,siluet,sigma,mrmsxl,nomnts
