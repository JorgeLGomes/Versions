c
c *** Namelist variables - post.
c
      integer*4 nfiles,incr,begin,ngrid,ncropx,ncropy
c
      character*8   post_type
c
      common /post/post_type
     .            ,nfiles,incr,begin,ngrid,ncropx,ncropy
c
c *** Namelist variables - avs.
c
      integer*4 maxfields,maxplots
      parameter (maxfields=25,maxplots=5)
c
      integer*4 corflg,nfields
c
      character*256 flddir,corfil
      character*8   cfield(maxfields)
c
      common /avs/corflg,nfields
     .           ,flddir,corfil,cfield
c
c *** Namelist variables - ncarg
c
      integer*4 level,ncplots,ncfields
c
      real*4 cint(maxfields,maxplots)
c
      character*8 ncfield(maxfields,maxplots)
c
      common /eta_ncarg/level,ncplots,ncfields,ncfield,cint
