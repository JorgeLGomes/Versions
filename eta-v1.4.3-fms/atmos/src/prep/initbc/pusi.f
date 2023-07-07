        subroutine pusi(spl,htll,qtll,utll,vtll,sfcgrid,ldm,idat)
!     *************************************************************     
!     *                                                           *
!     *  program for conversion from p to sigma system            *     
!     *  using cubic splines.                                     *     
!     *  programers - z.janjic, d. jovic and s. nickovic          *     
!     *         1981, 1999, fed. hydromet. inst., beograd, ncep   *     
!     *                                                           *     
!     *************************************************************     
!-----------------------------------------------------------------------
!      include "../include/deco.inc"
!-----------------------------------------------------------------------
!      include "../include/model.inc"
        include 'ecommons.h'
        include 'econstants.h'
!-----------------------------------------------------------------------
C      parameter(dtr=3.141592654/180.)
C      parameter(g=9.8,gor=9.8/287.04,nsmud=1,epsq=1.e-12)
      parameter(gor=9.8/287.04,nsmud=1,epsq=1.e-12)
!-----------------------------------------------------------------------
                             d i m e n s i o n
     & idat(5),dsg(lm),sgml(lm),sg(lm+1),dfl(lm+1)
                             d i m e n s i o n
     & zuv(lm),util(lm),vtil(lm)
     &,zh(lm+1),dg(lm+1)
     &,zqtil(lm),qtil(lm)
                             d i m e n s i o n
     & spl(ldm),zsl(ldm),rzsl(ldm),rhsl(ldm)
     &,y2(ldm+1),pp(ldm+1),qq(ldm+1)
     &,zslh(ldm+1),hsp(ldm+1)
     &,uij(ldm+1),vij(ldm+1),qij(ldm+1)
                             d i m e n s i o n
     & ihw(jm),ihe(jm),ivw(jm),ive(jm)
                             d i m e n s i o n
     & ptll(im,jm)
     &,sm(im,jm),hgt(im,jm),hgts(im,jm),fis(im,jm)
     &,pd(im,jm),ref(im,jm)
                             d i m e n s i o n
     & utll(im,jm,ldm),vtll(im,jm,ldm)
     &,htll(im,jm,ldm),qtll(im,jm,ldm)
                             d i m e n s i o n
     & u(im,jm,lm),v(im,jm,lm),sfcgrid(im,jm,12)
     &,t(im,jm,lm),q(im,jm,lm)
     &,z(im,jm,lm+1),alpi(im,jm,lm+1)
                             d i m e n s i o n
     & pdb(kb,2)
     &,ub(kb,lm,2),vb(kb,lm,2),tb(kb,lm,2),qb(kb,lm,2)

!-----------------------------------------------------------------------
      character*3 sfx
      character*64 infile,outfil 
!
Cmp      logical run,sigma,print
      logical run
!-----------------------------------------------------------------------
      data ntsd/0/     
C
Cmp  SPL is passed as an argument
C
Cmp      data spl/1000.,2000.,3000.,5000.,7000.,10000.,15000.,20000.
Cmp     &        ,25000.,30000.,35000.,40000.,45000.,50000.,55000.,60000.
Cmp     &        ,65000.,70000.,75000.,80000.,85000.,90000.,92500.,95000.
Cmp     &        ,97500.,100000./
Cmp      save spl

C	DATA PT/2500./
!-----------------------------------------------------------------------
	character*3 fname
                             d a t a
     & infile/'                                                         
     &       '/
     &,outfil/'                                                         
     &       '/
!-----------------------------------------------------------------------
 2300 format(' *** ',2(i2,'.'),i4,'.   ',i3,' gmt ***')  
 2500 format(' print=',l1)                                
 2001 format(' ',91f6.0)                                                
 2002 format(' ',3x,91f6.0)                                             
 2003 format(' sloj izmedju povrsina sigma=',f5.3,' I sigma=',f5.3,'  **
     &* ',3(i2,'.'),i3,' gmt ***')                                      
 2004 format(' temperatura')                                            
 2005 format(' zonalna komponenta brzine')                              
 2006 format(' meridionalna komponenta brzine')                         
 2007 format(' specificna vlaga')                         
 2008 format(' topografija/10')                         

CCCCCCCCC  READ FROM FILE METHOD  CCCCCCCCCCCCCCCCCCCCCCCCCCCC

	pt=ptinp
	write(6,*) 'set pt= ', pt

C        open(unit=16,file='deta',form='unformatted',access=
C     +                          'sequential')

C        REWIND 16
C        READ(16)dsg,LDUM

C        sg(1)=0.0
C	do l=1,lm
C	sg(l+1)=sg(l)+dsg(l)
C	enddo

C       sg(lm+1)=1.0

CCCCCCCCC EMPIRICAL FORMULA METHOD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	
          do l=1,lm-1
      x=float(l-1)/float(lm-1)
      sg(l)=0.8*x+1.03643*x**2-0.973*x**3+0.13673*x**4-0.00016*x**5
          end do
      dsgb=1.-sg(lm-1)
      sg(lm)=1.-dsgb/2

      sg(lm+1)=1.0

          do l=1,lm
          dsg(l)=sg(l+1)-sg(l)
          enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          do l=1,lm
      sgml(l)=0.5*(sg(l)+sg(l+1))
	write(6,*) 'L, sgml(L): ', L,sgml(L)
          enddo
!
C      open(unit=1,file='../output/deta.dat'
      open(unit=1,file='deta.dat'
     &    ,status='unknown',form='unformatted')
      write(1) dsg
      close(1)
!-----------------------------------------------------------------------
          do ld=1,ldm+1
      y2(ld)=0.
          enddo
!--------------pusi derived constants ----------------------------------
      print=.false.
!      print=.true.
      alpt=log(pt)
      ztop=alpt*alpt
!
          do ld=1,ldm
      alp=log(spl(ld))
      zsl(ld)=alp*alp
          enddo
!-----------------------------------------------------------------------
      zld2=zsl(ldm-2)                                             
      zld1=zsl(ldm-1)                                                
      zld0=zsl(ldm  )                                               
!--------------------------------------------------------------------- 
Cmp      open(unit=2,file='../output/mnts.bilin'
Cmp     &    ,status='unknown',form='unformatted')
Cmp      read(2) hgt,sm
Cmp      read(2) hgt,sm  ! hgt in second record not multiplied by hbm2
Cmp      close(2)
C
      open(1,file=topo_out,status='old',form='unformatted')
      rewind(1)
      read(1) hgt,sm
      close(1)
!-----------------------------------------------------------------------
          do j=1,jm
      ihw(j)=-mod(j,2)
      ihe(j)=ihw(j)+1
      ivw(j)=-mod(j+1,2)
      ive(j)=ivw(j)+1
          enddo
!--------------5-point smoothing of mountains---------------------------
                      if(nsmud.gt.0)    then
!-----------------------------------------------------------------------
                  do n=1,nsmud
!-----------------------------------------------------------------------
              do j=2,jm-1
          do i=1+mod(j,2),im-1+mod(j,2)
          if(sm(i,j).lt.0.5)    then
      hgts(i,j)=(hgt(i+ihw(j),j-1)+hgt(i+ihe(j),j-1)
     &          +hgt(i+ihw(j),j+1)+hgt(i+ihe(j),j+1)
     &          +hgt(i,j)*4.)*0.125
          else
      hgts(i,j)=hgt(i,j)
          endif
          enddo
              enddo
!
              do j=2,jm-1
          do i=2,im-1
      hgt(i,j)=hgts(i,j)
          enddo
              enddo
!-----------------------------------------------------------------------
                  enddo
!-----------------------------------------------------------------------
                      endif
!-----------------------------------------------------------------------
C	smooth topo on boundaries
C
C      call smdhld(hgt,sm,5,5)
	nlines=jm/20
      call smdhld(hgt,sm,nlines,10)
!-----------------------------------------------------------------------
!-------------4-point averaging of mountains along inner boundary-------

          do i=1,im-1
      hgt(i,2)=(hgt(i,1)+hgt(i+1,1)+hgt(i,3)+hgt(i+1,3))
     &  *0.25
!     &        *(1.-sm(i,2))*0.25
          enddo
          do i=1,im-1
      hgt(i,jm-1)=(hgt(i,jm-2)+hgt(i+1,jm-2)+
     +                          hgt(i,jm)+hgt(i+1,jm))
     &  *0.25
!     &        *(1.-sm(i,jm-1))*0.25
          enddo
          do j=4,jm-3,2
      hgt(1,j)=(hgt(1,j-1)+hgt(2,j-1)+
     +                          hgt(1,j+1)+hgt(2,j+1))
     &  *0.25
!     &        *(1.-sm(1,j))*0.25
          enddo
          do j=4,jm-3,2
      hgt(im-1,j)=(hgt(im-1,j-1)+hgt(im,j-1)+
     +                          hgt(im-1,j+1)+hgt(im,j+1))
     &  *0.25
!     &        *(1.-sm(im-1,j))*0.25
        enddo
C
!-----------------------------------------------------------------------
                  do l=1,lm
              do j=1,jm
          do i=1,im
      u(i,j,l)=0.
      v(i,j,l)=0.
      t(i,j,l)=0.
      q(i,j,l)=0.
          enddo
              enddo
                  enddo
!-----------------------------------------------------------------------
      print*,'*** Hi, this is pusi converting pressure data to sigma'
      ibc=tboco
C      ihr=0

        write(6,*) 'into pusi... idat(5), LDM=  ', idat(5), LDM
        ihr=idat(5)
!-----------------------------------------------------------------------
 1000 continue
!--------------read pressure data---------------------------------------

Cmp****** Data in WS Eta read in the the getdata subroutine which is
Cmp****** called before this routine

C      write(sfx,'(i3.3)') ihr
C      infile='../output/lltll.'//sfx
C      open(unit=3,file=infile,status='old',form='unformatted')
C      read(3) run,idat,ihrst,ihr,iend,ptll,htll,utll,vtll,qtll
C      close(3)
!
C      print*,'*** pusi read pressure data from ',infile
C      print*,'*** idat=',idat,' ',ihrst,' utc + ',ihr

C check incoming q data

        do l=1,ldm
        qmx=-9999.
        qmn=9999.
        do j=1,jm
        do i=1,im
        if (qtll(i,j,l) .gt. qmx) qmx=qtll(i,j,l)
        if (qtll(i,j,l) .lt. qmn) qmn=qtll(i,j,l)
        enddo
        enddo
C	write(6,*) 'extremes at L= ', l, 'are: ', qmn,qmx
	enddo

C end check
!--------------computation of pd, zeta and related variables------------
                  do j=1,jm
              do i=1,im                                                
!-----------------------------------------------------------------------
      hgtp=hgt(i,j)                                              
      hld0=htll(i,j,ldm)
!-----------------------------------------------------------------------
          if(hld0.gt.hgtp)    then
!--------------if sfc below lowest p sfc, extrapolate downwards---------
      hld2=htll(i,j,ldm-2)
      hld1=htll(i,j,ldm-1)
!
      d1=(zld0-zld1)/(hld0-hld1)
!      d2=(d1-(zld1-zld2)/(hld1-hld2))/(hld0-hld2)
      d2=0.
!
      x=hgtp-hld0
      zsp=d2*x*x+d1*x+zld0
!-----------------------------------------------------------------------
          else
!--------------otherwise, use spline to interpolate---------------------
          do ld=1,ldm
      rhsl(ld)=htll(i,j,ldm+1-ld)
      rzsl(ld)=zsl(ldm+1-ld)
      y2(ld)=0.
          enddo
!
      call spline_pus(ldm,rhsl,rzsl,y2,1,hgtp,zsp,pp,qq)
!-----------------------------------------------------------------------
          endif
!-----------------------------------------------------------------------
      z(i,j,lm+1)=zsp                                                
      plpi=sqrt(zsp)                                                
      alpi(i,j,lm+1)=plpi                                          
      pdp=exp(plpi)-pt                                            
      pd(i,j)=pdp                                                
!                                                               
      alpi(i,j,1)=alpt
      z(i,j,1)=ztop
!
          do l=2,lm                                             
      plpi=log(pt+pdp*sg(l))                                 
      alpi(i,j,l)=plpi                                      
      z(i,j,l)=plpi*plpi                                               
          enddo
!-----------------------------------------------------------------------
              enddo
                  enddo
	
	write(6,*) 'to velocity spline'
!--------------velocity spline interpolation inside the domain----------
                  do j=2,jm-1
              ivl=1+mod(j+1,2)
!
              do i=ivl,im-1
          do ld=1,ldm                                                   
      zslh(ld)=zsl(ld)
      uij(ld)=utll(i,j,ld)                                                 
      vij(ld)=vtll(i,j,ld)                                               
          enddo
!                                                                  
      zus=ztop                                                   
      zuw=ztop                                                  
      zue=ztop                                                        
      zun=ztop                                                       
!                                                                    
          do l=1,lm                              
      zls=z(i,j-1,l+1)
      zlw=z(i+ivw(j),j,l+1)
      zle=z(i+ive(j),j,l+1)
      zln=z(i,j+1,l+1)    
!                                                                    
      zuv(l)=.125*(zus+zls+zuw+zlw+zue+zle+zun+zln)
!                                                                     
      zus=zls                                                      
      zuw=zlw                                                     
      zue=zle                                                    
      zun=zln 
          enddo
!                                                                
      zslpu=(z(i,j-1,lm+1)+z(i+ivw(j),j,lm+1)
     &      +z(i+ive(j),j,lm+1)+z(i,j+1,lm+1))*0.25
!
          if(zslpu.le.zld0)    then
      lold=ldm
          else
      lold=ldm+1
      zslh(lold)=zslpu
      uij(lold)=uij(lold-1)
      vij(lold)=vij(lold-1)
          endif
!
          do ld=1,lold
      y2(ld)=0.
          enddo
!
!          if(i.eq.21.and.j.eq.41) then
!      print*,uij
!      print*,zslh
!      print*,vij
!          endif
      call spline_pus(lold,zslh,uij,y2,lm,zuv,util,pp,qq) 
      call spline_pus(lold,zslh,vij,y2,lm,zuv,vtil,pp,qq)
!          if(i.eq.21.and.j.eq.41) then
!      print*,util
!      print*,zuv
!      print*,vtil
!          endif
!               
!      if(i.eq.5.and.j.eq.36) then
!         print*,lold,lm
!         print*,zslh
!         print*,uij
!         print*,vij
!         print*,zuv
!         print*,util
!         print*,vtil
!         print*,y2
!         print*,pp
!         print*,qq
!      end if
!                                                                 
          do l=1,lm                                              
      u(i,j,l)=util(l)                      
      v(i,j,l)=vtil(l)                     
          enddo
              enddo
                  enddo
!-----------------------------------------------------------------------
                  do l=1,lm
              do j=1,jm
          do i=1,im
      if(abs(u(i,j,l)).gt.100.) print*,'pusi,i,j,l,u',i,j,l,u(i,j,l)
      if(abs(v(i,j,l)).gt.100.) print*,'pusi,i,j,l,v',i,j,l,v(i,j,l)
          enddo
              enddo
                  enddo
!--------------velocity spline interpolation at n and s boundaries------
                  do j=1,jm,jm-1
              do i=1,im-1
          do ld=1,ldm                                               
      zslh(ld)=zsl(ld)
      uij(ld)=utll(i,j,ld)                                                 
      vij(ld)=vtll(i,j,ld)                                                
          enddo
!                                                                   
      zuw=ztop                                                   
      zue=ztop 
!                                                                    
          do l=1,lm                                                 
      zlw=z(i+ivw(j),j,l+1) 
      zle=z(i+ive(j),j,l+1)                                           
!                                                               
      zuv(l)=.25*(zuw+zlw+zue+zle)                    
!                                                                    
      zuw=zlw  
      zue=zle                                                     
          enddo
!
      zslpu=(z(i+ivw(j),j,lm+1)+z(i+ive(j),j,lm+1))*0.5
!
          if(zslpu.le.zld0) then
      lold=ldm
          else
      lold=ldm+1
      zslh(lold)=zslpu
      uij(lold)=uij(lold-1)
      vij(lold)=vij(lold-1)
          endif
!
          do ld=1,lold
      y2(ld)=0.
          enddo
!
      call spline_pus(lold,zslh,uij,y2,lm,zuv,util,pp,qq)
      call spline_pus(lold,zslh,vij,y2,lm,zuv,vtil,pp,qq)
!
          do l=1,lm                                              
      u(i,j,l)=util(l)              
      v(i,j,l)=vtil(l)                     
          enddo
              enddo
                  enddo
!--------------velocity spline interpolation at e and w boundaries-----
                  do j=2,jm-1,2
              do i=1,im,im-1                  
          do ld=1,ldm                                                
      zslh(ld)=zsl(ld)
      uij(ld)=utll(i,j,ld)                                              
      vij(ld)=vtll(i,j,ld)                                             
          enddo
!                                                                
      zus=ztop                                                       
      zun=ztop                                                      
!                                                                  
          do l=1,lm                                               
      zls=z(i,j-1,l+1)
      zln=z(i,j+1,l+1)                                                   
!                                                                   
      zuv(l)=.25*(zus+zls+zun+zln)                               
! 
      zus=zls
      zun=zln 
          enddo
!                                                        
      zslpu=(z(i,j-1,lm+1)+z(i,j+1,lm+1))*0.5
!
          if(zslpu.le.zld0)    then
      lold=ldm
          else
      lold=ldm+1
      zslh(lold)=zslpu
      uij(lold)=uij(lold-1)
      vij(lold)=vij(lold-1)
          endif
!
          do ld=1,lold
      y2(ld)=0.
          enddo
!
      call spline_pus(lold,zslh,uij,y2,lm,zuv,util,pp,qq)
      call spline_pus(lold,zslh,vij,y2,lm,zuv,vtil,pp,qq)
!                                                                   
          do l=1,lm                                               
      u(i,j,l)=util(l)              
      v(i,j,l)=vtil(l)                  
          enddo
              enddo
                  enddo
!      if(ihr.eq.0) print*,u,v
!-----------------------------------------------------------------------
!      do l=1,lm
!      do j=1,jm
!      do i=1,im
!      if(abs(u(i,j,l)).gt.100.) print*,'i,j,l,u',i,j,l,u(i,j,l)
!      if(abs(v(i,j,l)).gt.100.) print*,'i,j,l,v',i,j,l,v(i,j,l)
!      enddo
!      enddo
!      enddo
!-----------------------------------------------------------------------
                      if(print)    then
                  do l=1,lm
      write(*,2003) sg(l),sg(l+1),idat,ihrst     
!
      write(*,2005)                                                
!                                                                 
              do j=1,jm
          if(mod(j,2).ne.0)   then                  
      write (*,2002) (u(i,j,l),i=1,im-1)                         
          else
      write (*,2001) (u(i,j,l),i=1,im)                               
          endif
              enddo
!                                                                    
      write(*,2006) 
!                                                                  
              do j=1,jm                                               
          if(mod(j,2).ne.0) then           
      write (*,2002) (v(i,j,l),i=1,im-1)
          else                      
      write (*,2001) (v(i,j,l),i=1,im)
          endif
              enddo
                  enddo 
                      endif
!--------------computation of sigma temperatures------------------------
                  do j=1,jm
              do i=1,im
!-----------------------------------------------------------------------
      zh(1)=ztop                                                      
!                                                                    
          do l=2,lm+1                                               
      zh(l)=z(i,j,l)                                               
          enddo
!
          do ld=1,ldm
      zslh(ld)=zsl(ld)
      hsp(ld)=htll(i,j,ld)
          enddo
!
          if(z(i,j,lm+1).le.zld0) then
      lold=ldm
          else
      lold=ldm+1
!
      zslh(lold)=z(i,j,lm+1)
      hsp(lold)=hgt(i,j)
          endif
!
          do ld=1,lold
      y2(ld)=0.
          enddo
!--------------temperatures inside the integration domain---------------
      call spline_pus(lold,zslh,hsp,y2,lm+1,zh,dg,pp,qq)                   
!                                                                  
          do l=1,lm
      t(i,j,l)=(dg(l)-dg(l+1))*gor/(alpi(i,j,l+1)-alpi(i,j,l))
          enddo                      
!-----------------------------------------------------------------------
              enddo
                  enddo
!-----------------------------------------------------------------------
                  do l=1,lm
              do j=1,jm
          do i=1,im
      if(abs(t(i,j,l)).lt.160.) print*,'COLD: i,j,l,t',i,j,l,t(i,j,l)
          enddo
              enddo
                  enddo
!-----------------------------------------------------------------------
              do j=1,jm
          do i=1,im
      fis(i,j)=hgt(i,j)*g
          enddo
              enddo
!-----------------------------------------------------------------------
                      if(print)    then
!                                                                    
                  do l=1,lm
!                                                                  
      write(*,2003) sg(l),sg(l+1),idat,ihrst 
      write(*,2004)                                              
!                                                                      
              do j=1,jm                        
          if(mod(i,2).ne.0)    then
      write (*,2001) (t(i,j,l)-273.,i=1,im)
          else
      write (*,2002) (t(i,j,l)-273.,i=1,im-1)
          endif
             enddo
                 enddo
                     endif                        
!--------------spec hum spline interpolation inside the domain----------
                  do j=1,jm
              do i=1,im
          do ld=1,ldm
      zslh(ld)=zsl(ld)
      qij(ld)=qtll(i,j,ld)
          enddo
!
      zu=ztop
          do l=1,lm
      zl=z(i,j,l+1)
      zqtil(l)=.5*(zu+zl)
      zu=zl
          enddo
!
          if(z(i,j,lm+1).le.zld0)    then
      lold=ldm
          else
      lold=ldm+1
!
      zss=z(i,j,lm+1)
      zslh(lold)=zss
!
      qld2=qij(ldm-2)
      qld1=qij(ldm-1)
      qld0=qij(ldm  )
!
      d1=(qld0-qld1)/(zld0-zld1)
!      d2=(d1-(qld1-qld2)/(zld1-zld2))/(zld0-zld2)
      d2=0.
      x=zss-zld0
!
      qij(lold)=d2*x*x+d1*x+qld0
          endif
!
          do ld=1,lold
      y2(ld)=0.
          enddo
!
      call spline_pus(lold,zslh,qij,y2,lm,zqtil,qtil,pp,qq)
!
          do l=1,lm
      qp=amax1(qtil(l),epsq)
      t(i,j,l)=t(i,j,l)/(qp*0.608+1.)
      q(i,j,l)=qp
	if(l .eq. 5 .and. qp .gt. 1.e10) then
	write(6,*) 'BIG Q: I,J,q(i,j,5) ', I,J,q(i,j,5)
	endif
          enddo
              enddo
                  enddo
!----------------------------------------------------------------------
!                  do l=1,lm
!              do j=1,jm
!          do i=1,im
!      if(abs(u(i,j,l)).gt.100.) print*,'i,j,l,u',i,j,l,u(i,j,l)
!      if(abs(v(i,j,l)).gt.100.) print*,'i,j,l,v',i,j,l,v(i,j,l)
!          enddo
!              enddo
!                  enddo
!!
                  do l=1,lm
!      umx=-99999.
!      umn=99999.
!      vmx=-99999.
!      vmn=99999.
!      tmx=-99999.
!      tmn=99999.
      qmx=-99999.
      qmn=99999.
              do j=1,jm
          do i=1,im
!      umx=amax1(u(i,j,l),umx)
!      umn=amin1(u(i,j,l),umn)
!      vmx=amax1(v(i,j,l),vmx)
!      vmn=amin1(v(i,j,l),vmn)
!      tmx=amax1(t(i,j,l),tmx)
!      tmn=amin1(t(i,j,l),tmn)
      qmx=amax1(q(i,j,l),qmx)
      qmn=amin1(q(i,j,l),qmn)
          enddo
              enddo
!!
!      print*,'l=',l,' umx=',umx,' umn=',umn,
!      print*,'l=',l,' vmx=',vmx,' vmn=',vmn
!      print*,'l=',l,' tmx=',tmx,' tmn=',tmn
!      print*,'l=',l,' qmx=',qmx,' qmn=',qmn
!! 
                  enddo
!-----------------------------------------------------------------------
                      if(print)    then
!
                  do l=1,lm
!
      write(*,2003) sg(l),sg(l+1),idat,ihrst
      write(*,2007)
!
              do j=1,jm
          if(mod(i,2).ne.0)    then
      write (*,2001) (q(i,j,l)*1000.,i=1,im)
          else
      write (*,2002) (q(i,j,l)*1000.,i=1,im-1)
          endif
             enddo
                 enddo
                     endif
!-----------------------------------------------------------------------
                      if(print)    then
!
      write(*,2008)
!
              do j=1,jm
          if(mod(i,2).ne.0)    then
      write (*,2001) (hgt(i,j)/100.,i=1,im)
          else
      write (*,2002) (hgt(i,j)/100.,i=1,im-1)
          endif
             enddo
                     endif
!-------------separation of boundary values-----------------------------
      n=1
              do i=1,im
          pdb(n,1)=pd(i,1)
          pdb(n,2)=0.
          do l=1,lm
      tb(n,l,1)=t(i,1,l)
      tb(n,l,2)=0.
      qb(n,l,1)=q(i,1,l)
      qb(n,l,2)=0.
          enddo
          n = n+1
              enddo
              do i=1,im
          pdb(n,1)=pd(i,jm)
          pdb(n,2)=0.
          do l=1,lm
      tb(n,l,1)=t(i,jm,l)
      tb(n,l,2)=0.
      qb(n,l,1)=q(i,jm,l)
      qb(n,l,2)=0.
          enddo
          n = n+1
              enddo
              do j=3,jm-2,2
          pdb(n,1)=pd(1,j)
          pdb(n,2)=0.
          do l=1,lm
      tb(n,l,1)=t(1,j,l)
      tb(n,l,2)=0.
      qb(n,l,1)=q(1,j,l)
      qb(n,l,2)=0.
          enddo
          n = n+1
              enddo
              do j=3,jm-2,2
          pdb(n,1)=pd(im,j)
          pdb(n,2)=0.
          do l=1,lm
      tb(n,l,1)=t(im,j,l)
      tb(n,l,2)=0.
      qb(n,l,1)=q(im,j,l)
      qb(n,l,2)=0.
          enddo
          n = n+1
              enddo
!-----------------------------------------------------------------------
      n=1
              do i=1,im-1
          do l=1,lm
      ub(n,l,1)=u(i,1,l)
      ub(n,l,2)=0.
      vb(n,l,1)=v(i,1,l)
      vb(n,l,2)=0.
          enddo
          n = n+1
              enddo
              do i=1,im-1
          do l=1,lm
      ub(n,l,1)=u(i,jm,l)
      ub(n,l,2)=0.
      vb(n,l,1)=v(i,jm,l)
      vb(n,l,2)=0.
          enddo
          n = n+1
              enddo
              do j=2,jm-1,2
          do l=1,lm
      ub(n,l,1)=u(1,j,l)
      ub(n,l,2)=0.
      vb(n,l,1)=v(1,j,l)
      vb(n,l,2)=0.
          enddo
          n = n+1
              enddo
              do j=2,jm-1,2
          do l=1,lm
      ub(n,l,1)=u(im,j,l)
      ub(n,l,2)=0.
      vb(n,l,1)=v(im,j,l)
      vb(n,l,2)=0.
          enddo
          n = n+1
              enddo
!-----------------------------------------------------------------------
      if(ihr.eq.0) then
!-----------------------------------------------------------------------
      run=.true.
      ntsd=0
          do l=1,lm+1
      dfl(l)=0.
          enddo
!
              do j=1,jm
          do i=1,im
      ref(i,j)=1.
          enddo
              enddo
!
      sigma=.true.
!-----------------------------------------------------------------------
C      open(unit=4,file='../output/nfcst'
C     &    ,status='unknown',form='unformatted')
!
C      write(4) run,idat,ihrst,ntsd,sigma,u,v
C      write(4) t,q,pd,fis,sm,ref,sg,dsg,sgml,dfl
C      close(4)
        NFCSTE=1

         l=index(init_out//' ',' ')-1
         if (init_out(l:l) .ne. '/') then
            l=l+1
            init_out(l:l)='/'
         endif

C       write(6,*) 'trying to open ', init_out(1:l)//'preproc.init'
      open(unit=NFCSTE,file=init_out(1:l)//'preproc.init'
     .       ,status='unknown',form='unformatted')
        write(6,*) 'writing to ',init_out(1:l)//'preproc.init'

      REWIND NFCSTE
C      write(NFCSTE) run,idat(1),idat(2),idat(3),idat(4),ntsd,sigma
      write(NFCSTE) run,idat(1),idat(2),idat(3),idat(4),ntsd,u,v
C		   
C     idat(5),dsg(lm),sgml(lm),sg(lm+1),dfl(lm+1)
	
Ctst
C	open(unit=4,file=init_out(1:l)//'sigma_temps.fil',
C     + status='unknown',form='unformatted')
C	write(4) t
Ctst

      write(NFCSTE) t,q,pd,fis,sm,ref,sg,pt,dsg,sgml,dfl
      write(NFCSTE) sfcgrid
      close(NFCSTE)


!-----------------------------------------------------------------------
      endif
!-----------------------------------------------------------------------
C      write(sfx,'(i3.3)') ihr
!
C      outfil='../output/bc.'//sfx
C      open (unit=5,file=outfil,status='unknown',form='unformatted')
C      write(5) run,idat,ihrst,pdb,tb,qb,ub,vb
C      close(5)
!
C      print*,'*** pusi wrote to ',outfil
!      print*,ub

      l=index(init_out//' ',' ')-1
      if (init_out(l:l) .ne. '/') then
         l=l+1
         init_out(l:l)='/'
      endif
      write(fname,'(i6.6)') idat(5)
        write(6,*) 'opening', init_out(1:l)//'preproc.bc.'//fname
      open (1,file=init_out(1:l)//'preproc.bc.'//fname
     .     ,status='unknown',form='unformatted')
      write(1) run,idat(1),idat(2),idat(3),idat(4),pdb,tb,qb,ub,vb

        if (idat(5) .eq. 0) then
      open (11,file=init_out(1:l)//'bndy.newstyle'
     .     ,status='unknown',form='unformatted')
        endif

	write(11) pdb,tb,qb,ub,vb

      close(1)

!-----------------------------------------------------------------------
C      ihr=ihr+ibc
C      if(ihr.le.iend) go to 1000
!-----------------------------------------------------------------------
	return
C      stop
      end 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine smdhld(h,s,lines,nsmud)
!     ******************************************************************
!     *                                                                *
!     *  routine for                                                   *
!     *                                                                *
!     *                                                                *
!     ******************************************************************
!-----------------------------------------------------------------------
!      include "../include/model.inc"
        include 'ecommons.h'
!-----------------------------------------------------------------------
      dimension ihw(jm),ihe(jm)
      dimension h(im,jm),s(im,jm),hbms(im,jm),hne(im,jm),hse(im,jm)
!-----------------------------------------------------------------------
          do j=1,jm
      ihw(j)=-mod(j,2)
      ihe(j)=ihw(j)+1
          enddo
!-----------------------------------------------------------------------
              do j=1,jm
          do i=1,im
      hbms(i,j)=1.-s(i,j)
          enddo
              enddo
!
      jmlin=jm-lines+1
      ibas=lines/2
      m2l=mod(lines,2)
!
              do j=lines,jmlin
          ihl=ibas+mod(j,2)+m2l*mod(j+1,2)
          ihh=im-ibas-m2l*mod(j+1,2)
!
          do i=ihl,ihh
      hbms(i,j)=0.
          enddo
              enddo
!-----------------------------------------------------------------------
                  do ks=1,nsmud
              do j=1,jm-1
          do i=1,im-1
      hne(i,j)=h(i+ihe(j),j+1)-h(i,j)
          enddo
              enddo
              do j=2,jm
          do i=1,im-1
Cmp??      hne(i,j)=h(i+ihe(j),j-1)-h(i,j)
      hse(i,j)=h(i+ihe(j),j-1)-h(i,j)
          enddo
              enddo
!
              do j=2,jm-1
          do i=1+mod(j,2),im-1
      h(i,j)=(hne(i,j)-hne(i+ihw(j),j-1)
     &       +hse(i,j)-hse(i+ihw(j),j+1))*hbms(i,j)*0.125+h(i,j)
          enddo
              enddo
                      enddo
!-----------------------------------------------------------------------
      return
      end

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine spline_pus(nold,xold,yold,y2,nnew,xnew,ynew,p,q)
!
!     ******************************************************************
!     *                                                                *
!     *  this is a one-dimensional cubic spline fitting routine        *
!     *  programed for a small scalar machine.                         *
!     *                                                                *
!     *  programer: z. janjic, yugoslav fed. Hydromet. Inst., beograd  *
!     *                                                                *
!     *                                                                *
!     *                                                                *
!     *  nold - number of given values of the function.  Must be ge 3. *
!     *  xold - locations of the points at which the values of the     *
!     *         function are given.  Must be in ascending order.       *
!     *  yold - the given values of the function at the points xold.   *
!     *  y2   - the second derivatives at the points xold.  If natural *
!     *         spline is fitted y2(1)=0. And y2(nold)=0. Must be      *
!     *         specified.                                             *
!     *  nnew - number of values of the function to be calculated.     *
!     *  xnew - locations of the points at which the values of the     *
!     *         function are calculated.  Xnew(k) must be ge xold(1)   *
!     *         and le xold(nold).                                     *
!     *  ynew - the values of the function to be calculated.           *
!     *  p, q - auxiliary vectors of the length nold-2.                *
!     *                                                                *
!     ******************************************************************
!
                             d i m e n s i o n
     2 xold(nold),yold(nold),y2(nold),p(nold),q(nold)
     3,xnew(nnew),ynew(nnew)
!-----------------------------------------------------------------------
      noldm1=nold-1
!
      dxl=xold(2)-xold(1)
      dxr=xold(3)-xold(2)
      dydxl=(yold(2)-yold(1))/dxl
      dydxr=(yold(3)-yold(2))/dxr
      rtdxc=.5/(dxl+dxr)
!
      p(1)= rtdxc*(6.*(dydxr-dydxl)-dxl*y2(1))
      q(1)=-rtdxc*dxr
!
      if(nold.eq.3) go to 700
!-----------------------------------------------------------------------
      k=3
!
 100  dxl=dxr

      dydxl=dydxr
      dxr=xold(k+1)-xold(k)
      dydxr=(yold(k+1)-yold(k))/dxr
      dxc=dxl+dxr
      den=1./(dxl*q(k-2)+dxc+dxc)
!
      p(k-1)= den*(6.*(dydxr-dydxl)-dxl*p(k-2))
      q(k-1)=-den*dxr
!
      k=k+1
      if(k.lt.nold) go to 100
!-----------------------------------------------------------------------
 700  k=noldm1
!
 200  y2(k)=p(k-1)+q(k-1)*y2(k+1)
!
      k=k-1
      if(k.gt.1) go to 200
!-----------------------------------------------------------------------
      k1=1
!
 300  xk=xnew(k1)
!
      do 400 k2=2,nold
      if(xold(k2).le.xk) go to 400
      kold=k2-1
      go to 450
 400  continue
      ynew(k1)=yold(nold)
      go to 600
!
 450  if(k1.eq.1)   go to 500
      if(k.eq.kold) go to 550
!
 500  k=kold
!
      y2k=y2(k)
      y2kp1=y2(k+1)
      dx=xold(k+1)-xold(k)
      rdx=1./dx
!
      ak=.1666667*rdx*(y2kp1-y2k)
      bk=.5*y2k
      ck=rdx*(yold(k+1)-yold(k))-.1666667*dx*(y2kp1+y2k+y2k)
!
 550  x=xk-xold(k)
      xsq=x*x
!
      ynew(k1)=ak*xsq*x+bk*xsq+ck*x+yold(k)
!
 600  k1=k1+1
      if(k1.le.nnew) go to 300
!-----------------------------------------------------------------------
      return
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

