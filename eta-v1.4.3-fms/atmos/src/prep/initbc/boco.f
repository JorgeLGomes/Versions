      subroutine eta_boco
c
c *** Creates the eta boundary condition file.
c
c *** Original code obtained from U. of Athens and modified at FSL.
c
      implicit none
c
      include 'ecommons.h'
c
      real*4 
     .    pdba(lb)  ,tba(lb,lm)  ,qba(lb,lm)  ,uba(lb,lm)  ,vba(lb,lm)
     .   ,pdbb(lb)  ,tbb(lb,lm)  ,qbb(lb,lm)  ,ubb(lb,lm)  ,vbb(lb,lm)  
     .   ,pdb (lb,2),tb (lb,lm,2),qb (lb,lm,2),ub (lb,lm,2),vb (lb,lm,2)
     .   ,q2b(lb,lm,2),cwmb(lb,lm,2)
     .   ,pdb2(lb,2),tb2(lb,lm,2),qb2(lb,lm,2),ub2(lb,lm,2),vb2(lb,lm,2)
     .   ,rtboco,run
c
      integer*4 itboco,ibc,nibc,idat(3),ihrst,i,j,k,l,len,n
c
      character*255 fname1,fname2,fname
      character*8   ftime
      character*8   ced
c
      real*4 epsq2
      parameter (epsq2=0.2)
c_______________________________________________________________________________
c
      print *,' '
      print *,'Creating eta boundary condition files...'
      print *, 'these are dimensioned as lb, which is= ', lb
      print *, 'kb,im,jm,imjm,imt,jmt= ',kb,im,jm,imjm,imt,jmt
c
      len=index(init_out,' ')-1
      if (init_out(len:len) .ne. '/') then
         len=len+1
         init_out(len:len)='/'
      endif
c
c *** Open boundary condition file.
c     Bndy conditions will be written sequentially into one file (NCEP style).
c
      write(ftime,'(4i2.2)') mod(iyear,100),imonth,idate,istrtim
Cmp      l=len+13
        l=len+9
Cmp      fname(1:l)=init_out(1:len)//'bndy.'//ftime
      fname(1:l)=init_out(1:len)//'BNDY.file'
      print *,'Open BC file  : ',fname(1:l)
      open(2,file=fname(1:l),status='unknown',form='unformatted')
c
c *** Loop through for all boco files.
c
      rtboco=1./(tboco*3600.)                                           
      itboco=nint(tboco)
      nibc=nhour/itboco
      do ibc=1,nibc
c
         if (ibc .eq. 1) then
            write(ced,'(i6.6)') (ibc-1)*itboco
            l=len+17
            fname1(1:l)=init_out(1:len)//'preproc.bc.'//ced
c
            print *,'Read from file: ',fname1(1:l)
            open(1,file=fname1(1:l),status='old',form='unformatted')
            read(1) run,idat,ihrst,pdb,tb,qb,ub,vb
            close(1)
        write(6,*) 'done reading from file'
         else
            do k=1,2
               do i=1,lb
                  pdb(i,k)=pdb2(i,k)
               enddo
               do j=1,lm
               do i=1,lb
                  tb(i,j,k)=tb2(i,j,k)
                  qb(i,j,k)=qb2(i,j,k)
                  ub(i,j,k)=ub2(i,j,k)
                  vb(i,j,k)=vb2(i,j,k)
               enddo
               enddo
            enddo
         endif
c
         write(ced,'(i6.6)') ibc*itboco
         l=len+17
         fname2(1:l)=init_out(1:len)//'preproc.bc.'//ced
c
         print *,'Read from file: ',fname2(1:l)
         open(1,file=fname2(1:l),status='old',form='unformatted')
         read(1) run,idat,ihrst,pdb2,tb2,qb2,ub2,vb2
         close(1)
c
         do n=1,lb
            pdbb(n)=pdb (n,1)
            pdba(n)=pdb2(n,1)
         enddo
c
         do l=1,lm
         do n=1,lb
            tbb(n,l)=tb (n,l,1)
            tba(n,l)=tb2(n,l,1)
            qbb(n,l)=qb (n,l,1)
            qba(n,l)=qb2(n,l,1)
            ubb(n,l)=ub (n,l,1)
            uba(n,l)=ub2(n,l,1)
            vbb(n,l)=vb (n,l,1)
            vba(n,l)=vb2(n,l,1)
         enddo
         enddo
c
         do n=1,lb                                                 
            pdb(n,1)=pdbb(n)                                           
            pdb(n,2)=(pdba(n)-pdbb(n))*rtboco                           
         enddo
c
         do l=1,lm                              
         do n=1,lb                                                 
            tb(n,l,1)=tbb(n,l)                                         
            qb(n,l,1)=qbb(n,l)                                          
            ub(n,l,1)=ubb(n,l)                                          
            vb(n,l,1)=vbb(n,l)                                          
c                                                                       
            tb(n,l,2)=(tba(n,l)-tbb(n,l))*rtboco                       
            qb(n,l,2)=(qba(n,l)-qbb(n,l))*rtboco                       
            ub(n,l,2)=(uba(n,l)-ubb(n,l))*rtboco                        
            vb(n,l,2)=(vba(n,l)-vbb(n,l))*rtboco                        
         enddo
         enddo
c
         write(ftime,'(4i2.2)') mod(iyear,100),imonth,idate,istrtim
Cmp         l=len+13
         l=len+9
Cmp         fname(1:l)=init_out(1:len)//'bndy.'//ftime
         fname(1:l)=init_out(1:len)//'BNDY.file'
         print *,'Write to file : ',fname(1:l)
         open(2,file=fname(1:l),status='unknown',form='unformatted')
         if (ibc .eq. 1) write(2) run,idat,ihrst,tboco
         write(2) float(ibc-1)*tboco
         write(2) pdb
         write(2) tb
         write(6,*) 'pdb(1) ', pdb(1,1),pdb(1,2)
         write(6,*) 'pdb(lb/2) ', pdb(lb/2,1),pdb(lb/2,2)
         write(6,*) 'pdb(last) ', pdb(lb,1),pdb(lb,2)
         write(6,*) 't(1,lm/2) ', tb(1,lm/2,1),tb(1,lm/2,2)
         write(6,*) 't(lb,lm/2) ', tb(lb,lm/2,1),tb(lb,lm/2,2)
         write(2) qb
         write(2) ub
         write(2) vb


        write(1000,*) 'T Q U V '
        do i=1,lb
        write(1000,633) i,tb(i,10,1),qb(i,10,1)*1.e6,ub(i,10,1),
     +             vb(i,10,1),(ub(i,10,1)**2.+ vb(i,10,1)**2.) ** 0.5
        enddo
  633   format(I10,1x,5(f6.1,1x))


Cmp         do k=1,2
         do j=1,lm
         do i=1,lb
            q2b(i,j,1)=epsq2
            q2b(i,j,2)=0.
            cwmb(i,j,1)=0.
            cwmb(i,j,2)=0.
         enddo
         enddo
Cmp         enddo
         write(2) q2b
         write(2) cwmb 
c
      enddo
      write(2) float(nibc)*tboco
Cmp
Cmp	duplicate last tendency
        write(6,*) 'duping last tendency with some modification'

C       tendencies dont matter here, but the value does.
c
C       want to dump the terminal value here, not the "before" value!
C
         do n=1,lb
            pdb(n,1)=pdba(n)
            pdb(n,2)=(pdba(n)-pdbb(n))*rtboco
         enddo
c
         do l=1,lm
         do n=1,lb
            tb(n,l,1)=tba(n,l)
            qb(n,l,1)=qba(n,l)
            ub(n,l,1)=uba(n,l)
            vb(n,l,1)=vba(n,l)
c
            tb(n,l,2)=(tba(n,l)-tbb(n,l))*rtboco
            qb(n,l,2)=(qba(n,l)-qbb(n,l))*rtboco
            ub(n,l,2)=(uba(n,l)-ubb(n,l))*rtboco
            vb(n,l,2)=(vba(n,l)-vbb(n,l))*rtboco
         enddo
         enddo
c
         write(2) pdb
         write(2) tb
         write(2) qb
         write(2) ub
         write(2) vb
         write(2) q2b
         write(2) cwmb


        write(6,*) 'pdb(1) ', pdb(1,1),pdb(1,2)
        write(6,*) 'pdb(last) ', pdb(lb,1),pdb(lb,2)
        write(6,*) 't(1,lm/2) ', tb(1,lm/2,1),tb(1,lm/2,2)
        write(6,*) 't(lb,lm/2) ', tb(lb,lm/2,1),tb(lb,lm/2,2)
C	write(6,*) 'u(20,lm/2) ', ub(20,lm/2,1),ub(20,lm/2,2)
C	write(6,*) 'v(20,lm/2) ', vb(20,lm/2,1),vb(20,lm/2,2)
cmp

      print *,'Close BC file : ',fname(1:l)
      close(2)
c
      return
      end                                        
