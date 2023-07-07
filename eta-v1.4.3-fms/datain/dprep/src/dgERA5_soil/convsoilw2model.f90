      subroutine conv2model(gribfile,outdir)
      use estab
      use constants
      implicit none
      character(len=255),intent(in)   :: gribfile
      character(len=255),intent(in)   :: outdir
!_____________________________________________________________________________________________________________
      integer                         :: nx,ny,nz
      integer                         :: dx,dy 
      integer                         :: proj,LatS,LonW,LatN,LonE
      integer                         :: i,j,k,l,n,it,nsfcfld,id
      integer                         :: n1,n2,n3,LenStr,iter
      integer                         :: icall,nrec
      real                            :: xe,mrsat,qs,pr_Pa
      integer, dimension(200)         :: kgds
!_____________________________________________________________________________________________________________
      real,    dimension(:,:,:) ,allocatable    :: slp            !slp(i,j,1)=LCMSK;slp(i,j,2)=PSLM
                                                                  !slp(i,j,3)=PSFC;slp(i,j,4)=ZSFC
                                                                  !slp(i,j,5 & 6) are 0-10 cm STC and SMC
                                                                  !slp(i,j,7 & 8) are 10-200 cm STC and SMC
!_____________________________________________________________________________________________________________
      real,    dimension(:,:)   ,allocatable    :: dummy 
!_____________________________________________________________________________________________________________
      character(len=255)              :: outfile,gdsfile
      character(len=20)               :: ModelDriver
      character(len=20)               :: atime
      character(len=8)                :: atime_r
      character(len=7)                :: model
      character(len=8)               :: timex
!_____________________________________________________________________________________________________________

      open(1,file='./InputModelInf.txt',form='formatted',status='old')
      rewind 1
      read(1,'(a20)')ModelDriver
      read(1,*)nx,ny,nz
      read(1,*)proj
      read(1,*)LatS
      read(1,*)LonW
      read(1,*)LatN
      read(1,*)LonE
      read(1,*)dx
      read(1,*)dy

      allocate(slp(nx,ny,12)) 
      allocate(dummy(nx,ny))  
      close(1)
!

! *** Fill pressure levels.
!
!mp      nsfcfld will be 12 although 4 of these are set to missing for...
      nsfcfld=12
!mp
        write(6,*) 'using hardwired stuff! ', nx,ny,nz
! *** Read in degrib data.
!
        n=index(outdir,' ')-1
        n2=index(ModelDriver,' ')-1
        write(6,*)ModelDriver(1:n2)

        open(12,file='./adate.txt')
        read(12,*) atime_r
        close(12)


      n=index(gribfile,' ')-1
       write(*,*) 'OPEN FILE ',gribfile(1:n)

      open(12,file=gribfile(1:n),form='unformatted',status='old',access='direct',recl=nx*ny*4)
      nrec=1
!Soil Temperature  layer 1   
      read(12,rec=nrec)((slp(i,j,5),i=1,nx),j=ny,1,-1) 
      nrec=nrec+1
!Soil Moisture	 layer 1     
      read(12,rec=nrec)((slp(i,j,6),i=1,nx),j=ny,1,-1) 
      nrec=nrec+1
!Soil Temperature layer 2
      read(12,rec=nrec)((slp(i,j,7),i=1,nx),j=ny,1,-1) 
      nrec=nrec+1
!Soil Moisture	layer 2
      read(12,rec=nrec)((slp(i,j,8),i=1,nx),j=ny,1,-1) 
      nrec=nrec+1
!Soil Temperature layer 3
      read(12,rec=nrec)((slp(i,j,9),i=1,nx),j=ny,1,-1) 
      nrec=nrec+1
!Soil Moisture layer 3
      read(12,rec=nrec)((slp(i,j,10),i=1,nx),j=ny,1,-1) 
      nrec=nrec+1
!Soil Temperature layer 4
      read(12,rec=nrec)((slp(i,j,11),i=1,nx),j=ny,1,-1)  
      nrec=nrec+1
!Soil Moisture layer 4
      read(12,rec=nrec)((slp(i,j,12),i=1,nx),j=ny,1,-1) 

      
      close(12)
      slp(:,:,4)=slp(:,:,4)/9.8

      where(slp(:,:,1)<0.1)
        slp(:,:,1)=0
      elsewhere(slp(:,:,1)>=0.1)
        slp(:,:,1)=1
      end where

      where(ht==9.9990003E+20)ht=-99999.
      where(uw==9.9990003E+20)uw=-99999.
      where(vw==9.9990003E+20)vw=-99999.
      where(tp==9.9990003E+20)tp=-99999.
      where(rh==9.9990003E+20)rh=-99999.
      where(q==9.9990003E+20)q=-99999.
      qnew=-99999.
      qnew2=-99999.
      where(slp==9.9990003E+20) slp=0.
!
! *** Check for any missing data.
      write(*,*) "Check for any missing data"
!
      do k=1,nz
          if (ht(1,1,k) .eq. -99999.) then
             print *,'Height data missing at level: ',pr(k)
             stop
          elseif (tp(1,1,k) .eq. -99999.) then
             print *,'Temperature data missing at level: ',pr(k)
             stop
          elseif (rh(1,1,k) .eq. -99999.) then
             print *,'RH data missing at level: ',pr(k)
             print *,'Calling RH patch.'
             call rh_fix(nx,ny,nz,rh,pr)
          elseif (uw(1,1,k) .eq. -99999.) then
             print *,'U-wind data missing at level: ',pr(k)
             stop
          elseif (vw(1,1,k) .eq. -99999.) then
             print *,'V-wind data missing at level: ',pr(k)
             stop
          elseif (q(1,1,k) .eq. -99999.) then
             print *,'Spec. Humid. data missing at level: ',pr(k)
             stop
          endif
      enddo

!mp
!mp     Handle missing surface data....
!mp

!
!      Just bogus to hold a place.
!
        write(6,*) 'going to use nsfcfld= ', nsfcfld

!	If rean put bogus data into slp
        where(slp(:,:, 6)>0.6)  slp(:,:, 6)=0.3
        where(slp(:,:, 8)>0.6)  slp(:,:, 8)=0.3
        where(slp(:,:,10)>0.6)  slp(:,:,10)=0.3
        where(slp(:,:,12)>0.6)  slp(:,:,12)=0.3

        do k=5,nsfcfld
          if (slp(1,1,k) .eq. -99999.) then
          write(6,*) 'SURFACE DATA MISSING... ', K       
          if (k.eq.9.or.k.eq.11) then
           write(6,*) 'filling soil temp data...'
           do j=1,ny
            do i=1,nx
             slp(i,j,k)=slp(i,j,7)
            enddo
           enddo
          elseif(k.eq.10.or.k.eq.12.) then
           write(6,*) 'filling soil moisture data...' 
           do j=1,ny
             do i=1,nx
              slp(i,j,k)=slp(i,j,8)
             enddo
           enddo
          endif
         endif
 
         do iter=1,15
          do j=1,ny
           do i=1,nx
            if (slp(i,j,k) .eq. 0 .and. slp(i+1,j,k) .ne. 0) then
             slp(i,j,k)=slp(i+1,j,k)
            endif
           enddo
          enddo
          do j=1,ny
           do i=2,nx
            if (slp(i,j,k) .eq. 0 .and. slp(i-1,j,k) .ne. 0) then
             slp(i,j,k)=slp(i-1,j,k)
            endif
           enddo
          enddo
         enddo
        enddo


!
! *** Convert 3d temp to theta.
! *** Compute Exner function.
! *** Convert 3d rh to mr.
!
      do k=1,nz
         pri(k)=1./pr(k)
      enddo
!
!JLG qnew nao usada - comentada em 20220824
!JLG motivo: acelerar o código
!      do k=1,nz
!      do j=1,ny
!      do i=1,nx
!        icall = 0
!        pr_Pa=pr(k)*1E+2
!        call TQadj (tp(i,j,k),qs,pr_Pa,icall)
!        qnew(i,j,k)=rh(i,j,k)*qs*1E-2
!      enddo
!      enddo
!      enddo
!JLG Fim 
      do k=1,nz
      do j=1,ny
      do i=1,nx
         th(i,j,k)=tp(i,j,k)*(1000.*pri(k))**kappa
!	 ex(i,j,k)=cp*tp(i,j,k)/th(i,j,k)
         it=tp(i,j,k)*100
         it=min(45000,max(15000,it))
         xe=esat(it)
         mrsat=0.00622*xe/(pr(k)-xe)
         qnew2(i,j,k)=rh(i,j,k)*mrsat
      enddo
      enddo
      enddo
!
! *** Create output file name.
!

      model='ETAwrk'
      
      n=index(gribfile,' ')-1
      outfile=gribfile(1:n)//'.'//model
      print*,outfile
      n1=index(outfile,' ')-1
      print *,model,' data ---> ',outfile(1:n1)
      open(1,file=outfile(1:n1),status='unknown',form='unformatted')
        
      gdsfile=outfile(1:n1)//'.gdsinfo'
      n2=index(gdsfile,' ')-1
      open(unit=14,file=gdsfile(1:n2),form='formatted',status='unknown')
!
!

!      open(unit=1010,file=outfile(1:n1)//'.qsin',form='formatted',status='unknown')
!      do k=1,nz
!         write(1010,*) "k ",k,"p ",pr(k),"qold ",q(nx/2,ny/2,k),"qnew ",qnew(nx/2,ny/2,k),"qnew2 ",qnew2(nx/2,ny/2,k)
!      enddo
!      close(1000)

!mp      nsfcfld will be 12 although 4 of these are set to missing for AVN...
      nsfcfld=12
!      write(1) nz

! new stuff for GDS file
!
        n=index(gribfile,' ')-1
        write(6,*) 'calling get_fullgds with ',gribfile(1:n)
!
        kgds(1)=proj
        kgds(2)=nx
        kgds(3)=ny
        kgds(4)=LatS
        kgds(5)=LonW
        kgds(6)=128
        kgds(7)=LatN
        kgds(8)=LonE
        kgds(9)=dx
        kgds(10)=dy
        kgds(11)=64
        kgds(20)=255

        write(6,*) 'back from get_fullgds'



!
!	MODIFY GDS TO REFLECT FLIPPING OF DATA IN N-S SENSE
!
      
!	kgds(4)= -kgds(4)
!	kgds(7)= -kgds(7)
        write(6,*) 'writing kgds ', (kgds(I),i=1,14)
        write(6,*) 'GDS written to ', gdsfile(1:n)
        write(14,*) kgds(1)
        write(14,*) kgds(2)
        write(14,*) kgds(3)
        write(14,*) kgds(4)
        write(14,*) kgds(5)
        write(14,*) kgds(6)
        write(14,*) kgds(7)
        write(14,*) kgds(8)
        write(14,*) kgds(9)
        write(14,*) kgds(10)
        write(14,*) kgds(11)
        write(14,*) kgds(20)
        close(14)
	write(6,*) 'past 14 close'
!
! end new stuff

!
! *** Write isobaric upper air data.
!
	write(6,*) 'starting writes to unit1'

      write(1) uw
      write(1) vw
      write(1) ht
!      write(1) rh
      write(1) qnew2
      write(1) tp
!     write(1) pr
      write(1) slp
!mp
	write(6,*) 'writing these sfc values to the ',model,' file...'
	do k=1,nsfcfld
	write(6,*) 'field, value ', k, slp(nx/2,ny/2,k)
	enddo
	do k=1,nsfcfld
	write(6,*) 'field, value ', k, slp(3*nx/4,3*ny/4,k)
	enddo

!
      close(1)
!
      return
      end
!
