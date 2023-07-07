      program degrib2model_driver
c
      implicit none
c
      integer nx,ny,nz          !Degrib grid dimensions
c
      character(LEN=255):: degrib_dir,degrib_file,outdir,gribfile
c
      real esat,es
      integer datsav,n,l,kpds(200)
c
      common /estab/esat(15000:45000),es(15000:45000)
      common /gdsinfo/datsav(11)
c_______________________________________________________________________________
c
c *** Initialize tables.
c
      call es_ini
c
c *** Read file names.
c
        read(5,'(a)') gribfile
        read(5,'(a)') outdir
	
      n=index(gribfile,' ')-1
      write(6,*) 'GRIBFILE is ', gribfile(1:n)
      CALL GET_GDS(gribfile(1:n),datsav,kpds)

c
C     X and Y dimensions now obtained from the degribbed GRIB file
c
       nx=datsav(1)
       ny=datsav(2)
       nz=22
c
       write(6,*) 'calling degrib2model ', nx,ny,nz
       call degrib2model(gribfile,outdir,nx,ny,nz)
c
      end
c
c===============================================================================
c
      subroutine degrib2model(gribfile,outdir,nx,ny,nz)
c
      implicit none
c
      real    cp,kappa
      parameter(cp=1004.,kappa=287./1004.)
c
      integer nx,ny,nz,i,j,k,l,n,n2,it,ip,jp,nsfcfld
     .         ,iyear,imonth,iday,ifcsthr
     .	       ,datsav,ii,jj
      integer kgds(200)

	real crot(nx*ny),srot(nx*ny),tbar,old,fact,biassum
	real z1000t,z975t
c
      real ht(nx,ny,nz)              !Isobaric heights (m)
     .      ,tp(nx,ny,nz)              !Isobaric temps (K)
     .      ,th(nx,ny,nz)              !Isobaric theta (K)
     .      ,uw(nx,ny,nz)              !Isobaric u-wind (m/s)
     .      ,vw(nx,ny,nz)              !Isobaric v-wind (m/s)
     .      ,rh(nx,ny,nz)              !Isobaric rh,mr (%,kg/kg)
     .      ,pr(nz),pri(nz)            !Isobaric pressures (mb)
     .      ,lat1,lat2,lon0,sw(2),ne(2)
     .      ,etalevs(22)
     .      ,xe,mrsat,esat,es
     .	    ,slp(nx,ny,12) 	  !slp(i,j,1)=EMSL;slp(i,j,2)=PMSL
				  !slp(i,j,3)=PSFC;slp(i,j,4)=ZSFC
c
      character(LEN=255):: degrib_dir,degrib_file
      character(LEN=255):: outdir,outfile,gribfile
      character(LEN=255):: gdsfile
      character(LEN=2)  :: gproj
      character(LEN=12) ::   atime
      character(LEN=6)  ::  model

	real urlat,urlon
c

       data ETALEVS/1020., 1000., 950., 925., 900., 
     +               850., 800., 750., 700., 650., 
     +               600., 550., 500., 450., 400., 
     +               350., 300., 250., 200., 150., 
     +               100., 50./

      common /estab/esat(15000:45000),es(15000:45000)
      common /gdsinfo/datsav(11)
c_______________________________________________________________________________
c
c *** Fill pressure levels.
c
Cmp
	nsfcfld=12

Cmp     pressure levels are hardwired here

      do k=1,nz
         pr(k)=ETALEVS(k)
	 write(6,*) 'PRESS ',pr(k)
      enddo

c 
c *** Read in degrib data.
c
      n=index(gribfile,' ')-1
	write(6,*) 'calling read_degrib'


      call read_degrib(gribfile(1:n)
     .                ,nx*ny,nz,pr,ht,tp,rh,uw,vw,slp,atime)

	write(6,*) 'BACK FROM read_degrib'

c *** Check for any missing data.
c
      do k=1,nz
	  if (ht(1,1,k) .eq. -99999.) then
	     print *,'Height data missing at level: ',pr(k)
	     stop
	  elseif (tp(1,1,k) .eq. -99999.) then
	     print *,'Temperature data missing at level: ',pr(k)
	     stop
	  elseif (rh(1,1,k) .eq. -99999.) then
	     print *,'RH data missing at level: ',k,pr(k)
	if (datsav(11) .ne. 104) then
	     print *,'Calling RH patch.'
	     call rh_fix(nx,ny,nz,rh,pr)
	else

Cmp     moisture every 50 for grid 104...handle in special way

        if (mod(k,2) .eq. 0 .or. k.eq.nz) then

        if (k.lt.(nz-1)) then
        write(6,*) 'ave value from ', pr(K+1),pr(K-1),k+1,k-1
        DO J=1,NY
        DO I=1,NX
        RH(I,J,K)=(RH(I,J,K+1)+RH(I,J,K-1))*0.5
        ENDDO
        ENDDO
        else
        write(6,*) 'copying value from ', pr(K-1),k-1
        DO J=1,NY
        DO I=1,NX
        RH(I,J,K)=RH(I,J,K-1)
        ENDDO
        ENDDO
        endif

        endif

	endif
	
	  elseif (uw(1,1,k) .eq. -99999.) then
	     print *,'U-wind data missing at level: ',pr(k)
	     stop
	  elseif (vw(1,1,k) .eq. -99999.) then
	     print *,'V-wind data missing at level: ',pr(k)
	     stop
	  endif
      enddo

C***
C***	Change 1000 hPa heights to match lowest temps
C***
	biassum=0.

	do J=1,NY
	do I=1,NX
	
C	old=ht(i,j,1)
	fact=287.*(pr(2)-pr(3))/(9.81*pr(2))
	z975t=(ht(i,j,3)-ht(i,j,2))/fact
	fact=287.*(pr(1)-pr(2))/(9.81*pr(1))
	z1000t=(ht(i,j,2)-ht(i,j,1))/fact

C
C this criteria empirically derived.  Was found that where 950-975 thick > 
C about 7 m greater than 975-1000 thick, lowest-lev temperatures became hosed.
C
	if (z975t-z1000t .gt. 2.5) then
C	tbar=(1./4.)*(tp(i,j,1)+tp(i,j,2)+tp(i,j,3)+tp(i,j,4))
	tbar=amax1(tp(i,j,1),tp(i,j,2),tp(i,j,3),tp(i,j,4))
	ht(i,j,1)=ht(i,j,2)-fact*tbar
	biassum=biassum+1
	endif
	
	
  633	format(7(f6.2,1x))
	enddo
	enddo

C	call smooth(ht,nx,ny,nz,2)

	write(6,*) 'modified Z1000 at ', biassum, ' points'

C***
C***
C***


Cmp
Cmp	Handle missing surface data....
Cmp
	do k=1,12
	  if (slp(1,1,k) .eq. -99999.) then
	    write(6,*) 'SURFACE DATA MISSING... ', K
	    if (datsav(11) .ne. 104) then
              if (k.eq.9.or.k.eq.11.) then
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
	    elseif (datsav(11) .eq. 104) then
	      if (k.eq.7 .or. k.eq.9 .or. k.eq.11.) then
                write(6,*) 'filling soil temp data...'
                do j=1,ny
                  do i=1,nx
                    slp(i,j,k)=slp(i,j,5)
                  enddo
                enddo
              elseif(k.eq.8 .or. k.eq.10.or.k.eq.12.) then
                write(6,*) 'filling soil moisture data...'
                do j=1,ny
                  do i=1,nx
                    slp(i,j,k)=slp(i,j,6)
                  enddo
                enddo
              endif  
	    endif
	  endif
	enddo
		
c
c *** Convert 3d temp to theta.
c *** Compute Exner function.
c *** Convert 3d rh to mr.
c
      do k=1,nz
        pri(k)=1./pr(k)
      enddo
c
      do k=1,nz
        do j=1,ny
          do i=1,nx
            th(i,j,k)=tp(i,j,k)*(1000.*pri(k))**kappa
C           ex(i,j,k)=cp*tp(i,j,k)/th(i,j,k)
            it=tp(i,j,k)*100
            it=min(45000,max(15000,it))
            xe=esat(it)
            mrsat=0.00622*xe/(pr(k)-xe)
            rh(i,j,k)=rh(i,j,k)*mrsat
          enddo
        enddo
      enddo
c
      print *,'ua :',ht(nx/2,ny/2,nz),th(nx/2,ny/2,nz),rh(nx/2,ny/2,nz)
     .       ,uw(nx/2,ny/2,nz),vw(nx/2,ny/2,nz)
      print *,'ua :',ht(nx/2,ny/2,1),th(nx/2,ny/2,1),rh(nx/2,ny/2,1)
     .       ,uw(nx/2,ny/2,1),vw(nx/2,ny/2,1)

Cmp
Cmp     Write out every 10th grid point in both directions?

Cnot    goto 1076
        write(6,*) 'level 10', nx,ny

        write(6,*) 'mixr * 1000.'
        do J=ny,1,-(ny/15)
          write(6,744) (rh(I,J,10)*1000.,I=1,nx,nx/9)
        enddo

        write(6,*) 'HT'
        do J=ny,1,-(ny/15)
          write(6,744) (ht(I,J,10),I=1,nx,nx/9)
        enddo

        write(6,*) 'UW'
        do J=ny,1,-(ny/15)
          write(6,744) (uw(I,J,10),I=1,nx,nx/9)
        enddo

        write(6,*) 'VW'
        do J=ny,1,-(ny/15)
          write(6,744) (vw(I,J,10),I=1,nx,nx/9)
        enddo

        write(6,*) 'TH'
        do J=ny,1,-(ny/15)
          write(6,744) (th(I,J,10),I=1,nx,nx/9)
        enddo

 1076   continue

  743   format(30(e9.3,1x))
  744   format(30(f7.1,1x))

	
c
c *** Create output file name.
c
CGSM      n=index(outdir,' ')-1
      n=index(gribfile,' ')-1
	if (datsav(11) .eq. 221) model='ETA221'
	if (datsav(11) .eq. 212) model='ETA212'
	if (datsav(11) .eq. 104) model='ETA104'
	if (datsav(11) .eq. 255) model='ETAwrk'
CGSM      outfile=outdir(1:n)//'/'//atime//'.'//model
      outfile=gribfile(1:n)//'.'//model
      n=index(outfile,' ')-1
      print *,model,' data ---> ',outfile(1:n)
      open(1,file=outfile(1:n),status='unknown',form='unformatted')
c
c *** Write header stuff.
c
C
      nsfcfld=12

	if (datsav(10) .eq. 3) gproj='LC'
	if (datsav(10) .eq. 5) gproj='PS'
	if (datsav(10) .eq. 0) gproj='CE'
	write(1) nz

        n2=index(gribfile,' ')-1
	write(6,*) 'calling get_fullgds with ',
     +    gribfile(1:n2)	 
	call get_fullgds(gribfile(1:n2),datsav(1),datsav(2),kgds)
	write(6,*) 'back from get_fullgds'

C new stuff for GDS file
C
CGSM        n=index(outdir,' ')-1
CGSM	gdsfile=outdir(1:n)//'/'//'gdsinfo.'//model
	gdsfile=outfile(1:n)//'.gdsinfo'
	n=index(gdsfile,' ')-1
	write(6,*) 'GDS written to ', gdsfile(1:n)
	open(unit=14,file=gdsfile(1:n),form='unformatted',
     +	access='sequential')
     
	write(6,*) 'writing kgds ', (kgds(I),i=1,14)
	write(14) kgds

	close(14)
C
C end new stuff

	write(6,*) '********************************************'
	write(6,*) 'writing out the grid with the following info'
	write(6,*) 'dimensions nx,ny,nz ', nx,ny,nz
	write(6,*) '********************************************'

c
c *** Write isobaric upper air data.
c
      write(1) uw
      write(1) vw
      write(1) ht
      write(1) rh
      write(1) pr
      write(1) slp

 	write(6,*) 'writing these sfc values to the ETA file...'
        write(6,*) 'at center of input grid'
        do k=1,nsfcfld
        write(6,*) 'field, value ', k, slp(nx/2,ny/2,k)
        enddo

c
      close(1)
c
      return
      end
c
c===============================================================================
c
      subroutine read_degrib(etafile,nxny,nz,pr
     .                      ,ht,tp,rh,uw,vw,slp,atime)
c
      implicit none
c
      integer nxny,nz,i,datsav
c
      real pr(nz),ht(nxny,nz),tp(nxny,nz)
     .      ,rh(nxny,nz),uw(nxny,nz),vw(nxny,nz)
     .      ,dummy,slp(nxny,12),tmp(nxny)
      real crot(nxny),srot(nxny),rmagb,rmagaft,ubef,vbef,urlat,
     .	    urlon
c
      integer kgds(200),kpds(50),len,jpds(200),jgds(200)
     .         ,lenpds,lenkgds,nwords,kpdsl
     .        ,j,k,nx,IRETO,KNUM,IRET1
c
      character*(*) etafile
      character(LEN=1):: pds(50)
      character(LEN=12):: atime
      character(LEN=2):: gproj
      logical bitmap(nxny)

      common /gdsinfo/datsav(11)
c_______________________________________________________________________________
c
      len=index(etafile//' ',' ')-1
c
c *** Check that the input dimensions match grib file dimensions. 
c
	write(6,*) 'inside read_degrib '

CCC
CCC  Should be able to derive all header info in the DEGRIB file
CCC from the GDS/PDS of the GRIB file.  
CCC
CCC
c
c *** Fill time stamp (assume all records are for same time).

	call get_gds(etafile(1:len),datsav,kpds)

       	write(6,*) 'KPDS vals ', kpds(10),kpds(9),kpds(8)

	if (kpds(8) .ge. 100) kpds(8)=kpds(8)-100

      write(atime,'(i2.2,i2.2,i2.2,i2.2,i4.4)') 
     .   kpds(8),kpds(9),kpds(10),kpds(11),kpds(14)
c
c
c *** Fill a missing value flag into first space of each variable.
c
      do k=1,nz
	 ht(1,k)=-99999.
	 tp(1,k)=-99999.
	 rh(1,k)=-99999.
	 uw(1,k)=-99999.
	 vw(1,k)=-99999.
      enddo

Cmp initialize surface fields to -99999.  so the interp code can handle
Cmp	appropriately 

	do k=1,12
	do j=1,nxny
	slp(j,k)=-99999.
	enddo
	enddo

Cmp
	
c
c *** Now put the data into the corresponding arrays.
c
Cmp
C  add something to read in surface fields in here
C

	call baopen(11,etafile,IRETO)
	if (IRETO .ne. 0) write(6,*) 'BAOPEN TROUBLE!!!! ', IRETO

	jpds=-1
	jgds=-1

C############# LAND/SEA ##############################
	jpds(5)=81
	jpds(6)=1
	jpds(7)=0

	call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,slp(1,1),IRET1)

	write(6,*) 'first GETGB, IRET1= ', IRET1
	write(6,*) 'LAND/SEA READ!!!!! '
	write(6,*) (slp(j,1),j=nwords/2,nwords/2+5)
C############# LAND/SEA ##############################


C############# PMSL ##############################
	jpds(5)=136
	jpds(6)=1
	jpds(7)=0

       call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,slp(1,2),IRET1)
        slp(:,2)=slp(:,2)
	write(6,*) 'PMSL READ!!!!! '
	write(6,*) (slp(j,2),j=nwords/2,nwords/2+8)

C############# PMSL ##############################


C############# PSFC ##############################
	jpds(5)=135
	jpds(6)=1
	jpds(7)=0

       call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,slp(1,3),IRET1)
        slp(:,3)=slp(:,3)
        write(6,*) 'PSFC READ!!!!! '
        write(6,*) (slp(j,3),j=nwords/2,nwords/2+8)

C############# PSFC ##############################

C############# ZSFC ##############################
	jpds(5)=132
	jpds(6)=1
	jpds(7)=0

       call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,slp(1,4),IRET1)

        write(6,*) 'ZSFC READ!!!!! '
        write(6,*) (slp(j,4),j=nwords/2,nwords/2+8)

C############# ZSFC ##############################

C
Cmp     SOIL FIELDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C

C	if (mod(kpds(7)-i,256).eq.0) then
Cmp	write(6,*) 'kpds(5)= ', kpds(5)
Cmp	write(6,*) 'lower is ', i
Cmp	write(6,*) 'upper is ', (kpds(7)-i)/256.
C	endif

	jpds(5)=191
	jpds(6)=1
	jpds(7)=0

       call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,slp(1,5),IRET1)

	if (IRET1. eq. 0) then
	write(6,*) 'found soil temp over ',jpds(7), 'layer!!'
        write(6,*) (slp(j,5),j=nwords/2,nwords/2+8)
	endif
        
	jpds(5)=191
	jpds(6)=1
	jpds(7)=0
       call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,slp(1,7),IRET1)
	if (IRET1. eq. 0) then
	write(6,*) 'found soil temp over ',jpds(7), 'layer!!', IRET1
        write(6,*) (slp(j,7),j=nwords/2,nwords/2+8)
	endif
	
	jpds(5)=175
	jpds(6)=1
	jpds(7)=0
       call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,slp(1,9),IRET1)
	if (IRET1. eq. 0) then
	write(6,*) 'found soil temp over ',jpds(7), 'layer!!',IRET1
        write(6,*) (slp(j,9),j=nwords/2,nwords/2+5)
	endif

	jpds(5)=175
        jpds(6)=1
        jpds(7)=0
       call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,slp(1,11),IRET1)
        if (IRET1. eq. 0) then
	write(6,*) 'found soil temp over ',jpds(7), 'layer!!',IRET1
        write(6,*) (slp(j,11),j=nwords/2,nwords/2+5)
	endif

C 	MOISTURE

	jpds(5)=182
	jpds(6)=1
	jpds(7)=0

       call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,slp(1,6),IRET1)

	if (IRET1. eq. 0) then
        write(6,*) 'found soil wet over ',jpds(7), 'layer!!', IRET1
        write(6,*) (slp(j,6),j=nwords/2,nwords/2+8)
	endif
        jpds(5)=182
	jpds(6)=1
        jpds(7)=0
       call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,slp(1,8),IRET1)
        if (IRET1. eq. 0) then
        write(6,*) 'found soil wet over ',jpds(7), 'layer!!', IRET1
        write(6,*) (slp(j,8),j=nwords/2,nwords/2+8)
	endif
        jpds(5)=183
	jpds(6)=1
        jpds(7)=0
       call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,slp(1,10),IRET1)
        if (IRET1. eq. 0) then
        write(6,*) 'found soil wet over ',jpds(7), 'layer!!',IRET1
        write(6,*) (slp(j,10),j=nwords/2,nwords/2+5)
	endif

        jpds(5)=183
	jpds(6)=1
	jpds(7)=0
       call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,slp(1,12),IRET1)
        if (IRET1. eq. 0) then
        write(6,*) 'found soil wet over ',jpds(7), 'layer!!',IRET1
        write(6,*) (slp(j,12),j=nwords/2,nwords/2+5)
	endif


C	jpds(5)=7
	jpds(6)=100

Cmp	doing this in proper order??????
Corig	do k=1,nz
	do k=nz,1,-1
	jpds(7)=nint(pr(k))
        write(6,*) 'LEVEL ', jpds(7)

C######################### RH ###################
	jpds(5)=52
        jpds(6)=100
      call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,rh(1,k),IRET1)
	if (rh(20,k).gt.1.e20) write(6,*) 'getgbd(1) rh of ',k,rh(20,k)
C######################### RH ###################
	
C######################### HT ###################
        jpds(5)=7
        jpds(6)=100
      call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,ht(1,k),IRET1)
	if (IRET1 .ne. 0) write(6,*) ' AT LEVEL ', jpds(7) , jpds(5)
C######################### HT ###################

C######################### TP ###################
	jpds(5)=11
        jpds(6)=100
      call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,tp(1,k),IRET1)
	if (IRET1 .ne. 0) write(6,*) ' AT LEVEL ', jpds(7) , jpds(5)
C######################### TP ###################

!	jpds(5)=52
!      call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
!     &     BITMAP,rh(1,k),IRET1)
!	if (rh(20,k).gt.1.e20) write(6,*) 'getgbd(2) rh of ',k,rh(20,k)
!	if (IRET1 .ne. 0) write(6,*) ' AT LEVEL ', jpds(7) , jpds(5)

C######################### UW ###################
	jpds(5)=33
        jpds(6)=100
      call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,uw(1,k),IRET1)
	if (IRET1 .ne. 0) write(6,*) ' AT LEVEL ', jpds(7) , jpds(5)
C######################### UW ###################

C######################### VW ###################
	jpds(5)=34
        jpds(6)=100
      call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,vw(1,k),IRET1)
	if (IRET1 .ne. 0) write(6,*) ' AT LEVEL ', jpds(7) , jpds(5)
C######################### VW ###################
	
	write(6,*) 'Z,T,Q,U,V ', ht(nxny/2+55,k),tp(nxny/2+55,k),
     +	rh(nxny/2+55,k),uw(nxny/2+55,k),vw(nxny/2+55,k)

	enddo


c
c *** Normal finish.
c
1000  continue


C	rotate the winds at this point
	nx=KGDS(2)

	write(6,*) 'using datsav(10)= ', datsav(10)
	if (datsav(10) .eq. 3) gproj='LC'
	if (datsav(10) .eq. 5) gproj='PS'

	if (gproj .eq. 'LC') then
	write(6,*) 'rotating a lambert projection'
	call rotate_lcc(kgds,crot,srot)
	elseif (gproj .eq. 'PS') then
	call rotate_str(kgds,crot,srot,urlat,urlon)
	else
	write(6,*) 'not doing a wind rotation...'
	goto 1075
	endif

	do K=1,nz
	do I=1,nxny
	rmagb=(uw(I,K)**2. + vw(I,K)**2.)**(0.5)
	ubef=uw(I,K)
	vbef=vw(I,K)

	uw(I,K)=crot(I)*ubef+srot(I)*vbef
	vw(I,K)=crot(I)*vbef-srot(I)*ubef

	rmagaft=(uw(I,K)**2. + vw(I,K)**2.)**(0.5)
	if (abs(rmagaft-rmagb).gt.3.) then
	write(6,*) 'MAG, I,K,old,new==> ',I,K,rmagb,rmagaft
	write(6,*) 'original components..', ubef,vbef
	write(6,*) 'new components..', uw(I,k),vw(I,K)
	write(6,*) 'rotation cosines ', crot(I),srot(I)
	write(6,*) '.................................'
	endif
	ENDDO
	ENDDO

 1075	continue

      return
c
c *** Premature end of file.
c
1100  continue
      print *,'Premature end of file.'
      print *,'Abort...'
      stop
c
      end
c
c===============================================================================
c
      subroutine rh_fix(nx,ny,nz,rh,pr)
c
      implicit none
c
      integer nx,ny,nz,i,j,k,kk
c
      real rh(nx,ny,nz),pr(nz)
c_______________________________________________________________________________
c
c *** Fix bottom levels if necessary.
c
      if (rh(1,1,1) .eq. -99999.) then
	 do k=2,nz
	    if (rh(1,1,k) .ne. -99999.) then

	       DO kk=k-1,1,-1
	       print *,'Copying',nint(pr(kk+1)),' mb to'
     .                , nint(pr(kk)),' mb.'
	       do j=1,ny
	       do i=1,nx
		  rh(i,j,kk)=rh(i,j,kk+1)
               enddo
               enddo

               ENDDO

	       goto 10
            endif
         enddo
	 print *,'RH patch did not work.'
	 stop
      endif
c
c *** Fix upper levels if necessary.
c
10    continue
      if (rh(1,1,nz) .eq. -99999.) then
	 do k=nz-1,1,-1
	write(6,*) 'checking RH at lev= ', pr(k)
	    if (rh(1,1,k) .ne. -99999.) then
	       do kk=k+1,nz
	       print *,'Copying',nint(pr(kk-1)),' mb to'
     .                , nint(pr(kk)),' mb.'
	       do j=1,ny
	       do i=1,nx
		  rh(i,j,kk)=rh(i,j,kk-1)
               enddo
               enddo
               enddo
	       goto 20
            endif
         enddo      
      endif
c
20    continue
      do k=1,nz
	 if (rh(1,1,k) .eq. -99999.) then
	    print *,'RH patch did not work.'
	    stop
	 endif
      enddo
c
      return
      end
c
c===============================================================================
c
      subroutine es_ini
c
      common /estab/esat(15000:45000),es(15000:45000)
c
c *** Create tables of the saturation vapour pressure with up to
c        two decimal figures of accuraccy:
c
      do it=15000,45000
         t=it*0.01
         p1 = 11.344-0.0303998*t
         p2 = 3.49149-1302.8844/t
         c1 = 23.832241-5.02808*alog10(t)
         esat(it) = 10.**(c1-1.3816E-7*10.**p1+
     .               8.1328E-3*10.**p2-2949.076/t)
         es(it) = 610.78*exp(17.269*(t-273.16)/(t-35.86))
      enddo
c
      return
      end
c
c===============================================================================

	subroutine rotate_lcc(kgds,crot,srot)

Cmp	stolen/adapted from iplib code gdswiz03
C
C SUBPROGRAM:  GDSWIZ03   GDS WIZARD FOR LAMBERT CONFORMAL CONICAL
C   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
C
C ABSTRACT: THIS SUBPROGRAM DECODES THE GRIB GRID DESCRIPTION SECTION
C           (PASSED IN INTEGER FORM AS DECODED BY SUBPROGRAM W3FI63)
C           AND RETURNS ONE OF THE FOLLOWING:
C             (IOPT=+1) EARTH COORDINATES OF SELECTED GRID COORDINATES
C             (IOPT=-1) GRID COORDINATES OF SELECTED EARTH COORDINATES
C           FOR LAMBERT CONFORMAL CONICAL PROJECTIONS.
C           IF THE SELECTED COORDINATES ARE MORE THAN ONE GRIDPOINT
C           BEYOND THE THE EDGES OF THE GRID DOMAIN, THEN THE RELEVANT
C           OUTPUT ELEMENTS ARE SET TO FILL VALUES.
C           THE ACTUAL NUMBER OF VALID POINTS COMPUTED IS RETURNED TOO.

C       LAMBERT CONFORMAL GRIDS
C          (2)   - NX NR POINTS ALONG X-AXIS
C          (3)   - NY NR POINTS ALONG Y-AXIS
C          (4)   - LA1 LAT OF ORIGIN (LOWER LEFT)
C          (5)   - LO1 LON OF ORIGIN (LOWER LEFT)
C          (6)   - RESOLUTION (RIGHT ADJ COPY OF OCTET 17)
C          (7)   - LOV - ORIENTATION OF GRID
C          (8)   - DX - X-DIR INCREMENT
C          (9)   - DY - Y-DIR INCREMENT
C          (10)  - PROJECTION CENTER FLAG
C          (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
C          (12)  - LATIN 1 - FIRST LAT FROM POLE OF SECANT CONE INTER
C          (13)  - LATIN 2 - SECOND LAT FROM POLE OF SECANT CONE INTER


	parameter (NPTS=100000)
      INTEGER KGDS(200)
      REAL RLON(NPTS),RLAT(NPTS)
      REAL CROT(NPTS),SROT(NPTS)
	real DLON,AN
	real  DE,DR
	REAL PI,DPR
      PARAMETER(RERTH=6.3712E6)
      PARAMETER(PI=3.14159265,DPR=180./PI)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


	FILL=-999

	LROT=1
	IROT=1
        IM=KGDS(2)
        JM=KGDS(3)
        RLAT1=KGDS(4)*1.E-3
        RLON1=KGDS(5)*1.E-3
        IROT=MOD(KGDS(6)/8,2)
        ORIENT=KGDS(7)*1.E-3
        DX=KGDS(8)
        DY=KGDS(9)
        IPROJ=MOD(KGDS(10)/128,2)
        ISCAN=MOD(KGDS(11)/128,2)
        JSCAN=MOD(KGDS(11)/64,2)
        NSCAN=MOD(KGDS(11)/32,2)
        RLATI1=KGDS(12)*1.E-3
        RLATI2=KGDS(13)*1.E-3
        H=(-1.)**IPROJ
        HI=(-1.)**ISCAN
        HJ=(-1.)**(1-JSCAN)
        DXS=DX*HI
        DYS=DY*HJ

        IF(RLATI1.EQ.RLATI2) THEN
          AN=SIN(H*RLATI1/DPR)
        ELSE
          AN=LOG(COS(RLATI1/DPR)/COS(RLATI2/DPR))/
     &       LOG(TAN((H*RLATI1+90)/2/DPR)/TAN((H*RLATI2+90)/2/DPR))
        ENDIF


        DE=RERTH*COS(RLATI1/DPR)*TAN((H*RLATI1+90)/2/DPR)**AN/AN
        IF(H*RLAT1.EQ.90) THEN
          XP=1
          YP=1
        ELSE
          DR=DE/TAN((H*RLAT1+90)/2/DPR)**AN
          DLON1=MOD(RLON1-ORIENT+180+3600,360.)-180
C	atmp=RLON1-ORIENT+180.+3600.
C	DLON1= atmp - INT (atmp/360.)*360. - 180.
          XP=1-H*SIN(AN*DLON1/DPR)*DR/DXS
          YP=1+COS(AN*DLON1/DPR)*DR/DYS
        ENDIF
        ANTR=1/(2*AN)
        DE2=DE**2
        XMIN=0
        XMAX=IM+1
        YMIN=0
        YMAX=JM+1
        NRET=0
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TRANSLATE GRID COORDINATES TO EARTH COORDINATES
C	XP=1
C	YP=1
        DO N=1,IM*JM
	J=INT((N-1)/IM)+1
	I=N-(J-1)*IM
            IF(I.GE.XMIN.AND.I.LE.XMAX.AND.
     &         J.GE.YMIN.AND.J.LE.YMAX) THEN
              DI=(I-XP)*DXS
              DJ=(J-YP)*DYS
              DR2=DI**2+DJ**2
              IF(DR2.LT.DE2*1.E-6) THEN
                RLON(N)=0.
                RLAT(N)=H*90.
              ELSE
                RLON(N)=MOD(ORIENT+H/AN*DPR*ATAN2(DI,-DJ)+3600,360.)
C	atmp=ORIENT+H/AN*DPR*ATAN2(DI,-DJ)+3600
C	RLON(N)= atmp - INT(atmp/360.) * 360.
                RLAT(N)=H*(2*DPR*ATAN((DE2/DR2)**ANTR)-90)
              ENDIF
              NRET=NRET+1
              IF(LROT.EQ.1) THEN
                IF(IROT.EQ.1) THEN
C	atmp=RLON(N)-ORIENT+180.+3600.
C	DLON=atmp - INT(atmp/360.)*360.-180.
                  DLON=MOD(RLON(N)-ORIENT+180+3600,360.)-180
                  CROT(N)=H*COS(AN*DLON/DPR)
                  SROT(N)=SIN(AN*DLON/DPR)
                ELSE
                  CROT(N)=1
                 SROT(N)=0
                ENDIF
              ENDIF
            ELSE
              RLON(N)=FILL
              RLAT(N)=FILL
            ENDIF
          ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       subroutine get_gds(etafile,gdsinfo,kpds)

      character*(*) etafile
      character(LEN=1):: pds(50)

      integer kgds(200),kpds(200),len,kerr
     .         ,lenpds,lenkgds,nwords,kpdsl
     .         ,j,k,gdsinfo(11)
     .         ,gdsav,IRETO,JGDS(200),JPDS(200)
      real tmp(100000)
      logical bitmap(100000)


	nxny=100000


	JPDS=-1
	JGDS=-1


        len=index(etafile//' ',' ')-1

	call baopen(11,etafile(1:len),IRETO)
	write(6,*) 'BAOPEN in get_gds: ', IRETO

        if (IRETO .ne. 0) then
         print *,'Error opening unit=11, file name = ',etafile(1:len)
     .          ,' iostat = ',kerr
         stop
        endif

	jpds(5)=11
	jpds(6)=100
	jpds(7)=500
	call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,tmp,IRET1)

	write(6,*) 'back from getgb in get_gds: ', IRET1
	
       gdsinfo(1)=KGDS(2)
       gdsinfo(2)=KGDS(3)
       gdsinfo(3)=KGDS(4)
       gdsinfo(4)=KGDS(5)
       gdsinfo(5)=KGDS(7)
       gdsinfo(6)=KGDS(8)
       gdsinfo(7)=KGDS(9)
       gdsinfo(8)=KGDS(12)
       gdsinfo(9)=KGDS(13)
	gdsinfo(10)=KGDS(1)
	gdsinfo(11)=KPDS(3)
       write(6,*) gdsinfo

	write(6,*) 'KPDSinfo ', (kpds(I),i=1,10)

        return
        print *,'GETGDS PROBLEM'
        stop

        end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                subroutine lambert(gdslatur,gdslonur)

C
C       Subroutine written 9 March 1999 by M. Pyle to support tiled 221 input.
C       Code adapted from GEMPAK routine gblamb.c.  Whole purpose is to get the
C       UR corner lat and lon for use in workstation Eta.

        integer latin1,latin2,nx,ny,la1,lo1,lov
        integer dx,dy,datsav

        real(kind=8):: earth_rad, const_cone, xll,yll,xur,yur,lat1,
     +  lon1,loncnt,angle1,angle2,x1,x2,y1,y2,alpha,rtemp

        real gdslatur,gdslonur
     +  ,gdslatll,gdslonll

        parameter(rpi=3.141592654)
        parameter(d2r=rpi/180.)
        parameter(r2d=180./rpi)
        parameter(radius=6370000.)

        common /gdsinfo/ datsav(11)

        latin1=datsav(8)
        latin2=datsav(9)
        la1=datsav(3)
        lo1=datsav(4)
        dx=datsav(6)
        dy=datsav(7)
        nx=datsav(1)
        ny=datsav(2)
        lov=datsav(5)

        write(6,*) 'values in lambert '
        write(6,*) latin1,latin2,la1,lo1,dx,dy,nx,ny,lov

        lat1= (la1/1000.0)*d2r
        if (lo1 .eq.  180000) lo1=lo1-360000
        if (lo1 .lt. -180000) lo1=lo1+360000
        lon1= (lo1/1000.0)*d2r

Cmp     now have LL corner in radians, W is negative

        if (lov .eq.  180000) lov=lov-360000
        if (lov .lt. -180000) lov=lov+360000
        loncnt= (lov/1000.0)*d2r

        angle1= (rpi/2.) - ( abs(latin1/1000.0) * d2r )
        angle2= (rpi/2.) - ( abs(latin2/1000.0) * d2r )

        if (latin1 .eq. latin2) then
        const_cone=cos(angle1)
        else
        const_cone= ( log ( sin(angle2) ) - log ( sin ( angle1 ) ) )/
     +  ( log ( tan ( angle2/2 ) ) - log ( tan ( angle1/2 ) ) )
        endif

        write(6,*) 'const_cone= ', const_cone

        earth_rad=radius/const_cone

cmp     assuming NH

        x1 = earth_rad * tan( (rpi/2.-lat1) / 2 )**(const_cone)*
     +  sin (const_cone * ( lon1 - loncnt ) )
        y1 = -earth_rad * tan( (rpi/2.-lat1) / 2 )**(const_cone)*
     +  cos (const_cone * ( lon1 - loncnt ) )

         alpha= (tan(angle1 / 2 )**const_cone)/sin (angle1)

        x2=x1 + ( nx - 1 ) * alpha * dx
        y2=y1 + ( ny - 1 ) * alpha * dy

        xll=min(x1,x2)
        xur=max(x1,x2)
        yll=min(y1,y2)
        yur=max(y1,y2)

        gdslatll= ( rpi/2. - 2 *
     + atan ( ( (abs(xll)**2.+abs(yll)**2.)**(0.5)/earth_rad)**
     + (1/const_cone) ) ) * r2d

	write(6,*) 'xll, yll, xll**2., yll**2. ', xll, yll, xll**2., 
     +	yll**2.
	write(6,*) 'lat pieces: ', rpi/2., (xll**2.+yll**2.), earth_rad,
     +	r2d


        rtemp= atan2 ( xll, -yll ) * ( 1 / const_cone ) + loncnt

        if ( rtemp .gt. rpi ) then
        gdslonll = ( rtemp - 2.*rpi ) * r2d
        else if ( rtemp .lt. -rpi ) then
        gdslonll = ( rtemp + 2.*rpi ) * r2d
        else
        gdslonll = rtemp * r2d
        endif


        gdslatur= ( rpi/2. - 2 *
     + atan ( ( (abs(xur)**2.+abs(yur)**2.)**(0.5)/earth_rad)**
     + (1/const_cone) ) ) * r2d

        rtemp= atan2 ( xur, -yur ) * ( 1 / const_cone ) + loncnt
        if ( rtemp .gt. rpi ) then
        gdslonur = ( rtemp - 2.*rpi ) * r2d
        else if ( rtemp .lt. -rpi ) then
        gdslonur = ( rtemp + 2.*rpi ) * r2d
        else
        gdslonur = rtemp * r2d
        endif

        write(6,*) 'output==> '
        write(6,*) 'LL points ', gdslatll,gdslonll
        write(6,*) 'UR points ', gdslatur,gdslonur

        return

        END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        subroutine rotate_str(kgds,crot,srot,urlat,urlon)

C
C Program adapted 26 March 99 for workstation Eta.  Original code from
C iplib program GDSWZD05
C
C

C SUBPROGRAM:  GDSWZD05   GDS WIZARD FOR POLAR STEREOGRAPHIC AZIMUTHAL
C   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
C
C ABSTRACT: THIS SUBPROGRAM DECODES THE GRIB GRID DESCRIPTION SECTION
C           (PASSED IN INTEGER FORM AS DECODED BY SUBPROGRAM W3FI63)
C           AND RETURNS ONE OF THE FOLLOWING:
C             (IOPT=+1) EARTH COORDINATES OF SELECTED GRID COORDINATES
C             (IOPT=-1) GRID COORDINATES OF SELECTED EARTH COORDINATES
C           FOR POLAR STEREOGRAPHIC AZIMUTHAL PROJECTIONS.
C           IF THE SELECTED COORDINATES ARE MORE THAN ONE GRIDPOINT
C           BEYOND THE THE EDGES OF THE GRID DOMAIN, THEN THE RELEVANT
C           OUTPUT ELEMENTS ARE SET TO FILL VALUES.
C           THE ACTUAL NUMBER OF VALID POINTS COMPUTED IS RETURNED TOO.
C           OPTIONALLY, THE VECTOR ROTATIONS AND THE MAP JACOBIANS
C           FOR THIS GRID MAY BE RETURNED AS WELL.

C       POLAR STEREOGRAPHIC GRIDS
C          (2)   - N(I) NR POINTS ALONG LAT CIRCLE
C          (3)   - N(J) NR POINTS ALONG LON CIRCLE
C          (4)   - LA(1) LATITUDE OF ORIGIN
C          (5)   - LO(1) LONGITUDE OF ORIGIN
C          (6)   - RESOLUTION FLAG  (RIGHT ADJ COPY OF OCTET 17)
C          (7)   - LOV GRID ORIENTATION
C          (8)   - DX - X DIRECTION INCREMENT
C          (9)   - DY - Y DIRECTION INCREMENT
C          (10)  - PROJECTION CENTER FLAG
C          (11)  - SCANNING MODE (RIGHT ADJ COPY OF OCTET 28)

        parameter (NPTS=80000)
      INTEGER KGDS(200)
      REAL RLON(NPTS),RLAT(NPTS)
      REAL CROT(NPTS),SROT(NPTS)
        real DLON,AN
        real  DE,DR,urlat,urlon
        REAL PI,DPR
      PARAMETER(RERTH=6.3712E6)
      PARAMETER(PI=3.14159265,DPR=180./PI)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        LROT=1
        IROT=1
        IM=KGDS(2)
        JM=KGDS(3)
        RLAT1=KGDS(4)*1.E-3
        RLON1=KGDS(5)*1.E-3
        IROT=MOD(KGDS(6)/8,2)
        ORIENT=KGDS(7)*1.E-3
        DX=KGDS(8)
        DY=KGDS(9)
        IPROJ=MOD(KGDS(10)/128,2)
        ISCAN=MOD(KGDS(11)/128,2)
        JSCAN=MOD(KGDS(11)/64,2)
        NSCAN=MOD(KGDS(11)/32,2)
        H=(-1.)**IPROJ
        HI=(-1.)**ISCAN
        HJ=(-1.)**(1-JSCAN)
        DXS=DX*HI
        DYS=DY*HJ

        DE=(1.+SIN(60./DPR))*RERTH
        DR=DE*COS(RLAT1/DPR)/(1+H*SIN(RLAT1/DPR))
        XP=1-H*SIN((RLON1-ORIENT)/DPR)*DR/DXS
        YP=1+COS((RLON1-ORIENT)/DPR)*DR/DYS

        DE2=DE**2
        XMIN=0
        XMAX=IM+1
        YMIN=0
        YMAX=JM+1
        NRET=0
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TRANSLATE GRID COORDINATES TO EARTH COORDINATES

        DO N=1,IM*JM
        J=INT((N-1)/IM)+1
        I=N-(J-1)*IM
            IF(I.GE.XMIN.AND.I.LE.XMAX.AND.
     &         J.GE.YMIN.AND.J.LE.YMAX) THEN
              DI=(I-XP)*DXS
              DJ=(J-YP)*DYS
              DR2=DI**2+DJ**2

              IF(DR2.LT.DE2*1.E-6) THEN
                RLON(N)=0.
                RLAT(N)=H*90.
              ELSE
                RLON(N)=MOD(ORIENT+H*DPR*ATAN2(DI,-DJ)+3600,360.)
                RLAT(N)=H*DPR*ASIN((DE2-DR2)/(DE2+DR2))
              ENDIF

              NRET=NRET+1

              IF(LROT.EQ.1) THEN
                IF(IROT.EQ.1) THEN
                 CROT(N)=H*COS((RLON(N)-ORIENT)/DPR)
                 SROT(N)=SIN((RLON(N)-ORIENT)/DPR)
                ELSE
                 CROT(N)=1
                 SROT(N)=0
                ENDIF
              ENDIF
	ENDIF

          ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	urlat=rlat(im*jm)
	urlon=rlon(im*jm)
	if (urlon .gt. 180) urlon=urlon-360.
	if (urlon .lt. -180) urlon=urlon+360.
	write(6,*) 'in stereo sub found lat,lon ', urlat,urlon
C	write(6,*) 'lat/lon in stereo '
	do N=1,IM*JM
        J=INT((N-1)/IM)+1
        I=N-(J-1)*IM
C	write(6,69) I,J,RLAT(N),RLON(N)
	enddo
   69   format(I3,1x,I3,1x,2(f7.3,1x))
        RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	        subroutine get_fullgds(etafile,nx,ny,kgds)

        character*(*) etafile
        character(LEN=1):: pds(50)

        integer kgds(200),kpds(200),len,jpds(200)
     .         ,lenpds,lenkgds,nwords,kpdsl,jgds(200)
     .         ,j,k,KNUM,nx,ny

	logical bitmap(nx*ny)
	real tmp(nx*ny)

	write(6,*) 'inside get_fullgds...', etafile


	nxny=nx*ny

        len=index(etafile//' ',' ')-1

	jpds=-1
	jgds=-1

	jpds(5)=11
	jpds(6)=100
	jpds(7)=500

	write(6,*) 'calling getgb '
C	write(6,*) 'jpds =  ', jpds
C	write(6,*) 'jgds =  ', jgds

	call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,tmp,IRET1)
     
        KGDS(4)=KGDS(4)*1.E+2
        KGDS(5)=KGDS(5)*1.E+2
        KGDS(7)=KGDS(7)*1.E+2
        KGDS(8)=KGDS(8)*1.E+2

	write(6,*) 'return from getgb ', IRET1

        if (IRET1 .ne. 0) then
         print *,'Error  getting GDS in get_fullgds ', IRET1
         stop
        endif

	write(6,*) 'leaving get_fullgds'
        return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine smooth(a,idim,jdim,ldim,npass)
	integer idim,jdim,ldim,npass
	real a(idim,jdim,ldim)

	DO N=1,NPASS


	do L=1,LDIM
	 do J=2,JDIM-1
	  do I=2,IDIM-1

	  a(i,j,l)=1./8.*(4.*a(i,j,l)+a(i+1,j,l)+a(i-1,j,l)
     +				     +a(i,j+1,l)+a(i,j-1,l))
	  enddo
         enddo
	enddo

	ENDDO
	
	return
	end
