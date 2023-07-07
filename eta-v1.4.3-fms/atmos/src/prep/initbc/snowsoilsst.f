
	SUBROUTINE SNOWSOILSST(im,jm,tlm0d,tph0d,dlmd,dphd,sfc_file,
     +	GGSTC1,GGSTC2,GGSMC1,GGSMC2,SSTRAW,SNOWRAW,GGLAND,GGICE,kgds)


	parameter (nx=192,ny=94)
	parameter(dtr=3.14159/180.,r2d=1./dtr)

CC	read and fill these 8 arrays...process rest in another sub

C	real,allocatable:: GGSTC1(:,:),GGSTC2(:,:)
C	real,allocatable:: GGSMC1(:,:),GGSMC2(:,:)
C	real,allocatable:: SSTRAW(:,:),SNOWRAW(:,:)
C	real,allocatable:: GGLAND(:,:),GGICE(:,:)

	REAL GGSTC1(nx,ny),GGSTC2(nx,ny),GGSMC1(nx,ny)
	REAL GGSMC2(nx,ny),SSTRAW(nx,ny),SNOWRAW(nx,ny)
	REAL GGLAND(nx,ny),GGICE(nx,ny)
	real,allocatable:: tmp(:)

	character*256 sfc_file

        integer ICOUNT,I,ixp,jyp

	      integer kgds(200),kpds(200),len
     .         ,nwords,IRETO,IRET1,JPDS(200),JGDS(200),KNUM

        logical BITMAP(nx*ny)

Cmp	eventually need a proper file name here
C
	len=index(sfc_file,' ')-1
C	write(6,*) 'trying to open ', sfc_file(1:len)

	call baopen(11,sfc_file(1:len),IRETO)
	if (IRETO .ne. 0) write(6,*) 'baopen error!'

        jpds=-1
        jgds=-1

        jpds(7)=10
        jpds(5)=11
        jpds(6)=112


CC	first call to get GDS of Gaussian grid
CC
	allocate(tmp(nx*ny))
        call getgb(11,0,nx*ny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,tmp,IRET1)

	if (IRET1 .ne. 0) write(6,*) 'TROUBLE AHEAD!!!'
!	write(6,*) 'first getgb, nwords: ' , nwords

!	write(6,*) 'raw vals for STC'
	do J=1,NX*NY,500
!	write(6,*) 'J, tmp(J): ', J, tmp(J)
	enddo

        jpds=-1
        jgds=-1

        jpds(7)=10
        jpds(5)=11
        jpds(6)=112

!	write(6,*) 'nx*ny= ', nx*ny

CC	0-10 STC
CC
!	write(6,*) 'calling getgb'

        call getgb(11,0,nx*ny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,tmp,IRET1)

C	allocate(GGSTC1(nx,ny))

	if (IRET1 .eq. 0) then
	call reord(nx,ny,tmp,GGSTC1)
	else
	write(6,*) 'IRET1= ', IRET1
	write(6,*) 'no GGSTC1'
	STOP 666
	endif

CC	0-10 SMC
CC
	jpds(5)=144

!	write(6,*) 'getgb for 0-10 SMC'
        call getgb(11,0,nx*ny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,tmp,IRET1)

	if (IRET1 .eq. 0) then
C	allocate (GGSMC1(nx,ny))
	call reord(nx,ny,tmp,GGSMC1)
!	write(6,*) 'past reord 0-10 SMC'
	else
	write(6,*) 'no GGSMC1'
	write(6,*) 'IRET1= ', IRET1
	STOP 666
	endif

CC	10-200 STC
CC
	jpds(5)=11
	jpds(7)=2760
!	write(6,*) 'calling getgb for 10-200 STC'
        call getgb(11,0,nx*ny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,tmp,IRET1)

	if (IRET1 .eq. 0) then
C	allocate (GGSTC2(nx,ny))
	call reord(nx,ny,tmp,GGSTC2)
	else
	write(6,*) 'no GGSTC2'
	endif

CC	10-200 SMC
CC

	jpds(5)=144
!	write(6,*) 'getgb for 10-200 SMC'
        call getgb(11,0,nx*ny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,tmp,IRET1)
	
	if (IRET1 .eq. 0) then
C	allocate (GGSMC2(nx,ny))
	call reord(nx,ny,tmp,GGSMC2)
	else
	write(6,*) 'no GGSMC2'
	endif

CC
CC 	SFC TEMP (for SST?)
CC	

	jpds=-1
	jpds(5)=11
	jpds(6)=1
	jpds(7)=0
!	write(6,*) 'getgb for  SFC TMP'
        call getgb(11,0,nx*ny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,tmp,IRET1)

	if (IRET1 .eq. 0) then
C	allocate (SSTRAW(nx,ny))
	call reord(nx,ny,tmp,SSTRAW)
	else
	write(6,*) 'no SSTRAW'
	endif
	

CC
CC	snow depth?
CC

	jpds(5)=65
	jpds(6)=1
	jpds(7)=0
!	write(6,*) 'getgb for SNOW DEPTH'
        call getgb(11,0,nx*ny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,tmp,IRET1)

	if (IRET1 .eq. 0) then
C	allocate (SNOWRAW(nx,ny))
	call reord(nx,ny,tmp,SNOWRAW)
	else
	write(6,*) 'no snowraw'
	endif

	write(6,*) 'snowraw after reordered'
	write(6,*) 'I=140-->165, J=80-->60'
!	do J=1,ny/2
	do J=13,35
	write(6,632) (snowraw(I,J),I=142,160)
	enddo
  632	format(50(f4.0))

	jpds(5)=81
!	write(6,*) 'getgb for land mask'
        call getgb(11,0,nx*ny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,tmp,IRET1)

!	write(6,*) 'NWORDS from ggland getgb: ', nwords

!	write(6,*) 'raw vals for GGLANd'
	do J=1,NX*NY,500
!	write(6,*) 'J, tmp(J): ', J, tmp(J)
	enddo
	
	if (IRET1 .eq. 0) then
C	allocate (GGLAND(nx,ny))
	call reord(nx,ny,tmp,GGLAND)
	else
	write(6,*) 'ggland problems ', IRET1
	endif

        jpds(5)=91
	write(6,*) 'getgb for ice mask'
        call getgb(11,0,nx*ny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS,
     &     BITMAP,tmp,IRET1)

	if (IRET1 .eq. 0) then
C	allocate (GGICE(nx,ny))
	call reord(nx,ny,tmp,GGICE)
	write(6,*) 'GGICE after reordered'
	do J=13,35
	write(6,632) (ggice(I,J),I=142,160)
	enddo
!	write(6,*) 'end GGICE'
	endif
CC
	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


	subroutine reord(nx,ny,input,output)

	real input(nx*ny),output(nx,ny)
C	write(6,*) 'inside reord ', nx,ny

	output=-9999.

	ICOUNT=0
	do J=1,ny
	do I=1,nx
	ICOUNT=ICOUNT+1
	output(I,J)=input(ICOUNT)
	enddo
	enddo

C	write(6,*) 'leaving reord'

	return
	end
