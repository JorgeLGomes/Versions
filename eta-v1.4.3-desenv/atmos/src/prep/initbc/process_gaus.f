
	SUBROUTINE process_gaus(im,jm,tlm0d,tph0d,dlmd,dphd,sm,stc,smc,
     +	sst,sice,si,GGSTC1,GGSTC2,GGSMC1,GGSMC2,SSTRAW,SNOWRAW,GGLAND,
     +	GGICE,KGDS)


	parameter (nx=192,ny=94)
	parameter (nsoil=4)
	parameter(dtr=3.14159/180.,r2d=1./dtr)


CC	read and fill these 8 arrays...process rest in another sub

        real GGSTC1(nx,ny),GGSTC2(nx,ny)
	real GGSMC1(nx,ny),GGSMC2(nx,ny)
	real SSTRAW(nx,ny),SNOWRAW(nx,ny)
	real GGLAND(nx,ny),GGICE(nx,ny)

	real, allocatable:: tmp(:)

	real sm(im,jm),si(im,jm),sice(im,jm)
	real stc(im,jm,nsoil),smc(im,jm,nsoil),sst(im,jm)

	real rlat,rlon,xpts(im,jm),ypts(im,jm)
C	real glat(im,jm),glon(im,jm),rout(im,jm),sfcgrid(im,jm,6)
	real glat(im,jm),glon(im,jm)

	character*256 sfc_file

C        logical*1 ltmp(nx,ny)
        integer ICOUNT,I,ixp,jyp

	      integer kgds(200),kpds(200),len,kerr
     .         ,nwords
     .         ,IRETO,IRET1,JPDS(200),JGDS(200),KNUM

        logical*1 BITMAP(nx*ny)


C	Define eta grid lat/lon
C
        WBD=-(float(IM)-1.)*DLMD
        SBD=(-(float(JM)-1.)/2.)*DPHD

	write(6,*) 'wbd, sbd= ', wbd,sbd

	 tph0=tph0d*dtr
      wb=wbd*dtr
      sb=sbd*dtr
      dlm=dlmd*dtr
      dph=dphd*dtr
      tdlm=dlm+dlm
      tdph=dph+dph
      tph=sb-dph
      ctph0=cos(tph0)
       stph0=sin(tph0)

      do j=1,jm
         tlm=wb-tdlm+mod(j+1,2)*dlm
         tph=tph+dph
         stph=sin(tph)
         ctph=cos(tph)
c
        do i=1,im
            tlm=tlm+tdlm
            sinphi=ctph0*stph+stph0*ctph*cos(tlm)
            glat(i,j)=asin(sinphi)
            coslam=ctph*cos(tlm)/(cos(glat(i,j))*ctph0)
     .            -tan(glat(i,j))*tan(tph0)
            coslam=min(coslam,1.)
            fact=1.
            if (tlm .gt. 0.0) fact=-1.
            glon(i,j)=-tlm0d*dtr+fact*acos(coslam)
        glat(i,j)=glat(i,j)*r2d
        glon(i,j)=glon(i,j)*r2d
C
        if (glon(i,j).lt.0) glon(i,j)=glon(i,j)+360.

	glon(i,j)=360.-glon(i,j)

         enddo
      enddo

CCC
CC	compute gaussian points corresponding to output e-grid	
CC
C
	call gauss_ijonce(kgds,-1,im*jm,-9999.,XPTS,YPTS,
     +			glon,glat,numret)
	write(6,*) 'numret from gauss_ijonce: ', numret


CC
CC	At this point have all of the pieces.  Begin by creating a SICE mask
CC

C 20010806  WHY IS SNOWRAW USED RATHER THAN GGICE TO GENERATE SICE????

	do J=1,NY
	do I=1,NX
	if (SNOWRAW(I,J) .gt. 0 .and. GGLAND(I,J) .eq. 0 .and. 
     +		GGICE(I,J) .eq. 0) then
	write(6,*) 'snow & water...possible ice ', I,J
	endif
	enddo
	enddo

	call reansnowint(SNOWRAW,nx,ny,im,jm,glat,glon,xpts,ypts,sm,si)
	call reaniceint(GGICE,nx,ny,im,jm,glat,glon,xpts,ypts,sm,sice)

        write(6,*) 'si from reaniceint'
        do J=jm,1,-JM/50
        write(6,372) (si(I,J),I=1,im,2)
        enddo

        write(6,*) 'sea ice from reaniceint'
        do J=jm,1,-JM/50
        write(6,372) (sice(I,J),I=1,im,2)
  372   format(50(f2.0))
        enddo

	DO J=1,JM
	DO I=1,IM
	 if (sm(I,J) .gt. 0.9 .and. SICE(I,J) .gt. 0.5) then
	  SI(I,J)=1.
	elseif (sm(I,J) .lt. 0.5) then

C	divide by 1000 to put in m, multiply by 5 to make snow depth
C	instead of liquid equiv
C
	si(I,J)=si(I,J)*(5./1000)

C
C	try to avoid bogus increase of ice mask!
C
	elseif (sm(I,J) .gt. 0.5 .and. SI(I,J) .gt. 0. .and. 
     +		SICE(I,J) .eq. 0) then
	write(6,*) 'water & snow ', i,j,si(I,J)

Ctoolittle	if (si(I,J) .eq. 25) then

	if (si(I,J) .ge. 25) then
	  SICE(I,J)=1.
	  SI(I,J)=1.
	else
	  SI(I,J)=0.
	endif

	endif

	ENDDO	
	ENDDO

	call reansoilint(GGSMC1,nx,ny,im,jm,glat,glon,xpts,ypts,
     +				ggland,sice,sm,smc(1,1,1))
	call reansoilint(GGSMC2,nx,ny,im,jm,glat,glon,xpts,ypts,
     +				ggland,sice,sm,smc(1,1,2))
	call reansoilint(GGSTC1,nx,ny,im,jm,glat,glon,xpts,ypts,
     +				ggland,sice,sm,stc(1,1,1))
	call reansoilint(GGSTC2,nx,ny,im,jm,glat,glon,xpts,ypts,
     +				ggland,sice,sm,stc(1,1,2))
!	call reansoilice(GGSTC1,nx,ny,im,jm,glat,glon,xpts,ypts,
!     +				ggice,sice,sm,stc(1,1,1))
!	call reansoilice(GGSTC2,nx,ny,im,jm,glat,glon,xpts,ypts,
!     +			ggice,sice,sm,stc(1,1,2))


	do N=1,nsoil
	do J=1,JM
	do I=1,IM
	if (STC(I,J,N) .eq. 0) stc(I,J,N)=273.15
!	if (SMC(I,J,N) .eq. 0) smc(I,J,N)=273.15
	enddo
	enddo
	enddo

	do KK=3,4
	 do J=1,jm
	  do I=1,im
		STC(I,J,KK)=STC(I,J,2)
		SMC(I,J,KK)=SMC(I,J,2)
	  enddo
	 enddo
	enddo

!	write(6,*) 'what is ggland?? ' , nx,ny
	do J=1,ny/2,2
!	write(6,611) (ggland(I,J),I=1,NX,8)
	enddo
  611	format(40f4.0)

	call reanseaint(SSTRAW,GGLAND,sm,sice,nx,ny,im,jm,glat,glon,
     +			xpts,ypts,sst)


	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine reansoilint(ginput,nx,ny,im,jm,glat,glon,xpts,ypts,
     +					rtest,sice,sm,output)

	real ginput(nx,ny),tmpo(nx,ny),glat(im,jm),glon(im,jm)
	real xpts(im,jm),ypts(im,jm),output(im,jm),sm(im,jm),
     +			sice(im,jm),rtest(nx,ny)
     


	tmpo=ginput

	output=-99.

  64	format(10(f5.1,1x))
  75	format(10(f6.1,1x))

	do J=1,jm
	do I=1,im

        if((sm(i,j).gt.0.5).or.(sice(i,j).eq.1.0)) then
c                           set non land-mass values
          output(i,j) = 0.0
        else
c                           set land-mass values
	y=ypts(i,j)
	x=xpts(i,j)
	  iy = y
          ix = x
          iyp1 = iy + 1
          ixp1 = ix + 1
c  First calculate the bilinear weights assuming a constant dx and dy
          wxy =     (ixp1-x) * (iyp1 - y)
          wxp1y =   (x-ix)   * (iyp1 - y)
          wxp1yp1 = (x-ix)   * (y-iy)
          wxyp1 = (ixp1 - x) * (y-iy)
c  Take care of the wrap around the 0 meridan when necessary
          if(ix.eq.nx) ixp1 = 1
c  Only use land mass points (ggsli=1) for the interpolation
C
C	here rtest is land
C
          if(rtest(ix,iy).eq.0.0) wxy = 0.0
          if(rtest(ixp1,iy).eq.0.0) wxp1y = 0.0
          if(rtest(ixp1,iyp1).eq.0.0) wxp1yp1 = 0.0
          if(rtest(ix,iyp1).eq.0.0) wxyp1 = 0.0
          sum = wxy + wxp1y + wxp1yp1 + wxyp1
          if(sum.eq.0.0) then
c  The eta land point is surrounded by gaussian grid non land
c  Look for nearest gaussian grid land point and assign it to the
c  Eta grid point. Find nearest at this latitude
            do jfnd = iy,nx
            do ifnd = ix,ny
            if(tmpo(ifnd,jfnd).gt.0.0) then
                wxy = 1.0
                sum = 1.0
                ix = ifnd
                iy = jfnd
                go to 3356
            end if
            end do
            do ifnd = ix,1,-1
            if(tmpo(ifnd,jfnd).gt.0.0) then
                wxy = 1.0
                sum = 1.0
                ix = ifnd
                iy = jfnd
                go to 3356
            end if
            end do
            end do
!       print *,"Oh Oh can't find a land point at:",glat(i,j),glon(i,j)
!       print *,((tmpo(ii,jj),ii=ix,ix+3),jj=iy,iy+2)
3356        continue
          end if
          wxy = wxy / sum
          wxp1y = wxp1y / sum
          wxp1yp1 = wxp1yp1 / sum
          wxyp1 = wxyp1 / sum
          output(i,j) = wxy     * tmpo(ix,iy)     +
     &                 wxp1y   * tmpo(ixp1,iy)   +
     &                 wxp1yp1 * tmpo(ixp1,iyp1) +
     &                 wxyp1   * tmpo(ix,iyp1)
	if (output(i,j) .eq. 0) 
     +		write(6,*) 'BAD OUTPUT in REANSOILINT', I,J
        end if




	enddo
	enddo

   65	format(30f6.1)

	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine reanseaint(tmpo,rland,sm,sice,nx,ny,im,jm,glat,glon,
     +				xpts,ypts,output)

	real ginput(nx*ny),tmpo(nx,ny),glat(im,jm),glon(im,jm)
	real xpts(im,jm),ypts(im,jm),output(im,jm),rland(nx,ny)
	real sm(im,jm),sice(im,jm)
     

	output=-99.


  64	format(10(f5.1,1x))
  75	format(10(f6.1,1x))

	do J=1,jm
	do I=1,im

C
        if((sm(i,j).lt.0.5).or.(sice(i,j).eq.1.0)) then
c                           set non land-mass values
Cmp          output(i,j) = 0.0
          output(i,j) = 273.15
        else
c                           set land-mass values
	y=ypts(i,j)
	x=xpts(i,j)
	  iy = y
          ix = x
          iyp1 = iy + 1
          ixp1 = ix + 1
c  First calculate the bilinear weights assuming a constant dx and dy
          wxy =     (ixp1-x) * (iyp1 - y)
          wxp1y =   (x-ix)   * (iyp1 - y)
          wxp1yp1 = (x-ix)   * (y-iy)
          wxyp1 = (ixp1 - x) * (y-iy)
c  Take care of the wrap around the 0 meridan when necessary
          if(ix.eq.nx) ixp1 = 1
c  Only use water points (rland=0) for the interpolation
          if(rland(ix,iy).eq.1.0) wxy = 0.0
          if(rland(ixp1,iy).eq.1.0) wxp1y = 0.0
          if(rland(ixp1,iyp1).eq.1.0) wxp1yp1 = 0.0
          if(rland(ix,iyp1).eq.1.0) wxyp1 = 0.0
          sum = wxy + wxp1y + wxp1yp1 + wxyp1
          if(sum.eq.0.0) then
!	write(6,*) 'dont want to be here'
c  The eta land point is surrounded by gaussian grid non land
c  Look for nearest gaussian grid land point and assign it to the
c  Eta grid point. Find nearest at this latitude
            do jfnd = iy,nx
            do ifnd = ix,ny
C            if(tmpo(ifnd,jfnd).gt.0.0) then
           if(rland(ifnd,jfnd).eq.0.0) then
                wxy = 1.0
                sum = 1.0
                ix = ifnd
                iy = jfnd
                go to 3356
            end if
            end do
            do ifnd = ix,1,-1
            if(rland(ifnd,jfnd).eq.0.0) then
                wxy = 1.0
                sum = 1.0
                ix = ifnd
                iy = jfnd
                go to 3356
            end if
            end do
            end do
       print *,"Oh Oh can't find a land point at:",glat(i,j),glon(i,j)
       print *,((tmpo(ii,jj),ii=ix,ix+3),jj=iy,iy+2)
3356        continue
          end if
          wxy = wxy / sum
          wxp1y = wxp1y / sum
          wxp1yp1 = wxp1yp1 / sum
          wxyp1 = wxyp1 / sum
          output(i,j) = wxy     * tmpo(ix,iy)     +
     &                 wxp1y   * tmpo(ixp1,iy)   +
     &                 wxp1yp1 * tmpo(ixp1,iyp1) +
     &                 wxyp1   * tmpo(ix,iyp1)
        end if

	enddo
	enddo

	write(6,*) 'OUTPUT SST'

	do J=jm,1,-jm/30
	write(6,65) (output(i,j),i=1,im,im/15)
	enddo
   65	format(30f6.1)


	return

	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine reansnowint(tmpo,nx,ny,im,jm,glat,glon,xpts,ypts,
     +				        sm,output)

	real ginput(nx*ny),tmpo(nx,ny),glat(im,jm),glon(im,jm)
	real xpts(im,jm),ypts(im,jm),output(im,jm),sm(im,jm)

	output=-99.

!	write(6,*) 'input snow field: '

	do J=1,ny/2
!	write(6,632) (tmpo(I,J),I=nx/2,nx,5)
	enddo
  632	format(50(f4.0,1x))

! 	write(6,*) 'end input snow:'

  64	format(10(f5.1,1x))
  75	format(10(f6.1,1x))

	do J=1,jm
	do I=1,im

C
C	interpolate blindly with regard to land/sea?
C
C        if((sm(i,j).gt.0.5).or.(sice(i,j).eq.1.0)) then
C        if((sm(i,j).gt.0.5)) then
c                           set non land-mass values
C          output(i,j) = 0.0
C        else
c                           set land-mass values
	y=ypts(i,j)
	x=xpts(i,j)
	  iy = y
          ix = x
          iyp1 = iy + 1
          ixp1 = ix + 1
c  First calculate the bilinear weights assuming a constant dx and dy
          wxy =     (ixp1-x) * (iyp1 - y)
          wxp1y =   (x-ix)   * (iyp1 - y)
          wxp1yp1 = (x-ix)   * (y-iy)
          wxyp1 = (ixp1 - x) * (y-iy)
c  Take care of the wrap around the 0 meridan when necessary
          if(ix.eq.nx) ixp1 = 1
c  Only use land mass points (ggsli=1) for the interpolation
C          if(tmpo(ix,iy).eq.0.0) wxy = 0.0
C          if(tmpo(ixp1,iy).eq.0.0) wxp1y = 0.0
C          if(tmpo(ixp1,iyp1).eq.0.0) wxp1yp1 = 0.0
C          if(tmpo(ix,iyp1).eq.0.0) wxyp1 = 0.0
          sum = wxy + wxp1y + wxp1yp1 + wxyp1
          if(sum.eq.0.0) then
!	write(6,*) 'dont want to be here'
c  The eta land point is surrounded by gaussian grid non land
c  Look for nearest gaussian grid land point and assign it to the
c  Eta grid point. Find nearest at this latitude
            do jfnd = iy,nx
            do ifnd = ix,ny
            if(tmpo(ifnd,jfnd).gt.0.0) then
                wxy = 1.0
                sum = 1.0
                ix = ifnd
                iy = jfnd
                go to 3356
            end if
            end do
            do ifnd = ix,1,-1
            if(tmpo(ifnd,jfnd).gt.0.0) then
                wxy = 1.0
                sum = 1.0
                ix = ifnd
                iy = jfnd
                go to 3356
            end if
            end do
            end do
       print *,"Oh Oh can't find a land point at:",glat(i,j),glon(i,j)
       print *,((tmpo(ii,jj),ii=ix,ix+3),jj=iy,iy+2)
3356        continue
          end if
          wxy = wxy / sum
          wxp1y = wxp1y / sum
          wxp1yp1 = wxp1yp1 / sum
          wxyp1 = wxyp1 / sum
          output(i,j) = wxy     * tmpo(ix,iy)     +
     &                 wxp1y   * tmpo(ixp1,iy)   +
     &                 wxp1yp1 * tmpo(ixp1,iyp1) +
     &                 wxyp1   * tmpo(ix,iyp1)
Cmp
!	if (output(i,j) .gt. 0 .and. sm(i,j) .eq. 0) then
C
C
!	if (tmpo(ix,iy) .ne. 100. .and. tmpo(ixp1,iy) .ne. 100 .and.
!     +	tmpo(ixp1,iyp1) .ne. 100 .and. tmpo(ix,iyp1) .ne. 100 
!     +	.AND. tmpo(ix,iy) .ne. 25. .and. tmpo(ixp1,iy) .ne. 25. .and.
!     +	tmpo(ixp1,iyp1) .ne. 25. .and. tmpo(ix,iyp1) .ne. 25.) then
!	write(6,*) 'would have created ice at ', I,J
!	write(6,*) 'leaving isolated sea ice at ', I,J
!	write(6,*) 'largest input: ', amax1(tmpo(ix,iy),tmpo(ixp1,iy),
!     +	tmpo(ixp1,iyp1),tmpo(ix,iyp1))
!	output(i,j)=0.
!	else
!	write(6,*) 'legitimate sea ice at ', I,J
!	endif
C
C
!	endif

C        end if

	enddo
	enddo

	write(6,*) 'output from reansnowint'

	do J=jm,1,-jm/30
	write(6,65) (output(i,j),i=1,im,im/15)
	enddo
   65	format(30f6.1)

	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine reansoilice(ginput,nx,ny,im,jm,glat,glon,xpts,ypts,
     +					rtest,sice,sm,output)

	real ginput(nx,ny),tmpo(nx,ny),glat(im,jm),glon(im,jm)
	real xpts(im,jm),ypts(im,jm),output(im,jm),sm(im,jm),
     +			sice(im,jm),rtest(nx,ny)
     


	tmpo=ginput

  64	format(10(f5.1,1x))
  75	format(10(f6.1,1x))

	do J=1,jm
	do I=1,im

C        if((sm(i,j).gt.0.5).or.(sice(i,j).eq.1.0)) then
c                           set non land-mass values
C          output(i,j) = 0.0
C        else
c                           set land-mass values

	if (sice(i,j) .eq. 1.0) then

	y=ypts(i,j)
	x=xpts(i,j)
	  iy = y
          ix = x
          iyp1 = iy + 1
          ixp1 = ix + 1
c  First calculate the bilinear weights assuming a constant dx and dy
          wxy =     (ixp1-x) * (iyp1 - y)
          wxp1y =   (x-ix)   * (iyp1 - y)
          wxp1yp1 = (x-ix)   * (y-iy)
          wxyp1 = (ixp1 - x) * (y-iy)
c  Take care of the wrap around the 0 meridan when necessary
          if(ix.eq.nx) ixp1 = 1
c  Only use land mass points (ggsli=1) for the interpolation
C
C	here rtest is ice
C	
          if(rtest(ix,iy).eq.0.0) wxy = 0.0
          if(rtest(ixp1,iy).eq.0.0) wxp1y = 0.0
          if(rtest(ixp1,iyp1).eq.0.0) wxp1yp1 = 0.0
          if(rtest(ix,iyp1).eq.0.0) wxyp1 = 0.0
          sum = wxy + wxp1y + wxp1yp1 + wxyp1
          if(sum.eq.0.0) then
c  The eta land point is surrounded by gaussian grid non land
c  Look for nearest gaussian grid land point and assign it to the
c  Eta grid point. Find nearest at this latitude
            do jfnd = iy,nx
            do ifnd = ix,ny
            if(rtest(ifnd,jfnd).eq.1.0) then
                wxy = 1.0
                sum = 1.0
                ix = ifnd
                iy = jfnd
                go to 3356
            end if
            end do
            do ifnd = ix,1,-1
            if(rtest(ifnd,jfnd).eq.1.0) then
                wxy = 1.0
                sum = 1.0
                ix = ifnd
                iy = jfnd
                go to 3356
            end if
            end do
            end do
       print *,"Oh Oh can't find a land point at:",glat(i,j),glon(i,j)
       print *,((tmpo(ii,jj),ii=ix,ix+3),jj=iy,iy+2)
                wxy=1.0
                sum=1.0
                tmpo(ix,iy)=260.0
3356        continue
          end if
          wxy = wxy / sum
          wxp1y = wxp1y / sum
          wxp1yp1 = wxp1yp1 / sum
          wxyp1 = wxyp1 / sum
          output(i,j) = wxy     * tmpo(ix,iy)     +
     &                 wxp1y   * tmpo(ixp1,iy)   +
     &                 wxp1yp1 * tmpo(ixp1,iyp1) +
     &                 wxyp1   * tmpo(ix,iyp1)
        end if

	enddo
	enddo

   65	format(30f6.1)

	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine reaniceint(tmpo,nx,ny,im,jm,glat,glon,xpts,ypts,
     +				        sm,output)

	real tmpo(nx,ny),glat(im,jm),glon(im,jm)
	real xpts(im,jm),ypts(im,jm),output(im,jm),sm(im,jm)

	output=-99.

!	write(6,*) 'input ice field: '

	do J=1,ny/2
!	write(6,632) (tmpo(I,J),I=nx/2,nx,5)
	enddo
  632	format(50(f4.0,1x))

!	write(6,*) 'end input snow:'

  64	format(10(f5.1,1x))
  75	format(10(f6.1,1x))

	do J=1,jm
	do I=1,im

C
C	interpolate blindly with regard to land/sea?
C
C        if((sm(i,j).gt.0.5).or.(sice(i,j).eq.1.0)) then
        if((sm(i,j).lt.0.5)) then
C                           set non land-mass values
          output(i,j) = 0.0
        else
c                           set land-mass values
	y=ypts(i,j)
	x=xpts(i,j)
	  iy = y
          ix = x
          iyp1 = iy + 1
          ixp1 = ix + 1
c  Take care of the wrap around the 0 meridan when necessary
          if(ix.eq.nx) ixp1 = 1

          output(i,j) =  tmpo(ix,iy)     
Cmp
	endif

	enddo
	enddo

	write(6,*) 'OUTPUT ICE MASK'

	do J=jm,1,-jm/30
	write(6,65) (output(i,j),i=1,im,im/15)
	enddo
   65	format(30f6.1)

	return
	end
