	subroutine coastalfix(egridval,rawval,sm,im,jm,nxin,nyin,
     +				dphd,dlmd,
     +				wbd,sbd,tph0d,tlm0d,gds)
C
C============================================================================
C
C Subroutine written 6 June 2000 to hopefully fix once and for all the
C coastal "smearing" problem seen in soil fields that are bilinearly
C interpolated onto the Eta model's e-grid.  
C
C This routine will make an effort to avoid interpolating soil values
C that are corrupted by being overly close to the coast.  If any of
C the neighbor input points are water, new values of the soil fields
C are calculated as an average of the surrounding non-water points.  
C Not perfect, but better than it was.
C
C This subroutine takes its inspiration (and the nearest neighbor code)
C from the SFCEDAS routine of M. Ek and K. Mitchell.
C
C============================================================================

        real egridval(im,jm,12),rawval(nxin,nyin,12),sm(im,jm)
        real hlat(im,jm),hlon(im,jm),glatr(im,jm),glonr(im,jm)
        real redo(im,jm),rsum(8)
        integer gds(200),icount(8)
        parameter (r2d=57.2957795)
        parameter (dtr=3.141592654/180.)
        parameter (nsoil=8)

       KRADUL=1
      print*,'KRADUL',KRADUL


        write(6,*) 'rawval land/sea mask '
        do J=NYIN,1,-10
        write(6,643) (rawval(I,J,1),I=1,NXIN,13)
        enddo
  643   format(30(f2.0))

C COMPUTE OUTPUT GRID lats AND lons

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
            glatr(i,j)=asin(sinphi)
            coslam=ctph*cos(tlm)/(cos(glatr(i,j))*ctph0)
     .            -tan(glatr(i,j))*tan(tph0)
            coslam=min(coslam,1.)
            fact=1.
            if (tlm .gt. 0.0) fact=-1.
            glonr(i,j)=-tlm0d*dtr+fact*acos(coslam)
        hlat(i,j)=glatr(i,j)*r2d
        hlon(i,j)=glonr(i,j)*r2d
C
        if (hlon(i,j).lt.0) hlon(i,j)=hlon(i,j)+360.
C
        if (mod(I,30).eq.0.and.mod(J,30).eq.0) then
C	write(6,*) 'hlat,hlon: ', hlat(i,j),hlon(i,j)
        endif
         enddo
      enddo

C
C  Begin looping through e-grid domain
C	
         redo=0.

        DO 200 J=1,JM
        DO 100 I=1,IM

C IF NON LAND-MASS POINT ON OUTPUT GRID, SKIP TO END OF LOOP (100)
C
          IF ((SM(I,J) .GT. 0.9)) GOTO 100

C DETERMINE LAT/LON OF TARGET (output) GRID POINT
C
        GLATD = HLAT(I,J)
        GLOND = 360. - HLON(I,J)
C
C DETERMINE NEAREST NEIGHBOR FROM INPUT GRID
C
        if (GDS(1).eq.5) then
          call str_ij(GLATD,glond,x,y,gds)
        elseif (GDS(1).eq.3) then
          call lcc_ij(GLATD,glond,x,y,gds)
        elseif(GDS(1).eq.0) then
          call ced_ij(GLATD,glond,x,y,gds)
        else
        write(6,*) 'BAD GDS(1) VALUE...NOT SUPPORTED...STOPPING!'
        write(6,*) 'GDS(1)= ', GDS(1)
        write(6,*) 'GDS(1)= ', GDS(1)
        STOP 666
        endif

        II=INT(X+0.5)
        JI=INT(Y+0.5)

C	rawval has land/sea mask opposite to that of sm (i.e., water=0.)
C
C SEARCH FOR ANY SURROUNDING WATER POINTS (TO A LIMIT OF KRAD=KRADUL)
C
        DO 50 KRAD=1,KRADUL 
          DO 40 JY = JI-KRAD,JI+KRAD 
            DO 30 IX = II-KRAD,II+KRAD
C
C-----------------------------------------------------------------------
C CHECK TO SEE IF THIS NEAREST NEIGHBOR (IX,JY) STILL WITHIN IMI,JMI
C DOMAIN
C
                  IF ( (IX .LT. 1) .OR. (IX .GT. NXIN) .OR.
     +                 (JY .LT. 1) .OR. (JY .GT. NYIN) )  THEN
                        GOTO 30
                  ENDIF

C IS THIS NEIGHBOR POINT WATER?  IF SO THE VALUE HERE IS SUSPICIOUS
C
C	write(6,*) 'NN , rawval ', IX,JY,rawval(ix,jy,1)
C
                  IF ( rawval(IX,JY,1) .lt. 0.4) THEN
			redo(i,j)=1.
			goto 52
		  ENDIF

  30            CONTINUE
  40          CONTINUE
  50        CONTINUE
C
C-----------------------------------------------------------------------
C
  52	continue

			IF (redo(i,j) .eq. 1) THEN

C	write(6,*) 'redoing at ', i,j

	icount=0
	rsum=0.

	DO 55 KRAD=1,KRADUL
	  DO 45 JY=JI-KRAD,JI+KRAD
	    DO 35 IX=II-KRAD,II+KRAD
C
                  IF ( (IX .LT. 1) .OR. (IX .GT. NXIN) .OR.
     +                 (JY .LT. 1) .OR. (JY .GT. NYIN) )  THEN
                        GOTO 35
                  ENDIF
C
C MAKE SURE THIS POSSIBLE REPLACEMENT ISN'T WATER
C
                  IF ( rawval(IX,JY,1) .lt. 0.4) THEN
                        goto 35
		  ELSE

		do IR=1,8
		rsum(IR)=rsum(IR)+rawval(IX,JY,IR+4)
		icount(IR)=icount(IR)+1
		enddo
			
                  ENDIF

C	
  35            CONTINUE
  45          CONTINUE
  55        CONTINUE

!GSM_JORGE          DO 70 K=1,NSOIL
          DO 70 K=1,NSOIL-4
	    tmphold=egridval(I,J,2*K+3)
Cmp
	if (icount(2*K-1) .gt. 0) then
            egridval(I,J,2*K+3)=rsum(2*K-1)/icount(2*K-1)
	endif

	if (icount(K) .eq. 0 .and. K .eq. 1) then
C	write(6,*) 'NO VALID POINTS....KEEP ORIGINAL VALUE!!!! '
	endif

	if (icount(2*K) .gt. 0) egridval(I,J,2*K+4)=rsum(2*K)/icount(2*K)
	    if (abs(tmphold-egridval(I,J,2*K+3)).gt. 5.) then
C	     write(6,*) 'soil temp: old,new ',tmphold,egridval(I,J,2*K+3)
	    endif
  70      CONTINUE

			ENDIF

 100	CONTINUE
 200	CONTINUE
	
	itot=0

	do J=1,jm
	do I=1,im
	itot=itot+redo(i,j)
	enddo
	enddo

	write(6,*) 'replaced soil values at : ', itot, ' points'

	RETURN
	END
