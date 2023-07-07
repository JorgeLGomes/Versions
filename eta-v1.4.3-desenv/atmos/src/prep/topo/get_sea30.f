	subroutine get_sea30(isea)

	include 'ecommons.h'

        parameter(ID1=1,ID2=43204,IRECL=ID2)
	parameter(dlat=1./120.)

	integer ounit,iunit,nunit
        parameter(OUNIT=50)
        character*12 fnameout
        INTEGER*2 LANDIN(ID1,ID2)
        INTEGER,ALLOCATABLE:: LAND2D(:,:)
        INTEGER GDS(200)
        character*2 uname
        CHARACTER*1 LANDIN_CHR(ID1,ID2)

        REAL, ALLOCATABLE:: GLATR(:,:),GLONR(:,:)
	real isea(imjm)
        real dum2d(im,jm)

	cshift=dlat/2.

C
C
C	details about target grid

C       compute lat lon on e-grid
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

        ALLOCATE(GLATR(IM,JM),GLONR(IM,JM))

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
C
         enddo
      enddo

	rnlat=glatr(im/2,jm)*r2d+3
	slat=glatr(1,1)*r2d-3
	west=glonr(1,jm)*r2d-3
	east=glonr(im,jm)*r2d+3

        IST=((west+180.-cshift)/dlat)+1
        IEND=((east+180.-cshift)/dlat)+1

        wtrue=(ist-1)*dlat+cshift-180.
        etrue=(iend-1)*dlat+cshift-180.

        JST=-(RNLAT-90)/dlat+1
        JEND=-(SLAT-90)/dlat+1

        ntrue=90.-(jst-1)*dlat-cshift
        strue=90.-(jend-1)*dlat-cshift

        ILIM=IEND-IST+1
        JLIM=JEND-JST+1
        allocate(land2d(ILIM,JLIM))

        JSUNIT=INT((JST-1)/600)+1
        JENDUNIT=INT((JEND-1)/600)+1

        do J=JSUNIT,JENDUNIT
        write(uname,'(i2.2)') INT(J)
        fnameout='smask.30s.'//uname

        IIUNIT=J+60
        open(unit=IIUNIT,file=fnameout,
     +                  access='direct',recl=43204)
        write(6,*) 'opened unit,file ', IIUNIT,fnameout
        enddo

        GDS=0

        GDS(1)=0
        GDS(2)=ILIM
        GDS(3)=JLIM
        GDS(4)=INT(NTRUE*1000)
        GDS(5)=INT(WTRUE*1000)
        GDS(6)=128
        GDS(7)=INT(STRUE*1000)
        GDS(8)=INT(ETRUE*1000)
        GDS(9)=INT(1./120.*1000)
        GDS(10)=INT(1./120.*1000)

        do 80 N=JST,JEND

        IUNIT=INT((N-1)/600)+61
        IREC=N-(IUNIT-61)*600

        READ(IUNIT,REC=IREC,ERR=100) LANDIN_CHR

          DO IID1 = 1,ID1
            DO IID2 = 1,ID2
              LANDIN(IID1,IID2) = ICHAR(LANDIN_CHR(IID1,IID2))
            END DO
          END DO

        do I=IST,IEND
        land2d(I-IST+1,N-JST+1)=LANDIN(1,I)
        enddo

   80   continue

        ISEA=1

        KHL22=2*IM+1
        KHH22=IMJM-2*IM

        do K=KHL22,KHH22
C       convert K value to get I,J
        IMT=2*IM-1
        JX=(K-1)/IMT+1
        IX=K-(JX-1)*IMT
        IF(IX.LE.IM)THEN
          I=IX
          J=2*JX-1
        ELSE
          I=IX-IM
          J=2*JX
        ENDIF

C DETERMINE LAT/LON OF TARGET (output) GRID POINT
C
        GLATD = GLATR(i,j)*r2d
        GLOND = GLONR(I,J)*r2d
        IF (GLOND .LT. 0) GLOND=GLOND+360
        GLOND = 360. - GLOND

        call ced_ij(glatd,glond,x,y,gds)
        II=INT(X+0.5)
        JI=INT(Y+0.5)

C       If input data is 1 --> land --> isea=0
C       If input data is 0 --> water --> isea=1
C


        if (land2d(II,JI) .eq. 1) then
        isea(k)=0
        elseif (land2d(II,JI) .eq. 0) then
        isea(k)=1
        else
        write(6,*) 'bad land2d values!!!! ', II,JI,land2d(II,JI)
        endif


        enddo

        call conh12(isea,imjm,1,dum2d,im,jm)
        do J=JM-1,1,-JM/40
        write(6,273)(dum2d(I,J),I=1,IM,IM/53)
        enddo
  273   format(71f2.0)

	RETURN

	STOP
  100	write(6,*) 'bad read!!!!'

	END
