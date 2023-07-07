
       program CORNERS
!
!
!     *  ROUTINE TO FIND EARTH LATITUDE/LONGITUDE FOR THE CORNER       *
!     *               POINTS OF AN ETA MODEL GRID                      *
!

!-----------------------------------------------------------------------
!                             D I M E N S I O N
!C     & KHL0  (JM),KHH0  (JM), GLAT(IMJM),GLON(IMJM)

	REAL,ALLOCATABLE::GLAT(:),GLON(:)
	INTEGER,ALLOCATABLE::KHL0(:),KHH0(:)

        DATA PI/3.141592654/

        REAL DLMD,DPHD,WBD,SBD,TPH0D,TLM0D,minlat,wlon,maxlat,elon

!C*******************************

       open(4,file= './namelist',status='old',action='read')
       read(4,*)im
       read(4,*)jm
       read(4,*)IMJM
       read(4,*)tph0d
       read(4,*)tlm0d
       read(4,*)dlmd
       read(4,*)dphd
       read(4,*)res
       close(4)

        WBD=-(float(IM)-1.)*DLMD
        SBD=(-(float(JM)-1.)/2.)*DPHD

      DTR=PI/180.
      TPH0=TPH0D*DTR
      WB=WBD*DTR
      SB=SBD*DTR
      DLM=DLMD*DTR
      DPH=DPHD*DTR
      TDLM=DLM+DLM
      TDPH=DPH+DPH
!C
      STPH0=SIN(TPH0)
      CTPH0=COS(TPH0)
!C
	ALLOCATE(KHL0(JM),KHH0(JM))
	ALLOCATE(GLAT(IMJM),GLON(IMJM))

      DO J=1,JM
        KHL0(J)=IM*(J-1)-(J-1)/2+1
        KHH0(J)=IM*J-J/2
        !WRITE(6,9999) J, KHL0(J), KHH0(J)
        !FORMAT(2X,3(I10,1X))
      enddo
!C--------------GEOGRAPHIC LAT AND LONG OF TLL GRID POINTS---------------
              TPH=SB-DPH
        maxlat=-999.
	minlat=99.
      doout: DO J=1,JM
             KHL=KHL0(J)
             KHH=KHH0(J)
!C
             TLM=WB-TDLM+MOD(J+1,2)*DLM
             TPH=TPH+DPH
             STPH=SIN(TPH)
             CTPH=COS(TPH)
             doin: DO K=KHL,KHH
                   TLM=TLM+TDLM
                   SPH=CTPH0*STPH+STPH0*CTPH*COS(TLM)
                   GLAT(K)=ASIN(SPH)
                   CLM=CTPH*COS(TLM)/(COS(GLAT(K))*CTPH0)-TAN(GLAT(K))*TAN(TPH0)
                   IF(CLM.GT.1.)      CLM=1.
                   FACT=1.
                   IF(TLM.GT.0.)      FACT=-1.
                   GLON(K)=(-TLM0D*DTR+FACT*ACOS(CLM))/DTR

                   !Cmp  at this point GLON is in DEGREES WEST
                   if (GLON(K) .lt. 0) GLON(K)=GLON(K)+360.
                   if (GLON(K) .gt. 360.) GLON(K)=GLON(K)-360.
                   if (GLON(K) .lt. 180) GLON(K)=-GLON(K)         ! make WH negative
                   if (GLON(K) .gt. 180) GLON(K)=360.-GLON(K)     ! make EH

                   GLAT(K)=GLAT(K)/DTR

                   if (glat(k) .gt. maxlat) maxlat=glat(k)
	           if (glat(k) .lt. minlat) minlat=glat(k)

             enddo doin
      enddo doout

	DEALLOCATE(KHL0,KHH0)
	
	if (TPH0D .ge. 0) then
        wlon=glon(imjm-im+1)
        elon=glon(imjm)
	else
	wlon=glon(1)
	elon=glon(im)
	endif

!	write(6,*) 'raw lon values (w,e) ', wlon, elon
        if (tlm0d .lt. 0 .and. wlon .gt. 0) wlon=wlon-360.
        if (tlm0d .gt. 0 .and. elon .lt. 0) elon=elon+360.

        rnlat=maxlat+0.1
        slat=minlat-0.1
        west=wlon-0.1
        east=elon+0.1
        open(5,file= './corners.out',status='unknown')
        write(5,*) 'west= ', west
        write(5,*) 'east= ', east
        write(5,*) 'north= ', rnlat
        write(5,*) 'south= ', slat
        write(5,*) 'res=  ',res
        write(5,*) int((abs(west-east))/res)
	write(5,*) int((abs(rnlat-slat))/res)
        close(5)
	DEALLOCATE(GLAT,GLON)
      STOP
      END

