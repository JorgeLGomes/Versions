      program select

      character*256 init_in(250),init_gdsdir,init_out,
     .              froot,sname,fname
      character*1 londir,latdir
      character*2 alat
      character*3 alon
      namelist/model_grids/tlm0d,tph0d,im,jm,lm,ptinp,dlmd,dphd
     .                    ,dt,idtad,imonth,idate,iyear,istrtim
     .                    ,nsoil
     .                    ,ninit,tboco,init_in,
     .                     init_gdsdir,init_out
      namelist/lsmkdata/ froot,res

      real minlat,maxlat
      real minlon,maxlon	!isr: incluído
      parameter(ntilesize=1)

      open(1,file='ETAIN',form='formatted',status='old')
      read(1,model_grids)

      IMJM=IM*JM-JM/2
      
      open(30,file='seamaskdata.dat',status='old')
      read(30,lsmkdata)
      dlat=1./res
      npt=int(res)
      close(30)

      call corners(im,jm,imjm,tph0d,tlm0d,dlmd,dphd, !isr: incluído
     .	           minlat,minlon,maxlat,maxlon)      !isr: incluído
     
!     --------------------------------------------------------------

      west=minlon-0.1
      east=maxlon+0.1 
      slat=minlat-0.1
      rnlat=maxlat+0.1

      write(*,*) 'west:  ', west
      write(*,*) 'east:  ', east
      write(*,*) 'north: ', rnlat
      write(*,*) 'south: ', slat
      write(*,*) res,dlat
      
      cshift=dlat/2.     

      IST=((west+180.-cshift)/dlat)+1
      IEND=((east+180.-cshift)/dlat)+1

      JST=-(RNLAT-60)/dlat+1
      JEND=-(SLAT-60)/dlat+1

      JSUNIT=INT((JST-1)/npt)+1
      JENDUNIT=INT((JEND-1)/npt)+1

      ISUNIT=INT((IST-1)/npt)+1
      IENDUNIT=INT((IEND-1)/npt)+1
      
!     -------------------------------------------------------------- 	
      
      open(11,file='tmp.smask',form='formatted')
      
!     --------------------------------------------------------------

      write(*,*)ISUNIT,JENDUNIT
      write(*,*)JSUNIT,JENDUNIT

      do J=JSUNIT,JENDUNIT
      lat=60-J*ntilesize
      if(lat.lt.0) latdir='s'
      if(lat.ge.0) latdir='n'

      do I=ISUNIT,IENDUNIT
      lon=-180+(I-1)*ntilesize
      if(lon.lt.0) londir='w'
      if(lon.ge.0) londir='e'
	
      write(alon,'(i3.3)') abs(lon)
      write(alat,'(i2.2)') abs(lat)
      sname=londir//alon//latdir//alat
      write(11,*) trim(sname)
      enddo                                
      enddo

!     ---------------------------------------------------------------		

      close(11)

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!isr: inclui a corners contida no programa etatopo.F

       subroutine CORNERS(im,jm,IMJM,tph0d,tlm0d,dlmd,dphd,minlat,wlon,
     +                                                maxlat,elon)
C
C
C     *  ROUTINE TO FIND EARTH LATITUDE/LONGITUDE FOR THE CORNER       *
C     *               POINTS OF AN ETA MODEL GRID                      *
C

C-----------------------------------------------------------------------
C                             D I M E N S I O N
C     & KHL0  (JM),KHH0  (JM), GLAT(IMJM),GLON(IMJM)

	REAL,ALLOCATABLE::GLAT(:),GLON(:)
	INTEGER,ALLOCATABLE::KHL0(:),KHH0(:)

                             D A T A
     & PI/3.141592654/

        REAL DLMD,DPHD,WBD,SBD,TPH0D,TLM0D,minlat,wlon,maxlat,elon

C*******************************
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
C
      STPH0=SIN(TPH0)
      CTPH0=COS(TPH0)
C
	ALLOCATE(KHL0(JM),KHH0(JM))
	ALLOCATE(GLAT(IMJM),GLON(IMJM))

         DO 100 J=1,JM
      KHL0(J)=IM*(J-1)-(J-1)/2+1
      KHH0(J)=IM*J-J/2
C     WRITE(6,9999) J, KHL0(J), KHH0(J)
C9999 FORMAT(2X,3(I10,1X))
  100 CONTINUE
C--------------GEOGRAPHIC LAT AND LONG OF TLL GRID POINTS---------------
              TPH=SB-DPH
        maxlat=-999.
	minlat=99.
       outa:DO J=1,JM
              KHL=KHL0(J)
              KHH=KHH0(J)
C
              TLM=WB-TDLM+MOD(J+1,2)*DLM
              TPH=TPH+DPH
              STPH=SIN(TPH)
              CTPH=COS(TPH)
      inna:DO K=KHL,KHH
      TLM=TLM+TDLM
      SPH=CTPH0*STPH+STPH0*CTPH*COS(TLM)
      GLAT(K)=ASIN(SPH)
      CLM=CTPH*COS(TLM)/(COS(GLAT(K))*CTPH0)-TAN(GLAT(K))*TAN(TPH0)
          IF(CLM.GT.1.)      CLM=1.
      FACT=1.
          IF(TLM.GT.0.)      FACT=-1.
      GLON(K)=(-TLM0D*DTR+FACT*ACOS(CLM))/DTR

Cmp     at this point GLON is in DEGREES WEST
        if (GLON(K) .lt. 0) GLON(K)=GLON(K)+360.
        if (GLON(K) .gt. 360.) GLON(K)=GLON(K)-360.
        if (GLON(K) .lt. 180) GLON(K)=-GLON(K)         ! make WH negative
        if (GLON(K) .gt. 180) GLON(K)=360.-GLON(K)     ! make EH

        GLAT(K)=GLAT(K)/DTR

        if (glat(k) .gt. maxlat) maxlat=glat(k)
	if (glat(k) .lt. minlat) minlat=glat(k)

       enddo inna
       enddo outa

	DEALLOCATE(KHL0,KHH0)
	
	
	if (TPH0D .ge. 0) then
        wlon=glon(imjm-im+1)
        elon=glon(imjm)
	else
	wlon=glon(1)
	elon=glon(im)
	endif

C	write(6,*) 'raw lon values (w,e) ', wlon, elon
        if (tlm0d .lt. 0 .and. wlon .gt. 0) wlon=wlon-360.
        if (tlm0d .gt. 0 .and. elon .lt. 0) elon=elon+360.

	DEALLOCATE(GLAT,GLON)

        RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
