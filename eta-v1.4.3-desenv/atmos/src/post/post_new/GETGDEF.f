	SUBROUTINE GETGDEF

C*	SUBROUTINE WRITTEN AS FIRST STEP IN A "FLEXIBLE" POST
C*	
C*	USED IN ASSOCIATION WITH THE WORKSTATION ETA, IT READS
C*	THE NECESSARY GRID INFORMATION FROM THE NAMELIST FILE 
C*	ETAIN, AND CALCULATES THE NECESSARY INFORMATION ABOUT THE
C*	OUTPUT GRID.  
C*
C*	98-12-21 Pyle--Initial version which only handles the native
C*	e-grid.  Based largely on the e-grid boundary finding code
C*	of E. Rogers.
C*
C*	99-04-06 PYLE--Changed to reflect new model_grids namelist
C*
C*	99-04-06 PYLE--Changed so copygb input fields would be generated
C*	here instead of in another script as was done previously
C*
C*	00-03-03 PYLE--Changed so will interpolate properly across the DL
C*	A bit of a fudge for lmbc grids since copygb seems troubled by
C*	the DL.  Also made some fixes for the southern hemisphere
C*

	INCLUDE "GDS.com"
	include 'ecommons.h'
        INCLUDE "mpif.h"

      namelist/model_grids/tlm0d,tph0d,im,jm,lm,ptinp,dlmd,dphd
     .                    ,dt,w,idtad,imonth,idate,iyear,istrtim
     .                    ,nsoil
     .                    ,ninit,init_in,rean_sfc,rean,
     .			   init_gdsdir,init_out
     .                    ,tboco,nhour
     .                    ,fcst_out
c
C
                             P A R A M E T E R
     & (IMAX=2500,JMAX=4500,IMJMMAX=IMAX*JMAX-JMAX/2)
                          D I M E N S I O N
     & KHL0  (JMAX),KHH0  (JMAX), GLAT(IMJMMAX),GLON(IMJMMAX)
                             D A T A
     & PI/3.141592654/

	character*4 tag
	character*3 proj

C-----------------------------------------------------------------------
c *** Read ETA namelist.
c
      open(1,file='ETAIN',form='formatted',status='old',err=900)
      read(1,model_grids,end=901)
      
      IF ( ME .EQ. 0 ) write(6,*) 'past read of ETAIN'
      close(1)
	WBD=-(IM-1.)*DLMD
	SBD=(-(JM-1.)/2.)*DPHD
	IF ( ME .EQ. 0 )write(6,*) 'IM= ', IM
 	IF ( ME .EQ. 0 )write(6,*) 'JM= ', JM	
	IMJM=IM*JM-JM/2
C
      IF ( ME .EQ. 0 ) 
     & write(6,*) 'grid centered at ',TPH0D, TLM0D

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
      DO 100 J=1,JM
      KHL0(J)=IM*(J-1)-(J-1)/2+1
      KHH0(J)=IM*J-J/2
  100 CONTINUE

C--------------GEOGRAPHIC LAT AND LONG OF TLL GRID POINTS---------------
              TPH=SB-DPH
              doout200: DO J=1,JM
              KHL=KHL0(J)
              KHH=KHH0(J)
C
              TLM=WB-TDLM+MOD(J+1,2)*DLM
              TPH=TPH+DPH
              STPH=SIN(TPH)
              CTPH=COS(TPH)
C
          doin200: DO K=KHL,KHH
      TLM=TLM+TDLM
      SPH=CTPH0*STPH+STPH0*CTPH*COS(TLM)
      GLAT(K)=ASIN(SPH)
      CLM=CTPH*COS(TLM)/(COS(GLAT(K))*CTPH0)-TAN(GLAT(K))*TAN(TPH0)
C      print*,CTPH,COS(TLM),COS(GLAT(K)),CTPH0,TAN(GLAT(K)),TAN(TPH0)
          IF(CLM.GT.1.)      CLM=1.
      FACT=1.
          IF(TLM.GT.0.)      FACT=-1.
      GLON(K)=(-TLM0D*DTR+FACT*ACOS(CLM))/DTR

Cmp	at this point GLON is in DEGREES WEST
	if (GLON(K) .lt. 0) GLON(K)=GLON(K)+360.
	if (GLON(K) .gt. 360.) GLON(K)=GLON(K)-360.
	if (GLON(K) .lt. 180) GLON(K)=-GLON(K)         ! make WH negative
	if (GLON(K) .gt. 180) GLON(K)=360.-GLON(K)     ! make EH 

	GLAT(K)=GLAT(K)/DTR

       END DO doin200
       END DO doout200
C
	write(6,*) ' '
	ilat=(glat(1)*1000.)
	ilon=(glon(1)*1000.)
	iclat=(TPH0D*1000.)
	iclon=(TLM0D*1000.)
	idelx=(DLMD*1000.)
	idely=(DPHD*1000.)
Cmp	write(6,*) ilat,ilon,iclat,iclon,idelx,idely
 1020 format('LOWER LEFT  POINT= ',f12.5,x,f12.5)
 1021 FORMAT('UPPER RIGHT POINT= ',F12.5,X,F12.5)
C
Cmp	define the grid in w3fi71 format
	IGDS(1)=0
	IGDS(2)=255
	IGDS(3)=203
	IGDS(4)=IM
	IGDS(5)=JM
	IGDS(6)=ilat
	IGDS(7)=ilon
	IGDS(8)=136
	IGDS(9)=iclat
	IGDS(10)=iclon
	IGDS(11)=idelx
	IGDS(12)=idely
	IGDS(13)=64
	IGDS(14)=0

	IF ( ME .EQ. 0 )write(6,*) 'IGDS= ',(IGDS(I),I=1,18)

	tag='ETWS'
	proj='CED'
	ilat2=0
	IMT=IM*2-1
	JMT=JM
	frac=0.5
	open(unit=63,file='grdnav.tbl',form='formatted')
	write(63,207) tag,IGDS(2),proj,int(TPH0D),int(TLM0D),ilat2,
     +	glat(1),glon(1),glat(imjm),glon(imjm),IMT,JMT,frac,ilat2
  207	format(A4,x,I3,x,A3,x,I3,x,I4,x,I3,x,f7.3,x,f8.3,x,f7.3,x,f8.3,x,
     +	I3,x,I3,x,f4.2,x,I1)


Cmp	get the lat-lon info for copygb

Cmp	set dlat a bit coarser than the filled e-grid
Cmp	dlat=DPHD*1.25
	dlattmp=(1./(DPHD*1.25))
	dlat=1./dlattmp
	if (TPH0D .ge. 0) wlon=(glon(1))+1
	if (TPH0D .lt. 0) wlon=(glon(imjm-im+1))+1
        slat=(TPH0D+SBD)+1
Cmp
	if (wlon .gt. 0 .and. TLM0D .lt. 0) then
	wloncalc=wlon-360.
	else
	wloncalc=wlon
	endif

        delta=abs(TLM0D-wloncalc)
        elon=(TLM0D+delta)-1
        nlat=(GLAT(IMJM))-1

Cmp     search the lat values corresponding to wlon

C	if (tph0d .ge. 0) then
	
	IF ( ME .EQ. 0 ) 
     &  write(6,*) 'target longitude is ' , wlon

        do K=1,imjm
        dif=abs(wlon-glon(k))
        if (dif .lt. 1 .and. dif .gt. 0.) then
C        IF ( ME .EQ. 0 )
C    &  write(6,*) 'changing tnlat to ', tnlat
        tnlat=glat(k)
        endif
        enddo

C	elseif (tph0d .lt. 0) then

C	do K=imjm,1,-1
C	dif=abs(wlon-glon(k))
C        if (dif .lt. 1 .and. dif .gt. 0.) then
C        tnlat=glat(k)
C        endif
C        enddo

C	endif

Cmp	dateline add
       IF ( ME .EQ. 0 )
     & write(6,*) 'going into check ', wlon, elon
C	extend from West into East
	if (TLM0D .lt. 0 .and. wlon .gt. 0) then
	wloncalc=wlon-360.
	eloncalc=elon
C	extend from East into West
	elseif (TLM0D .gt. 0 .and. elon .lt. 0) then
	wloncalc=wlon
	eloncalc=elon+360.
	else
C	normal
	eloncalc=elon
	wloncalc=wlon
	endif
       IF ( ME .EQ. 0 )
     & write(6,*) 'found eloncalc,wloncalc ', eloncalc,wloncalc
C	STOP

Cmp end dateline add
        nlat=(tnlat)-1
        jdim=((nlat-slat)/dlat)+1
Cmp        idim=((elon-wlon)/dlat)+1
        idim=((eloncalc-wloncalc)/dlat)+1
        islat=slat*1000
        iwlon=wlon*1000
        nlat=nlat*1000
Cmp
Cmp	if (elon .gt. 180) elon=elon-360.
        ielon=elon*1000
        idlat=(dlat*1000)
        open (unit=23,file='outjob_input_lat',form='formatted',
     +  access='sequential',status='unknown')
        write(23,299) IDIM,JDIM,islat,iwlon,nlat,ielon,idlat,idlat
  299   format('255 0 ',2(I4,x),I6,x,I7,x,'128 ',I6,x,I7,x,2(I3,x),'64')
	close(23)

Cmp     now compute the values needed for a lambert conic conformal
Cmp     projection.  Use dimensions of filled e-grid

        IDIM=IM*2-1
        JDIM=JM

        idx=INT(DPHD*111*1000)
        idy=idx

Cmp        LATONE=INT(GLAT(1)*1000)
Cmp        LONONE=INT(GLON(1)*1000)
        LATONE=INT(GLAT(2)*1000)
        LONONE=INT(GLON(2)*1000)

Cmp
	IF ( ME .EQ. 0 )write(6,*) 'LONONE= ', LONONE
	if (LONONE .gt. 0 .and. TLM0D .lt. 0) then
	tlmuse=179.
C	tlmuse=TLM0D
	else
	tlmuse=TLM0D
	endif
	
        open (unit=24,file='outjob_input_lmbc',form='formatted',
     +  access='sequential',status='unknown')
        write(24,301) IDIM,JDIM,latone,lonone,int(tlmuse*1000),idx,idy,
     +  int(TPH0D*1000),INT(TPH0D*1000)

  301   format('255 3 ',2(I3,x),I6,x,I7,x,'8 ',I7,x,2(I6,x),'0 64',
     +  2(x,I6))

  900	IF ( ME .EQ. 0 )print*, 'end of namelist file'
  901	IF ( ME .EQ. 0 )print*, 'model list trouble'
	RETURN
      END

C***************************************************************************

	BLOCK DATA GDSVAR
	include "GDS.com"
	
	DATA IGDS/18*0/
	end

C***********************************************************************

