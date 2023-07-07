	program select


	character*256 init_in(250),init_gdsdir,init_out
      namelist/model_grids/tlm0d,tph0d,im,jm,lm,ptinp,dlmd,dphd
     .                    ,dt,idtad,imonth,idate,iyear,istrtim
     .                    ,nsoil
     .                    ,ninit,tboco,init_in,
     .                     init_gdsdir,init_out

	real minlat,maxlat

	open(1,file='ETAIN',form='formatted',status='old')
        read(1,model_grids)

	IMJM=IM*JM-JM/2

	call corners(im,jm,imjm,tph0d,tlm0d,dlmd,dphd,
     +			minlat,maxlat)

C       write(6,*) 'latitude extremes are ', minlat, maxlat

	costrt=90-(maxlat+3)
	coend=90-(minlat-3)

C       write(6,*) 'costrt, coend: ', costrt, coend

	istrt=costrt/5+1
	iend=coend/5+1

C       write(6,*) 'istrt, iend ', istrt, iend

	open(unit=11,file='tmp.smask',form='formatted')
	write(11,273) istrt
	write(11,273) iend
  273	format(I2)

	END program select

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       subroutine CORNERS(im,jm,IMJM,tph0d,tlm0d,dlmd,dphd,minlat,
     +                                                maxlat)

        REAL,ALLOCATABLE::GLAT(:)
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
        ALLOCATE(GLAT(IMJM))

         DO 100 J=1,JM
      KHL0(J)=IM*(J-1)-(J-1)/2+1
      KHH0(J)=IM*J-J/2
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
      GLAT(K)=GLAT(K)/DTR

        if (glat(k) .gt. maxlat) maxlat=glat(k)
        if (glat(k) .lt. minlat) minlat=glat(k)

        enddo inna
        enddo outa

        DEALLOCATE(KHL0,KHH0)
        DEALLOCATE(GLAT)

        RETURN
      END

