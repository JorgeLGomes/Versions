	program staids_prep


C	program written 18 June 2001 for workstation Eta.
C
C	It will read the information from the ETAIN file (IM,JM,
C	DLMD,DPHD,TLM0D,TPH0D), compute the needed WBD and SBD, then
C	write a couple of lines to the STALST.f code


	REAL:: tlm0d,tph0d,ptinp,dlmd,dphd,dt,tboco
	INTEGER:: im,jm,lm,idtad,imonth,idate,iyear,istrtim,
     +		  nsoil,ninit,nhour
	CHARACTER*256:: init_in(50),init_gdsdir,init_out,rean_sfc
	LOGICAL:: rean

      namelist/model_grids/tlm0d,tph0d,im,jm,lm,ptinp,dlmd,dphd
     .                    ,dt,idtad,imonth,idate,iyear,istrtim
     .                    ,nsoil
     .                    ,ninit,tboco,init_in,init_gdsdir
     .                    ,init_out


	open(1,file='ETAIN',form='formatted',status='old')
	read(1,model_grids)

	NHOUR=(NINIT-1)*TBOCO

	write(6,*) 'IM= ', IM
	write(6,*) 'JM= ', JM
	write(6,*) 'DLMD= ', DLMD
	write(6,*) 'DPHD= ', DPHD
	write(6,*) 'TLM0D= ', TLM0D
	write(6,*) 'TPH0D= ', TPH0D

C	compute WBD, SBD

	WBD = -(IM-1)*DLMD
	SBD = -((JM-1)/2.)*DPHD
	write(6,*) 'computed WBD, SDB= ', WBD, SBD

	open(unit=63,file='stalst_add.txt',form='formatted')
	write(63,644) WBD,SBD,TPH0D,TLM0D
	write(63,645) DPHD,DLMD
	write(63,646)

	open(unit=64,file='parmsndp',form='formatted')
	write(64,744) LM,NHOUR+1

  644	format(7x,'DATA WBD/',f9.5,'/,SBD/',f9.5,'/,TPH0D/',
     +	       f6.2,'/,TLM0D/',f7.2,'/')
  645	format(7x,'DATA DPHD/',f8.6,'/,DLMD/',f8.6,'/')
  646	format(7x,'END BLOCK DATA STN')
  744	format(7x,'PARAMETER (LM=',I3,', NPNT= ',I2,')')


	end
