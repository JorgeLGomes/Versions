
! --------------------------------------------------------------------------
!  PROGRAMA gera_sondagens
!
!           Autor: ERNANI DE LIMA NASCIMENTO
!           Data: OUTUBRO de 2006
!
!  Este programa lê arquivo binário do grads gerado pelo script SOUND3.EXEC
!  e gera um arquivo de saída ASCII contendo as 3264 sondagens do ETA para o caso
!  do tornado de Indaiatuba. A formatação desta saída satisfaz a leitura para
!  o programa indices_severos5.f90
! ---------------------------------------------------------------------------
       program gera_sondagens
       IMPLICIT NONE
!   Parameter    
       integer, parameter   :: Max_Sond_Levs=200
       real,    parameter   :: T0=273.15
 
  ! Parametros do Modelo Passados via namelist
       integer              :: im	      ! Numero de pontos em x
       integer              :: jm	      ! Numero de pontos em y
       integer              :: lm 	      ! Numero de niveis
       integer              :: FreqModOut     ! Frequencia de saida 
       real                 :: loni	      ! Longitude a oeste 
       real                 :: lati	      ! Latitude ao sul
       real                 :: res	      ! Resolucao do modelo
       integer              :: idia	      ! Dia da (passado via namelist)
       integer              :: imes           ! Mes    (passado via namelist)
       integer              :: iano           ! Ano    (passado via namelist)
       integer              :: ihora          ! Hora   (passado via namelist)

!____________________________________________________________________________  

! Variaveis 
       integer              :: imjm
       integer              :: lmp1

  ! Variáveis de trabalho
       integer              :: i,j,k,l,n,irec,geop,lnlev,umrl,strm_opt,ind,l500hPa
       real                 :: est, estd,wvmr2m,wvmr,vap
       real                 :: mlevs(Max_Sond_Levs)

       integer, allocatable   :: station(:),istart(:),nlev(:)
! Variáveis de entrada
       integer            :: datsav(11),kpds(200)
       character(len=60) 	   :: ARQIN,ARQOUT
       character(LEN=12):: atime
       real, allocatable  :: lev(:,:)
       real, allocatable  :: umes(:,:),zgeo(:,:),temp(:,:),td(:,:),uvel(:,:),vvel(:,:)
       real, allocatable  :: theta(:,:),rh(:,:),wvmrkg(:,:)

       namelist/model_grids/im,jm,lm,loni,lati,res,FreqModOut,        &
     &                     idia,imes,iano,ihora,ARQIN,ARQOUT,strm_opt 
                       
       namelist/model_levs/mlevs
 
       open(1,file='namelist_grib',form='formatted',status='old')
       rewind(1)
       read(1,model_grids)
       rewind(1)
       read(1,model_levs)
       close(1)

       imjm=im*jm
       lmp1=lm+1       

      n=index(ARQIN,' ')-1
      write(6,*) 'GRIBFILE is ', ARQIN(1:n)
      CALL GET_GDS(ARQIN(1:n),datsav,kpds)

!
!     X and Y dimensions now obtained from the degribbed GRIB file
!
       im=datsav(1)
       jm=datsav(2)
       imjm=im*jm
 
       allocate(lev(lmp1,imjm))
       ALLOCATE(station(imjm+1),istart(imjm),nlev(imjm))
       allocate(zgeo(lmp1,imjm),temp(lmp1,imjm),td(lmp1,imjm),uvel(lmp1,imjm))
       allocate(vvel(lmp1,imjm),theta(lmp1,imjm),rh(lmp1,imjm),wvmrkg(lmp1,imjm))
       allocate(umes(lmp1,imjm))
      
       write(*,*)idia,imes,iano,ihora &
      &                    ,ARQIN,ARQOUT,strm_opt
       write(*,*)im,jm,lm,loni,lati,res,FreqModOut,lmp1,imjm
       lev(1,1:imjm)=9999.
       do k=1, imjm
        lev(2:lmp1,k)=mlevs(1:lm)
       enddo      
       write(*,*)lev(1:lmp1,1)
       write(*,*) 'Utilizar deslocamento observado da tempestade? (1=sim 5=não)'
       write(*,*) 'strm_opt ',strm_opt
       if ((strm_opt /= 1).AND.(strm_opt /= 5)) then
        write(*,*) 'Opção não disponível. Programa interrompido.'
        stop
       endif
       write(*,*)'call read_degrib'
       call read_degrib(ARQIN(1:n),imjm,lmp1,lev,zgeo,temp,td,uvel,vvel,rh,umes,atime)

!    Tranforma Temp (K) => Temp (C)
!       variab(1:lmp1,3,1:imjm)=variab(1:lmp1,3,1:imjm)-T0
!       variab(1:lmp1,4,1:imjm)=variab(1:lmp1,4,1:imjm)-T0
        temp(1:lmp1,1:imjm)=temp(1:lmp1,1:imjm)-T0
        td(1:lmp1,1:imjm)=td(1:lmp1,1:imjm)-T0

       do k=1, imjm
        theta(1,k)=exp(alog(temp(1,k)+T0)+0.286*alog(1000/lev(1,k)))
        est=6.112*exp((17.67*temp(1,k))/(temp(1,k)+243.5))
        estd=6.112*exp((17.67*td(1,k))/(td(1,k)+243.5))
        rh(1,k)=100*(estd/est)
        wvmr2m=(0.62197*estd)/(lev(1,k)-estd)
        wvmrkg(1,k)=wvmr2m*1000
       enddo
       do k=1, imjm
        do l=2, lmp1
         theta(l,k)=exp(alog(temp(l,k)+T0)+0.286*alog(1000/lev(l,imjm)))
         wvmr=(umes(l,k)/(1-umes(l,k)))
         vap=(lev(l,imjm)*wvmr)/(0.62197+wvmr)
         td(l,k)=(243.5*alog(vap/6.112))/(17.67-alog(vap/6.112))
	 wvmrkg(l,k)=wvmr*1000
        enddo
       enddo

! i= 20 níveis (19 níveis de pressão mais o nível de superfície)
! j= 9 variáveis
! k= 3264 pontos da matriz (64 x 51)
!
! Agora checamos o número de níveis de cada perfil. Isto vai depender de quantos
! níveis estão "abaixo" do solo em cada ponto de grade. A sondagem vai sempre
! começar do nível de superfície, não do nível de 1000hPa.
!
! Define o indice do nivel de 500 hPa  
        write(*,*)lev(1:lmp1,1)
        do l=2,lmp1
         if (lev(l,1).eq.500.) then
          l500hPa=l
         endif
        enddo 
        write(*,*)'l500hPa',l500hPa
       

       do k=1,imjm
        nlev(k)=0
        lnlev=11
        do l=l500hPa,2,-1
         if (lev(1,k) > lev(l,k)) then
          nlev(k)=lnlev
          istart(k)=l
         endif
         lnlev=lnlev+1
        enddo 
       enddo
       geop=0
       umrl=0
!      Definindo po arquivo de saída contendo as sondagens para serem
       open (3, file='eta.sound.txt', FORM='FORMATTED',status='UNKNOWN')
!
       station(1)=75000
       write(3,299)im,jm,loni,lati,res,FreqModOut,lmp1,strm_opt
       ind=index(ARQOUT,' ')-1
       write(3,*)ARQOUT(1:ind)
       do k=1,imjm   ! Loop pelo número de pontos de grade.
        write(3,300) station(k),nlev(k),ihora,idia,imes,iano
        if (nlev(k) == 0) then
         print*,' Ponto sem sondagem: ', k
        else
         geop=nint(zgeo(1,k))
         umrl=nint(rh(1,k))
!        Escreve o primeiro nível do ponto de grade k (sempre o nível de superfície).
         write(3,301) lev(1,k),geop,temp(1,k),td(1,k),  &
         umrl,wvmrkg(1,k),uvel(1,k),vvel(1,k),theta(1,k)
!        Loop pelo número de níveis do ponto de grade k.
         do i=istart(k), lmp1
          geop=nint(zgeo(i,k))
          umrl=nint(rh(i,k))
          write(3,301) lev(i,k),geop,temp(i,k),td(i,k),  &
          umrl,wvmrkg(i,k),uvel(i,k),vvel(i,k),theta(i,k)
         enddo
         station(k+1)=station(k)+1
        endif
       enddo
299  format (2(1x,I3),3(1x,f5.1),1x,I2,1x,I3,1x,I2)
300  format (1X,I6,1x,I4,3(1x,I2),1X,I4)
301  format (1X,F7.1,I7,2F7.1,I7,F7.2,2F8.2,F9.3)
305  format (1X,I3)

     stop
!    END OF PROGRAM GERA SONDAGENS
     end
!===============================================================================
!
      subroutine read_degrib(etafile,imjm,lmp1,lev,zgeo,temp,td,uvel,vvel,rh,umes,atime)
!
      implicit none
!
      integer  :: imjm,lmp1,i,datsav,n,l
!
      real   ::  lev(lmp1,imjm)
      real   ::  zgeo(lmp1,imjm),temp(lmp1,imjm),td(lmp1,imjm)
      real   ::  uvel(lmp1,imjm),vvel(lmp1,imjm),rh(lmp1,imjm),umes(lmp1,imjm)
      real   ::  tmp(imjm,1)
!
      integer :: kgds(200),kpds(50),len,kerr
      integer :: jpds(200),jgds(200) 
      integer :: lenpds,lenkgds,nwords,kpdsl 
      integer :: j,k,nx,IRETO,KNUM,IRET1
!
      character*(*) etafile
      character(LEN=1):: pds(50)
      character(LEN=12):: atime
      character(LEN=2):: gproj
      logical bitmap(imjm)
      common /gdsinfo/datsav(11)


!_______________________________________________________________________________
!
      len=index(etafile//' ',' ')-1
!
! *** Check that the input dimensions match grib file dimensions. 
!
	write(6,*) 'inside read_degrib '

!  Should be able to derive all header info in the DEGRIB file
! from the GDS/PDS of the GRIB file.  
!
!
!
! *** Fill time stamp (assume all records are for same time).

	call get_gds(etafile(1:len),datsav,kpds)

	if (kpds(8).ge.100) kpds(8)=kpds(8)-100
!      write(atime,'(i2.2,i2.2,i2.2,i2.2,i4.4)') kpds(8),kpds(9),kpds(10),kpds(11),kpds(14)
!
!
! *** Fill a missing value flag into first space of each variable.
!
	   zgeo=-9999.
!	   
!
! *** Now put the data into the corresponding arrays.
!
!Cmp
!  add something to read in surface fields in here
!

	call baopen(11,etafile,IRETO)
	if (IRETO .ne. 0) write(6,*) 'BAOPEN TROUBLE!!!! ', IRETO

	jpds=-1
	jgds=-1

	jpds(5)=132
	jpds(6)=1
	jpds(7)=0

	call getgb(11,0,imjm,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS, &
     &     BITMAP,tmp(1,1),IRET1)
        zgeo(1,1:imjm)=tmp(1:imjm,1) 
	write(6,*) 'first GETGB, IRET1= ', IRET1
	write(6,*) 'Topo READ!!!!! '
	write(6,*) (zgeo(1,j),j=nwords/2,nwords/2+5)

	jpds=-1
	jgds=-1

	jpds(5)=135
	jpds(6)=1
	jpds(7)=0

	call getgb(11,0,imjm,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS, &
     &     BITMAP,tmp(1,1),IRET1)
        lev(1,1:imjm)=tmp(1:imjm,1) 
	write(6,*) 'first GETGB, IRET1= ', IRET1
	write(6,*) 'PSLC READ!!!!! '
	write(6,*) (lev(1,j),j=nwords/2,nwords/2+5)

	jpds=-1
	jgds=-1

	jpds(5)=128
	jpds(6)=1
	jpds(7)=0

	call getgb(11,0,imjm,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS, &
     &     BITMAP,tmp(1,1),IRET1)
        temp(1,1:imjm)=tmp(1:imjm,1) 
	write(6,*) 'first GETGB, IRET1= ', IRET1
	write(6,*) 'TMP2m READ!!!!! '
	write(6,*) (temp(1,j),j=nwords/2,nwords/2+5)

	jpds=-1
	jgds=-1

	jpds(5)=129
	jpds(6)=1
	jpds(7)=0

	call getgb(11,0,imjm,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS, &
     &     BITMAP,tmp(1,1),IRET1)
        td(1,1:imjm)=tmp(1:imjm,1) 
	write(6,*) 'first GETGB, IRET1= ', IRET1
	write(6,*) 'DP2m READ!!!!! '
	write(6,*) (td(1,j),j=nwords/2,nwords/2+5)


	jpds(5)=130
	jpds(6)=1
	jpds(7)=0

	call getgb(11,0,imjm,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS, &
     &     BITMAP,tmp(1,1),IRET1)
        uvel(1,1:imjm)=tmp(1:imjm,1) 
	write(6,*) 'first GETGB, IRET1= ', IRET1
	write(6,*) 'U10m READ!!!!! '
	write(6,*) (uvel(1,j),j=nwords/2,nwords/2+5)

	jpds=-1
	jgds=-1

	jpds(5)=131
	jpds(6)=1
	jpds(7)=0

	call getgb(11,0,imjm,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS, &
     &     BITMAP,tmp(1,1),IRET1)
        vvel(1,1:imjm)=tmp(1:imjm,1) 
	write(6,*) 'first GETGB, IRET1= ', IRET1
	write(6,*) 'V10mm READ!!!!! '
	write(6,*) (vvel(1,j),j=nwords/2,nwords/2+5)

	jpds(6)=100
       	do l=2,lmp1
	 jpds(7)=nint(lev(l,1))
	 jpds(5)=7
         write(6,*) 'LEVEL ', jpds(7), jpds(5)
         call getgb(11,0,imjm,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS, &
     &     BITMAP,tmp(1,1),IRET1)
         if (IRET1 .ne. 0) write(6,*) 'ZGEO AT LEVEL ', jpds(7) , jpds(5)
         zgeo(l,1:imjm)=tmp(1:imjm,1)
	 write(6,*) (zgeo(l,j),j=nwords/2,nwords/2+5)

  	 jpds(5)=11
         call getgb(11,0,imjm,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS, &
     &     BITMAP,tmp(1,1),IRET1)
         if (IRET1 .ne. 0) write(6,*) 'TEMP AT LEVEL ', jpds(7) , jpds(5)
         temp(l,1:imjm)=tmp(1:imjm,1)
	 write(6,*) (temp(l,j),j=nwords/2,nwords/2+5)
	 
	 jpds(5)=33
         call getgb(11,0,imjm,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS, &
     &     BITMAP,tmp(1,1),IRET1)
         if (IRET1 .ne. 0) write(6,*) 'UVEL AT LEVEL ', jpds(7) , jpds(5)
         uvel(l,1:imjm)=tmp(1:imjm,1)
 	 jpds(5)=34
         call getgb(11,0,imjm,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS, &
     &       BITMAP,tmp(1,1),IRET1)
         if (IRET1 .ne. 0) write(6,*) 'VVEL AT LEVEL ', jpds(7) , jpds(5)
         vvel(l,1:imjm)=tmp(1:imjm,1)

	 jpds(5)=52
         call getgb(11,0,imjm,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS, &
     &      BITMAP,tmp(1,1),IRET1)
         if (IRET1 .ne. 0) write(6,*) ' AT LEVEL ', jpds(7) , jpds(5)
         rh(l,1:imjm)=tmp(1:imjm,1)
	 jpds(5)=51
         call getgb(11,0,imjm,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS, &
     &     BITMAP,tmp(1,1),IRET1)
         if (IRET1 .ne. 0) write(6,*) ' AT LEVEL ', jpds(7) , jpds(5)
         umes(l,1:imjm)=tmp(1:imjm,1)
	enddo

!
! *** Normal finish.
!
1000  continue

      return
!
! *** Premature end of file.
!
1100  continue
      print *,'Premature end of file.'
      print *,'Abort...'
      stop
!
      end
!
!===============================================================================

       subroutine get_gds(etafile,gdsinfo,kpds)

      character*(*) etafile
      character(LEN=1):: pds(50)

      integer kgds(200),kpds(200),len,kerr &
     &         ,lenpds,lenkgds,nwords,kpdsl &
     &         ,j,k,gdsinfo(11) &
     &         ,gdsav,IRETO,JGDS(200),JPDS(200)
      real tmp(100000)
      logical bitmap(100000)


	nxny=100000


	JPDS=-1
	JGDS=-1


        len=index(etafile//' ',' ')-1

	call baopen(11,etafile(1:len),IRETO)
	write(6,*) 'BAOPEN in get_gds: ', IRETO

        if (IRETO .ne. 0) then
         print *,'Error opening unit=11, file name = ',etafile(1:len) &
     &          ,' iostat = ',kerr
         stop
        endif

	jpds(5)=11
	jpds(6)=100
	jpds(7)=500
	call getgb(11,0,nxny,0,JPDS,JGDS,nwords,KNUM,KPDS,KGDS, &
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

