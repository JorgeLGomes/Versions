
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
       integer              :: i,j,k,l,irec,geop,lnlev,umrl,strm_opt,ind
       real                 :: est, estd,wvmr2m,wvmr,vap
       real                 :: mlevs(Max_Sond_Levs)

       integer, allocatable   :: station(:),istart(:),nlev(:)
! Variáveis de entrada
       character(len=60) 	   :: ARQIN,ARQOUT
       real, allocatable  :: variab(:,:,:)
       real, allocatable  :: lev(:)
       real, allocatable  :: umes(:,:)

       namelist/model_grids/im,jm,lm,loni,lati,res,FreqModOut,        &
     &                     idia,imes,iano,ihora,ARQIN,ARQOUT,strm_opt 
                       
       namelist/model_levs/mlevs
 
       open(1,file='namelist',form='formatted',status='old')
       rewind(1)
       read(1,model_grids)
       rewind(1)
       read(1,model_levs)
       close(1)

       imjm=im*jm
       lmp1=lm+1
       allocate(lev(lmp1))
       
       write(*,*)idia,imes,iano,ihora &
      &                    ,ARQIN,ARQOUT,strm_opt
       write(*,*)im,jm,lm,loni,lati,res,FreqModOut,lmp1,imjm
       lev(1)=9999.
       lev(2:lmp1)=mlevs(1:lm)
       write(*,*)lev
       ALLOCATE(station(imjm),istart(imjm),nlev(imjm))
       ALLOCATE(variab(lmp1,9,imjm))
       ALLOCATE(umes(lm,imjm))
       variab=-9999.
  
!       write(*,*) 'Utilizar deslocamento observado da tempestade? (1=sim 5=não)'
!       write(*,*) 'strm_opt ',strm_opt
       if ((strm_opt /= 1).AND.(strm_opt /= 5)) then
        write(*,*) 'Opção não disponível. Programa interrompido.'
        stop
       endif
       open (10,file=ARQIN,form='unformatted', & 
      &      access='direct',recl=imjm*4,status='old')
 
       read(10,rec=3) variab(1,1,:) !topo 
       write(31,*) variab(1,1,1:imjm) 
       read(10,rec=2) variab(1,2,:) !pslc
       read(10,rec=5) variab(1,3,:) !tp2m
       write(30,*) variab(1,3,1:imjm) 
       read(10,rec=6) variab(1,4,:) !dp2m
       read(10,rec=7) variab(1,7,:) !u10m
       read(10,rec=8) variab(1,8,:) !v10m
       irec=51
       do l=2,lmp1
        read(10,rec=irec) variab(l,1,:)  ! zgeo
        irec=irec+1
       enddo
       irec=111
       do l=2,lmp1
        read(10,rec=irec) variab(l,3,:) !temp(l,:)
        irec=irec+1
       enddo 

       irec=131
       do l=2,lmp1
        read(10,rec=irec) variab(l,5,:) !umrl(:)
        irec=irec+1
       enddo
       irec=71
       do l=2,lmp1
        read(10,rec=irec) variab(l,7,:) !uvel(l,:)
        irec=irec+1
       enddo
       irec=91
       do l=2, lmp1
        read(10,rec=irec) variab(l,8,:) !vvel(l,:)
        irec=irec+1
       enddo
       irec=171
       do l=2,lmp1
        read(10,rec=irec) umes(l,:)
        irec=irec+1
       enddo 

!    Tranforma Temp (K) => Temp (C)
       variab(1:lmp1,3,1:imjm)=variab(1:lmp1,3,1:imjm)-T0
       variab(1:lmp1,4,1:imjm)=variab(1:lmp1,4,1:imjm)-T0

       do k=1, imjm
        variab(1,9,k)=exp(alog(variab(1,3,k)+T0)+0.286*alog(1000/variab(1,2,k)))
        est=6.112*exp((17.67*variab(1,3,k))/(variab(1,3,k)+243.5))
        estd=6.112*exp((17.67*variab(1,4,k))/(variab(1,4,k)+243.5))
        variab(1,5,k)=100*(estd/est)
        wvmr2m=(0.62197*estd)/(variab(1,2,k)-estd)
        variab(1,6,k)=wvmr2m*1000
       enddo
       do k=1, imjm
        do l=2, lmp1
	 variab(l,2,k)=lev(l)
         variab(l,9,k)=exp(alog(variab(l,3,k)+T0)+0.286*alog(1000/lev(l)))
         wvmr=(umes(l,k)/(1-umes(l,k)))
         vap=(lev(l)*wvmr)/(0.62197+wvmr)
         variab(l,4,k)=(243.5*alog(vap/6.112))/(17.67-alog(vap/6.112))
	 variab(l,6,k)=wvmr*1000
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
       do k=1,imjm
        nlev(k)=0
        lnlev=11
        do l=12,2,-1
         if (variab(1,2,k) > lev(l)) then
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
         geop=nint(variab(1,1,k))
         umrl=nint(variab(1,5,k))
!        Escreve o primeiro nível do ponto de grade k (sempre o nível de superfície).
         write(3,301) variab(1,2,k),geop,variab(1,3,k),variab(1,4,k),  &
         umrl,variab(1,6,k),variab(1,7,k),variab(1,8,k),variab(1,9,k)
!        Loop pelo número de níveis do ponto de grade k.
         do i=istart(k), lmp1
          geop=nint(variab(i,1,k))
          umrl=nint(variab(i,5,k))
          write(3,301) lev(i),geop,variab(i,3,k),variab(i,4,k),  &
          umrl,variab(i,6,k),variab(i,7,k),variab(i,8,k),variab(i,9,k)
         enddo
         station(k+1)=station(k)+1
        endif
       enddo
299  format (2(1x,I3),3(1x,f5.1),1x,I2,1x,I3,1x,I2)
300  format (1X,I5,1x,I4,3(1x,I2),1X,I4)
301  format (1X,F7.1,I7,2F7.1,I7,F7.2,2F8.2,F9.3)
305  format (1X,I3)

     stop
!    END OF PROGRAM GERA SONDAGENS
     end
