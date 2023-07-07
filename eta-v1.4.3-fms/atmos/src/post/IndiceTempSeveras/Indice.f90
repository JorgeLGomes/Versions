
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
! Variáveis de trabalho
  integer i,j,k,irec,geop,umrl,station,istart,nlev
  dimension station(3270),istart(3270),nlev(3270)
! Variáveis de entrada
  real :: variab,lev
  dimension  variab(20,9,3270),lev(20)
!
! Início do programa executável...
!
! Lendo o arquivo de entrada contendo as saída binária do GrADS:
! O arquivo binário gerado pelo SOUND3.EXEC para o caso Indaiatuba tem 9 variáveis
! à supercície mais 9 variáveis em 19 níveis de pressão em matrizes de 64 X 51 pontos.
!
  open (2, file='sonds_2400_2418.le.bin',FORM='UNFORMATTED',status=   &
        'UNKNOWN',ACCESS='DIRECT',recl=3264*9*20*4)
 irec=1
 read (2,rec=irec) (((variab(i,j,k),k=1,3264),j=1,9),i=1,20)
! i= 20 níveis (19 níveis de pressão mais o nível de superfície)
! j= 9 variáveis
! k= 3264 pontos da matriz (64 x 51)
!
! Agora checamos o número de níveis de cada perfil. Isto vai depender de quantos
! níveis estão "abaixo" do solo em cada ponto de grade. A sondagem vai sempre
! começar do nível de superfície, não do nível de 1000hPa.
!
 do k=1,3264
  if (variab(1,2,k) > 1000.0) then
   nlev(k)=20
   istart(k)=2
  else
   if (variab(1,2,k) > 925.0) then
    nlev(k)=19
    istart(k)=3
   else
    if (variab(1,2,k) > 900.0) then
     nlev(k)=18
     istart(k)=4
    else
     if (variab(1,2,k) > 850.0) then
      nlev(k)=17
      istart(k)=5
     else
      if (variab(1,2,k) > 800.0) then
       nlev(k)=16
       istart(k)=6
      else
       if (variab(1,2,k) > 750.0) then
        nlev(k)=15
        istart(k)=7
       else
        if (variab(1,2,k) > 700.0) then
         nlev(k)=14
         istart(k)=8
        else
         if (variab(1,2,k) > 650.0) then
          nlev(k)=13
          istart(k)=9
         else
          if (variab(1,2,k) > 600.0) then
           nlev(k)=12
           istart(k)=10
          else
           if (variab(1,2,k) > 550.0) then
            nlev(k)=11
            istart(k)=11
           else
            nlev(k)=0
           endif
          endif
         endif
        endif
       endif
      endif
     endif
    endif
   endif
  endif
 enddo
 lev(2)=1000.0
 lev(3)=925.0
 lev(4)=900.0
 do i=5,20
  lev(i)=lev(i-1)-50.0
 enddo
 geop=0
 umrl=0
! Definindo po arquivo de saída contendo as sondagens para serem
  open (3, file='eta.sound_2400_2418.txt', FORM='FORMATTED',status='UNKNOWN')
! A LINHA ACIMA TEM QUE SER EDITADA TODA VEZ QUE MUDARMOS O ARQUIVO DE ENTRADA.
!
  station(1)=75000
  write(3,*)' 3264'
  do k=1,3264   ! Loop pelo número de pontos de grade.
   write(3,*) station(k),' ',nlev(k),' 18 24 05 2005'
! A LINHA ACIMA TEM QUE SER EDITADA TODA VEZ QUE MUDARMOS O ARQUIVO DE SONDAGEM.
   if (nlev(k) == 0) then
    print*,' Ponto sem sondagem: ', k
   else
    geop=nint(variab(1,1,k))
    umrl=nint(variab(1,5,k))
! Escreve o primeiro nível do ponto de grade k (sempre o nível de superfície).
    write(3,301) variab(1,2,k),geop,variab(1,3,k),variab(1,4,k),  &
           umrl,variab(1,6,k),variab(1,7,k),variab(1,8,k),variab(1,9,k)
! Loop pelo número de níveis do ponto de grade k.
    do i=istart(k),20
     geop=nint(variab(i,1,k))
     umrl=nint(variab(i,5,k))
     write(3,301) lev(i),geop,variab(i,3,k),variab(i,4,k),  &
           umrl,variab(i,6,k),variab(i,7,k),variab(i,8,k),variab(i,9,k)
    enddo
    station(k+1)=station(k)+1
   endif
  enddo
!300  format (1X,I5,3I4,I5)
301  format (1X,F7.1,I7,2F7.1,I7,F7.2,2F8.2,F9.3)
305  format (1X,I3)
  stop
! END OF PROGRAM GERA SONDAGENS
  end



! --------------------------------------------------------------------------
!  PROGRAMA índices_severos5
!
!           Autor: ERNANI DE LIMA NASCIMENTO
!           Data: Outubro de 2006
!
!  Este programa calcula índices de tempo severo a partir de perfis
!  verticais atmosféricos gerados por modelos numéricos.
!  NOTA 1: ESTA VERSÃO CALCULA ÍNDICES SEVEROS TERMODINÂMICOS E CINEMÁTICOS.
!  NOTA 2: ESTA VERSÃO É IDÊNTICA À VERSÃO índices_severos4, EXCETO QUE LÊ
!          UM ARQUIVO DE ENTRADA COM CONTEÚDO DIFERENTE E GERA UM ARQUIVO
!          BINÁRIO PARA O GRADS, ALÉM DO ARQUIVO ASCII PARA O ARPS.
!  NOTA 3: ESTA VERSÃO GERA SAÍDA DE SONDAGENS EM FORMATO PARA EXECUÇÃO NO
!          ARPS. POR EXEMPLO, LÊ O ARQUIVO eta.sound_0212_0318.txt E GERA
!          O ARQUIVO arps.eta.sound_0212_0318.txt. DEPOIS ESTE ARQUIVO
!          DEVE SER EDITADO PARA INCLUIR O CABEÇALHO DESEJADO E ENTÃO
!          SER SALVO COM EXTENSÃO .snd PARA SER UTILIZADO COM O ARPS
!          (LEIA O ARQUIVO LEIA-ME PARA MAIORES DETALHES).
!
!
!  Exemplo de arquivo de entrada necessário para este programa (este formato
!  é aquele gerado pelo programa gera_sond_2.f90):
!
! 78122          20 18 24 05 2005
!  1000.0    125   25.4   19.7     70  14.65   -1.24   -7.64  298.581
!   925.0    804   20.3   16.4     79  12.83   -1.26   -7.91  300.023
!   900.0   1040   18.3   15.0     81  11.97   -1.31   -8.07  300.422
!   850.0   1528   15.4   11.6     78  10.18   -2.30  -10.65  302.267
!   800.0   2041   13.0    7.5     69   8.17   -2.44  -10.54  305.057
!   750.0   2581    9.8    3.8     67   6.72   -1.68   -7.88  307.192
!   700.0   3151    6.3    0.2     65   5.54   -0.90   -5.52  309.413
!   650.0   3757    3.1   -4.2     59   4.30   -0.38   -4.89  312.487
!
!   ...    ....    ...   ...     ...   ...    ...     ...   ....   ....  ....
!  etc.....
!
! O que significa cada número no exemplo acima:
! 78122 => identificador da "estação" (irrelevante neste programa)
! 20 => número de níveis do perfil
! 18 => horário UTC
! 24 => dia
! 05 => mês
! 2005 => ano
! Colunas: pressão (hPa), altitude (m), temperatura (C), temperatura do ponto
! de orvalho (C), umidade relativa (%), razão de mistura (g/kg), componente zonal
! do vento (m/s), componente meridional do vento (m/s), temperatura potencial (K).
!
!
! Índices calculados: CAPE de superfície, CAPE com correção de water loading,
! indice de levantamento, nível de condensação por levantamento, nível de convecção
! espontânea, inibição convectiva, componentes u e v do movimento esperado de tempestades,
! escoamento relativo à tempestade em níveis baixos e médios, helicidade relativa à
! tempestade, número de Richardson volumétrico, denominador do número de Richardson
! volumétrico, índice de energia-helicidade.
!
! Com exceção da rotina de cálculo do índice de energia-helicidade (interiramente nova),
! todas as rotinas foram modificadas a partir do ARPS5.0_IHOP6 e ADAS.
!
!
! ---------------------------------------------------------------------------
  program indices_severos
! Variáveis miscelâneas
  integer i,j,k
! Variáveis de entrada
  integer :: nrad,estacao,hora,dia,mes,ano
  integer :: altm,urporc,dddeg,vvkt,nlev,strm_opt
  dimension :: estacao(3270),hora(3270),dia(3270),mes(3270),ano(3270),  &
               nlev(3270)
  dimension  altm(102,3270),urporc(102,3270)
!  dddeg(100,3270),vvkt(100,3270) vvmps(100,3270)
  real :: presshpa,tc,tdc,wvgkg,thetae,thetav,presspa,tk,wvkgkg, &
          u,v,arg,u2dd,v2dd,theta,specif
  dimension  presshpa(102,3270),tc(102,3270),tdc(102,3270),wvgkg(102,3270), &
             presspa(102,3270),tk(102,3270),wvkgkg(102,3270),u(102,3270),   &
             v(102,3270),theta(102,3270),specif(102,3270)
! Variáveis para os índices de tempo severo
  real :: cape,mcape,li,ncl,nce,cin,ustrm,vstrm,llsrm,mlsrm,heli,nrv,  &
          dnrv,ieh,sup,tcap,thetab,bl,blcon,shr37,ffstrm,ddstrm,dnrv2km
  dimension  cape(3270),mcape(3270),li(3270),ncl(3270),nce(3270),cin(3270),el(3270), &
             ustrm(3270),vstrm(3270),llsrm(3270),mlsrm(3270),heli(3270),nrv(3270), &
             dnrv(3270),ieh(3270),tcap(3270),thetab(3270),bl(3270),shr37(3270),    &
             blcon(3270),sup(3270),ffstrm(3270),ddstrm(3270),dnrv2km(3270)
! Obs: tcap = lid strength; thetab = twdf on arpsplt.f (max wet bulb potential
! temperature difference)
!
! Variáveis para cálculos em coluna
  real :: p1d,ht1d,t1d,tv1d,td1d,wmr1d,partem,buoy,wload,mbuoy,pbesnd, &
          mbesnd
  dimension  p1d(102),ht1d(102),t1d(102),tv1d(102),td1d(102),          &
             wmr1d(102),partem(102),buoy(102),wload(102),mbuoy(102),   &
             pbesnd(102),mbesnd(102)
! Variável para nome do arquivo de entrada
  character(60) inname,outname1,outname2
!
! Início do programa executável...
!
! Lendo o arquivo de entrada contendo as sondagens:
  write(*,*) 'Este programa calcula índices de tempo severo.'
  write(*,*) 'Ernani L. Nascimento (OUT/2006) (modificado do ARPS). '
  write(*,*) 'Entre com o nome do arquivo de entrada (até 60 caracteres):'
  read(*,*) inname
!  write(*,*) 'Utilizar deslocamento observado da tempestade? (1=sim 5=não)'
  read(*,*) strm_opt
  if ((strm_opt /= 1).AND.(strm_opt /= 5)) then
   write(*,*) 'Opção não disponível. Programa interrompido.'
   stop
  endif
  outname1='grads.'//inname(1:42)
  outname2='arps.'//inname(1:42)
  open (10,file=inname,status='old')
  read(10,*) nrad
  do k=1,nrad
   read (10,*) estacao(k),nlev(k),hora(k),dia(k),mes(k),ano(k)
   if (nlev(k) /= 0) then
    do i=1,nlev(k)
     read (10,*) presshpa(i,k),altm(i,k),tc(i,k),tdc(i,k),urporc(i,k),   &
                wvgkg(i,k),u(i,k),v(i,k),theta(i,k)
     presspa(i,k)=presshpa(i,k)*100.
     tk(i,k)=tc(i,k)+273.15
     specif(i,k)=(wvgkg(i,k)/1000.)/(1+(wvgkg(i,k)/1000.))
!      print *,'specif(1,k)= ',specif(1,k), k
     arg=0.0
    end do
!     write(*,*) k,' ',presspa(1,k)

! As linhas abaixo geram 2 níveis adicionais na sondagem. O ETA possui 19 níveis, mas
! para as rodadas do ARPS tornou-se conveniente acrescentar dois níveis a mais para
! garantir um topo de domínio mais alto. Estas linhas podem ser convenientemente comentadas.
    altm(nlev(k)+1,k)=altm(nlev(k),k)+(altm(nlev(k),k)-altm(nlev(k)-1,k))
    altm(nlev(k)+2,k)=altm(nlev(k)+1,k)+(altm(nlev(k),k)-altm(nlev(k)-1,k))
    specif(nlev(k)+1,k)=specif(nlev(k),k)
    specif(nlev(k)+2,k)=specif(nlev(k),k)
    u(nlev(k)+1,k)=u(nlev(k),k)
    u(nlev(k)+2,k)=u(nlev(k),k)
    v(nlev(k)+1,k)=v(nlev(k),k)
    v(nlev(k)+2,k)=v(nlev(k),k)
   endif
!   print *,'presshpa(19,k)= ',presshpa(19,k),' ', k
  end do
  close(10)
!
!     if (j==1) then
!      write(*,*) 'dddeg(1,k)= ',dddeg(j,k),' para k= ',k
!      write(*,*) 'vvmps(1,k)= ',vvmps(j,k),' para k= ',k
!     endif
!     call ddff2uv(dddeg,vvmps,uu,vv,j)
!
! Transformando graus e magnitude do vento em componentes zonal e meridional
!     arg = (dddeg(j,k) * (3.141592654/180.))
!     uu = -vvmps(j,k) * sin(arg)
!     vv = -vvmps(j,k) * cos(arg)
!     u(j,k)=uu
!     v(j,k)=vv
!     if (j==1) then
!      write(*,*) 'u(1,k)= ',u(j,k),' para k= ',k
!      write(*,*) 'v(1,k)= ',v(j,k),' para k= ',k
!     endif

!
!   Chama rotina de cálculo de índices termodinâmicos
!
!   print *,'hgt(m) antes: ', altm(1,1)
!   do k=1,nrad
!     print *,'presshpa(19,k)= ',presshpa(19,k),' ', k
!   enddo
!   write(*,*) k,' PASSEI AQUI 1.'
   call arps_be(nrad,presshpa,altm,tc,tdc,wvgkg,ncl,nce,el,thetab,li,cape,   &
               mcape,cin,tcap,p1d,ht1d,t1d,tv1d,td1d,wmr1d,partem,          &
               buoy,wload,mbuoy,pbesnd,mbesnd,nlev,estacao,hora,dia,        &
               mes,ano)
!
!  Cálculo de índices cinemáticos
!
!
!   write(*,*) k,' PASSEI AQUI 2'
   call calcshr(nrad,altm,bl,presspa,tk,u,v,cape,shr37,ustrm,vstrm,llsrm,  &
               mlsrm,heli,nrv,dnrv,dnrv2km,blcon,u1,v1,nlev,estacao,hora,dia,  &
               mes,ano,strm_opt)
   call calcehi_sup(nrad,nlev,presshpa,cape,heli,dnrv,ieh,sup)

!  Arquivo de saída: binário para o grads
  open(12,file=outname1,FORM='UNFORMATTED',status='unknown',ACCESS='DIRECT',   &
       recl=3264*17*4)
!   write(12,*) ' EST UTC DIA MÊS ANO SHR37 ustrm vstrm ddstrm ffstrm helic BRN BRNSHR B2km IEH SUP'
  irec=1
  do k=1,nrad
   u2dd=ustrm(k)
   v2dd=vstrm(k)
   call uv2ddff(u2dd,v2dd,dd,ff)
   ddstrm(k)=dd
   ffstrm(k)=ff
   write(12,rec=irec) ncl(k),nce(k),el(k),thetab(k),li(k),cape(k),mcape(k),cin(k),shr37(k),       &
                ustrm(k),vstrm(k),heli(k),nrv(k),dnrv(k),dnrv2km(k),ieh(k),sup(k)
!    if ((ieh(k)<-1.0).and.(sup(k)<-1.0)) then
!      write(*,*) 'EHI AND SUP LESS THAN -1 (k EHI SUP): ',k,ieh(k),sup(k)
!    endif
  enddo
  close(12)

!  Arquivo de saída: sondagem com formato para o ARPS
!   open(13,file=outname2,status='unknown')
!   write(13,304) nrad
!   do k=1,nrad
!    write(13,300) estacao(k),hora(k),dia(k),mes(k),ano(k)
!    write (13,*) ' ZSND   THSND   QVSND    USND    VSND '
!  As duas linhas abaixo acrescentam dois níveis para theta.
!    theta(nlev(k)+1,k)=theta(nlev(k),k)+(theta(nlev(k),k)-theta(nlev(k)-1,k))
!    theta(nlev(k)+2,k)=theta(nlev(k)+1,k)+(theta(nlev(k),k)-theta(nlev(k)-1,k))
!  O loop abaixo acrescenta os dois níveis adicionais para o ARPS.
!    do i=1,nlev(k)+2
!     write(13,303) float(altm(nlev(k)+3-i,k)),theta(nlev(k)+3-i,k),specif(nlev(k)+3-i,k),u(nlev(k)+3-i,k), &
!                  v(nlev(k)+3-i,k)
!    enddo
!   enddo
!   close(13)
300  format (1X,I5,3I4,I5)
301  format (1X,I5,3I4,I5,F8.1,1X,2F8.1,2F7.1,3F8.1)
302  format (1X,I5,3I4,I5,F8.5,2F7.3,F8.1,F7.1,F9.2,F6.1,2F7.1,2F6.1)
303  format (1X,F9.3,F12.5,F12.8,2F12.5)
304  format (1X,I4)
  stop
! END OF PROGRAM ÍNDICES SEVEROS
  end
!
!
! BELOW: FUNCTIONS AND SUBROUTINES USED BY PROGRAM ÍNDICES SEVEROS
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ARPS_BE                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######    and modified by Ernani L. Nascimento at Simepar   ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE arps_be(nrad,pres,hgt,tc,tdc,wmr,lcl,lfc,el,twdf,li,cape,mcape,   &
          cin,tcap,p1d,ht1d,t1d,tv1d,td1d,wmr1d,partem,buoy,wload,      &
          mbuoy,pbesnd,mbesnd,nlev,estacao,hora,dia,mes,ano)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the lifting condensation level (lcl), level of free
!  convection (lfc), equilibrium level (el), max wet-bulb potential
!  temperature difference (twdf), lifted index (LI), Convective
!  Available Potential Energy (CAPE), Moist Convective Potential
!  Energy (MCAPE, includes water loading), convective inhibition
!  (CIN) and lid strength (tcap)  over the ARPS grid.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  (Keith Brewster)
!  Cleaned-up, removed OLAPS artifacts.
!
!  FEB 13-14 2004 (ERNANI L. NASCIMENTO)
!  Modified from ARPS to compute indices from rawinsonde data.
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
 INTEGER :: nrad
!
!-----------------------------------------------------------------------
!
!  "1-D" input variables (Ernani L. Nascimento)
!
!-----------------------------------------------------------------------
!
 INTEGER :: nlev(nrad),estacao(nrad),hora(nrad),dia(nrad),mes(nrad),   &
            ano(nrad)
!
!-----------------------------------------------------------------------
!
!  "2-D" input variables (Ernani L. Nascimento)
!
!-----------------------------------------------------------------------
!
 REAL :: pres(102,nrad)
! hgt is integer! Different from the original code: Ernani L. Nascimento
 INTEGER :: hgt(102,nrad)
 REAL :: tc(102,nrad)
 REAL :: wmr(102,nrad)
 REAL :: tdc(102,nrad)
!
!-----------------------------------------------------------------------
!
!  Output variables ("1-D") (Ernani L. Nascimento)
!
!-----------------------------------------------------------------------
!
 REAL :: lcl(nrad)
 REAL :: lfc(nrad)
 REAL :: el(nrad)
 REAL :: twdf(nrad)
 REAL :: li(nrad)
 REAL :: cape(nrad)
 REAL :: mcape(nrad)
 REAL :: cin(nrad)
 REAL :: tcap(nrad)
!
!-----------------------------------------------------------------------
!
!  Scratch space for calculations
!
!-----------------------------------------------------------------------
!
 REAL :: p1d(102),ht1d(102),tv1d(102),td1d(102)
 REAL :: partem(102),buoy(102),wload(102)
 REAL :: mbuoy(102),pbesnd(102),mbesnd(102)
 REAL :: wmr1d(102),t1d(102)
!
!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------
!
 REAL :: wmr2td,tctotv
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
 INTEGER :: i,j,k,n,nlevel,bla
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!  Loop over all soundings (Ernani L. Nascimento)
!-----------------------------------------------------------------------
!  print *, 'hgt(m) depois: ',hgt(1,1)
 DO k=1,nrad
  if (nlev(k)==0) then
   lcl(k)=-999.9
   lfc(k)=-999.9
   el(k)=-999.9
   twdf(k)=-999.9
   li(k)=-999.9
   cape(k)=-999.9
   mcape(k)=-999.9
   cin(k)=-999.9
   tcap(k)=-999.9
  else
   j=nlev(k)
   if (pres(j,k)>300.0) then
    write(*,*) 'Sondagem abaixo não chegou aos 300hPa (estacao,hora,dia,mes,ano):'
    write(*,*) estacao(k),hora(k),dia(k),mes(k),ano(k)
    write(*,*) pres(j,k),nlev(k)
    write(*,*) 'Buscando próxima sondagem.'
    cycle
   endif
!-----------------------------------------------------------------------
! Loop over all levels of the sounding
!-----------------------------------------------------------------------
   DO i=1,nlev(k)
     p1d(i) = pres(i,k)
     ht1d(i) = float(hgt(i,k))
     wmr1d(i) = wmr(i,k)
     t1d(i) = tc(i,k)
!      td1d(i) = wmr2td(p1d(i),wmr1d(i))
     td1d(i) = tdc(i,k)
     tv1d(i) = tctotv(t1d(i),wmr1d(i))
   END DO
!
!    print *, 'ht1d(1): ', ht1d(1)
   nlevel = nlev(k)
!
   CALL sindex(nlevel,p1d,ht1d,t1d,tv1d,td1d,wmr1d,partem,buoy,wload,    &
               mbuoy,pbesnd,mbesnd,lcl(k),lfc(k),el(k),twdf(k),li(k),    &
               cape(k),mcape(k),cin(k),tcap(k))
  endif
 END DO
RETURN
END SUBROUTINE arps_be
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SINDEX                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######    and modified by Ernank L. Nascimento at Simepar   ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE sindex(nlevel,p,ht,t,tv,td,w,                                &
          partem,buoy,wload,mbuoy,pbesnd,mbesnd,                       &
          lcl_pbe,lfc,el,twdf,li,cape,mcape,cin,tcap)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate Convective Available Potential Energy (CAPE)
!  in each column of ARPS grid.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  Cleaned-up, removed OLAPS artifacts.
!  5/13/1996  Added cap strength.
!  Feb 13 2004 ERNANI L. NASCIMENTO: removed variable maxlev
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
!
 INTEGER :: nlevel
!
 REAL :: p(nlevel+2),ht(nlevel+2),t(nlevel+2),tv(nlevel+2),td(nlevel+2),w(nlevel+2)
 REAL :: partem(nlevel+2),buoy(nlevel+2),wload(nlevel+2),mbuoy(nlevel+2)
 REAL :: pbesnd(nlevel+2),mbesnd(nlevel+2)
!
!  Returned from sindex
!
 REAL :: lfc,el,twdf,li,cape,mcape,cin,mcin,tcap
!
!  Potbe variables
!
 REAL :: plcl_pbe,tlcl_pbe,lcl_pbe,thepcl
 REAL :: velneg,mvelneg
!
!  Functions
!
 REAL :: oe
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
 thepcl=oe(t(1),td(1),p(1))
!
!  print *, ' theta-e of parcel: ',thepcl
!
 CALL ptlcl(p(1),t(1),td(1),plcl_pbe,tlcl_pbe)
!  print *, ' press and temp at LCL: ',plcl_pbe,tlcl_pbe
!
!  Find height of LCL
!

 CALL intrpr(nlevel,p,ht,plcl_pbe,lcl_pbe)
!  print *, ' NCL: ', lcl_pbe
!
!  Calculate the CAPE and such
!
 CALL potbe(nlevel,p(1),t(1),w(1),                                     &
            thepcl,plcl_pbe,tlcl_pbe,lcl_pbe,                          &
            p,ht,t,tv,td,w,                                            &
            partem,buoy,wload,mbuoy,pbesnd,mbesnd,                     &
            cin,velneg,cape,mcin,mvelneg,mcape,lfc,el,tcap)
!
!  Calculate Lifted Index
!
 CALL calcli(nlevel,thepcl,p,t,li)
!
!  Calculate max and min wet bulb potential temperature
!
 CALL thwxn(nlevel,p,ht,t,td,ht(1),twdf)
!
 RETURN
END SUBROUTINE sindex
!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE ptlcl(p,t,td,pc,tc)
!
!   this subroutine estimates the pressure pc (mb) and the temperature
!   tc (celsius) at the lifted condensation level (lcl), given the
!   initial pressure p (mb), temperature t (celsius) and dew point
!   (celsius) of the parcel.  the approximation is that lines of
!   constant potential temperature and constant mixing ratio are
!   straight on the skew t/log p chart.
!
!    baker,schlatter   17-may-1982   original version
!
!   teten's formula for saturation vapor pressure as a function of
!   pressure was used in the derivation of the formula below.  for
!   additional details, see math notes by t. schlatter dated 8 sep 81.
!   t. schlatter, noaa/erl/profs program office, boulder, colorado,
!   wrote this subroutine.
!
!   akap = (gas constant for dry air) / (specific heat at constant
!       pressure for dry air)
!   cta = difference between kelvin and celsius temperatures
!
 DATA akap,cta/0.28541,273.16/
 c1 = 4098.026/(td+237.3)**2
 c2 = 1./(akap*(t+cta))
 pc = p*EXP(c1*c2*(t-td)/(c2-c1))
 tc = t+c1*(t-td)/(c2-c1)
 RETURN
END SUBROUTINE ptlcl
!
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INTRPR                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######   and modified by Ernani L. Nascimento at Simepar    ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE intrpr(nlev,p,var,plvl,varatp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate variable "var" linearly in log-pressure (p).
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  Cleaned-up, removed OLAPS artifacts.
!
!  Feb 13 2004  Ernani L. Nascimento
!  Removed variable maxlev
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
 INTEGER :: nlev
 REAL :: p(nlev+2),var(nlev+2)
 REAL :: plvl
 REAL :: varatp
!
 INTEGER :: k
 REAL :: w1
!
 DO k=2,nlev
   IF(p(k) < plvl) EXIT
 END DO
!
 w1=ALOG(p(k)/plvl)/ALOG(p(k)/p(k-1))
 varatp = w1*var(k-1) + (1.-w1)*var(k)
!
 RETURN
END SUBROUTINE intrpr
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE POTBE                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######   and modified by Ernani L. Nascimento at Simepar    ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE potbe(nlevel,pmean,tmean,wmean,                              &
          blthte,plcl,tlcl,lcl,                                        &
          p,ht,t,tv,td,w,                                              &
          partem,buoy,wload,mbuoy,pbesnd,mbesnd,                       &
          pbeneg,velneg,pos_max,                                       &
          mbeneg,mvelneg,mpos_max,lfc,el,tcap)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate Convective Available Potential Energy (CAPE)
!  in each column of ARPS grid.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  February, 1994  Based on OLAPS, hence LAPS, version of same.
!                  from FSL, by Steve Albers 1991
!
!  MODIFICATION HISTORY:
!  3/11/1996  Cleaned-up, removed OLAPS artifacts.
!
!  Feb/14/2004  Ernani L. Nascimento: removed variable maxlev
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
!
 INTEGER :: nlevel
 REAL :: pmean,tmean,wmean,blthte,plcl,tlcl,lcl
 REAL :: p(nlevel+2),ht(nlevel+2)
 REAL :: t(nlevel+2),tv(nlevel+2),td(nlevel+2),w(nlevel+2)
 REAL :: partem(nlevel+2),buoy(nlevel+2),wload(nlevel+2),mbuoy(nlevel+2)
 REAL :: pbesnd(nlevel+2),mbesnd(nlevel+2)
 REAL :: pbeneg,velneg,pos_max
 REAL :: mbeneg,mvelneg,mpos_max,lfc,el,tcap
!
!  Parameters
!
 REAL :: g,gamma
 PARAMETER (g=9.80665,                                                 &
            gamma = .009760)   ! Dry Adiabatic Lapse Rate Deg/m
!
!  Functions
!
 REAL :: tsa_fast,tctotv
!
!  Misc internal variables
!
 INTEGER :: n,nel
 REAL :: deltah,delta_ht_dry,delta_ht_wet
 REAL :: sntlcl,buoy_lcl,wsat,partv
 REAL :: nbe_min,pbe_wet,pbe_dry,pos_area
 REAL :: wlow,htzero,adjeng
!
!-----------------------------------------------------------------------
!
!  Function f_mrsat and inline directive for Cray PVP
!
!  Note from Ernani L. Nascimento: this directive is ignored whenever
!  this code is compiled in a machine different from Cray. Thus, no need to
!  delete the piece of code !fpp$, !dir$, !*$* below.
!
!-----------------------------------------------------------------------
!
 REAL :: f_mrsat

!fpp$ expand (f_mrsat)
!dir$ inline always f_mrsat
!*$*  inline routine (f_mrsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!  Reset output variables.
!
!  These should be the same as what is assigned
!  no positive area found.  They are limited in
!  range to allow for contouring when positive
!  areas exist in some columns of a domain and
!  not in others.
!
 pbeneg=-400.
 velneg=20.
 pos_max=0.
 mbeneg=-400.
 mvelneg=20.
 mpos_max=0.
 lfc=10000.
 el=0.
 tcap=0.
!
!  Initialize parcel path arrays
!
 partem(1) = t(1)
 buoy(1) = 0.
 wload(1) = 0.
 mbuoy(1) = 0.
 pbesnd(1) = 0.
 mbesnd(1) = 0.

!  WRITE(6,810)pmean,tmean,wmean,plcl,tlcl,lcl
! 810 format(' pmean,tmean,wmean,plcl,tlcl,lcl',2F10.2,F10.5,2F10.2
!    +   ,F5.1)

 DO n=2,nlevel
   deltah = ht(n) - ht(n-1)
   IF(plcl < p(n-1))THEN ! lower level is below LCL
     IF(plcl < p(n))THEN ! upper level is below LCL
!        WRITE(6,*)' DRY CASE'
       partem(n)=partem(n-1)-gamma*deltah
       partv=tctotv(partem(n),w(1))
       buoy(n)=(partv-tv(n))/tv(n)
       pbesnd(n)=pbesnd(n-1)+g*0.5*(buoy(n)+buoy(n-1))*deltah
       wload(n)=0.
       mbuoy(n)=buoy(n)
       mbesnd(n)=pbesnd(n)
       IF((p(1)-p(n)) < 300.) tcap=AMAX1(tcap,(tv(n)-partv))

     ELSE ! Upper level is above LCL
!
!  BRACKETING CASE AROUND lcl - DRY ADIABATIC PART
!
!        WRITE(6,*)' DRY ADIABATIC PART'
       delta_ht_dry = lcl - ht(n-1)
!        WRITE(6,307)tlcl
!307        format(' PARCEL TEMP AT lcl= ',F10.3)
       CALL intrpr(nlevel,p,tv,plcl,sntlcl)
       partv=tctotv(tlcl,w(1))
       buoy_lcl=(partv-sntlcl)/sntlcl
       pbe_dry=g*0.5*(buoy_lcl+buoy(n-1))*delta_ht_dry
       IF((p(1)-plcl) < 300.) tcap=AMAX1(tcap,(sntlcl-partv))
!        WRITE(6,777)N,P(N),tlcl,sntlcl,buoy_lcl
!#          ,buoy(N-1),delta_ht_dry,HT(N),pbesnd(N-1)+pbe_dry
!
!        MOIST ADIABATIC PART
!
!        WRITE(6,*)' MOIST ADIABATIC PART'
       delta_ht_wet=deltah-delta_ht_dry

       partem(n) = tsa_fast(blthte,p(n))
       wsat=1000.*f_mrsat( p(n)*100., partem(n)+273.15 )
       partv=tctotv(partem(n),wsat)
       buoy(n)=(partv-tv(n))/tv(n)
       pbe_wet = g*0.5*(buoy(n)+buoy_lcl)*delta_ht_wet
       pbesnd(n)=pbesnd(n-1) + pbe_dry + pbe_wet
!
       wload(n)=0.001*(w(1)-wsat)
       mbuoy(n)=buoy(n) - wload(n)
       pbe_wet = g*0.5*(mbuoy(n)+buoy_lcl)*delta_ht_wet
       mbesnd(n)=mbesnd(n-1) + pbe_dry + pbe_wet
       IF((p(1)-plcl) < 300.) tcap=AMAX1(tcap,(tv(n)-partv))

     END IF ! Upper level below LCL (Dry or bracket)
   ELSE ! Lower Level is above LCL
!      WRITE(6,*)' GETTING PARCEL TEMPERATURE FOR MOIST CASE'
     partem(n) = tsa_fast(blthte,p(n))

     wsat=1000.*f_mrsat( p(n)*100., partem(n)+273.15 )
     partv=tctotv(partem(n),wsat)
     buoy(n)=(partv-tv(n))/tv(n)
     pbesnd(n)=pbesnd(n-1)+g*0.5*(buoy(n)+buoy(n-1))*deltah
!
     wload(n)=0.001*(w(1)-wsat)
     mbuoy(n)=buoy(n) - wload(n)
     mbesnd(n)=mbesnd(n-1)+g*0.5*(mbuoy(n)+mbuoy(n-1))*deltah
     IF((p(1)-p(n)) < 300.) tcap=AMAX1(tcap,(tv(n)-partv))

   END IF

!    WRITE(6,777)N,P(N),partem(N),T(N),(buoy(n)*1000.),pbesnd(n)
!777    format(' PBE: P,partem,t,b,pbe=',I3,F6.1,4F8.2)
 END DO
!
!  DETERMINE ENERGY EXTREMA
!  Find heights with nuetral buoyancy
!
 pos_area=0.
 nbe_min=0.
 DO n=2,nlevel
!    WRITE(6,940)N
!940    format(
!    :' LOOKING FOR NEUTRAL BUOYANCY - ENERGY EXTREMUM, LEVEL',I3)

   IF((buoy(n)*buoy(n-1)) < 0.)THEN
     wlow=buoy(n)/(buoy(n)-buoy(n-1))
     htzero=ht(n)*(1.-wlow) + wlow*ht(n-1)
     deltah=htzero-ht(n-1)
     adjeng=pbesnd(n-1)+g*0.5*buoy(n-1)*deltah
!
     IF (p(n) >= 500.)  THEN
       nbe_min=AMIN1(adjeng,nbe_min)
     END IF
!
     pos_area=adjeng-nbe_min
     pos_max=AMAX1(pos_area,pos_max)
   END IF
 END DO

!  WRITE(6,464)ICP,ICT,N1,NLEVEL
!464  format(' ICP,ICT,N1,NLEVEL',4I5)
!
!  Case when equlibrium level is above top of domain
!
 pos_area=pbesnd(nlevel)-nbe_min
 pos_max=AMAX1(pos_area,pos_max)
!
!  At least one region of positive area in sounding
!  Make sure there is at least 1 J/kg to avoid some
!  round-off errors esp near LCL.
!
 IF(pos_max > 1.0)THEN
   pbeneg=AMAX1(nbe_min,-400.)
   velneg=SQRT(2.0*ABS(pbeneg))
   velneg=AMIN1(velneg,20.)
 ELSE ! Case when no positive area exists anywhere in sounding
   pos_max=0.0
   pbeneg =-400.
   velneg = 20.
 END IF
!  WRITE(6,485)pos_max,PBENEG,VELNEG
!485  format(' pos_max',F10.1,' PBENEG',F10.1,' VELNEG',F10.1)

!
!  DETERMINE ENERGY EXTREMA FOR MOIST BUOYANCY
!  Find heights with nuetral buoyancy
!
 pos_area=0.
 nbe_min=0.
 DO n=2,nlevel
!    WRITE(6,940)N

   IF((mbuoy(n)*mbuoy(n-1)) < 0.)THEN
     wlow=mbuoy(n)/(mbuoy(n)-mbuoy(n-1))
     htzero=ht(n)*(1.-wlow) + wlow*ht(n-1)
     deltah=htzero-ht(n-1)
     adjeng=mbesnd(n-1)+g*0.5*mbuoy(n-1)*deltah
!
     IF (p(n) >= 500.)  THEN
       nbe_min=AMIN1(adjeng,nbe_min)
     END IF
!
     pos_area=adjeng-nbe_min
     mpos_max=AMAX1(pos_area,mpos_max)
   END IF
 END DO

!  WRITE(6,464)ICP,ICT,N1,NLEVEL
!
!  Case when equlibrium level is above top of domain
!
 pos_area=mbesnd(nlevel)-nbe_min
 mpos_max=AMAX1(pos_area,mpos_max)
!
!  At least one region of positive area in sounding
!  Make sure there is at least 1 J/kg to
!  spurious pos energy due to round off.
!
 IF(mpos_max > 1.0)THEN
   mbeneg=AMAX1(nbe_min,-400.)
   mvelneg=SQRT(2.0*ABS(pbeneg))
   mvelneg=AMIN1(mvelneg,20.)
 ELSE ! Case when no positive area exists anywhere in sounding
   mpos_max=0.0
   mbeneg =-400.
   mvelneg = 20.
 END IF
!  WRITE(6,486)mpos_max,PBENEG,VELNEG
!486  format(' Mpos_max',F10.1,' MBENEG',F10.1,' mVELNEG',F10.1)
!
!    Case when equlibrium level is above top of domain
!
 mpos_max = MAX(mpos_max,(mbesnd(nlevel) - nbe_min))
!
!  Find EL and LFC
!  Unxts are set to km ASL
!
 IF(pos_max > 1.0) THEN
   IF(buoy(nlevel) > 0.) THEN
     nel=nlevel
     el=0.001*ht(nlevel)
   ELSE
     DO  n=nlevel-1,2,-1
       IF(buoy(n) > 0.) EXIT
     END DO
!      1201     CONTINUE
     nel=n
     wlow=buoy(n+1)/(buoy(n+1)-buoy(n))
     el=0.001 * (ht(n+1)*(1.-wlow) + ht(n)*wlow)
   END IF
!
   DO n=nel,1,-1
     IF(buoy(n) < 0.) EXIT
   END DO
!    1301   CONTINUE
   IF(n > 0) THEN
     wlow=buoy(n+1)/(buoy(n+1)-buoy(n))
     lfc=ht(n+1)*(1.-wlow) + ht(n)*wlow
   ELSE
     lfc=ht(1)
   END IF
 ELSE
   el=0.
   lfc=10000.
 END IF
 lfc=AMIN1(lfc,10000.)
 RETURN
END SUBROUTINE potbe
!
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CALCLI                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######  and modified by Ernani L. Nascimento at Simepar     ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE calcli(nlevel,thepcl,p,t,li)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate Convective Available Potential Energy (CAPE)
!  in each column of ARPS grid.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  Cleaned-up, removed OLAPS artifacts.
!
!  Fev/14/2004  Ernani L. Nascimento: removed variable maxlev
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
!
!  Input variables
!
 INTEGER :: nlevel
 REAL :: thepcl
 REAL :: p(nlevel+2),t(nlevel+2)
!
!  Output variable
!
 REAL :: li
!
!  Functions
!
 REAL :: tsa_fast
!
!  Misc internal variables
!
 INTEGER :: n
 REAL :: dp,wlow,t500,par500
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
 DO n=2,nlevel-1
   IF( p(n) <= 500.) EXIT
 END DO
!  101 CONTINUE
 dp=ALOG(p(n-1)/p(n))
 wlow=ALOG(500./p(n))/dp
 t500=t(n)*(1.-wlow) + t(n-1)*wlow
 par500=tsa_fast(thepcl,500.)
!
 li=t500-par500
!
 RETURN
END SUBROUTINE calcli
!
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE THWXN                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######   and modified by Ernani L. Nascimento at Simepar    ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE thwxn(nlevel,p,ht,t,td,elev,twdf)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  Cleaned-up, removed OLAPS artifacts.
!
!  Feb/14/2004 Ernani L. Nascimento: removed variable maxlev
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
!
!  Input variables
!
 INTEGER :: nlevel
 REAL :: p(nlevel+2),ht(nlevel+2),t(nlevel+2),td(nlevel+2)
 REAL :: elev
!
!  Output variables
!
 REAL :: twdf
!
!  Functions
!
 REAL :: oe,tsa_fast
!
!  Misc internal variables
!
 INTEGER :: n
 REAL :: h3km,thaec,thw,twx,twn
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
 twx=-999.
 twn=999.
 h3km=elev+3000.
 DO n=1,nlevel
   IF(ht(n) >= elev) THEN
     IF(ht(n) > h3km) EXIT
     thaec=oe(t(n),td(n),p(n))
     thw=tsa_fast(thaec,1000.)
     twx=AMAX1(twx,thw)
     twn=AMIN1(twn,thw)
   END IF
 END DO
!  101 CONTINUE
!
!  Find difference between max and min
!
 twdf=twx-twn
 RETURN
END SUBROUTINE thwxn
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION WMR2TD                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
 FUNCTION wmr2td(pres,wmr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate Convective Available Potential Energy (CAPE)
!  in each column of ARPS grid.  => WRONG!!
!
!  CONVERTS WATER VAPOR MIXING RATIO INTO DEW POINT TEMPERATURE
!  (Ernani L. Nascimento)
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on GEMPAK routine of same name.
!
!  MODIFICATION HISTORY:
!
!
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
 REAL :: wmr2td
 REAL :: pres
 REAL :: wmr
 REAL :: wkgkg,e,evap
!
 wkgkg = 0.001 * wmr
 wkgkg = AMAX1(wmr,0.00005)
 e= (pres*wkgkg) / (0.62197 + wkgkg)
 evap = e /(1.001 + (( pres - 100.) /900.) * 0.0034)
 wmr2td = ALOG(evap/6.112) * 243.5 /( 17.67 - ALOG (evap/6.112))

 RETURN
 END FUNCTION wmr2td
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION TCTOTV                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
 FUNCTION tctotv(tt,ww)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Virtual Temperature
!
!  Given T in Celcius and mixing ratio in g/kg
!  find the virtual temperature.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
!
 REAL :: tctotv,tt,ww
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
 tctotv=(tt+273.15)*(1.+0.0006*ww)
 RETURN
 END FUNCTION tctotv
!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION oe(t,td,p)
!
!    g.s. stipanuk     1973          original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982

!   this function returns equivalent potential temperature oe (celsius)
!   of a parcel of air given its temperature t (celsius), dew point
!   td (celsius) and pressure p (millibars).
!   find the wet bulb temperature of the parcel.

 atw = tw(t,td,p)

!   find the equivalent potential temperature.

 oe = os(atw,p)
 RETURN
 END FUNCTION oe
!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION os(t,p)
!
!    g.s. stipanuk     1973          original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982

!   this function returns the equivalent potential temperature os
!   (celsius) for a parcel of air saturated at temperature t (celsius)
!   and pressure p (millibars).
 DATA b/2.6518986/
!   b is an empirical constant approximately equal to the latent heat
!   of vaporization for water divided by the specific heat at constant
!   pressure for dry air.

 tk = t+273.15
 osk= tk*((1000./p)**.286)*(EXP(b*w(t,p)/tk))
 os= osk-273.15
 RETURN
 END FUNCTION os
!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION w(t,p)
!
!    g.s. stipanuk     1973              original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982
!
!  this function returns the mixing ratio (grams of water vapor per
!  kilogram of dry air) given the dew point (celsius) and pressure
!  (millibars). if the temperture  is input instead of the
!  dew point, then saturation mixing ratio (same units) is returned.
!  the formula is found in most meteorological texts.

 x= esat(t)
 w= 622.*x/(p-x)
 RETURN
 END FUNCTION w
!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION esat(t)
!
!    g.s. stipanuk     1973           original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982
!
!   this function returns the saturation vapor pressure over
!   water (mb) given the temperature (celsius).
!   the algorithm is due to nordquist, w.s.,1973: "numerical approxima-
!   tions of selected meteorlolgical parameters for cloud physics prob-
!   lems," ecom-5475, atmospheric sciences laboratory, u.s. army
!   electronics command, white sands missile range, new mexico 88002.

 tk = t+273.15
 p1 = 11.344-0.0303998*tk
 p2 = 3.49149-1302.8844/tk
 c1 = 23.832241-5.02808*ALOG10(tk)
 esat = 10.**(c1-1.3816E-7*10.**p1+8.1328E-3*10.**p2-2949.076/tk)
 RETURN
 END FUNCTION esat
!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION tw(t,td,p)

!    g.s. stipanuk     1973           original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982

!   this function returns the wet-bulb temperature tw (celsius)
!   given the temperature t (celsius), dew point td (celsius)
!   and pressure p (mb).  see p.13 in stipanuk (1973), referenced
!   above, for a description of the technique.
!
!
!   determine the mixing ratio line thru td and p.

 aw = w(td,p)
!
!   determine the dry adiabat thru t and p.

 ao = o(t,p)
 pi = p

!   iterate to locate pressure pi at the intersection of the two
!   curves .  pi has been set to p for the initial guess.

 DO i= 1,10
   x= .02*(tmr(aw,pi)-tda(ao,pi))
   IF (ABS(x) < 0.01) EXIT
   pi= pi*(2.**(x))
 END DO

!   find the temperature on the dry adiabat ao at pressure pi.

 ti= tda(ao,pi)

!   the intersection has been located...now, find a saturation
!   adiabat thru this point. function os returns the equivalent
!   potential temperature (c) of a parcel saturated at temperature
!   ti and pressure pi.

 aos= os(ti,pi)

!   function tsa returns the wet-bulb temperature (c) of a parcel at
!   pressure p whose equivalent potential temperature is aos.

 tw = tsa(aos,p)
 RETURN
 END FUNCTION tw
!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION o(t,p)
!
!    g.s. stipanuk     1973          original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982

!   this function returns potential temperature (celsius) given
!   temperature t (celsius) and pressure p (mb) by solving the poisson
!   equation.

 tk= t+273.15
 ok= tk*((1000./p)**.286)
 o= ok-273.15
 RETURN
 END FUNCTION o
!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION tmr(w,p)
!
!    g.s. stipanuk     1973           original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982

!   this function returns the temperature (celsius) on a mixing
!   ratio line w (g/kg) at pressure p (mb). the formula is given in
!   table 1 on page 7 of stipanuk (1973).
!
!   initialize constants

 DATA c1/.0498646455/,c2/2.4082965/,c3/7.07475/
 DATA c4/38.9114/,c5/.0915/,c6/1.2035/

 x= ALOG10(w*p/(622.+w))
 tmrk= 10.**(c1*x+c2)-c3+c4*((10.**(c5*x)-c6)**2.)
 tmr= tmrk-273.15
 RETURN
 END FUNCTION tmr
!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION tda(o,p)
!
!    g.s. stipanuk     1973           original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982

!   this function returns the temperature tda (celsius) on a dry adiabat
!   at pressure p (millibars). the dry adiabat is given by
!   potential temperature o (celsius). the computation is based on
!   poisson's equation.

 ok= o+273.15
 tdak= ok*((p*.001)**.286)
 tda= tdak-273.15
 RETURN
 END FUNCTION tda
!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 FUNCTION tsa(os,p)
!
!    g.s. stipanuk     1973           original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982

!   this function returns the temperature tsa (celsius) on a saturation
!   adiabat at pressure p (millibars). os is the equivalent potential
!   temperature of the parcel (celsius). sign(a,b) replaces the
!   algebraic sign of a with that of b.
!   b is an empirical constant approximately equal to 0.001 of the latent
!   heat of vaporization for water divided by the specific heat at constant
!   pressure for dry air.

 DATA b/2.6518986/
 a= os+273.15

!   tq is the first guess for tsa.

 tq= 253.15

!   d is an initial value used in the iteration below.

 d= 120.

!   iterate to obtain sufficient accuracy....see table 1, p.8
!   of stipanuk (1973) for equation used in iteration.

 DO i= 1,12
   tqk= tq-273.15
   d= d/2.
   x= a*EXP(-b*w(tqk,p)/tq)-tq*((1000./p)**.286)
   IF (ABS(x) < 1E-7) GOTO 2
   tq= tq+SIGN(d,x)
 END DO
2 tsa= tq-273.15
 RETURN
 END FUNCTION tsa
!
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                FUNCTION TSA_FAST                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
FUNCTION tsa_fast(os,p)
!
!   THIS FUNCTION RETURNS THE TEMPERATURE TSA (CELSIUS) ON A SATURATION
!   ADIABAT AT PRESSURE P (MILLIBARS). OS IS THE EQUIVALENT POTENTIAL
!   TEMPERATURE OF THE PARCEL (CELSIUS). SIGN(A,B) REPLACES THE
!   ALGEBRAIC SIGN OF A WITH THAT OF B.
!
!    BAKER,SCHLATTER 17-MAY-1982     Original version
!    Modification for better convergence, Keith Brewster, Feb 1994.
!
!   B IS AN EMPIRICAL CONSTANT APPROXIMATELY EQUAL TO THE LATENT HEAT
!   OF VAPORIZATION FOR WATER DIVIDED BY THE SPECIFIC HEAT AT CONSTANT
!   PRESSURE FOR DRY AIR.
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate Convective Available Potential Energy (CAPE)
!  in each column of ARPS grid.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  Cleaned-up, removed OLAPS artifacts.
!
!-----------------------------------------------------------------------
!
 REAL :: b
!  PARAMETER (B=2.6518986)
 PARAMETER (b=2651.8986)
 a= os+273.15
!
!   Above 200 mb figure all the moisture is wrung-out, so
!   the temperature is that which has potential temp of theta-e.
!   Otherwise iterate to find combo of moisture and temp corresponding
!   to thetae.
!
 IF( p < 200.) THEN
   tq=a*((p/1000.)**.286)
 ELSE
!   D IS AN INITIAL VALUE USED IN THE ITERATION BELOW.
   d= 120.
!   TQ IS THE FIRST GUESS FOR TSA.
   tq= 253.15
   x = 0.
!
!   ITERATE TO OBTAIN SUFFICIENT ACCURACY....SEE TABLE 1, P.8
!   OF STIPANUK (1973) FOR EQUATION USED IN ITERATION.
   DO i= 1,25
     d= 0.5*d

     x_last = x

     x= a*EXP(-b*f_mrsat(p*100.,tq)/tq)-tq*((1000./p)**.286)

     IF (ABS(x) < 1E-3) GO TO 2
!
     IF (x_last * x < 0.) THEN
       slope = (x-x_last) / (tq - tq_last)
       delta = - x / slope
       ad = AMIN1(ABS(delta),d)
       tq_last = tq
       tq = tq + SIGN(ad,delta)
     ELSE
       tq_last = tq
       tq= tq+SIGN(d,x)
     END IF

   END DO
 END IF
2 tsa_fast = tq-273.15
 RETURN
END FUNCTION tsa_fast
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION F_MRSAT                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
FUNCTION f_mrsat( p, t )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the saturation water vapor mixing ratio using enhanced
!  Teten's formula.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  01/08/1998
!
!  MODIFICATION HISTORY:
!
!  16 FEB 2004
!  Ernani L. Nascimento: modified to include constant rddrv without
!  INCLUDE command.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    p        Pressure (Pascal)
!    t        Temperature (K)
!
!  OUTPUT:
!
!    f_mrsat  Saturation water vapor mixing ratio (kg/kg).
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE

 REAL :: p         ! Pressure (Pascal)
 REAL :: t         ! Temperature (K)
 REAL :: f_mrsat   ! Saturation water vapor mixing ratio (kg/kg)

 REAL :: rd        ! Gas constant for dry air  (m**2/(s**2*K))
 PARAMETER( rd     = 287.0 )

 REAL :: rv        ! Gas constant for water vapor  (m**2/(s**2*K)).
 PARAMETER( rv     = 461.0 )

 REAL :: rddrv
 PARAMETER( rddrv  = rd/rv )
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
 REAL :: fes
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
! Removed line below. ERNANI L. NASCIMENTO
!  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Function f_es and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
 REAL :: f_es
!fpp$ expand (f_es)
!dir$ inline always f_es
!*$*  inline routine (f_es)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
 fes = f_es( p,t )
 f_mrsat = rddrv * fes / (p-fes)

 RETURN
END FUNCTION f_mrsat
!
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION F_ES                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
FUNCTION f_es( p, t )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the saturation specific humidity using enhanced Teten's
!  formula.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  01/08/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    p        Pressure (Pascal)
!    t        Temperature (K)
!
!  OUTPUT:
!
!    f_es     Saturation water vapor pressure (Pa)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE

 REAL :: p         ! Pressure (Pascal)
 REAL :: t         ! Temperature (K)
 REAL :: f_es      ! Saturation water vapor pressure (Pa)
!
!-----------------------------------------------------------------------
!
!  Function f_es and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
 REAL :: f_esl, f_esi

!fpp$ expand (f_esl)
!fpp$ expand (f_esi)
!dir$ inline always f_esl, f_esi
!*$*  inline routine (f_esl, f_esi)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
 IF ( t >= 273.15 ) THEN      ! for water
   f_es = f_esl( p,t )
 ELSE                            ! for ice
   f_es = f_esi( p,t )
 END IF

 RETURN
END FUNCTION f_es
!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
FUNCTION f_esl( p, t )
!-----------------------------------------------------------------------
!  Calculate the saturation water vapor over liquid water using
!  enhanced Teten's formula.
!
!  Feb 16 2004:
!  Modified by Ernani L. Nascimento to include constants satfwa,satfwb,
!  satewa,satewb,satewc without INCLUDE command.
!-----------------------------------------------------------------------
 IMPLICIT NONE

 REAL :: p         ! Pressure (Pascal)
 REAL :: t         ! Temperature (K)
 REAL :: f_esl     ! Saturation water vapor pressure over liquid water

 REAL :: f

 REAL :: satfwa, satfwb
 PARAMETER ( satfwa = 1.0007 )
 PARAMETER ( satfwb = 3.46E-8 )  ! for p in Pa

 REAL :: satewa, satewb, satewc
 PARAMETER ( satewa = 611.21 )   ! es in Pa
 PARAMETER ( satewb = 17.502 )
 PARAMETER ( satewc = 32.18 )

!  Removed line below. ERNANI L. NASCIMENTO
!  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
 f = satfwa + satfwb * p
 f_esl = f * satewa * EXP( satewb*(t-273.15)/(t-satewc) )

 RETURN
END FUNCTION f_esl
!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
FUNCTION f_esi( p, t )
!-----------------------------------------------------------------------
!  Calculate the saturation water vapor over ice using enhanced
!  Teten's formula.
!
!  Feb 16 2004:
!  Modified by Ernani L. Nascimento to include constants satfia,satfib,
!  sateia,sateib,sateic without INCLUDE command.
!-----------------------------------------------------------------------
 IMPLICIT NONE

 REAL :: p         ! Pressure (Pascal)
 REAL :: t         ! Temperature (K)
 REAL :: f_esi     ! Saturation water vapor pressure over ice (Pa)

 REAL :: f

 REAL :: satfia, satfib
 PARAMETER ( satfia = 1.0003 )
 PARAMETER ( satfib = 4.18E-8 )  ! for p in Pa

 REAL :: sateia, sateib, sateic
 PARAMETER ( sateia = 611.15 )   ! es in Pa
 PARAMETER ( sateib = 22.452 )
 PARAMETER ( sateic = 0.6 )

!  Removed line below. ERNANI L. NASCIMENTO
!  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
 f = satfia + satfib * p
 f_esi = f * sateia * EXP( sateib*(t-273.15)/(t-sateic) )

 RETURN
END FUNCTION f_esi

!
!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                SUBROUTINE DDFF2UV                    ######
!     ######                                                      ######
!     ######                     Developed by                     ######
!     ######     Center for Analysis and Prediction of Storms     ######
!     ######                University of Oklahoma                ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!
     SUBROUTINE ddff2uv(dd,ff,u,v)
!
!#######################################################################
!
!     PURPOSE:
!
!     Calculate u and v wind components from direction and speed.
!
!#######################################################################
!
!
!     AUTHOR: Keith Brewster
!     3/11/1996
!
!     09/feb/2004 (ERNANI L. NASCIMENTO)
!       THIS SUBROUTINE IS EXACTLY LIKE SUBROUTINE  ddff2uv IN
!     /src/adas/thermo3d.f  I´VE JUST CHANGED THE NAME.
!
!#######################################################################
!
  implicit none
  integer j,k
  real :: u,v,dd,ff
!   dimension  u(100,35),v(100,35),dd(100,35),ff(100,35)
  real :: arg
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!   arg = (dd(j,k) * (3.141592654/180.))
  arg = (dd * (3.141592654/180.))
!   if (m==1) then
!    write(*,*) 'Dentro do ddff2uv: ',dd,ff,arg
!   endif
!   u(j,k) = -ff(j,k) * sin(arg)
  u = -ff * sin(arg)
!   v(j,k) = -ff(j,k) * cos(arg)
  v = -ff * cos(arg)
  RETURN
END SUBROUTINE ddff2uv
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE UV2DDFF                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE uv2ddff(u,v,dd,ff)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate direction and speed from u and v.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  3/11/1996
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
 REAL :: u,v,dd,ff
 REAL :: dlon
 REAL :: r2deg
 PARAMETER (r2deg=180./3.141592654)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
 ff = SQRT(u*u + v*v)

 IF(v > 0.) THEN
   dlon=r2deg*ATAN(u/v)
 ELSE IF(v < 0.) THEN
   dlon=180. + r2deg*ATAN(u/v)
 ELSE IF(u >= 0.) THEN
   dlon=90.
 ELSE
   dlon=-90.
 END IF

 dd= dlon + 180.
 dd= dd-360.*(nint(dd)/360)
 RETURN
END SUBROUTINE uv2ddff
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CALCSHR                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE calcshr(nrad,zp,sigma,                                       &
          p_pa,t3d,u3d,v3d,cape,                                       &
          shr37,ustrm,vstrm,srlfl,srmfl,helicity,brn,brnu,brnu2km,blcon, &
          tem2,tem3,nlev,estacao,hora,dia,mes,ano,strm_opt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate various wind shear parameters useful for gauging
!  the potential for severe storms.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  February, 1994   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  Keith Brewster
!  Added storm-relative flows, general clean-up to
!  meet ARPS coding standards.
!
!  3/22/1996  (Keith Brewster)
!  Fixed some bugs, added smoothing at the end.
!
!  06/20/2000 (Eric Kemp and Keith Brewster)
!  Changed BRN Shear to be the denominator of BRN, instead of wind
!  speed, now has units of speed squared.
!
!  04/19/2004 (Ernani L. Nascimento)
!  Modified the code to fit in índices_severos3.f90, and adapted
!  the computation of expected storm-motion for the Southern Hemisphere.
!
!-----------------------------------------------------------------------
!
!  Calculates some of the shear related variables from the
!  ARPS 3D wind field.
!
!  shr37         Magnitude of wind shear between 3 and 7 km AGL
!  ustrm,vstrm   Estimated storm motion (modified from Bob Johns)
!  srlfl         Low-level storm-relative wind
!  srmfl         Mid-level storm-relative wind
!  helicity      Helicity, storm relative
!  brn           Bulk Richardson Number (Weisman and Klemp)
!  brnu          Shear parameter of BRN, "U"
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
 INTEGER :: nrad,strm_opt
!
!-----------------------------------------------------------------------
!
!  "1-D" input variables (Ernani L. Nascimento)
!
!-----------------------------------------------------------------------
!
 INTEGER :: nlev(nrad),estacao(nrad),hora(nrad),dia(nrad),mes(nrad),   &
            ano(nrad)
 REAL :: sigma(nrad)
 REAL :: cape(nrad)
!
!-----------------------------------------------------------------------
!
!  "2-D" input variables (Ernani L. Nascimento)
!
!-----------------------------------------------------------------------
!
 INTEGER :: zp(102,nrad)
 REAL :: p_pa(102,nrad)   ! Pressure in Pascals
 REAL :: t3d(102,nrad)    ! Temperature in Kelvin
 REAL :: u3d(102,nrad)
 REAL :: v3d(102,nrad)

!
!-----------------------------------------------------------------------
!
!  Output variables
!
!-----------------------------------------------------------------------
!
 REAL :: shr37(nrad)      ! 7km - 3km wind shear
 REAL :: ustrm(nrad)
 REAL :: vstrm(nrad)  ! Estimated storm motion (Bob Johns)
 REAL :: srlfl(nrad)
 REAL :: srmfl(nrad)
 REAL :: helicity(nrad)   ! Helicity, storm relative
 REAL :: brn(nrad)        ! Bulk Richardson Number (Weisman and Klemp)
 REAL :: brnu(nrad)       ! Shear parameter of BRN, "U"
 REAL :: brnu2km(nrad)    ! 2km shear parameter of BRN
 REAL :: blcon(nrad)
!
!-----------------------------------------------------------------------
!
!  Temporary variables
!
!-----------------------------------------------------------------------
!

 REAL :: tem2(nrad)
 REAL :: tem3(nrad)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
!  INTEGER :: imid,jmid
 REAL :: elev(nrad),hgt3d(102,nrad)
 INTEGER :: i,j,k,ksfc,k2km,k3km
 REAL :: h3km,u3km,v3km,h7km,u7km,v7km,u2,v2
 REAL :: p2km,t2km,h2km,u2km,v2km,p9km,t9km,h9km,u9km,v9km
 REAL :: p500m,t500m,h500m,u500m,v500m,p6km,t6km,h6km,u6km,v6km
 REAL :: sumu,sumv,sump,sumh,wlow,whigh,dx,dy,dz,dp,dx2,dy2
 REAL :: rhohi,rholo,rhoinv,arg,new_ustrm,new_vstrm
 REAL :: dirmean,spmean,ushr,vshr,ddir,perc,obs_ddstrm,obs_ffstrm
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!  Loop over all soundings (Ernani L. Nascimento)
!-----------------------------------------------------------------------
!
 write(*,*) 'Entrada da subrotina calcshr'
 DO k=1,nrad
  if (nlev(k)==0) then
   shr37(k)=-999.9
   ustrm(k)=-999.9
   vstrm(k)=-999.9
   srlfl(k)=-999.9
   srmfl(k)=-999.9
   helicity(k)=-999.9
   brn(k)=-999.9
   brnu(k)=-999.9
   brnu2km(k)=-999.9
   blcon(k)=-999.9
  else
   print*, estacao(k),hora(k),dia(k),mes(k),ano(k)
   j=nlev(k)
   if ((p_pa(j,k)*0.01)>300.0) then
    cycle
   endif

   arg=0.0
   elev(k)=real(zp(1,k))
!    write(*,*) 'elev(k)= ',elev(k),' para k= ',k
   DO j=1,nlev(k)
     hgt3d(j,k)=real(zp(j,k))
   ENDDO
!   write(*,*) 'PASSEI AQUI 3, com k= ',k
!
!-----------------------------------------------------------------------
!
!  Find mass weighted mean wind in first 500m
!
!-----------------------------------------------------------------------
!
     sumu=0.
     sumv=0.
     sump=0.
     h500m=elev(k)+500.
     DO j=2,nlev(k)
       IF( hgt3d(j,k) < h500m ) THEN
         dp=p_pa(j-1,k)-p_pa(j,k)
         rhohi=p_pa(j,k)/t3d(j,k)
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u3d(j,k))
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v3d(j,k))
         sump=sump+dp
       ELSE
         dz=hgt3d(j,k)-hgt3d(j-1,k)
         wlow=(hgt3d(j,k)-h500m)/dz
         u500m=u3d(j,k)*(1.-wlow) + u3d(j-1,k)*wlow
         v500m=v3d(j,k)*(1.-wlow) + v3d(j-1,k)*wlow
         p500m=p_pa(j,k)*(1.-wlow) + p_pa(j-1,k)*wlow
         t500m=t3d(j,k)*(1.-wlow) + t3d(j-1,k)*wlow
         dp=p_pa(j-1,k)-p500m
         rhohi=p500m/t500m
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u500m)
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v500m)
         sump=sump+dp
!          print *, ' sumu,sumv,sump = ',sumu,sumv,sump,'para k= ',k
         EXIT
       END IF
     END DO
!      121    CONTINUE
     u500m=sumu/sump
     v500m=sumv/sump
!     write(*,*) 'PASSEI AQUI 4, com k= ',k
!      write(*,*) '(u500m,v500m)= ',u500m,v500m,' para k= ',k
!
!-----------------------------------------------------------------------
!
!  Find mass weighted mean wind sfc-2km AGL
!
!-----------------------------------------------------------------------
!
     sumu=0.
     sumv=0.
     sump=0.
     h2km=elev(k)+2000.
     DO j=2,nlev(k)
       IF( hgt3d(j,k) < h2km ) THEN
         dp=p_pa(j-1,k)-p_pa(j,k)
         rhohi=p_pa(j,k)/t3d(j,k)
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u3d(j,k))
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v3d(j,k))
         sump=sump+dp
       ELSE
         dz=hgt3d(j,k)-hgt3d(j-1,k)
         wlow=(hgt3d(j,k)-h2km)/dz
         u2km=u3d(j,k)*(1.-wlow) + u3d(j-1,k)*wlow
         v2km=v3d(j,k)*(1.-wlow) + v3d(j-1,k)*wlow
         p2km=p_pa(j,k)*(1.-wlow) + p_pa(j-1,k)*wlow
         t2km=t3d(j,k)*(1.-wlow) + t3d(j-1,k)*wlow
         dp=p_pa(j-1,k)-p2km
         rhohi=p2km/t2km
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u2km)
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v2km)
         sump=sump+dp
         EXIT
       END IF
     END DO
!      141    CONTINUE
!     write(*,*) 'PASSEI AQUI 4.5, com k= ',k
     u2km=sumu/sump
     v2km=sumv/sump
!     write(*,*) 'u2km,v2km: ',u2km,v2km
     tem2(k)=u2km
     tem3(k)=v2km
!     write(*,*) 'tem2(k),tem3(k): ',tem2(k),tem3(k)
!     write(*,*) 'PASSEI AQUI 5, com k= ',k
!
!-----------------------------------------------------------------------
!
!  Find mass weighted mean wind sfc-6km AGL
!
!-----------------------------------------------------------------------
!
     sumu=0.
     sumv=0.
     sump=0.
     h6km=elev(k)+6000.
     DO j=2,nlev(k)
       IF( hgt3d(j,k) < h6km ) THEN
         dp=p_pa(j-1,k)-p_pa(j,k)
         rhohi=p_pa(j,k)/t3d(j,k)
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u3d(j,k))
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v3d(j,k))
         sump=sump+dp
       ELSE
         dz=hgt3d(j,k)-hgt3d(j-1,k)
         wlow=(hgt3d(j,k)-h6km)/dz
         u6km=u3d(j,k)*(1.-wlow) + u3d(j-1,k)*wlow
         v6km=v3d(j,k)*(1.-wlow) + v3d(j-1,k)*wlow
         p6km=p_pa(j,k)*(1.-wlow) + p_pa(j-1,k)*wlow
         t6km=t3d(j,k)*(1.-wlow) + t3d(j-1,k)*wlow
         dp=p_pa(j-1,k)-p6km
         rhohi=p6km/t6km
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u6km)
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v6km)
         sump=sump+dp
         EXIT
       END IF
     END DO
     u6km=sumu/sump
     v6km=sumv/sump
!     write(*,*) 'PASSEI AQUI 6, com k= ',k,' e (u6km,v6km)= ',u6km,v6km
!
!-----------------------------------------------------------------------
!
!  Find mass weighted mean wind 2km-9km AGL
!
!-----------------------------------------------------------------------
!
     sumu=0.
     sumv=0.
     sump=0.
     h9km=elev(k)+9000.
     DO j=2,nlev(k)
       IF( hgt3d(j,k) > h2km ) EXIT
     END DO
!      181   CONTINUE
     k2km=j
     dz=hgt3d(j,k)-hgt3d(j-1,k)
     wlow=(hgt3d(j,k)-h2km)/dz
     u2=u3d(j,k)*(1.-wlow) + u3d(j-1,k)*wlow
     v2=v3d(j,k)*(1.-wlow) + v3d(j-1,k)*wlow
     p2km=p_pa(j,k)*(1.-wlow) + p_pa(j-1,k)*wlow
     t2km=t3d(j,k)*(1.-wlow) + t3d(j-1,k)*wlow
     dp=p2km-p_pa(j,k)
     rholo=p2km/t2km
     rhohi=p_pa(j,k)/t3d(j,k)
     rhoinv=1./(rhohi+rholo)
     sumu=sumu+dp*rhoinv*(rholo*u2+rhohi*u3d(j,k))
     sumv=sumv+dp*rhoinv*(rholo*v2+rhohi*v3d(j,k))
     sump=sump+dp
     DO j=k2km+1,nlev(k)
       IF( hgt3d(j,k) < h9km ) THEN
         dp=p_pa(j-1,k)-p_pa(j,k)
         rhohi=p_pa(j,k)/t3d(j,k)
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u3d(j,k))
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v3d(j,k))
         sump=sump+dp
       ELSE
         dz=hgt3d(j,k)-hgt3d(j-1,k)
         wlow=(hgt3d(j,k)-h9km)/dz
         u9km=u3d(j,k)*(1.-wlow) + u3d(j-1,k)*wlow
         v9km=v3d(j,k)*(1.-wlow) + v3d(j-1,k)*wlow
         p9km=p_pa(j,k)*(1.-wlow) + p_pa(j-1,k)*wlow
         t9km=t3d(j,k)*(1.-wlow) + t3d(j-1,k)*wlow
         dp=p_pa(j-1,k)-p9km
         rhohi=p9km/t9km
         rholo=p_pa(j-1,k)/t3d(j-1,k)
         rhoinv=1./(rhohi+rholo)
         sumu=sumu+dp*rhoinv*(rholo*u3d(j-1,k)+rhohi*u9km)
         sumv=sumv+dp*rhoinv*(rholo*v3d(j-1,k)+rhohi*v9km)
         sump=sump+dp
         EXIT
       END IF
     END DO
!      191    CONTINUE
     u9km=sumu/sump
     v9km=sumv/sump
!     write(*,*) 'PASSEI AQUI 7, com k= ',k
!
!-----------------------------------------------------------------------
!
!  Storm motion estimation
!  From Davies and Johns, 1993
!  "Some wind and instability parameters associated With
!  strong and violent tornadoes."
!  AGU Monograph 79, The Tornado...(Page 575)
!
!  Becuase of the discontinuity produced by that method
!  at the 15.5 m/s cutoff, their rules have been modified
!  to provide a gradual transition, and accomodate all the
!  data they mention in the article.
!
!  (04/19/2004) Modified by Ernani L. Nascimento. Adapting for the
!  Southern Hemisphere.
!
!-----------------------------------------------------------------------
!
     CALL uv2ddff(u6km,v6km,dirmean,spmean)
!      write(*,*) 'PASSAGEM 5.1, com k= ',k
!      write(*,*) 'PASSAGEM 5.1, com k= ',k,' e spmean= ',spmean
!      write(*,*) 'PASSAGEM 5.1, com k= ',k,' e dirmean= ',dirmean
     IF(spmean >= 20.0) THEN
!        write(*,*) 'Hei, passei aqui!'
       dirmean=dirmean-18.
       IF(dirmean <= 0.) dirmean=360.+dirmean
       spmean=spmean*0.89
     ELSE IF (spmean > 8.0) THEN
!        write(*,*) 'Não! Passei foi aqui!'
       whigh=(spmean - 8.0)/12.
       wlow =1.-whigh
       ddir=wlow*32.0 + whigh*18.0
       perc=wlow*0.75 + whigh*0.89
       dirmean=dirmean-ddir
       IF(dirmean <= 0.) dirmean=360.+dirmean
       spmean=spmean*perc
     ELSE
!        write(*,*) 'Na verdade, passei foi aqui!'
       dirmean=dirmean-32.
       IF(dirmean <= 0.) dirmean=360.+dirmean
       spmean=spmean*0.75
     END IF
!      write(*,*) 'PASSAGEM 5.2, com k= ',k
     arg = (dirmean * (3.141592654/180.))
     ustrm(k) = -spmean * sin(arg)
     vstrm(k) = -spmean * cos(arg)
! Utilizando movimento estimado da célula da esquerda em 9 de outubro de 2003
     IF (strm_opt == 1) then
      IF ((estacao(k)==83827).and.(hora(k)==00).and.(dia(k)==09).and. &
          (mes(k)==10).and.(ano(k)==2003)) THEN
       obs_ddstrm=270.0
       obs_ffstrm=10.7
       CALL ddff2uv(obs_ddstrm,obs_ffstrm,new_ustrm,new_vstrm)
       ustrm(k)=new_ustrm
       vstrm(k)=new_vstrm
      ENDIF
! Utilizando movimento estimado da célula da célula de 24 de maio de 2005
!       IF ((estacao(k)==99999).and.(hora(k)==18).and.(dia(k)==24).and. &
!           (mes(k)==05).and.(ano(k)==2005)) THEN
!        obs_ddstrm=297.0
!        obs_ffstrm=17.8
!        CALL ddff2uv(obs_ddstrm,obs_ffstrm,new_ustrm,new_vstrm)
!        ustrm(k)=new_ustrm
!        vstrm(k)=new_vstrm
!       ENDIF
     ENDIF
!      CALL ddff2uv(dirmean,spmean,ustrm(k),vstrm(k),5)
!      write(*,*) 'PASSAGEM 6, com k= ',k
!
!-----------------------------------------------------------------------
!
!  Storm-relative low-level flow
!
!-----------------------------------------------------------------------
!
     srlfl(k)=SQRT((ustrm(k)-u2km)*(ustrm(k)-u2km) +             &
                     (vstrm(k)-v2km)*(vstrm(k)-v2km))
!
!-----------------------------------------------------------------------
!
!  Storm relative mid-level flow
!
!-----------------------------------------------------------------------
!
     srmfl(k)=SQRT((ustrm(k)-u9km)*(ustrm(k)-u9km) +             &
                     (vstrm(k)-v9km)*(vstrm(k)-v9km))
!
!-----------------------------------------------------------------------
!
!  Shear parameter for Bulk Richardson number
!
!-----------------------------------------------------------------------
!
!    print *, ' density-weight mean 0-500 m ',u500m,v500m
!    print *, ' density-weight mean 0-6  km ',u6km,v6km
!
     brnu(k)=0.5*( (u6km-u500m)*(u6km-u500m) +                       &
                     (v6km-v500m)*(v6km-v500m) )
! Added the 2km BRNSHR (Ernani L. Nascimento, 03/feb/2005
     brnu2km(k)=0.5*( (u2km-u500m)*(u2km-u500m) +                    &
                     (v2km-v500m)*(v2km-v500m) )
!     write(*,*) 'PASSAGEM 7, com k= ',k
!
!-----------------------------------------------------------------------
!
!  Bulk Richardson number
!  A limit of 200 is imposed, since this could
!  go to inifinity.
!
!-----------------------------------------------------------------------
!
     IF(brnu(k) > 0.) THEN
       brn(k)=cape(k)/brnu(k)
       brn(k)=AMIN1(brn(k),200.)
     ELSE
       brn(k)=200.
     END IF
!    write(*,*) 'PASSEI AQUI 8, com k= ',k
!
!
!-----------------------------------------------------------------------
!
!  Calculate Helicity and 3km to 7km shear
!  since both involve the 3km wind.
!
!  For more efficient computation the Helicity is
!  computed for zero storm motion and the storm
!  motion is accounted for by adding a term at the end.
!  This is mathematically equivalent to accounting
!  for the storm motion at each level.
!
!-----------------------------------------------------------------------
!
     h3km=elev(k)+3000.
     h7km=elev(k)+7000.
!
!-----------------------------------------------------------------------
!
!  Find level just above 3 km AGL
!  Note, it is assumed here that there is at least
!  one level between the sfc and 3 km.
!
!-----------------------------------------------------------------------
!
     sumh=0.
     DO j=2,nlev(k)
       IF(hgt3d(j,k) > h3km) EXIT
       sumh=sumh +                                                     &
           ( u3d(j,k)*v3d(j-1,k) ) -                                   &
           ( v3d(j,k)*u3d(j-1,k) )
     END DO
!      240    CONTINUE
     k3km=j
     dz=hgt3d(j,k)-hgt3d(j-1,k)
     wlow=(hgt3d(j,k)-h3km)/dz
     u3km=u3d(j,k)*(1.-wlow) + u3d(j-1,k)*wlow
     v3km=v3d(j,k)*(1.-wlow) + v3d(j-1,k)*wlow
     sumh=sumh +                                                       &
           ( u3km*v3d(j-1,k) ) -                                       &
           ( v3km*u3d(j-1,k) )
     ushr=u3km-u3d(1,k)
     vshr=v3km-v3d(1,k)
     helicity(k)=sumh + vshr*ustrm(k) - ushr*vstrm(k)
!      write(*,*) 'PASSAGEM 9, com k= ',k
!
!-----------------------------------------------------------------------
!
!  Now Find 7km wind for 3-to-7km shear
!
!-----------------------------------------------------------------------
!
     DO j=k3km,nlev(k)
       IF(hgt3d(j,k) > h7km) EXIT
     END DO
!      260   CONTINUE
     dz=hgt3d(j,k)-hgt3d(j-1,k)
     wlow=(hgt3d(j,k)-h7km)/dz
     u7km=u3d(j,k)*(1.-wlow) + u3d(j-1,k)*wlow
     v7km=v3d(j,k)*(1.-wlow) + v3d(j-1,k)*wlow
     shr37(k)=(SQRT( (u7km-u3km)*(u7km-u3km) +                       &
                       (v7km-v3km)*(v7km-v3km) ))/4000.
  endif
 END DO
!  write(*,*) 'Saída da subrotina calcshr.'

!-----------------------------------------------------------------------
!
!  Calculate the low-level convergence
!
!-----------------------------------------------------------------------
!
!  dx=(x(2)-x(1))
!  dy=(y(2)-y(1))
!  dx2=2.*dx
!  dy2=2.*dy
!  DO j=2,ny-2
!    DO i=2,nx-2
!      blcon(i,j)=-1000.*(                                               &
!                 (tem2(i+1,j)-tem2(i-1,j))/(sigma(i,j)*dx2)+            &
!                 (tem3(i,j+1)-tem3(i,j-1))/(sigma(i,j)*dy2) )
!    END DO
!  END DO
!
!  DO j=2,ny-2
!    blcon(1,j)=-1000.*(                                                 &
!                (tem2(2,j)-tem2(1,j))/(sigma(1,j)*dx)+                  &
!                (tem3(1,j+1)-tem3(1,j-1))/(sigma(1,j)*dy2) )
!    blcon(nx-1,j)=-1000.*(                                              &
!           (tem2(nx-1,j)-tem2(nx-2,j))/(sigma(nx-1,j)*dx)+              &
!           (tem3(nx-1,j+1)-tem3(nx-1,j-1))/(sigma(nx-1,j)*dy2) )
!  END DO
!  DO i=2,nx-2
!    blcon(i,1)=-1000.*(                                                 &
!                (tem2(i+1,1)-tem2(i-1,1))/(sigma(i,1)*dx2)+             &
!                (tem3(i,2)-tem3(i,1))/(sigma(i,1)*dy) )
!    blcon(i,ny-1)=-1000.*(                                              &
!            (tem2(i+1,ny-1)-tem2(i-1,ny-1))/(sigma(i,ny-1)*dx2)+        &
!            (tem3(i,ny-1)-tem3(i,ny-2))/(sigma(i,ny-1)*dy) )
!  END DO
!  blcon(1,1)=blcon(2,2)
!  blcon(nx-1,ny-1)=blcon(nx-2,ny-2)
!
!-----------------------------------------------------------------------
!
!  Smooth the output arrays
!
!-----------------------------------------------------------------------
!
!  CALL smooth9p(   shr37,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(   ustrm,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(   vstrm,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(   srlfl,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(   srmfl,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(helicity,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(     brn,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(    brnu,nx,ny,1,nx-1,1,ny-1,tem1)
!  CALL smooth9p(   blcon,nx,ny,1,nx-1,1,ny-1,tem1)
!
 RETURN
END SUBROUTINE calcshr
!
!
!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                SUBROUTINE CALCEHI_SUP                ######
!     ######                                                      ######
!     ######                     Developed by                     ######
!     ######     Center for Analysis and Prediction of Storms     ######
!     ######                University of Oklahoma                ######
!     ######        and by Ernani L. Nascimento at SIMEPAR        ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!
SUBROUTINE calcehi_sup(nrad,nlev,presshpa,cape,heli,dnrv,ehi,sup)
!
!#######################################################################
!
!     PURPOSE:
!
!     Calculate the energy helicity-index (EHI) (DAVIES 1993; RASMUSSEN AND
!     BLANCHARD 1998), and the supercell composite index (THOMPSON et al
!     2003).
!
!#######################################################################
!
!     AUTHOR: Ernani L. Nascimento, adapting code from arpspltderive.f
!     09 February, 2004
!
!     MODIFICATION HISTORY:
!     THIS IS A NEW SUBROUTINE, NOT AVAILABLE IN THE ORIGINAL ARPS CODE
!
!     (04/19/2004) (Ernani L. Nascimento)
!     Code was modified to fit in índices_severos3.f90
!
!     (08/02/2004) (Ernani L. Nascimento)
!     Code was modified to include the computation of the supercell
!     composite index.
!
!#######################################################################
!
!     Calculates EHI and SUP from cape, helicity and bulk Richardson
!     shear.
!
!
!     nrad          Number of soundings for which EHI will be computed
!     cape          Convective available potential energy (J/kg)
!     heli          Helicity, storm relative (m²/s²)
!     dnrv          Bulk Richardson number shear (m²/s²)
!     ehi           Energy-helicity index (non-dimensional)
!     sup           Supercell composite index
!
!#######################################################################
!
   implicit none
!
!#######################################################################
!
!     Input variables
!
!#######################################################################
!
   integer :: nrad,k,j
   integer :: nlev(nrad)
   real :: cape(nrad), heli(nrad), dnrv(nrad)
   real :: presshpa(102,nrad)
!
!#######################################################################
!
!     Output variables
!
!#######################################################################
!
   real :: ehi(nrad)      ! Energy-helicity index
   real :: sup(nrad)      ! Supercell composite index

!#######################################################################
!     Work variable (constant)
!#######################################################################

   real :: denomi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
   denomi=160000.0
   DO  k=1,nrad
    if (nlev(k)==0) then
     ehi(k)=-999.9
     sup(k)=-999.9
    else
     j=nlev(k)
     if (presshpa(j,k)>300.0) then
      cycle
     endif
     ehi(k) = ( cape(k)*heli(k) )/denomi
     sup(k) = ( (cape(k)/1000.)*(heli(k)/150.)*(dnrv(k)/40.) )
    endif
   ENDDO

 RETURN
END SUBROUTINE calcehi_sup
