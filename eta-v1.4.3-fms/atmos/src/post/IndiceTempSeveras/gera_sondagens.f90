
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
