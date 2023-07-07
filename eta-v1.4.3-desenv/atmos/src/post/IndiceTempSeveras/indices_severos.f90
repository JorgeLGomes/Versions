! --------------------------------------------------------------------------
!  PROGRAMA �dices_severos5
!
!           Autor: ERNANI DE LIMA NASCIMENTO
!           Data: Outubro de 2006
!
!  Este programa calcula �dices de tempo severo a partir de perfis
!  verticais atmosf�icos gerados por modelos num�icos.
!  NOTA 1: ESTA VERS� CALCULA �DICES SEVEROS TERMODIN�ICOS E CINEM�ICOS.
!  NOTA 2: ESTA VERS� �ID�TICA �VERS� �dices_severos4, EXCETO QUE L�!          UM ARQUIVO DE ENTRADA COM CONTE�O DIFERENTE E GERA UM ARQUIVO
!          BIN�IO PARA O GRADS, AL� DO ARQUIVO ASCII PARA O ARPS.
!  NOTA 3: ESTA VERS� GERA SA�A DE SONDAGENS EM FORMATO PARA EXECU�O NO
!          ARPS. POR EXEMPLO, L�O ARQUIVO eta.sound_0212_0318.txt E GERA
!          O ARQUIVO arps.eta.sound_0212_0318.txt. DEPOIS ESTE ARQUIVO
!          DEVE SER EDITADO PARA INCLUIR O CABE�LHO DESEJADO E ENT�
!          SER SALVO COM EXTENS� .snd PARA SER UTILIZADO COM O ARPS
!          (LEIA O ARQUIVO LEIA-ME PARA MAIORES DETALHES).
!
!
!  Exemplo de arquivo de entrada necess�io para este programa (este formato
!  �aquele gerado pelo programa gera_sond_2.f90):
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
! O que significa cada nmero no exemplo acima:
! 78122 => identificador da "esta�o" (irrelevante neste programa)
! 20 => nmero de n�eis do perfil
! 18 => hor�io UTC
! 24 => dia
! 05 => m�
! 2005 => ano
! Colunas: press� (hPa), altitude (m), temperatura (C), temperatura do ponto
! de orvalho (C), umidade relativa (%), raz� de mistura (g/kg), componente zonal
! do vento (m/s), componente meridional do vento (m/s), temperatura potencial (K).
!
!
! �dices calculados: CAPE de superf�ie, CAPE com corre�o de water loading,
! indice de levantamento, n�el de condensa�o por levantamento, n�el de convec�o
! espont�ea, inibi�o convectiva, componentes u e v do movimento esperado de tempestades,
! escoamento relativo �tempestade em n�eis baixos e m�ios, helicidade relativa �! tempestade, nmero de Richardson volum�rico, denominador do nmero de Richardson
! volum�rico, �dice de energia-helicidade.
!
! Com exce�o da rotina de c�culo do �dice de energia-helicidade (interiramente nova),
! todas as rotinas foram modificadas a partir do ARPS5.0_IHOP6 e ADAS.
!
!
! ---------------------------------------------------------------------------
       program indices_severos
       implicit none
  
  
! Vari�eis miscel�eas
       integer                 :: i, j, k, nrad, irec,im,jm,ind
       integer                 :: lm,m_lm
       integer*2               :: imes   
       real                    :: arg, u2dd, v2dd, dd, ff
       integer                 :: Freq_Mod_Out
       real                    :: loni, lati, res
       character*2             :: cdia,chora             
! Vari�eis de entrada
       integer                 :: strm_opt
       integer, allocatable    :: estacao(:), hora(:), dia(:), mes(:), ano(:), nlev(:)
       integer, allocatable    :: altm(:,:), urporc(:,:), dddeg(:,:), vvkt(:,:)   
       real,    allocatable    :: presshpa(:,:), tc(:,:), tdc(:,:), wvgkg(:,:), presspa(:,:), tk(:,:)
       real,    allocatable    :: wvkgkg(:,:), u(:,:), v(:,:), theta(:,:), specif(:,:)
! Vari�eis para os �dices de tempo severo
       real,    allocatable    :: cape(:), mcape(:), li(:), ncl(:), nce(:), cin(:), el(:)
       real,    allocatable    :: ustrm(:), vstrm(:), llsrm(:), mlsrm(:), heli(:)
       real,    allocatable    :: nrv(:), dnrv(:), ieh(:), tcap(:), thetab(:) 
       real,    allocatable    :: bl(:), shr37(:), blcon(:), sup(:), ffstrm(:), ddstrm(:), dnrv2km(:)

       character(len=3)        :: cmes(12)

! Obs: tcap = lid strength; thetab = twdf on arpsplt.f (max wet bulb potential
! temperature difference)
!
! Vari�eis para c�culos em coluna
       real,    ALLOCATABLE    :: p1d(:), ht1d(:), t1d(:), tv1d(:), td1d(:)
       real,    ALLOCATABLE    :: wmr1d(:), partem(:), buoy(:), wload(:), mbuoy(:)
       real,    ALLOCATABLE    :: pbesnd(:), mbesnd(:)
       real,    ALLOCATABLE    :: u1(:), v1(:)

! Vari�el para nome do arquivo de entrada
       character(len=60)       :: inname, outname1, outname2
!
       DATA cmes/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/
! In�io do programa execut�el...
!      
! Lendo o arquivo de entrada contendo as sondagens:
!       write(*,*) 'Este programa calcula �dices de tempo severo.'
!       write(*,*) 'Ernani L. Nascimento (OUT/2006) (modificado do ARPS). '
!      write(*,*) 'Entre com o nome do arquivo de entrada (at�60 caracteres):'
       open (10,file='eta.sound.txt',status='old')
       read(10,*) im,jm,loni,lati,res, Freq_Mod_Out, m_lm, strm_opt
       read(10,*) outname1
 !      read(*,*) inname
!       write(*,*) 'Utilizar deslocamento observado da tempestade? (1=sim 5=n�)'
!       read(*,*) strm_opt
!       write(*,*) 'strm_opt=',strm_opt
       if ((strm_opt /= 1).AND.(strm_opt /= 5)) then
        write(*,*) 'Op�o n� dispon�el. Programa interrompido.'
        stop
       endif
       nrad=im*jm
       lm=m_lm+4
!       outname2='arps.'//inname(1:42)
!       write(*,*) "Numero de Pontos do Modelo ",nrad
       ALLOCATE(estacao(nrad),hora(nrad),dia(nrad),mes(nrad),ano(nrad),nlev(nrad))
       ALLOCATE(altm(lm,nrad),urporc(lm,nrad),dddeg(lm,nrad),vvkt(lm,nrad))   
       ALLOCATE(presshpa(lm,nrad),tc(lm,nrad),tdc(lm,nrad),wvgkg(lm,nrad),presspa(lm,nrad), tk(lm,nrad))
       ALLOCATE(wvkgkg(lm,nrad),u(lm,nrad),v(lm,nrad),theta(lm,nrad),specif(lm,nrad))
! Vari�eis para os �dices de tempo severo
       ALLOCATE(cape(nrad),mcape(nrad),li(nrad),ncl(nrad),nce(nrad),cin(nrad),el(nrad))
       ALLOCATE(ustrm(nrad),vstrm(nrad),llsrm(nrad),mlsrm(nrad),heli(nrad))
       ALLOCATE(nrv(nrad),dnrv(nrad),ieh(nrad),tcap(nrad),thetab(nrad)) 
       ALLOCATE(bl(nrad),shr37(nrad),blcon(nrad),sup(nrad),ffstrm(nrad),ddstrm(nrad),dnrv2km(nrad))
       ALLOCATE(u1(nrad),v1(nrad))
       ALLOCATE (p1d(lm),ht1d(lm),t1d(lm),tv1d(lm),td1d(lm))
       ALLOCATE (wmr1d(lm),partem(lm),buoy(lm),wload(lm),mbuoy(lm))
       ALLOCATE (pbesnd(lm),mbesnd(lm))
       do k=1,nrad
        read (10,*) estacao(k),nlev(k),hora(k),dia(k),mes(k),ano(k)
        if (nlev(k) /= 0) then
         do i=1,nlev(k)
          read (10,*) presshpa(i,k),altm(i,k),tc(i,k),tdc(i,k),urporc(i,k),   &
                      wvgkg(i,k),u(i,k),v(i,k),theta(i,k)
          presspa(i,k)=presshpa(i,k)*100.
          tk(i,k)=tc(i,k)+273.15
          specif(i,k)=(wvgkg(i,k)/1000.)/(1+(wvgkg(i,k)/1000.))
!          print *,'specif(1,k)= ',specif(1,k), k
         enddo
!         write(*,*) k,' ',presspa(1,k)

! As linhas abaixo geram 2 n�eis adicionais na sondagem. O ETA possui 19 n�eis, mas
! para as rodadas do ARPS tornou-se conveniente acrescentar dois n�eis a mais para
! garantir um topo de dom�io mais alto. Estas linhas podem ser convenientemente comentadas.
         altm(nlev(k)+1,k)=altm(nlev(k),k)+(altm(nlev(k),k)-altm(nlev(k)-1,k))
         altm(nlev(k)+2,k)=altm(nlev(k)+1,k)+(altm(nlev(k),k)-altm(nlev(k)-1,k))
         specif(nlev(k)+1,k)=specif(nlev(k),k)
         specif(nlev(k)+2,k)=specif(nlev(k),k)
         u(nlev(k)+1,k)=u(nlev(k),k)
         u(nlev(k)+2,k)=u(nlev(k),k)
         v(nlev(k)+1,k)=v(nlev(k),k)
         v(nlev(k)+2,k)=v(nlev(k),k)
        endif
!        print *,'presshpa(19,k)= ',presshpa(19,k),' ', k
       enddo
       close(10)
!
!       if (j==1) then
!        write(*,*) 'dddeg(1,k)= ',dddeg(j,k),' para k= ',k
!        write(*,*) 'vvmps(1,k)= ',vvmps(j,k),' para k= ',k
!       endif
!       call ddff2uv(dddeg,vvmps,uu,vv,j)
!
! Transformando graus e magnitude do vento em componentes zonal e meridional
!       arg = (dddeg(j,k) * (3.141592654/180.))
!       uu = -vvmps(j,k) * sin(arg)
!       vv = -vvmps(j,k) * cos(arg)
!       u(j,k)=uu
!       v(j,k)=vv
!       if (j==1) then
!        write(*,*) 'u(1,k)= ',u(j,k),' para k= ',k
!        write(*,*) 'v(1,k)= ',v(j,k),' para k= ',k
!       endif

!
!   Chama rotina de c�culo de �dices termodin�icos
!
!       print *,'hgt(m) antes: ', altm(1,1)
!       do k=1,nrad
!        print *,'presshpa(19,k)= ',presshpa(19,k),' ', k
!       enddo
!        write(*,*) k,' PASSEI AQUI 1.'
       call arps_be(nrad,lm,presshpa,altm,tc,tdc,wvgkg,ncl,nce,el,thetab,li,cape,   &
                    mcape,cin,tcap,p1d,ht1d,t1d,tv1d,td1d,wmr1d,partem,          &
                    buoy,wload,mbuoy,pbesnd,mbesnd,nlev,estacao,hora,dia,        &
                    mes,ano)
!
!  C�culo de �dices cinem�icos
!
!
!       write(*,*) k,' PASSEI AQUI 2'
       call calcshr(nrad,lm,altm,bl,presspa,tk,u,v,cape,shr37,ustrm,vstrm,llsrm,  &
                    mlsrm,heli,nrv,dnrv,dnrv2km,blcon,u1,v1,nlev,estacao,hora,dia,  &
                    mes,ano,strm_opt)
       call calcehi_sup(nrad,lm,nlev,presshpa,cape,heli,dnrv,ieh,sup)

!  Arquivo de sa�a: bin�io para o grads
       ind=index(outname1,' ')-1
       open(12,file=outname1,FORM='UNFORMATTED',status='unknown',ACCESS='DIRECT',   &
            recl=nrad*4)
!       write(12,*) ' EST UTC DIA M� ANO SHR37 ustrm vstrm ddstrm ffstrm helic BRN BRNSHR B2km IEH SUP'
       irec=1
       do k=1,nrad
        u2dd=ustrm(k)
        v2dd=vstrm(k)
        call uv2ddff(u2dd,v2dd,dd,ff)
        ddstrm(k)=dd
        ffstrm(k)=ff
       enddo
       irec=1
       write(12,rec=irec) (ncl(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (nce(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (el(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (thetab(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (li(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (cape(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (mcape(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (cin(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (shr37(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (ustrm(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (vstrm(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (heli(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (nrv(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (dnrv(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (dnrv2km(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (ieh(k),k=1,nrad)
       irec=irec+1
       write(12,rec=irec) (sup(k),k=1,nrad) 
   
!       if ((ieh(k)<-1.0).and.(sup(k)<-1.0)) then
!        write(*,*) 'EHI AND SUP LESS THAN -1 (k EHI SUP): ',k,ieh(k),sup(k)
!       endif
       close(12)
  
       write(*,*) 'Fim do Processamento'
!  Arquivo de sa�a: sondagem com formato para o ARPS
!       open(13,file=outname2,status='unknown')
!       write(13,304) nrad
!       do k=1,nrad
!        write(13,300) estacao(k),hora(k),dia(k),mes(k),ano(k)
!        write (13,*) ' ZSND   THSND   QVSND    USND    VSND '
!  As duas linhas abaixo acrescentam dois n�eis para theta.
!        theta(nlev(k)+1,k)=theta(nlev(k),k)+(theta(nlev(k),k)-theta(nlev(k)-1,k))
!        theta(nlev(k)+2,k)=theta(nlev(k)+1,k)+(theta(nlev(k),k)-theta(nlev(k)-1,k))
!  O loop abaixo acrescenta os dois n�eis adicionais para o ARPS.
!        do i=1,nlev(k)+2
!         write(13,303) float(altm(nlev(k)+3-i,k)),theta(nlev(k)+3-i,k), &
!     &                 specif(nlev(k)+3-i,k),u(nlev(k)+3-i,k), &
!     &                 v(nlev(k)+3-i,k)
!        enddo
!       enddo
!       close(13)
       write(cdia(1:2),'(I2)')dia(1)
       if (dia(1).lt.10) cdia(1:1)='0'
       write (chora(1:2),'(I2)')hora(1)   
       if (hora(1).lt.10) chora(1:1)='0'
       open(15,file=outname1(1:ind-4)//'.ctl',status='unknown')
       write(15,400)outname1
       write(15,401)
       write(15,402)
       write(15,403)im,loni,res
       write(15,404)jm,lati,res
       write(15,405)
       imes=mes(1)
       write(15,406)chora(1:2),cdia(1:2),cmes(imes),ano(1),Freq_Mod_Out
       write(15,407)
       write(15,408)
       write(15,409)
       write(15,410)
       write(15,411)
       write(15,412)
       write(15,413)
       write(15,414)
       write(15,415)
       write(15,416)
       write(15,417)
       write(15,418)
       write(15,419)
       write(15,420)
       write(15,421)
       write(15,422)
       write(15,423)
       write(15,424)
       write(15,425)
       close(15)
 300   format (1X,I5,3I4,I5)
 301   format (1X,I5,3I4,I5,F8.1,1X,2F8.1,2F7.1,3F8.1)
 302   format (1X,I5,3I4,I5,F8.5,2F7.3,F8.1,F7.1,F9.2,F6.1,2F7.1,2F6.1)
 303   format (1X,F9.3,F12.5,F12.8,2F12.5)
 304   format (1X,I4)
 
 400   format('DSET ^',A60)
 401   format('UNDEF -9999.           ')
 402   format('TITLE Indices de Tempestades Severas') 
 403   format('XDEF  ',I4,' LINEAR  ',f7.3,'  ',f8.4)
 404   format('YDEF  ',I4,' LINEAR  ',f7.3,'  ',f8.4)
 405   format('ZDEF   1 LEVELS 1000  ')
 406   format('TDEF   1 LINEAR ',a2,'Z',a2,A3,I4,' ',I2'hr') 
 407   format('VARS  17') 
 408   format('ncl      0  99     Nivel de Cond. Levantamento         (m)') 
 409   format('nce      0  99     Nivel de Conv. Espontanea           (m)') 
 410   format('ne       0  99     Nivel de Equilibrio                 (m)') 
 411   format('thetab   0  99     Max dif. Temp. Bulbo umido          (K)') 
 412   format('li       0  99     Lifted Index                   (degree)') 
 413   format('cape     0  99     CAPE                             (J/kg)') 
 414   format('mcape    0  99     CAPE                             (J/kg)') 
 415   format('cin      0  99     CINE                             (J/kg)')
 416   format('shr37    0  99     Total  6h Precip.                   (m)')
 417   format('ustrm    0  99     Desl. Esp. da Tempestade          (m/s)')
 418   format('vstrm    0  99     Desl. Esp. da Tempestade          (m/s)')
 419   format('heli     0  99     Helicidaade Rel. a Tempestade (m^2/s^2)')
 420   format('nrv      0  99     Numero de Richardson Volumetrico     ()')
 421   format('dnrv     0  99     Den. de Num.de R. Volumetrico (m^2/s^2)')
 422   format('dnrv2km  0  99     Den. de Num.de R. V. Cam. 2km (m^2/s^2)')
 423   format('ieh      0  99     Indice de Energia-Helicidade         ()')
 424   format('sup      0  99     Parametro de Super Celula            ()')
 425   format('ENDVARS')

!       stop
! END OF PROGRAM �DICES SEVEROS
       end
