* ERNANI DE LIMA NASCIMENTO 24/OUT/2005
*
* Este "script" gera, a partir do GRADS, um arquivo (binário) de saída contendo
* os perfis verticais do modelo ETA-CPTEC para diversas variáveis para a região
* em torno da cidade de Indaiatuba no dia 24 de Maio de 2005. São um total
* de 3200 pontos.
*
* Note que, como eu executo este script em um PC, eu gero o binário em
* formato LITTLE ENDIAN (se executar este script em um máquina de grande
* porte, então utilizar a opção BIG ENDIAN; i.e., set fwrite -be).
*
* O objetivo final é ler este arquivo binário com o código gera_sond_2.f90 (executável
* chamado de gera2) que gera as sondagens no formato para futura execução do
* código índices_severos5.f90.
*
reinit
open LL40GANL40km2006111300+2006111400.ctl
set gxout fwrite
set fwrite -le sonds_2400_2418.le.bin

set lat -35.0 -15.0
set lon 298.2 323.4

*Surface variables

define logthe2m=log(tp2m)+0.286*log(1000/pslc)
define theta2m=exp(logthe2m)
define dp2mC=dp2m-273.15
define tp2mC=tp2m-273.15
define est=6.112*exp((17.67*tp2mC)/(tp2mC+243.5))
define estd=6.112*exp((17.67*dp2mC)/(dp2mC+243.5))
define umrl2m=100*(estd/est)
define wvmr2m=(0.62197*estd)/(pslc-estd)
d topo
d pslc
d tp2mC
d dp2mc
d umrl2m
d wvmr2m*1000
d u10m
d v10m
d theta2m

*Pressure level variables
set lev 1000
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 925
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 900
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 850
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 800
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 750
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 700
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 650
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 600
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 550
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 500
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 450
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 400
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 350
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 300
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 250
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 200
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 150
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

set lev 100
define logthe=log(temp)+0.286*log(1000/lev)
define theta=exp(logthe)
define tempC=temp-273.15
define wvmr=umes/(1-umes)
define vap=(lev*wvmr)/(0.62197+wvmr)
define td=(243.5*log(vap/6.112))/(17.67-log(vap/6.112))
d zgeo
d lev
d tempC
d td
d umrl
d wvmr*1000
d uvel
d vvel
d theta

disable fwrite
