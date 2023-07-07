'reinit'

'open /stornext/online16/etamc/DiegoChagas/Gera_fig_est_Eta/TEMPO/Eta_teste1_8km_BMJ_FER_2021021400.ctl.lonmod'
'open /stornext/online16/etamc/DiegoChagas/Gera_fig_est_Eta/TEMPO/sst_20210214.ctl'

'q dims'
lonl=sublin(result,2)
latl=sublin(result,3)
lati=subwrd(latl,6)
latf=subwrd(latl,8)
loni=subwrd(lonl,6)
lonf=subwrd(lonl,8)

'set display color white'
'c'

* CTL do Eta foi alterado XDEF para ficar com longitude-360
'set lat 'lati' 'latf''
'set lon 'loni' 'lonf''
say 'set lat 'lati' 'latf''
say 'set lon 'loni' 'lonf''

'q file 2'
vsstlin=sublin(result,7)
varsst=subwrd(vsstlin,1)


#faz a comparacao da sst inicial
if (264 = 000)

'/stornext/online16/etamc/DiegoChagas/RECURSOS/grads_scripts/mul.gs 2 1 1 1'
'set grads off'
*intervalo de controle entre as linhas (default 1)
_PARcint=1

'set cint '_PARcint''
'set mpdraw off'
'd 'varsst'.2(t=1)-273.15'
'!cp /stornext/online16/etamc/DiegoChagas/RECURSOS/grads_scripts/lpoly_lowres.asc /stornext/online16/etamc/DiegoChagas/Gera_fig_est_Eta/TEMPO'
if (1 = 1)
'set rgb 50 255 255 255 0'
'set gxout shaded'
'set cmin 0.9'
'set clevs 0.9'
'set ccols 50 15'
'd lsmk.1(t=1)'
else
'set shpopts 15'
'draw shp /stornext/online16/etamc/DiegoChagas/RECURSOS/SHAPES/Continentes/level1'
*'/stornext/online16/etamc/DiegoChagas/RECURSOS/grads_scripts/basemap L 15 15'
endif
'set gxout contour'

'draw title SST'

'/stornext/online16/etamc/DiegoChagas/RECURSOS/grads_scripts/mul.gs 2 1 2 1'
'set grads off'
'set cint '_PARcint''
'set mpdraw off'
'd tsfc.1(t=1)-273.15'

if (1 = 1)
'set rgb 50 255 255 255 0'
'set gxout shaded'
'set cmin 0.9'
'set clevs 0.9'
'set ccols 50 15'
'd lsmk.1(t=1)'
else
'set shpopts 15'
'draw shp /stornext/online16/etamc/DiegoChagas/RECURSOS/SHAPES/Continentes/level1'
*'/stornext/online16/etamc/DiegoChagas/RECURSOS/grads_scripts/basemap L 15 15'
endif

'draw title Eta_teste1_8km_BMJ_FER_'

'set string 1 tc 4 0'
'set strsiz 0.17'
*'draw string 5.5 8 SSTxTSFC t=1'

'gxprint /stornext/online16/etamc/DiegoChagas/Gera_fig_est_Eta/TEMPO/figuras8km/sst/sst_inicial_compare.png'
'c'

'close 2'
'close 1'

else

***********************************************
*      TSFC MEDIA PERIODO TODO                *
***********************************************
'reset'
'open /scratchout/grupos/grpeta/home/gustavo.sueiro/compila/tempo/eta/Etaoper8km/binctl_XC50/2D/Eta_teste1_8km_BMJ_FER_2021021400_2D.ctl'

say 'tsfcfin=ave(tsfc.1,t=1,t=265)-273.15'
'define tsfcfin=ave(tsfc.1,t=1,t=265)-273.15'

'set display color white'
'c'
'set grads off'
'set grid off'
'set mpdraw off'
'set gxout contour'
*intervalo de controle entre as linhas (default 2)
_PARcint=1
'set cint '_PARcint''
'set xlopts 1 3 0.13'
'set ylopts 1 3 0.13'
'd tsfcfin'

'draw title Media da Rodada \ TSFC (C) Eta_teste1_8km_BMJ_FER_' 
'gxprint /stornext/online16/etamc/DiegoChagas/Gera_fig_est_Eta/TEMPO/figuras8km/sst/tsfc_Eta_MediaDaRodada.png'
'!mv /stornext/online16/etamc/DiegoChagas/Gera_fig_est_Eta/TEMPO/figuras8km/sst/tsfc_Eta_MediaDaRodada.png /stornext/online16/etamc/DiegoChagas/Gera_fig_est_Eta/TEMPO/figuras8km/sst/tsfc_Eta_teste1_8km_BMJ_FER__MediaDaRodada.png'
endif

'quit'
