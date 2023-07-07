#! /bin/bash -x

#########################################################################
## Script que gera figuras de diagnostico para modelo Eta versão de Tempo
#########################################################################
#OBS: É possivel alterar alguns parametros dentros dos scripts de acordo 
# com a resolucao e o dominio, pesquise dentro do script por _PAR

#VARIAVEIS


_CALDATE=${Eta_util}/caldate.3.0
#GRADS LANDSCAPE (l) ou PORTRAIT (p)
_grads=p
_gradsr=${gradsPath}/grads
_RECURSOS="${gradsPath}/RECURSOS"


#######################################################
##               CONFIGURACAO DA RODADA
#######################################################
_RunDate=${Run_Date}
_Fct=${Fct}
_exp=${LabRod}
_res=${Res}
_dirin2D=${Eta_binctl}/2D
_dirin3D=${Eta_binctl}/3D
_sstfile=${_RECURSOS}/template_SST.ctl
_freqMODEL=1
_ctl2D=${_dirin2D}/${_exp}${_RunDate}_2D.ctl
_ctl3D=${_dirin3D}/${_exp}${_RunDate}_3D.ctl
_dataf=$(${_CALDATE} ${_RunDate} + ${_Fct}h 'yyyymmddhh')

###Gera media prec ao longo to tempo (dominio)
loni=-80
lonf=-30
lati=-35
latf=10

#######################################################
##                CONFIGURACAO DA SAIDA
#######################################################
#Diretorio onde serao criadas as figuras
_dirFIG=${Eta_postfig}

#Diretorio 
_dirtrab=${Eta_run}/Figuras

#########plota shape MUNICIPIOS BRASIL (ON/OFF)########
##OFF plota mpdset brmap_hires
_PARshape=OFF

##### Scripts GS (ON/OFF) ######
_gsacum24h=ON
_gsacum6h=ON
_gsvariaveis=ON
_gssst=ON
_gsmeteog=ON


if [[ ! -f "${_ctl2D}" ]] ; then 
echo "ATENCAO:"
echo "O arquivo ${_ctl2D} não existe ou ainda não foi criado !!!!"
fi
echo "${_ctl2D} ENCONTRADO !!!"

#######################################################
##             LOOP DE HORAS
#######################################################
##ALTERAR AQUI SE QUISER EXECUTAR UM INTERVALO ESPECIFICO 
#ex: (FHRi=024 FHRf=036)

FHRi=000
FHRi=`printf "%03d" "${FHRi}"`
FHRf=`printf "%03d" "${_Fct}"`


mkdir -p  ${_dirFIG}/{meteogramas,templ_4x4,templ_prec_acum,sst,ciclo_diurno} 

_DATAT=$(${_CALDATE} ${_RunDate} + 000h 'dd/mm/yyyy hh')

if [  $_freqMODEL -eq 1 ] ; then
let tf=${_Fct}+1
#limita ate 7 dias o meteograma
if [[ ${tf} -gt 169  ]];then tf=169 ;fi
else
let tf=${_Fct}/3
let tf=${tf}+1
#limita ate 7 dias o meteograma
if [[ ${tf} -gt 56  ]];then tf=56 ;fi
fi

 
 FHR=`printf "%03d" "${FHRi}"`
 let TT=${FHR}+1
 let TT2=${FHR}+1

#VERIFICA SE EXISTE TOPO NO CTL
 _vertopo=$(grep -iF topo ${_ctl2D} | wc -l)


_templ4=${FHRi}
if [[ ${FHR} -eq 000  ]];then _geracum6=ON  ;fi 
#####################################
if [  $_freqMODEL -eq 1 ] ; then
_acum6=6
_acum6i=2
_acum6f=7
_acum6inc=6
_acum6inc2=6
else
_acum6=1
_acum6i=2
_acum6f=3
_acum6inc=2
_acum6inc2=1
fi


HPi=0
HPf=6
HPi=`printf "%03d" "${HPi}"`
HPf=`printf "%03d" "${HPf}"`
#####################################

while [ $((10#${FHR})) -le  $((10#${FHRf})) ]  ; do


#######################################################
## VERIFICA SE EXISTE O PROXIMO BIN OU GRB DA RODADA
#######################################################

_dataarq=$(${_CALDATE} ${_RunDate} + ${FHR}h 'yyyymmddhh')
_dataarq_6h=$(${_CALDATE} ${_RunDate} + ${HPf}h 'yyyymmddhh')


_proxdata=$(${_CALDATE} ${_RunDate} + ${_templ4}h 'yyyymmddhh')
_proximo2D=${_exp}${_RunDate}+${_proxdata}_2D.bin
_proximo3D=${_exp}${_RunDate}+${_proxdata}_3D.bin
_verifica=$(ls ${_dirin2D}/${_proximo2D} 2>/dev/null | wc -l )


_tentativa=0
while [[ ${_verifica} -eq "0" ]];do
echo "ARQUIVO ${_proximo2D} NAO ENCONTRADO"
echo "PROCURANDO ... ... ${_proximo2D}"
echo "TENTATIVA ${_tentativa}"
sleep 300
_verifica=$(ls ${_dirin2D}/${_proximo2D} | wc -l)
let _tentativa=${_tentativa}+1

if [  $_tentativa -eq 20 ] ; then
echo "ARQUIVO ${_proximo2D} NAO ENCONTRADO, ${_tentativa} tentativas"
exit
fi
done
echo "ARQUIVO ${_dirin2D}/${_proximo2D} ENCONTRADO"
echo "EXECUTANDO FHR ${FHR}"

########################################################
# FAZ FIGURAS ACUMULADAS DE PREC EM 24H
########################################################

RD=OFF
#acumulos da precipitacao na frequencia de 1h
if [  $_freqMODEL -eq 1 ] ; then
if [[ $((10#${FHR})) -eq 24  ]];then RD=ON HP=024 TII=2   TFI=25  ;fi
if [[ $((10#${FHR})) -eq 36  ]];then RD=ON HP=036 TII=14  TFI=37  ;fi
if [[ $((10#${FHR})) -eq 48  ]];then RD=ON HP=048 TII=26  TFI=49  ;fi
if [[ $((10#${FHR})) -eq 60  ]];then RD=ON HP=060 TII=38  TFI=61  ;fi
if [[ $((10#${FHR})) -eq 72  ]];then RD=ON HP=072 TII=50  TFI=73  ;fi
if [[ $((10#${FHR})) -eq 84  ]];then RD=ON HP=084 TII=62  TFI=85  ;fi
if [[ $((10#${FHR})) -eq 96  ]];then RD=ON HP=096 TII=74  TFI=97  ;fi
if [[ $((10#${FHR})) -eq 108 ]];then RD=ON HP=108 TII=86  TFI=109 ;fi
if [[ $((10#${FHR})) -eq 120 ]];then RD=ON HP=120 TII=98  TFI=121 ;fi
if [[ $((10#${FHR})) -eq 132 ]];then RD=ON HP=132 TII=110 TFI=133 ;fi
if [[ $((10#${FHR})) -eq 144 ]];then RD=ON HP=144 TII=122 TFI=145 ;fi
if [[ $((10#${FHR})) -eq 156 ]];then RD=ON HP=156 TII=134 TFI=157 ;fi
if [[ $((10#${FHR})) -eq 168 ]];then RD=ON HP=168 TII=146 TFI=169 ;fi
if [[ $((10#${FHR})) -eq 180 ]];then RD=ON HP=180 TII=158 TFI=181 ;fi
if [[ $((10#${FHR})) -eq 192 ]];then RD=ON HP=192 TII=170 TFI=193 ;fi
if [[ $((10#${FHR})) -eq 204 ]];then RD=ON HP=204 TII=182 TFI=205 ;fi
if [[ $((10#${FHR})) -eq 216 ]];then RD=ON HP=216 TII=194 TFI=217 ;fi
if [[ $((10#${FHR})) -eq 228 ]];then RD=ON HP=228 TII=206 TFI=229 ;fi
if [[ $((10#${FHR})) -eq 240 ]];then RD=ON HP=240 TII=218 TFI=241 ;fi
if [[ $((10#${FHR})) -eq 252 ]];then RD=ON HP=252 TII=230 TFI=252 ;fi
if [[ $((10#${FHR})) -eq 264 ]];then RD=ON HP=264 TII=242 TFI=265 ;fi

else
#acumulos da precipitacao na frequencia de 3h
if [[ $((10#${FHR})) -eq  24  ]];then RD=ON HP=024 TII=2  TFI=9   ;fi
if [[ $((10#${FHR})) -eq  36  ]];then RD=ON HP=036 TII=6  TFI=13  ;fi
if [[ $((10#${FHR})) -eq  48  ]];then RD=ON HP=048 TII=10 TFI=17  ;fi
if [[ $((10#${FHR})) -eq  60  ]];then RD=ON HP=060 TII=14 TFI=21  ;fi
if [[ $((10#${FHR})) -eq  72  ]];then RD=ON HP=072 TII=18 TFI=25  ;fi
if [[ $((10#${FHR})) -eq  84  ]];then RD=ON HP=084 TII=22 TFI=29  ;fi
if [[ $((10#${FHR})) -eq  96  ]];then RD=ON HP=096 TII=26 TFI=33  ;fi
if [[ $((10#${FHR})) -eq  108 ]];then RD=ON HP=108 TII=30 TFI=37  ;fi
if [[ $((10#${FHR})) -eq  120 ]];then RD=ON HP=120 TII=34 TFI=41  ;fi
if [[ $((10#${FHR})) -eq  132 ]];then RD=ON HP=132 TII=38 TFI=45  ;fi
if [[ $((10#${FHR})) -eq  144 ]];then RD=ON HP=144 TII=42 TFI=49  ;fi
if [[ $((10#${FHR})) -eq  156 ]];then RD=ON HP=156 TII=46 TFI=53  ;fi
if [[ $((10#${FHR})) -eq  168 ]];then RD=ON HP=168 TII=50 TFI=57  ;fi
if [[ $((10#${FHR})) -eq  180 ]];then RD=ON HP=180 TII=54 TFI=61  ;fi
if [[ $((10#${FHR})) -eq  192 ]];then RD=ON HP=192 TII=58 TFI=65  ;fi
if [[ $((10#${FHR})) -eq  204 ]];then RD=ON HP=204 TII=62 TFI=69  ;fi
if [[ $((10#${FHR})) -eq  216 ]];then RD=ON HP=216 TII=66 TFI=73  ;fi
if [[ $((10#${FHR})) -eq  228 ]];then RD=ON HP=228 TII=70 TFI=77  ;fi
if [[ $((10#${FHR})) -eq  240 ]];then RD=ON HP=240 TII=74 TFI=81  ;fi
if [[ $((10#${FHR})) -eq  252 ]];then RD=ON HP=252 TII=78 TFI=85  ;fi
if [[ $((10#${FHR})) -eq  264 ]];then RD=ON HP=264 TII=82 TFI=89  ;fi
fi

dtemp=$(${_CALDATE} ${_RunDate} + ${HP}h 'yyyymmddhh')
DT_ANT=$(${_CALDATE} ${dtemp} - 1d 'dd/mm/yyyy, hh')
DTLE=$(${_CALDATE} ${_RunDate} + ${HP}h 'dd/mm/yyyy, hh')


if [[ $RD = "ON" ]] ; then
cat << EOF > ${_dirtrab}/templ_prec_acum24h.gs
'q gxinfo'
tmp=sublin(result,2)
_tela=subwrd(tmp,4)
if ( _tela = 11 ) ; _grads=landscape ; endif
if ( _tela = 8.5 ) ; _grads=portrait ; endif

'open ${_ctl2D}'

********* + HP precip  **********
'set t 1'
'set gxout shaded'
'set display color white'
'c'
'set string 1 c 5'

if ( _grads = landscape )
'set strsiz .15 .14'
'draw string 5.6 8.25 DIMNT/CGCT/INPE -  Model ${_exp} - Daily Precipitation (mm)'
'draw string 5.6 8.0  Fct ${_DATAT}+${HP}h, from ${DT_ANT}UTC to ${DTLE}UTC'
'set parea 0.6 10.8 0.8 7.8'
barra="0.8 10.5 0.25 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'set strsiz .12 .11'
'draw string 4.3  10.65 DIMNT/CGCT/INPE -  Model ${_exp} - Daily Precipitation (mm)'
'draw string 4.3  10.4  Fct ${_DATAT}+${HP}h, from ${DT_ANT}UTC to ${DTLE}UTC'         
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

rc = RGB()

'set ccols 0 82 83 84 85 18 21 23 24  26  28 30 31 32 33 34 35 36 37 38 80 79 78 77 75  68 65 64 63 62 61 60 59 58 '
'set clevs 0.2 0.4 0.6 0.8 1 2 3 4  5  6  7  8 9 10  12  14 16  18  20  25  30 40 50  60  70  80 90 100 125 150'

ti=${TII}
tf=${TFI}

'set mproj scaled'
'set grads off'
if (${_PARshape}=ON)
'set mpdraw off'
'd sum(prec,t='ti',t='tf')*1000'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd sum(prec,t='ti',t='tf')*1000'
endif
'run ${_RECURSOS}/xcbar.gs '%barra
'printim ${_dirFIG}/templ_prec_acum/precip_24acum_${dtemp}_${HP}.png png white'

********* + HP precip CONV **********
'set t 1'
'set gxout shaded'
'set display color white'
'c'
'set string 1 c 5'

if ( _grads = landscape )
'set strsiz .15 .14'
'draw string 5.6 8.25 DIMNT/CGCT/INPE -  Model ${_exp} - Daily Precipitation CV (mm)'
'draw string 5.6 8.0  Fct ${_DATAT}+${HP}h, from ${DT_ANT}UTC to ${DTLE}UTC'
'set parea 0.6 10.8 0.8 7.8'
barra="0.8 10.5 0.25 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'set strsiz .12 .11'
'draw string 4.3  10.65 DIMNT/CGCT/INPE -  Model ${_exp} - Daily Precipitation CV (mm)'
'draw string 4.3  10.4  Fct ${_DATAT}+${HP}h, from ${DT_ANT}UTC to ${DTLE}UTC'         
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

rc = RGB()

'set ccols 0 82 83 84 85 18 21 23 24  26  28 30 31 32 33 34 35 36 37 38 80 79 78 77 75  68 65 64 63 62 61 60 59 58 '
'set clevs 0.2 0.4 0.6 0.8 1 2 3 4  5  6  7  8 9 10  12  14 16  18  20  25  30 40 50  60  70  80 90 100 125 150'

ti=${TII}
tf=${TFI}

'set mproj scaled'
'set grads off'
if (${_PARshape}=ON)
'set mpdraw off'
'd sum(prcv,t='ti',t='tf')*1000'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd sum(prcv,t='ti',t='tf')*1000'
endif
'run ${_RECURSOS}/xcbar.gs '%barra
'printim ${_dirFIG}/templ_prec_acum/precip_cnvc_24acum_${dtemp}_${HP}.png png white'


********* + HP precip GE  **********
'set t 1'
'set gxout shaded'
'set display color white'
'c'
'set string 1 c 5'

if ( _grads = landscape )
'set strsiz .15 .14'
'draw string 5.6 8.25 DIMNT/CGCT/INPE -  Model ${_exp} - Daily Precipitation GE (mm)'
'draw string 5.6 8.0  Fct ${_DATAT}+${HP}h, from ${DT_ANT}UTC to ${DTLE}UTC'
'set parea 0.6 10.8 0.8 7.8'
barra="0.8 10.5 0.25 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'set strsiz .12 .11'
'draw string 4.3  10.65 DIMNT/CGCT/INPE -  Model ${_exp} - Daily Precipitation GE (mm)'
'draw string 4.3  10.4  Fct ${_DATAT}+${HP}h, from ${DT_ANT}UTC to ${DTLE}UTC'         
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

rc = RGB()

'set ccols 0 82 83 84 85 18 21 23 24  26  28 30 31 32 33 34 35 36 37 38 80 79 78 77 75  68 65 64 63 62 61 60 59 58 '
'set clevs 0.2 0.4 0.6 0.8 1 2 3 4  5  6  7  8 9 10  12  14 16  18  20  25  30 40 50  60  70  80 90 100 125 150'

ti=${TII}
tf=${TFI}

'set mproj scaled'
'set grads off'
if (${_PARshape}=ON)
'set mpdraw off'
'd sum(prge,t='ti',t='tf')*1000'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd sum(prge,t='ti',t='tf')*1000'
endif
'run ${_RECURSOS}/xcbar.gs '%barra
'printim ${_dirFIG}/templ_prec_acum/precip_prge_24acum_${dtemp}_${HP}.png png white'


'quit'

function RGB()
'set rgb 99 140 140 140'

'set rgb  16    0    0   74'
'set rgb  17    0    0   90'
'set rgb  18    0    0  122'
'set rgb  19    0    0  150'
'set rgb  20    0    0  175'
'set rgb  21    0    0  200'
'set rgb  22    0    0  225'
'set rgb  23    0   10  255'
'set rgb  24    0   91  255'
'set rgb  25    0  120  255'
'set rgb  26  100  150  255'
'set rgb  27  120  160  255'
'set rgb  28  162  189  255'


'set rgb  30    0   80    0'
'set rgb  31    0  100    0'
'set rgb  32    0  116    0'
'set rgb  33    0  136    0'
'set rgb  34    0  150    0'
'set rgb  35    0  174    0'
'set rgb  36    0  200    0'
'set rgb  37    0  220    0'
'set rgb  38    0  255    0'

'set rgb  58  30     0    0'
'set rgb  59  60     0    0'
'set rgb  60  90     0    0'
'set rgb  61  124    0    0'
'set rgb  62  180    0    0'
'set rgb  63  200    0    0'
'set rgb  64  220   10    0'
'set rgb  65  255   20   0'
'set rgb  66  255   45   15'
'set rgb  68  255   80   30'
'set rgb  75  255  110   60'
'set rgb  77  255  150   75'
'set rgb  78  255  180   99'
'set rgb  79  255  210   99'
'set rgb  80  255  255   80'

*** GRAY ***
'set rgb 81 230 230 230'
'set rgb 82 204 204 204'
'set rgb 83 179 179 179'
'set rgb 84 153 153 153'
'set rgb 85 128 128 128'
return

EOF
if [[ ${_gsacum24h} = "ON" ]] ; then
${_gradsr} -b${_grads}c "run ${_dirtrab}/templ_prec_acum24h.gs"
fi

fi


if [[ ${_geracum6} = "ON" ]] ; then
########################################################
# FAZ FIGURAS ACUMULADAS DE PREC EM 06H
########################################################

RD=OFF

if [  $_freqMODEL -eq 1 ] ; then
if [[ $((10#${FHR})) -eq ${_acum6}  ]];then RD=ON TII=${_acum6i} TFI=${_acum6f} ;fi
else
let FHR2=${FHR}/3
let fctf=${_Fct}/6
if [ ${FHR2} -eq ${_acum6} ] && [ ${FHR2} -le ${fctf} ] ; then RD=ON TII=${_acum6i} TFI=${_acum6f};fi
fi

DT_ANT=$(${_CALDATE} ${_RunDate} + ${HPi}h 'dd/mm/yyyy, hh')
DTLE=$(${_CALDATE}   ${_RunDate} + ${HPf}h 'dd/mm/yyyy, hh')

if [[ $RD = "ON" ]] ; then
cat << EOF > ${_dirtrab}/templ_prec_acum06h.gs
'q gxinfo'
tmp=sublin(result,2)
_tela=subwrd(tmp,4)
if ( _tela = 11 ) ; _grads=landscape ; endif
if ( _tela = 8.5 ) ; _grads=portrait ; endif

'open ${_ctl2D}'

********* + HP precip  **********
'set t 1'
'set gxout shaded'
'set display color white'
'c'
'set string 1 c 5'

if ( _grads = landscape )
'set strsiz .15 .14'
'draw string 5.6 8.25 DIMNT/CGCT/INPE -  Model ${_exp} - 6-hour Precipitation (mm)'
'draw string 5.6 8.0  Fct ${_DATAT}+${HPf}h, from ${DT_ANT}UTC to ${DTLE}UTC'
'set parea 0.6 10.8 0.8 7.8'
barra="0.8 10.5 0.25 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'set strsiz .12 .11'
'draw string 4.3  10.65 DIMNT/CGCT/INPE -  Model ${_exp} - 6-hour Precipitation (mm)'
'draw string 4.3  10.4  Fct ${_DATAT}+${HPf}h, from ${DT_ANT}UTC to ${DTLE}UTC'         
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

rc = RGB()

'set ccols 0 82 83 84 85 18 21 23 24  26  28 30 31 32 33 34 35 36 37 38 80 79 78 77 75  68 65 64 63 62 61 60 59 58 '
'set clevs 0.2 0.4 0.6 0.8 1 2 3 4  5  6  7  8 9 10  12  14 16  18  20  25  30 40 50  60  70  80 90 100 125 150'

ti=${TII}
tf=${TFI}

'set mproj scaled'
'set grads off'
if (${_PARshape}=ON)
'set mpdraw off'
'd sum(prec,t='ti',t='tf')*1000'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd sum(prec,t='ti',t='tf')*1000'
endif
'run ${_RECURSOS}/xcbar.gs '%barra
'printim ${_dirFIG}/templ_prec_acum/precip_6acum_${_dataarq_6h}_${HPf}.png png white'


********* + HP precip CONCV **********
'set t 1'
'set gxout shaded'
'set display color white'
'c'
'set string 1 c 5'

if ( _grads = landscape )
'set strsiz .15 .14'
'draw string 5.6 8.25 DIMNT/CGCT/INPE -  Model ${_exp} - 6-hour Precipitation CV (mm)'
'draw string 5.6 8.0  Fct ${_DATAT}+${HPf}h, from ${DT_ANT}UTC to ${DTLE}UTC'
'set parea 0.6 10.8 0.8 7.8'
barra="0.8 10.5 0.25 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'set strsiz .12 .11'
'draw string 4.3  10.65 DIMNT/CGCT/INPE -  Model ${_exp} - 6-hour Precipitation CV (mm)'
'draw string 4.3  10.4  Fct ${_DATAT}+${HPf}h, from ${DT_ANT}UTC to ${DTLE}UTC'         
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

rc = RGB()

'set ccols 0 82 83 84 85 18 21 23 24  26  28 30 31 32 33 34 35 36 37 38 80 79 78 77 75  68 65 64 63 62 61 60 59 58 '
'set clevs 0.2 0.4 0.6 0.8 1 2 3 4  5  6  7  8 9 10  12  14 16  18  20  25  30 40 50  60  70  80 90 100 125 150'

ti=${TII}
tf=${TFI}

'set mproj scaled'
'set grads off'
if (${_PARshape}=ON)
'set mpdraw off'
'd sum(prcv,t='ti',t='tf')*1000'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd sum(prcv,t='ti',t='tf')*1000'
endif
'run ${_RECURSOS}/xcbar.gs '%barra
'printim ${_dirFIG}/templ_prec_acum/precip_prcv_6acum_${_dataarq_6h}_${HPf}.png png white'


********* + HP precip GE **********
'set t 1'
'set gxout shaded'
'set display color white'
'c'
'set string 1 c 5'

if ( _grads = landscape )
'set strsiz .15 .14'
'draw string 5.6 8.25 DIMNT/CGCT/INPE -  Model ${_exp} - 6-hour Precipitation GE (mm)'
'draw string 5.6 8.0  Fct ${_DATAT}+${HPf}h, from ${DT_ANT}UTC to ${DTLE}UTC'
'set parea 0.6 10.8 0.8 7.8'
barra="0.8 10.5 0.25 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'set strsiz .12 .11'
'draw string 4.3  10.65 DIMNT/CGCT/INPE -  Model ${_exp} - 6-hour Precipitation GE (mm)'
'draw string 4.3  10.4  Fct ${_DATAT}+${HPf}h, from ${DT_ANT}UTC to ${DTLE}UTC'         
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

rc = RGB()

'set ccols 0 82 83 84 85 18 21 23 24  26  28 30 31 32 33 34 35 36 37 38 80 79 78 77 75  68 65 64 63 62 61 60 59 58 '
'set clevs 0.2 0.4 0.6 0.8 1 2 3 4  5  6  7  8 9 10  12  14 16  18  20  25  30 40 50  60  70  80 90 100 125 150'

ti=${TII}
tf=${TFI}

'set mproj scaled'
'set grads off'
if (${_PARshape}=ON)
'set mpdraw off'
'd sum(prge,t='ti',t='tf')*1000'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd sum(prge,t='ti',t='tf')*1000'
endif
'run ${_RECURSOS}/xcbar.gs '%barra
'printim ${_dirFIG}/templ_prec_acum/precip_prge_6acum_${_dataarq_6h}_${HPf}.png png white'

'quit'

function RGB()
'set rgb 99 140 140 140'

'set rgb  16    0    0   74'
'set rgb  17    0    0   90'
'set rgb  18    0    0  122'
'set rgb  19    0    0  150'
'set rgb  20    0    0  175'
'set rgb  21    0    0  200'
'set rgb  22    0    0  225'
'set rgb  23    0   10  255'
'set rgb  24    0   91  255'
'set rgb  25    0  120  255'
'set rgb  26  100  150  255'
'set rgb  27  120  160  255'
'set rgb  28  162  189  255'

'set rgb  30    0   80    0'
'set rgb  31    0  100    0'
'set rgb  32    0  116    0'
'set rgb  33    0  136    0'
'set rgb  34    0  150    0'
'set rgb  35    0  174    0'
'set rgb  36    0  200    0'
'set rgb  37    0  220    0'
'set rgb  38    0  255    0'

'set rgb  58  30     0    0'
'set rgb  59  60     0    0'
'set rgb  60  90     0    0'
'set rgb  61  124    0    0'
'set rgb  62  180    0    0'
'set rgb  63  200    0    0'
'set rgb  64  220   10    0'
'set rgb  65  255   20   0'
'set rgb  66  255   45   15'
'set rgb  68  255   80   30'
'set rgb  75  255  110   60'
'set rgb  77  255  150   75'
'set rgb  78  255  180   99'
'set rgb  79  255  210   99'
'set rgb  80  255  255   80'

*** GRAY ***
'set rgb 81 230 230 230'
'set rgb 82 204 204 204'
'set rgb 83 179 179 179'
'set rgb 84 153 153 153'
'set rgb 85 128 128 128'
return

EOF
if [[ ${_gsacum6h} = "ON" ]] ; then
${_gradsr} -b${_grads}c "run ${_dirtrab}/templ_prec_acum06h.gs"
fi

let HPi=10#${HPi}+6
let HPf=10#${HPf}+6
HPi=`printf "%03d" "${HPi}"`
HPf=`printf "%03d" "${HPf}"`
let _acum6=${_acum6}+${_acum6inc2}
let _acum6i=${_acum6i}+${_acum6inc}
let _acum6f=${_acum6f}+${_acum6inc}
fi
fi



RD4=OFF
if [[ $((10#${FHR})) -eq ${_templ4}  ]];then RD4=ON ;fi
if [[ $RD4 = "ON" ]] ; then

gera3D=OFF
_verifica2=$(ls ${_dirin3D}/${_proximo3D} 2>/dev/null | wc -l )
if [[ ${_verifica2} -eq "1" ]]; then gera3D=ON ; fi

cat << EOF > ${_dirtrab}/templ_4x4_zoom.gs
'q gxinfo'
tmp=sublin(result,2)
_tela=subwrd(tmp,4)
if ( _tela = 11 ) ; _grads=landscape ; endif
if ( _tela = 8.5 ) ; _grads=portrait ; endif

'open ${_ctl2D}'


********************  PREC  ********************
'reset'

'run ${_RECURSOS}/rgb_right.gs'

'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0 Precipitacao (mm) - ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 10.8 0.8 7.8'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'draw string 4.4  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.4  10.4  Precipitacao (mm) - ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.1 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

'set mproj scaled'
'set grads off'
'set display color white'

'set cthick 6'

'set gxout shaded'
'set t ${TT2}'

rc = RGB_PREC()
'set ccols 0 82 83 84 85 18 21 23 24  26  28 30 31 32 33 34 35 36 37 38 80 79 78 77 75  68 65 64 63 62 61 60 59 58'
'set clevs 0.2 0.4 0.6 0.8 1 2 3 4  5  6  7  8 9 10  12  14 16  18  20  25  30 40 50  60  70  80 90 100 125 150'

if (${_PARshape}=ON)
'set mpdraw off'
'd prec.1*1000'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd prec.1*1000'
endif

'run ${_RECURSOS}/xcbar.gs '%barra
'printim ${_dirFIG}/templ_4x4/prec_${_dataarq}_${FHR}.png png  white'

********************  PRCV  ********************
'reset'

'run ${_RECURSOS}/rgb_right.gs'

'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0 Precipitacao Conv (mm) - ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 10.8 0.8 7.8'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'draw string 4.4  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.4  10.4  Precipitacao Conv (mm) - ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.1 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

'set mproj scaled'
'set grads off'
'set display color white'

'set cthick 6'

'set gxout shaded'
'set t ${TT2}'

rc = RGB_PREC()
'set ccols 0 82 83 84 85 18 21 23 24  26  28 30 31 32 33 34 35 36 37 38 80 79 78 77 75  68 65 64 63 62 61 60 59 58'
'set clevs 0.2 0.4 0.6 0.8 1 2 3 4  5  6  7  8 9 10  12  14 16  18  20  25  30 40 50  60  70  80 90 100 125 150'

if (${_PARshape}=ON)
'set mpdraw off'
'd prcv.1*1000'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd prcv.1*1000'
endif

'run ${_RECURSOS}/xcbar.gs '%barra
'printim ${_dirFIG}/templ_4x4/prcv_${_dataarq}_${FHR}.png png  white'

********************  PRGE  ********************
'reset'

'run ${_RECURSOS}/rgb_right.gs'

'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0 Precipitacao GE (mm) - ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 10.8 0.8 7.8'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'draw string 4.4  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.4  10.4  Precipitacao GE (mm) - ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.1 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

'set mproj scaled'
'set grads off'
'set display color white'

'set cthick 6'

'set gxout shaded'
'set t ${TT2}'

rc = RGB_PREC()
'set ccols 0 82 83 84 85 18 21 23 24  26  28 30 31 32 33 34 35 36 37 38 80 79 78 77 75  68 65 64 63 62 61 60 59 58'
'set clevs 0.2 0.4 0.6 0.8 1 2 3 4  5  6  7  8 9 10  12  14 16  18  20  25  30 40 50  60  70  80 90 100 125 150'

if (${_PARshape}=ON)
'set mpdraw off'
'd prge.1*1000'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd prge.1*1000'
endif

'run ${_RECURSOS}/xcbar.gs '%barra
'printim ${_dirFIG}/templ_4x4/prge_${_dataarq}_${FHR}.png png  white'

rc = RGB()
********************  PSLM  ********************
'reset'

*intervalo de controle entre as linhas (default 2)
_PARcint=2

'run ${_RECURSOS}/rgb_right.gs'

'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0 Sea Mean Level Pressure (hPa) - ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 10.8 0.8 7.8'
else
'draw string 4.4  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.4  10.4  Sea Mean Level Pressure (hPa) - ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 8.3 0.8 10.2'
endif

'set mproj scaled'
'set grads off'
'set display color white'

'set cthick 6'

'set t ${TT2}'


if (${_PARshape}=ON)
'${_RECURSOS}/color -gxout contour 978 1028 2 -kind purple->blue->green->orange->red->magenta'
'set mpdraw off'
'd pslm.1'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
rc = RGB()
'set map 1'
'set mpdset brmap_hires'
'${_RECURSOS}/color -gxout contour 978 1028 2 -kind purple->blue->green->orange->red->magenta'
'd pslm.1'
endif

'printim ${_dirFIG}/templ_4x4/pslm_${_dataarq}_${FHR}.png png  white'



********************  OCIS  ********************
'reset'

'run ${_RECURSOS}/rgb_right.gs'

'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  Ocis (W/m2) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 10.8 0.8 7.8'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'draw string 4.4  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.4  10.4  Ocis (W/m2) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.1 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

rc = RGB_TP2M()

'set mproj scaled'
'set grads off'
'set display color white'
*'set clevs 200 250 300 350 400 450 500 550 600 650 700 750 800'
*'set ccols 72 71 70 82 37 26 27 28 29 45 30 31 32 33'

'set clevs 400 500 600 700 800 900 1000 1100 1200'
'set ccols 37  26  27  28  29  45  30   31   32 33'

'set gxout shaded'
'set t ${TT2}'

if (${_PARshape}=ON)
'set mpdraw off'
'd ocis.1'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd ocis.1'
endif

'run ${_RECURSOS}/xcbar.gs '%barra

'printim ${_dirFIG}/templ_4x4/ocis_${_dataarq}_${FHR}.png png white'


********************  OLIS  ********************
'reset'

'run ${_RECURSOS}/rgb_right.gs'

'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  Olis (W/m2) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 10.8 0.8 7.8'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'draw string 4.4  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.4  10.4  Olis (W/m2) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.1 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

rc = RGB_TP2M()

'set mproj scaled'
'set grads off'
'set display color white'

'set clevs  50 100 200 300 400 500 600 700 800 900 1000 1100'
'set ccols 72 71 70 82 37 26 27 28 29 45 30 31 32 33'

'${_RECURSOS}/color -gxout shaded 0 1200 100 -kind white->deepskyblue->lightskyblue->lightyellow->khaki->gold->orange->orangered->red'

'set t ${TT2}'

if (${_PARshape}=ON)
'set mpdraw off'
'd olis.1'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd olis.1'
endif

'run ${_RECURSOS}/xcbar.gs '%barra

'printim ${_dirFIG}/templ_4x4/olis_${_dataarq}_${FHR}.png png white'


********************  ROLE  ********************
'reset'

'run ${_RECURSOS}/rgb_right.gs'

'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  Role (W/m2) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 10.8 0.8 7.8'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'draw string 4.4  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.4  10.4  Role (W/m2) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.1 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

rc = RGB_TP2M()

'set mproj scaled'
'set grads off'
'set display color white'

'set clevs 50 100 200 300 400 500 600 700 800'
'set ccols 72 71 70 82 37 26 27 28 29 45 30 31 32 33'

'${_RECURSOS}/color -gxout shaded 50 400 50 -kind mediumblue->green->yellow->orange->red'

'set t ${TT2}'

if (${_PARshape}=ON)
'set mpdraw off'
'd role.1'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd role.1'
endif

'run ${_RECURSOS}/xcbar.gs '%barra

'printim ${_dirFIG}/templ_4x4/role_${_dataarq}_${FHR}.png png white'


********************  TP2M  ********************
'reset'

'run ${_RECURSOS}/rgb_right.gs'

'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  2 Metre Temperature (C) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 10.8 0.8 7.8'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'draw string 4.4  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.4  10.4  2 Metre Temperature (C) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.1 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

rc = RGB_TP2M()

'set mproj scaled'
'set grads off'
'set display color white'
'set clevs 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38'
'set ccols 77 76 75 74 73 72 71 70 82 37 26 27 28 29 45 30 31 32 33'

*'${_RECURSOS}/color -gxout shaded 4 38 2 -kind mediumblue->deepskyblue->lightskyblue->lightyellow->khaki->orange->red'
'set gxout shaded'
'set t ${TT2}'

if (${_PARshape}=ON)
'set mpdraw off'
'd tp2m.1-273.15'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd tp2m.1-273.15'
endif

'run ${_RECURSOS}/xcbar.gs '%barra

'printim ${_dirFIG}/templ_4x4/tp2m_${_dataarq}_${FHR}.png png white'


********************  TMIN  ********************
'reset'

'run ${_RECURSOS}/rgb_right.gs'

'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  Min Temperature (C) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 10.8 0.8 7.8'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'draw string 4.4  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.4  10.4  Min Temperature (C) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.1 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

rc = RGB_TP2M()

'set mproj scaled'
'set grads off'
'set display color white'
'set clevs 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38'
'set ccols 77 76 75 74 73 72 71 70 82 37 26 27 28 29 45 30 31 32 33'

*'${_RECURSOS}/color -gxout shaded 4 38 2 -kind mediumblue->deepskyblue->lightskyblue->lightyellow->khaki->orange->red'
'set gxout shaded'
'set t ${TT2}'

if (${_PARshape}=ON)
'set mpdraw off'
'd mntp.1-273.15'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd mntp.1-273.15'
endif

'run ${_RECURSOS}/xcbar.gs '%barra

'printim ${_dirFIG}/templ_4x4/mntp_${_dataarq}_${FHR}.png png white'


********************  TMAX  ********************
'reset'

'run ${_RECURSOS}/rgb_right.gs'

'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  Max Temperature (C) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 10.8 0.8 7.8'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'draw string 4.4  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.4  10.4  Max Temperature (C) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.1 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

rc = RGB_TP2M()

'set mproj scaled'
'set grads off'
'set display color white'
'set clevs 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38'
'set ccols 77 76 75 74 73 72 71 70 82 37 26 27 28 29 45 30 31 32 33'
'set gxout shaded'

*'${_RECURSOS}/color -gxout shaded 4 38 2 -kind mediumblue->deepskyblue->lightskyblue->lightyellow->khaki->orange->red'
'set gxout shaded'
'set t ${TT2}'

if (${_PARshape}=ON)
'set mpdraw off'
'd mxtp.1-273.15'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd mxtp.1-273.15'
endif

'run ${_RECURSOS}/xcbar.gs '%barra

'printim ${_dirFIG}/templ_4x4/mxtp_${_dataarq}_${FHR}.png png white'


********************  CSSF  ********************
'reset'

'run ${_RECURSOS}/rgb_right.gs'

'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  Sensible heat (W/m2) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 10.8 0.8 7.8'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'draw string 4.4  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.4  10.4  Sensible heat (W/m2) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.1 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

rc = RGB_TP2M()

'set mproj scaled'
'set grads off'
'set display color white'
'set clevs 0 50 100 150 200 250 300 350 400 450 500 550'
'set ccols 77 76 75 74 73 72 71 70 82 37 26 27 28 29 45 30 31 32 33'
'set gxout shaded'

'${_RECURSOS}/color -gxout shaded 0 560 50 -kind grainbow'

'set t ${TT2}'

if (${_PARshape}=ON)
'set mpdraw off'
'd cssf.1*(-1)'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd cssf.1*(-1)'
endif

'run ${_RECURSOS}/xcbar.gs '%barra

'printim ${_dirFIG}/templ_4x4/cssf_${_dataarq}_${FHR}.png png white'


********************  CLSF  ********************
'reset'

'run ${_RECURSOS}/rgb_right.gs'

'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  Latent heat (W/m2) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 10.8 0.8 7.8'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'draw string 4.4  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.4  10.4  Latent heat (W/m2) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 8.3 0.8 10.2'
barra="0.1 8.4 0.25 0.5 -fw 0.1 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

rc = RGB_TP2M()

'set mproj scaled'
'set grads off'
'set display color white'
'set clevs 0 50 100 150 200 250 300 350 400 450 500 550'
'set ccols 77 76 75 74 73 72 71 70 82 37 26 27 28 29 45 30 31 32 33'
'set gxout shaded'

'${_RECURSOS}/color -gxout shaded 0 560 50 -kind grainbow'

'set t ${TT2}'

if (${_PARshape}=ON)
'set mpdraw off'
'd clsf.1*(-1)'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd clsf.1*(-1)'
endif

'run ${_RECURSOS}/xcbar.gs '%barra

'printim ${_dirFIG}/templ_4x4/clsf_${_dataarq}_${FHR}.png png white'


***************  V10M  ***************
rc = RGB()
'reset'

'set string 1 c 5'
'set strsiz .13 .12'

if ( _grads = landscape )
'set parea 0.6 10.8 0.8 7.8'
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  10m Wind Speed (m/s) -  ${_DATAT}UTC fct='${TT}-1'h'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'set parea 0.6 8.3 0.8 10.2'
'draw string 4.2  10.65  DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.2  10.4   10m Wind Speed (m/s) -  ${_DATAT}UTC fct='${TT}-1'h'
barra="0.1 8.4 0.3 0.5  -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

'set display color white'
'set mproj scaled'
'set grads off'
'set clevs 2 4 6 8 10 12 14'
'set ccols 0 46 48 50 52 54 56 58'
'set t ${TT2}'
'set gxout shaded'

if (${_PARshape}=ON)
'set mpdraw off'
'd mag(u10m.1,v10m.1)'

*VECTOR OU STREAM (defaut vector)
'set gxout vector'
*'set gxout stream'
'set cthick 1'
'set strmden -4 0.5 0.06 1'
'set ccolor 1'
'set arrowhead 0.1'
'set arrlab off'
'd skip(u10m.1,25);v10m.1'
*USE ABAIXO SE FOR STREAM
*'d u10m.1;v10m.1'

'set line 99 1 3'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else

'set map 1'
'set mpdset brmap_hires'
'd mag(u10m.1,v10m.1)'

*VECTOR OU STREAM (defaut vector)
'set gxout vector'
*'set gxout stream'
'set cthick 1'
'set strmden -4 0.5 0.06 1'
'set ccolor 1'
'set arrowhead 0.1'
'set arrlab off'
'd skip(u10m.1,20);v10m.1'
*USE ABAIXO SE FOR STREAM
*'d u10m.1;v10m.1'

endif

'run ${_RECURSOS}/xcbar.gs '%barra

'printim ${_dirFIG}/templ_4x4/vento10m_${_dataarq}_${FHR}.png png white'


***************  V100M  ***************
rc = RGB()
'reset'

'set string 1 c 5'
'set strsiz .13 .12'

if ( _grads = landscape )
'set parea 0.6 10.8 0.8 7.8'
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  100m Wind Speed (m/s) -  ${_DATAT}UTC fct='${TT}-1'h'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'set parea 0.6 8.3 0.8 10.2'
'draw string 4.2  10.65  DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.2  10.4   100m Wind Speed (m/s) -  ${_DATAT}UTC fct='${TT}-1'h'
barra="0.1 8.4 0.3 0.5  -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

'set display color white'
'set mproj scaled'
'set grads off'
'set clevs 2 4 6 8 10 12 14'
'set ccols 0 46 48 50 52 54 56 58'
'set t ${TT2}'
'set gxout shaded'

if (${_PARshape}=ON)
'set mpdraw off'
'd mag(u100.1,v100.1)'

*VECTOR OU STREAM (defaut vector)
'set gxout vector'
*'set gxout stream'
'set cthick 1'
'set strmden -4 0.5 0.06 1'
'set ccolor 1'
'set arrowhead 0.1'
'set arrlab off'
'd skip(u100.1,20);v100.1'
*USE ABAIXO SE FOR STREAM
*'d u100.1;v100.1'

'set line 99 1 3'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else

'set map 1'
'set mpdset brmap_hires'
'd mag(u100.1,v100.1)'

*VECTOR OU STREAM (defaut vector)
'set gxout vector'
*'set gxout stream'
'set cthick 1'
'set strmden -4 0.5 0.06 1'
'set ccolor 1'
'set arrowhead 0.1'
'set arrlab off'
'd skip(u100.1,20);v100.1'
*USE ABAIXO SE FOR STREAM
*'d u100.1;v100.1'

endif

'run ${_RECURSOS}/xcbar.gs '%barra

'printim ${_dirFIG}/templ_4x4/vento100m_${_dataarq}_${FHR}.png png white'


*faz essa figura se existir o arquivo de topografia
***********************************************
if (${_vertopo} = 1)
**********  V10M e OROGRAPHY  ****************
rc = RGB()
'reset'

'run ${_RECURSOS}/rgb_topo.gs'

'set string 1 c 5'
'set strsiz .13 .12'

if ( _grads = landscape )
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  Orography (m) and 10 Metre V-Wind (mph) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 10.2 0.8 7.8'
barra="10.3 10.5 0.6 7.7 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir v -line on -lc 15 -edge triangle"
else
'draw string 4.2  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.2  10.4  Orography (m) and 10 Metre V-Wind (mph) -  ${_DATAT}UTC fct='${TT}-1'h'
'set parea 0.6 7.8 0.8 10.2'
barra="7.9 8.1 1.8 9 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir v -line on -lc 15 -edge triangle"
endif

'set mproj scaled'
'set grads off'
'set display color white'
'set gxout shaded'
'set lev 950'

'set clevs 0 100 200 300 400 500 600 800 1000 1200 1400'
'set ccols 0 0 31 32 34 35 36 37 45 46 47 49'
'set grads off'

if (${_PARshape}=ON)
'set mpdraw off'
'd topo'
'run ${_RECURSOS}/xcbar.gs '%barra
rc = RGB()
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd topo'
endif



*********VENTO EM BARBELA***************
'set lev 1000'
'set gxout barb'
'set ccolor 20'
'set cthick 8'
'set t ${TT2}'
'd skip(u10m.1*2.23694,30);v10m.1*2.23694;mag(u10m.1,v10m.1)*2.23694'
****************************************
if ( _grads = landscape )
'run ${_RECURSOS}/xcbar.gs 1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle'
else
'run ${_RECURSOS}/xcbar.gs  0.1 8.4 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle'
endif

'printim ${_dirFIG}/templ_4x4/v10m_${_dataarq}_${FHR}.png png  white '
endif

if ( ${gera3D} = ON )

'open ${_ctl3D}'

***************  UVEL e VVEL  ***************
rc = RGB()
'reset'

'set string 1 c 5'
'set strsiz .13 .12'

if ( _grads = landscape )
'set parea 0.6 10.8 0.8 7.8'
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  Stream Lines/Wind Speed 250hPa(m/s) -  ${_DATAT}UTC fct='${TT}-1'h'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'set parea 0.6 8.3 0.8 10.2'
'draw string 4.2  10.65  DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.2  10.4   Stream Lines/Wind Speed 250hPa(m/s) -  ${_DATAT}UTC fct='${TT}-1'h'
barra="0.1 8.4 0.3 0.5  -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

'set display color white'
'set strmden -5'
'set arrscl 0.4 5'
'set arrowhead 0.3'
'set cthick 5'

'set mproj scaled'
'set gxout stream'
'set grads off'
*'set clevs 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58'
'set clevs 2 4 6 8 10 12 14 16 18 20 24 26 28 30 32 34'
'set lev 250'
'set t ${TT2}'

if (${_PARshape}=ON)
'set mpdraw off'
'd uvel.2;vvel.2;mag(uvel.2,vvel.2)'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd uvel.2;vvel.2;mag(uvel.2,vvel.2)'
endif

'run ${_RECURSOS}/xcbar.gs '%barra

'printim ${_dirFIG}/templ_4x4/ST2_${_dataarq}_${FHR}.png png white'


***************  UVEL e VVEL 850 ***************
rc = RGB()
'reset'

'set string 1 c 5'
'set strsiz .13 .12'

if ( _grads = landscape )
'set parea 0.6 10.8 0.8 7.8'
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  Stream Lines/Wind Speed 250hPa(m/s) -  ${_DATAT}UTC fct='${TT}-1'h'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
else
'set parea 0.6 8.3 0.8 10.2'
'draw string 4.2  10.65  DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.2  10.4   Stream Lines/Wind Speed 250hPa(m/s) -  ${_DATAT}UTC fct='${TT}-1'h'
barra="0.1 8.4 0.3 0.5  -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

'set display color white'
'set strmden -5'
'set arrscl 0.4 5'
'set arrowhead 0.3'
'set cthick 5'

'set mproj scaled'
'set gxout stream'
'set grads off'
*'set clevs 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58'
'set clevs 2 4 6 8 10 12 14 16 18 20 24 26 28 30 32 34'
'set lev 850'
'set t ${TT2}'

if (${_PARshape}=ON)
'set mpdraw off'
'd uvel.2;vvel.2;mag(uvel.2,vvel.2)'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd uvel.2;vvel.2;mag(uvel.2,vvel.2)'
endif

'run ${_RECURSOS}/xcbar.gs '%barra

'printim ${_dirFIG}/templ_4x4/ST850_${_dataarq}_${FHR}.png png white'


****************  ZGEO  *****************
rc = RGB()
'reset'

#intervalo de controle entre as linhas (default 2)
_PARcint=2

'run ${_RECURSOS}/rgb_right.gs'

'set mproj scaled'
'set grads off'
'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'set parea 0.6 10.8 0.8 7.8'
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  Geopotential Height 500hPa -  ${_DATAT}UTC fct='${TT}-1'h'
else
'set parea 0.6 8.3 0.8 10.2'
'draw string 4.2  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.2  10.4  Geopotential Height 500hPa -  ${_DATAT}UTC fct='${TT}-1'h'
endif

'set display color white'
'set cint '_PARcint''
'set cthick 7'
'set lev 500'
'set t ${TT2}'

'${_RECURSOS}/color -gxout contour 5000 6000 30 -kind purple->blue->green->orange->red->magenta'

if (${_PARshape}=ON)
'set mpdraw off'
'd zgeo.2'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'${_RECURSOS}/color -gxout contour 5000 6000 30 -kind purple->blue->green->orange->red->magenta'
'set mpdset brmap_hires'
'd zgeo.2'
endif

'printim ${_dirFIG}/templ_4x4/zgeo500_${_dataarq}_${FHR}.png png white'


****************  ZGEO  200 *****************
rc = RGB()
'reset'

#intervalo de controle entre as linhas (default 2)
_PARcint=2

'run ${_RECURSOS}/rgb_right.gs'

'set mproj scaled'
'set grads off'
'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'set parea 0.6 10.8 0.8 7.8'
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  Geopotential Height 200hPa -  ${_DATAT}UTC fct='${TT}-1'h'
else
'set parea 0.6 8.3 0.8 10.2'
'draw string 4.2  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.2  10.4  Geopotential Height 200hPa -  ${_DATAT}UTC fct='${TT}-1'h'
endif

'set display color white'
'set cint '_PARcint''
'set cthick 7'
'set lev 200'
'set t ${TT2}'

if (${_PARshape}=ON)
'set mpdraw off'
'${_RECURSOS}/color -gxout contour 11000 13000 50 -kind purple->blue->green->orange->red->magenta'
'd zgeo.2'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'${_RECURSOS}/color -gxout contour 11000 13000 50 -kind purple->blue->green->orange->red->magenta'
'set mpdset brmap_hires'
'd zgeo.2'
endif

'printim ${_dirFIG}/templ_4x4/zgeo200_${_dataarq}_${FHR}.png png white'


****************  ZGEO  850 *****************
rc = RGB()
'reset'

#intervalo de controle entre as linhas (default 2)
_PARcint=2

'run ${_RECURSOS}/rgb_right.gs'

'set mproj scaled'
'set grads off'
'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'set parea 0.6 10.8 0.8 7.8'
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0  Geopotential Height 850hPa -  ${_DATAT}UTC fct='${TT}-1'h'
else
'set parea 0.6 8.3 0.8 10.2'
'draw string 4.2  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.2  10.4  Geopotential Height 850hPa -  ${_DATAT}UTC fct='${TT}-1'h'
endif

'set display color white'
'set cint '_PARcint''
'set cthick 7'
'set lev 850'
'set t ${TT2}'

if (${_PARshape}=ON)
'set mpdraw off'
'${_RECURSOS}/color -gxout contour 1000 1600 20 -kind purple->blue->green->orange->red->magenta'
'd zgeo.2'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
say "PASSOU COLOR"
'${_RECURSOS}/color -gxout contour 1000 1600 20 -kind purple->blue->green->orange->red->magenta'
'd zgeo.2'
endif

'printim ${_dirFIG}/templ_4x4/zgeo850_${_dataarq}_${FHR}.png png white'


****************  TEMP  *****************
rc = RGB()
'reset'

'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'set parea 0.6 10.8 0.8 7.8'
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0 Absolute Temperature (C) 850hPa -  ${_DATAT}UTC fct='${TT}-1'h'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle "
else
'set parea 0.6 8.3 0.8 10.2'
'draw string 4.2  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.2  10.4 Absolute Temperature (C) 850hPa -  ${_DATAT}UTC fct='${TT}-1'h'
barra="0.1 8.4 0.25 0.5 -fw 0.1 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif



'set mproj scaled'
'set grads off'
'set display color white'
'set gxout shaded'
'set lev 850'
'set clevs 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23'
'set ccols 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 62 64 68 70'

'${_RECURSOS}/color -gxout shaded -10 26 2 -kind mediumblue->lightyellow->orange->red->darkred'

'set grads off'
'set t ${TT2}'

if (${_PARshape}=ON)
'set mpdraw off'
'd temp.2-273.15'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd temp.2-273.15'
endif

'run ${_RECURSOS}/xcbar.gs '%barra

'printim ${_dirFIG}/templ_4x4/temp850_${_dataarq}_${FHR}.png png white'


****************  TEMP 500 *****************
rc = RGB()
'reset'

'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'set parea 0.6 10.8 0.8 7.8'
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0 Absolute Temperature (C) 500hPa -  ${_DATAT}UTC fct='${TT}-1'h'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle "
else
'set parea 0.6 8.3 0.8 10.2'
'draw string 4.2  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.2  10.4 Absolute Temperature (C) 500hPa -  ${_DATAT}UTC fct='${TT}-1'h'
barra="0.1 8.4 0.25 0.5 -fw 0.1 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

'set mproj scaled'
'set grads off'
'set display color white'
'set gxout shaded'
'set lev 500'
'set clevs 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23'
'set ccols 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 62 64 68 70'

'${_RECURSOS}/color -gxout shaded -34 2 2 -kind mediumblue->lightyellow->orange->red->darkred'

'set grads off'
'set t ${TT2}'

if (${_PARshape}=ON)
'set mpdraw off'
'd temp.2-273.15'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd temp.2-273.15'
endif

'run ${_RECURSOS}/xcbar.gs '%barra

'printim ${_dirFIG}/templ_4x4/temp500_${_dataarq}_${FHR}.png png white'

****************  TEMP 200 *****************
rc = RGB()
'reset'


'set string 1 c 5'
'set strsiz .15 .14'

if ( _grads = landscape )
'set parea 0.6 10.8 0.8 7.8'
'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 5.5 8.0 Absolute Temperature (C) 200hPa -  ${_DATAT}UTC fct='${TT}-1'h'
barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle "
else
'set parea 0.6 8.3 0.8 10.2'
'draw string 4.2  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
'draw string 4.2  10.4 Absolute Temperature (C) 200hPa -  ${_DATAT}UTC fct='${TT}-1'h'
barra="0.1 8.4 0.25 0.5 -fw 0.1 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
endif

'set mproj scaled'
'set grads off'
'set display color white'
'set gxout shaded'
'set lev 200'
'set clevs -72 -70 -68 -66 -64 -62 -60 -58 -56 -54 -52 -50 -48 -46 -44 -42 -40 -38 -36 -34 -32 -30'
'set ccols 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 62 64 68 70'

'${_RECURSOS}/color -gxout shaded -64 -40 2 -kind mediumblue->lightyellow->orange->red'
'set grads off'
'set t ${TT2}'

if (${_PARshape}=ON)
'set mpdraw off'
'd temp.2-273.15'
'set line 99 1 1'
'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
else
'set map 1'
'set mpdset brmap_hires'
'd temp.2-273.15'
endif

'run ${_RECURSOS}/xcbar.gs '%barra

'printim ${_dirFIG}/templ_4x4/temp200_${_dataarq}_${FHR}.png png white'


******  AGUA DE NUVEM  *******
#'reset'

#'set string 1 c 5'
#'set strsiz .15 .14'

#if ( _grads = landscape )
#'set parea 0.6 10.8 0.8 7.8'
#'draw string 5.5 8.25 DIMNT/CGCT/INPE -  Model ${_exp}'
#'draw string 5.5 8.0 Cloud Water (mm) -  ${_DATAT}UTC fct='${TT}-1'h'
#barra="1 10 0.3 0.5 -fw 0.09 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle "
#else
#'set parea 0.6 8.3 0.8 10.2'
#'draw string 4.2  10.65 DIMNT/CGCT/INPE -  Model ${_exp}'
#'draw string 4.2  10.4 Cloud Water (mm) -  ${_DATAT}UTC fct='${TT}-1'h'
#barra="0.1 8.4 0.25 0.5 -fw 0.1 -fh 0.14 -ft 4 -fo center -dir h -line on -lc 15 -edge triangle"
#endif

#'set mproj scaled'
#'set grads off'
#'set display color white'
#'set gxout shaded'
#'set grads off'
#'set t ${TT2}'

#rc = RGB_PREC()
#'set ccols 0 82 83 84 85 18 21 23 24  26  28 30 31 32 33 34 35 36 37 38 80 79 78 77 75  68 65 64 63 62 61 60 59 58'
#'set clevs 0.2 0.4 0.6 0.8 1 2 3 4  5  6  7'

#'define agnu=sum(wtnv.2,z=1,z=22)*1000'


#if (${_PARshape}=ON)
#'set mpdraw off'
#'d agnu'
#'set line 99 1 1'
#'draw shp ${_RECURSOS}/SHAPES/BR/municipios_2010'
#else
#'set map 1'
#'set mpdset brmap_hires'
#'d agnu'
#endif

#'run ${_RECURSOS}/xcbar.gs '%barra

#'printim ${dirFIG}/templ_4x4/wtnv${_dataarq}_${FHR}.png png white'

endif
#fim do if do gera3D


#########################################################
##CICLO DIURNO
#########################################################

if (${FHR} = 000)

'reset'

'sdfopen /lustre_xc50/chou_sinchan/20210301-20210321_ERA5_SFCVARS.nc'

'q files'
'set display color white'
'c'
'set t 1 last'
'set x 1'
'set y 1'

r=1
while(r<=4)

if (r=1) ; reg='amz' ; lon1=295 ; lon2=300 ; lon1eta=lon1-360 ; lon2eta=lon2-360 ; lat1=-5  ; lat2=0   ; endif
if (r=2) ; reg='neb' ; lon1=315 ; lon2=320 ; lon1eta=lon1-360 ; lon2eta=lon2-360 ; lat1=-10 ; lat2=-5  ; endif
if (r=3) ; reg='sud' ; lon1=310 ; lon2=315 ; lon1eta=lon1-360 ; lon2eta=lon2-360 ; lat1=-24 ; lat2=-20 ; endif
if (r=4) ; reg='sul' ; lon1=305 ; lon2=310 ; lon1eta=lon1-360 ; lon2eta=lon2-360 ; lat1=-30 ; lat2=-25 ; endif

'c'
'set grads off'
'set grid off'

*'set vrange -20 400'
'set xlopts 1 6 0.16'
'set ylopts 1 6 0.16'

'set dfile 3'
'set x 1'
'set y 1'
*'set time '

'set ccolor 1'
'set cthick 8'
'set cmark 2'
'set cstyle 1'
'd tloop(aave(sshf.3(z=1)*0.00027777777778*-1,lon='lon1',lon='lon2',lat='lat1',lat='lat2'))'
'set dfile 1'
'set x 1'
'set y 1'
'set ccolor 4'
'set cthick 4'
'set cmark 6'
'set cstyle 1'
'd tloop(aave(cssf*-1,lon='lon1eta',lon='lon2eta',lat='lat1',lat='lat2'))'
'draw title DIMNT/CGCT/INPE - CSSF-'reg'-Eta_${Res}km'
'gxprint ${_dirFIG}/ciclo_diurno/ciclo_cssf_${_RunDate}_'reg'.png x1024 y768'

'c'
'set grads off'
'set grid off'

*'set vrange -20 600'
'set xlopts 1 6 0.16'
'set ylopts 1 6 0.16'

'set dfile 3'
'set ccolor 1'
'set cthick 8'
'set cmark 2'
'set cstyle 3'
'd tloop(aave(slhf.3(z=1)*0.00027777777778*-1,lon='lon1',lon='lon2',lat='lat1',lat='lat2'))'
'set dfile 1'
'set ccolor 4'
'set cthick 4'
'set cmark 6'
'set cstyle 3'
'd tloop(aave(clsf*-1,lon='lon1eta',lon='lon2eta',lat='lat1',lat='lat2'))'
'draw title DIMNT/CGCT/INPE - CLSF-'reg' -Eta_${Res}km'

'gxprint ${_dirFIG}/ciclo_diurno/ciclo_clsf_${_RunDate}_'reg'.png x1024 y768'

'c'
'set grads off'
'set grid off'

*'set vrange 10 40'
'set xlopts 1 6 0.16'
'set ylopts 1 6 0.16'

'set dfile 3'
'set ccolor 1'
'set cthick 8'
'set cmark 2'
'set cstyle 3'
'd tloop(aave((t2m.3(z=1)-273.15),lon='lon1',lon='lon2',lat='lat1',lat='lat2'))'
'set dfile 1'
'set ccolor 4'
'set cthick 4'
'set cmark 6'
'set cstyle 3'
'd tloop(aave(tp2m-273.15,lon='lon1eta',lon='lon2eta',lat='lat1',lat='lat2'))'
'draw title DIMNT/CGCT/INPE - TP2M-'reg' -Eta_${Res}km'
'gxprint ${_dirFIG}/ciclo_diurno/ciclo_tp2m_${_RunDate}_'reg'.png x1024 y768'

'c'
'set grads off'
'set grid off'

*'set vrange 0 5'
'set xlopts 1 6 0.16'
'set ylopts 1 6 0.16'

'set dfile 3'
'set ccolor 1'
'set cthick 8'
'set cmark 2'
'set cstyle 3'
'd tloop(aave((mag(u10(z=1),v10(z=1))),lon='lon1',lon='lon2',lat='lat1',lat='lat2'))'
'set dfile 1'
'set ccolor 4'
'set cthick 4'
'set cmark 6'
'set cstyle 3'
'd tloop(aave((mag(u10m,v10m)),lon='lon1eta',lon='lon2eta',lat='lat1',lat='lat2'))'
'draw title DIMNT/CGCT/INPE - Magv10m-'reg'-Eta_${Res}km'

'gxprint ${_dirFIG}/ciclo_diurno/ciclo_vento10m_${_RunDate}_'reg'.png x1024 y768'

r=r+1
endwhile

'close 3'
endif

'quit'

function RGB()

'set rgb 17  50  50  50'

* Dark Blue to Light Blue

'set rgb 18   0   0 102'
'set rgb 20   9  23 116'
'set rgb 22  19  46 130'
'set rgb 24  28  70 144'
'set rgb 26  37  93 158'
'set rgb 28  46 116 172'
'set rgb 30  56 139 185'
'set rgb 32  65 162 199'
'set rgb 34  74 185 213'
'set rgb 36  83 209 227'
'set rgb 38  93 232 241'
'set rgb 40 102 255 255'

* Light Yellow to Dark Brown

'set rgb 42 209 246 255'
'set rgb 44 255 248 228'
'set rgb 46 255 252 161'
'set rgb 48 255 235 106'
'set rgb 50 247 213  72'
'set rgb 52 238 190  43'
'set rgb 54 242 171   0'
'set rgb 56 229 137  20'
'set rgb 58 217  93   0'
'set rgb 60 197  65   0'
'set rgb 62 163  44   0'
'set rgb 64 123  29  13'
'set rgb 68 111   0  32'
'set rgb 70  67  15  20'

'set rgb 71 204 0 0'
'set rgb 72 204 51 0'
'set rgb 73 204 102 0'
'set rgb 74 204 153 0'
'set rgb 75 204 204 0'
'set rgb 76 153 204 0'
'set rgb 77 102 204 0'
'set rgb 78 51 204 0'
'set rgb 79 0 204 0'
'set rgb 80 0 204 51'
'set rgb 81 0 204 102'
'set rgb 82 0 204 153'
'set rgb 83 0 204 204'
'set rgb 84 0 153 204'
'set rgb 85 0 102 204'
'set rgb 86 0 51 204	'
'set rgb 87 0 0 204'
'set rgb 88 51 0 204'
'set rgb 89 102 0 204'
'set rgb 90 153 0 204'
'set rgb 91 204 0 204'
'set rgb 92 204 0 153'
'set rgb 93 204 0 102'
'set rgb 94 204 0 51'

'set rgb 98 255 255 255'
'set rgb 99 140 140 140'
return

function RGB_TP2M()
'set rgb 26 255 252 161'
'set rgb 27 255 235 106'
'set rgb 28 247 213  72'
'set rgb 29 238 190  43'
'set rgb 30 229 137  20'
'set rgb 31 217  93   0'
'set rgb 32 197  65   0'
'set rgb 33 163  44   0'
'set rgb 34 123  29  13'
'set rgb 35  37 176 250'
'set rgb 36   5  68   2'
'set rgb 37 255 248 228'
'set rgb 38 111   0  32'
'set rgb 39  67  15  20'
'set rgb 40 237 124  94'
'set rgb 41 103  11  75'
'set rgb 42  45  25 126'
'set rgb 43  27   8  71'
'set rgb 44  20  16  27'
'set rgb 45 242 171   0'

'set rgb 70 102 255 255'
'set rgb 71  93 232 241'
'set rgb 72  83 209 227'
'set rgb 73  74 185 213'
'set rgb 74  65 162 199'
'set rgb 75  56 139 185'
'set rgb 76  46 116 172'
'set rgb 77  37  93 158'
'set rgb 78  28  70 144'
'set rgb 79  19  46 130'
'set rgb 80   9  23 116'
'set rgb 81   0   0 102'
'set rgb 82 209 246 255'
return

function RGB_PREC()
'set rgb  99 140 140 140'
'set rgb  16    0    0   74'
'set rgb  17    0    0   90'
'set rgb  18    0    0  122'
'set rgb  19    0    0  150'
'set rgb  20    0    0  175'
'set rgb  21    0    0  200'
'set rgb  22    0    0  225'
'set rgb  23    0   10  255'
'set rgb  24    0   91  255'
'set rgb  25    0  120  255'
'set rgb  26  100  150  255'
'set rgb  27  120  160  255'
'set rgb  28  162  189  255'
'set rgb  30    0   80    0'
'set rgb  31    0  100    0'
'set rgb  32    0  116    0'
'set rgb  33    0  136    0'
'set rgb  34    0  150    0'
'set rgb  35    0  174    0'
'set rgb  36    0  200    0'
'set rgb  37    0  220    0'
'set rgb  38    0  255    0'
'set rgb  58  30     0    0'
'set rgb  59  60     0    0'
'set rgb  60  90     0    0'
'set rgb  61  124    0    0'
'set rgb  62  180    0    0'
'set rgb  63  200    0    0'
'set rgb  64  220   10    0'
'set rgb  65  255   20   0'
'set rgb  66  255   45   15'
'set rgb  68  255   80   30'
'set rgb  75  255  110   60'
'set rgb  77  255  150   75'
'set rgb  78  255  180   99'
'set rgb  79  255  210   99'
'set rgb  80  255  255   80'
'set rgb 81 230 230 230'
'set rgb 82 204 204 204'
'set rgb 83 179 179 179'
'set rgb 84 153 153 153'
'set rgb 85 128 128 128'
return

EOF
if [[ ${_gsvariaveis} = "ON" ]] ; then
${_gradsr} -b${_grads}c "run ${_dirtrab}/templ_4x4_zoom.gs"
fi

if [[ ${_gssst} = "ON" ]] ; then
#se for primeiro horario verifica sst inicial
#se for o ultimo faz a media de toda rodada
if [[ ${FHR} = "000" || ${FHR} = ${FHRf} ]] ; then
#################################################
### GERA COMPARACAO SST     ###
#################################################
edata=$(echo $_RunDate | cut -c1-8)
#modifica o template
sed  "s:DATESST:${_RunDate}:g" ${_sstfile} > ${_dirtrab}/sst_${edata}.ctl
sed -i "s:CAMINHO:${Eta_run}:g" ${_dirtrab}/sst_${edata}.ctl

_xlon2=$(cat ${_ctl2D} | grep -iF XDEF |awk '{print $4}')
_xlon=$(echo "360+${_xlon2}" | bc)

sed "s:${_xlon2}:${_xlon}:g" ${_ctl2D} > ${_dirtrab}/${_exp}_${_RunDate}.ctl.lonmod
sed -i "s:[Dd][Ss][Ee][Tt] ^:DSET ${_dirin2D}/:g" ${_dirtrab}/${_exp}_${_RunDate}.ctl.lonmod

_tip=$(ls ${_dirin}/$_proximo 2>/dev/null | tail -1  | cut -d+ -f2 | cut -d. -f2)

if [[ $_tip = "grb" ]] ; then
gribmap -i ${_exp}_${_RunDate}.ctl.lonmod
fi


#VERIFICA SE EXISTE A MASCARA
 _verlsmk=$(grep -iF lsmk ${_ctl2D} | wc -l)

tdef=$(cat ${_ctl2D} | grep -iF tdef | awk '{print $2}')

cat << EOF > ${_dirtrab}/verifica_sst.gs
'reinit'

'open ${_dirtrab}/${_exp}_${_RunDate}.ctl.lonmod'
'open ${_dirtrab}/sst_${edata}.ctl'

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
if (${FHR} = 000)

'${_RECURSOS}/mul.gs 2 1 1 1'
'set grads off'
*intervalo de controle entre as linhas (default 1)
_PARcint=1

'set cint '_PARcint''
'set mpdraw off'
'd 'varsst'.2(t=1)-273.15'
'!cp ${_RECURSOS}/lpoly_lowres.asc ${_dirtrab}'
if (${_verlsmk} = 1)
'set rgb 50 255 255 255 0'
'set gxout shaded'
'set cmin 0.9'
'set clevs 0.9'
'set ccols 50 15'
'd lsmk.1(t=1)'
else
'set shpopts 15'
'draw shp ${_RECURSOS}/SHAPES/Continentes/level1'
*'${_RECURSOS}/basemap L 15 15'
endif
'set gxout contour'

*GSM'draw title SST'

'${_RECURSOS}/mul.gs 2 1 2 1'
'set grads off'
'set cint '_PARcint''
'set mpdraw off'
'd tsfc.1(t=1)-273.15'

if (${_verlsmk} = 1)
'set rgb 50 255 255 255 0'
'set gxout shaded'
'set cmin 0.9'
'set clevs 0.9'
'set ccols 50 15'
'd lsmk.1(t=1)'
else
'set shpopts 15'
'draw shp ${_RECURSOS}/SHAPES/Continentes/level1'
*'${_RECURSOS}/basemap L 15 15'
endif

'draw title ${_exp}'

'set string 1 tc 4 0'
'set strsiz 0.17'
*'draw string 5.5 8 SSTxTSFC t=1'

'gxprint ${_dirFIG}/sst/sst_inicial_compare.png'
'c'

'close 2'
'close 1'

else

***********************************************
*      TSFC MEDIA PERIODO TODO                *
***********************************************
'reset'
'open ${_ctl2D}'

say 'tsfcfin=ave(tsfc.1,t=1,t=${tdef})-273.15'
'define tsfcfin=ave(tsfc.1,t=1,t=${tdef})-273.15'
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

'draw title Media da Rodada \ TSFC (C) ${_exp}' 
'gxprint ${_dirFIG}/sst/tsfc_Eta_MediaDaRodada.png'
'!mv ${_dirFIG}/sst/tsfc_Eta_MediaDaRodada.png ${_dirFIG}/sst/tsfc_${_exp}_MediaDaRodada.png'
endif

'quit'
EOF
if [[ $((10#${FHR})) -le 0 ]] ; then
${_gradsr} -blc "run ${_dirtrab}/verifica_sst.gs"
else
${_gradsr} -b${_grads}c "run ${_dirtrab}/verifica_sst.gs"
fi
fi
fi

let _templ4=${_templ4}+${_freqMODEL}
let TT2=${TT2}+1
fi

let TT=${TT}+1
let FHR=10#${FHR}+1
FHR=`printf "%03d" "${FHR}"`
done


#########################################################
##GERA PREC MEDIA AO LONGO DO TEMPO
#########################################################
cat << EOF > ${_dirtrab}/gera_prec_medtempo.gs
'reinit'
'open ${_ctl2D}'
'set display color white'
'c'
'set gxout stat'
'set x 1'
'set y 1'
'set t 1 last'

'd tloop(aave((prec*1000),lon=${loni},lon=${lonf},lat=${lati},lat=${latf}))'
tmp=sublin(result,9)
vmax=subwrd(tmp,6)
vmin=subwrd(tmp,5)

'set x 1'
'set y 1'
'set t 1 last'
'set ccolor 2'
'set cthick 4'
'set cstyle 1'
'set line 2 1 4'
'set strsiz 0.15 0.15'
'set vrange 0 'vmax''
'set gxout contour'
'set cint 0.2'
'set cmark 0'
'set line 2 1 4'
'd tloop(aave((prec*1000),lon=${loni},lon=${lonf},lat=${lati},lat=${latf}))'
'draw title Prec (mm) Media Tempo ${_exp}'
'gxprint ${_dirFIG}/templ_prec_acum/prec_media_tempo${_exp}.png'
'quit'
EOF
${_gradsr} -blc "run ${_dirtrab}/gera_prec_medtempo.gs"


#########################################################
##METEOGRAMA
#########################################################
cat << EOF > ${_dirtrab}/gera_meteograma.gs
'reinit'
* Abre o ctl da rodada
'open ${_ctl2D}'

linha=sublin(result,3)
lonOdominio=subwrd(linha,4)
lonLdominio=subwrd(linha,5)
linha=sublin(result,4)
latSdominio=subwrd(linha,4)
latNdominio=subwrd(linha,5)


Say 'Dominio: lat ('latSdominio' 'latNdominio') lon ('lonOdominio' 'lonLdominio')' 


* Escolhe o numero de tempos, numero de dias (73 para 3 dias)
_tf=${tf}
t=1
_c=1
while(t<=_tf)
  'set t  't
  'q dim'
  tmp=sublin(result,5)
  tmp=subwrd(tmp,6)
  _hh._c=substr(tmp,1,3)
  _ddmm._c=substr(tmp,4,5)
  t=t+12
  _c=_c+1
endwhile


'set display grey white'
'set display color'
*************************************************
* LOCALIDADES
cid.1 =  ' RIO_DE_JANEIRO, RJ, BR         ';
cid.2 =  ' ARACAJU, SE, BR                ';
cid.3 =  ' BELO_HORIZONTE, MG, BR         ';
cid.4 =  ' BOA_VISTA, RR, BR              ';
cid.5 =  ' BRASILIA, DF, BR               ';
cid.6 =  ' BUENOS_AIRES, ARGENTINA, AG    ';
cid.7 =  ' BOGOTA, COLOMBIA, CO           ';
cid.8 =  ' BELEM, PA, BR                  ';
cid.9 =  ' CACHOEIRA_PAULISTA, SP, BR     ';
cid.10 = ' CUIABA, MT, BR                 ';
cid.11 = ' CARACAS, VENEZUELA, VZ         ';
cid.12 = ' CAYENE, GUIANA_FRANCESA, GF    ';
cid.13 = ' CURITIBA, PR, BR               ';
cid.14 = ' CAMPO_GRANDE, MS, BR           ';
cid.15 = ' FERNANDO_DE_NORONHA, PE, BR    ';
cid.16 = ' FORTALEZA, CE, BR              ';
cid.17 = ' FLORIANOPOLIS, SC, BR          ';
cid.18 = ' GUYANA, GY                     ';
cid.19 = ' GOIANIA, GO, BR                ';
cid.20 = ' JOAO_PESSOA, PB, BR            ';
cid.21 = ' LA_PAZ, BOLIVIA, BO            ';
cid.22 = ' LIMA, PERU, PU                 ';
cid.23 = ' MACEIO, AL, BR                 ';
cid.24 = ' MANAUS, AM, BR                 ';
cid.25 = ' MACAPA, AP, BR                 ';
cid.26 = ' MONTEVIDEO, URUGUAY, UR        ';
cid.27 = ' NATAL, RN, BR                  ';
cid.28 = ' PORTO_ALEGRE, RS, BR           ';
cid.29 = ' PALMAS, TO, BR                 ';
cid.30 = ' PARAMARIBO, SURINAME, SR       ';
cid.31 = ' PORTO_VELHO, RO, BR            ';
cid.32 = ' PASSO_FUNDO, RS, BR            ';
cid.33 = ' ASUNCION, PARAGUAY, PY         ';
cid.34 = ' RECIFE, PE, BR                 ';
cid.35 = ' RIO_BRANCO, AC, BR             ';
cid.36 = ' SAO_PAULO, SP, BR              ';
cid.37 = ' SALVADOR, BA, BR               ';
cid.38 = ' SAO_LUIZ, MA, BR               ';
cid.39 = ' SANTIAGO, CHILE, CH            ';
cid.40 = ' TERESINA, PI, BR               ';
cid.41 = ' VITORIA, ES, BR                ';
cid.42 = ' SAO_JOSE_DOS_CAMPOS, SP, BR    ';
cid.43 = ' ITABIRA, MG, BR                ';
cid.44 = ' CAMPINAS, SP, BR               ';
cid.45 = ' CANTAREIRA, SP, BR             ';
cid.46 = ' FOZ_DO_IGUACU, PR, BR          ';
cid.47 = ' COTIA, SP, BR                  ';
cid.48 = ' DIAMANTINA, MG, BR             ';
cid.49 = ' ALTO_TIETE, SP, BR             ';
cid.50 = ' ILHEUS, BA, BR                 ';
cid.51 = ' NEPOMUCENO, MG, BR             ';
cid.52 = ' POUSO_ALEGRE, MG, BR           ';
cid.53 = ' SANTOS, SP, BR                 ';
cid.54 = ' SAO_FRANCISCO_DO_SUL, SC, BR   ';
cid.55 = ' PARANAGUA, PR, BR              ';
cid.56 = ' RIO_GRANDE, RS, BR             ';
cid.57 = ' BACABAL, MA, BR                ';
cid.58 = ' IMPERATRIZ, MA, BR             ';
cid.59 = ' CAXIAS, MA, BR                 ';
cid.60 = ' BALSAS, MA, BR                 ';
cid.61 = ' QUITO, EQUADOR, EQ             ';
cid.62 = ' LENCOIS_PAULISTA, SP, BR       ';
cid.63 = ' QUATA, SP, BR                  ';
cid.64 = ' BENTO_GONCALVES, RS, BR        ';
cid.65 = ' VACARIA, RS, BR                ';
cid.66 = ' URUGUAIANA, RS, BR             ';
cid.67 = ' SANTANA_DO_LIVRAMENTO, RS, BR  ';
cid.68 = ' BAGE, RS, BR                   ';
cid.69 = ' SAO_JOAQUIM, SC, BR            ';
cid.70 = ' TELEMACO_BORBA, PR, BR         ';
cid.71 = ' CASCAVEL, PR, BR               ';
cid.72 = ' GUARAPUAVA, PR, BR             ';
cid.73 = ' UMUARAMA, PR, BR               ';
cid.74 = ' LONDRINA, PR, BR               ';
cid.75 = ' PRESIDENTE_PRUDENTE, SP, BR    ';
cid.76 = ' ITAPEVA, SP, BR                ';
cid.77 = ' DOURADOS, MS, BR               ';
cid.78 = ' CAMPOS_DO_JORDAO, SP, BR       ';
cid.79 = ' POCOS_DE_CALDAS, MG, BR        ';
cid.80 = ' UBERABA, MG, BR                ';
cid.81 = ' UBERLANDIA, MG, BR             ';
cid.82 = ' LAGES, SC, BR                  ';
cid.83 = ' SAO_MIGUEL_D_OESTE, SC, BR     ';
cid.84 = ' SAO_JOSE_DO_RIO_PRETO, SP, BR  ';
cid.85 = ' IJUI, RS, BR                   ';
cid.86 = ' PALMAS, PR, BR                 ';
cid.87 = ' VIDEIRA, SC, BR                ';
cid.88 = ' SANTA_MARIA, RS, BR            ';
cid.89 = ' MENDONZA, ARGENTINA, AG        ';
cid.90 = ' ANDRADINA, SP, BR              ';
cid.91 = ' CANANEIA, SP, BR               ';
cid.92 = ' RIO_CLARO, SP, BR              ';
cid.93 = ' TATUI, SP, BR                  ';
cid.94 = ' UBATUBA, SP, BR                ';
cid.95 = ' VOTUPORANGA, SP, BR            ';
cid.96 = ' TRAMANDAI, RS, BR              ';
cid.97 = ' SAO_SEBASTIAO, SP, BR          ';
cid.98 = ' ANGRA_DOS_REIS-ILHA_GRANDE, RJ, BR ';
cid.99 = ' REGENCIA, ES, BR               ';
cid.100 =' GUAMARE, RN, BR                ';
cid.101 =' SOLIMOES, AM, BR               ';
cid.102 =' E&P-BCSUL, RJ, BR              ';
cid.103 =' EP-RNCE, CE, BR                ';
cid.104 =' JUIZ_DE_FORA, MG, BR           ';
cid.105 =' VICOSA, MG, BR                 ';
cid.106 =' MONTES_CLAROS, MG, BR          ';
cid.107 =' GOVERNADOR_VALADARES, MG, BR   ';
cid.108 =' OURO_PRETO, MG, BR             ';
cid.109 =' SAO_JOAO_DEL-REI, MG, BR       ';
cid.110 =' SAO_LOURENCO, MG, BR           ';
cid.111 =' ARAXA, MG, BR                  ';
cid.112 =' CACHOEIRO_ITAPEMIRIM, ES, BR   ';
cid.113 =' CAMPOS, RJ, BR                 ';
cid.114 =' MACAE, RJ, BR                  ';
cid.115 =' ITATIAIA, RJ, BR               ';
cid.116 =' RIBEIRAO_PRETO, SP, BR         ';
cid.117 =' IVAIPORA, PR, BR               ';
cid.118 =' RIO_VERDE, GO, BR              ';
cid.119 =' ITUMBIARA, GO, BR              ';
cid.120 =' ARAGARCAS, GO, BR              ';
cid.121 =' PRAIA_RICA, GO, BR             ';
cid.122 =' ACORIZAL, MT, BR               ';
cid.123 =' RONDONOPOLIS, MT, BR           ';
cid.124 =' FZ_CAIABI_ALT_FLRST, MT, BR    ';
cid.125 =' PASSO_DO_LONTRA, MS, BR        ';
cid.126 =' TEMUCO, CHILE, CH              ';
cid.127 =' PUERTO_MONTT, CHILE, CH        ';
cid.128 =' BALMACEDA, CHILE, CH           ';
cid.129 =' GUAJARA_MIRIM, RO, BR          ';
cid.130 =' VILHENA, RO, BR                ';
cid.131 =' REBIO_JARU, RO, BR             ';
cid.132 =' JI_PARANA, RO, BR              ';
cid.133 =' FAZENDA_NSENHORA, RO, BR       ';
cid.134 =' CALDAS_NOVAS, GO, BR           ';
cid.135 =' SANTAREM, PA, BR               ';
cid.136 =' CORDOBA, ARGENTINA, AG         ';
cid.137 =' CHAMICAL, ARGENTINA, AG        ';
cid.138 =' CERES, ARGENTINA, AG           ';
cid.139 =' SANTIAGO_ESTERO, ARGENTINA,AG  ';
cid.140 =' PAMPA_GUANACO, ARGENTINA, AG   ';
cid.141 =' SALTA, ARGENTINA, AG           ';
cid.142 =' M_ESTIGARRIBIA, PARAGUAI, PA   ';
cid.143 =' ROBORE, BOLIVIA, BO            ';
cid.144 =' TRINIDAD, BOLIVIA, BO          ';
cid.145 =' YAQUIBA, BOLIVIA, BO           ';
cid.146 =' P_MALDONADO, PERU, PE          ';
cid.147 =' CRUZEIRO_SUL, AC, BR           ';
cid.148 =' CAXIUANA, PA, BR               ';
cid.149 =' CHUI, RS, BR                   ';
cid.150 =' PELOTAS, RS, BR                ';
cid.151 =' BAURU, SP, BR                  ';
cid.152 =' PORTO_ACRE, AC, BR             ';
cid.153 =' SENA_MADUREIRA, AC, BR         ';
cid.154 =' SANTA_ROSA_DOS_PURUS, AC, BR   ';
cid.155 =' MAZAGAO, AP, BR                ';
cid.156 =' OIAPOQUE, AP, BR               ';
cid.157 =' FERREIRA_GOMES, AP, BR         ';
cid.158 =' SERRA_DO_NAVIO, AP, BR         ';
cid.159 =' BARCELOS, AM, BR               ';
cid.160 =' ITACOATIRA, AM, BR             ';
cid.161 =' BARRA_DA_CORDA, MA, BR         ';
cid.162 =' NOVA_XAVANTINA, MT, BR         ';
cid.163 =' SINOP, MT, BR                  ';
cid.164 =' ALTAMIRA, PA, BR               ';
cid.165 =' CAPANEMA, PA, BR               ';
cid.166 =' SANTANA_DO_ARAGUAIA, PA, BR    ';
cid.167 =' COSTA_MARQUES, RO, BR          ';
cid.168 =' JARU, RO, BR                   ';
cid.169 =' PRESIDENTE_MEDICI, RO, BR      ';
cid.170 =' CARACARAI, RR, BR              ';
cid.171 =' NORMANDIA, RR, BR              ';
cid.172 =' BONFIM, RR, BR                 ';
cid.173 =' ARAGUAINA, TO, BR              ';
cid.174 =' CONCEICAO_DO_TOCANTINS, TO,BR  ';
cid.175 =' NATIVIDADE, TO, BR             ';
cid.176 =' TOCANTINOPOLIS, TO ,BR         ';
cid.177 =' PACARAIMA, RR, BR              ';
cid.178 =' CATANDUVA, PR, BR              ';
cid.179 =' CORBELIA, PR, BR               ';
cid.180 =' MARINGA, PR, BR                ';
cid.181 =' PATO_BRANCO, PR, BR            ';
cid.182 =' TOLEDO, PR, BR                 ';
cid.183 =' XINGU_FUNAI, MT, BR            ';
cid.184 =' GAVIAO_PEIXOTO, SP, BR         ';
cid.185 =' CANDIOTA, RS, BR               ';
cid.186 =' LAGUNA,   SC, BR               ';
cid.187= 'ASUNCION  , PY, PY              ';
cid.188= 'LA_PAZ    , PY, PY              ';
cid.189= 'COOPERATIVA_FRIESLAND, PY,PY    ';
cid.190= 'SAN_JUAN_BAUTISTA, PY, PY       ';
cid.191= 'JHECHAPYRA, PY, PY              ';
cid.192= 'UNION_CURUPAYTY, PY, PY         ';
cid.193= 'RAUL_PENA, PY, PY               ';
cid.194= 'COLONIAS_UNIDAS, PY, PY         ';
cid.195= 'PINDO, PY, PY                   ';
cid.196= 'NEULAND(PIRIZAL), PY, PY        ';
cid.197= 'NEULAND(CENTRO), PY, PY         ';
cid.198= 'SOMMERFELD, PY, PY              ';
cid.199= 'FERNHEIM(CAMPO-I), PY, PY       ';
cid.200= 'CHORTITZER(PARATODO), PY, PY    ';
cid.201= 'FERNHEIM(RIBEIRA), PY, PY       ';
cid.202= 'VOLENDAM(HURON), PY, PY         ';
cid.203= 'CHORTITZER(CAMPO_LEON),PY,PY    ';
cid.204= 'COOPASAM(MINGA_PORA),PY,PY      ';
cid.205= 'COLONIA_LA_ESTANZUELA, URUGUAY, UR     ';
cid.206= 'CANELONES_LAS_BRUJAS, URUGUAY, UR      ';
cid.207= 'SALTO_GRANDES, URUGUAY, UR             ';
cid.208= 'PAYSANDU_GLENCOE, URUGUAY, UR          ';
cid.209= 'TREINTA_Y_TRES_PALO_PIQUE, URUGUAY, UR ';
cid.210= 'TACUAREMBO_LA_MAGNOLIA, URUGUAY, UR    ';
cid.211= 'RIO_GRANDE, BR, RS                     ';

*************************************************
* LATITUDES, LONGITUDES E ALTITUDES  
lonlat.1 =  '  43.16W - 22.90S      0 '
lonlat.2 =  '  37.05W - 10.92S     20  '
lonlat.3 =  '  43.44W - 19.91S    845  '
lonlat.4 =  '  60.60W - 02.73N    115  '
lonlat.5 =  '  47.80W - 15.88S   1038  '
lonlat.6 =  '  58.48W - 34.58S     20  '
lonlat.7 =  '  74.13W - 04.70N   2644  '
lonlat.8 =  '  48.50W - 01.45S      0  '
lonlat.9 =  '  45.00W - 22.67S    845  '
lonlat.10 = '  56.10W - 15.65S    186  '
lonlat.11 = '  66.88W - 10.50N    845  '
lonlat.12 = '  52.30W - 04.93N     20  '
lonlat.13 = '  49.17W - 25.52S    845  '
lonlat.14 = '  54.67W - 20.47S    520  '
lonlat.15 = '  32.42W - 03.85S      0  '
lonlat.16 = '  38.53W - 03.78S     20  '
lonlat.17 = '  48.61W - 27.79S      0  '
lonlat.18 = '  58.15W - 06.80N      0  '
lonlat.19 = '  49.25W - 16.67S    845  '
lonlat.20 = '  34.86W - 07.10S     20  '
lonlat.21 = '  68.17W - 16.50S   4188  '
lonlat.22 = '  77.12W - 12.00S      0  '
lonlat.23 = '  35.93W - 09.50S    115  '
lonlat.24 = '  60.06W - 03.03S     63  '
lonlat.25 = '  51.07W - 00.05N      0  '
lonlat.26 = '  56.20W - 34.75S     63  '
lonlat.27 = '  35.40W - 06.05S     63  '
lonlat.28 = '  51.18W - 30.00S      0  '
lonlat.29 = '  48.36W - 10.21S    186  '
lonlat.30 = '  55.20W - 05.83N      0  '
lonlat.31 = '  63.91W - 08.62S    115  '
lonlat.32 = '  52.40W - 28.26S    672  '
lonlat.33 = '  57.63W - 25.27S    115  '
lonlat.34 = '  34.88W - 08.05S      0  '
lonlat.35 = '  67.81W - 09.97S    186  '
lonlat.36 = '  46.48W - 23.71S    845  '
lonlat.37 = '  38.37W - 12.89S      0  '
lonlat.38 = '  44.15W - 02.55S      0  '
lonlat.39 = '  70.70W - 33.45S    672  '
lonlat.40 = '  42.87W - 04.85S     63  '
lonlat.41 = '  40.32W - 20.20S     63  '
lonlat.42 = '  45.89W - 23.18S    600  '
lonlat.43 = '  43.22W - 19.62S    845  '
lonlat.44 = '  47.13W - 23.00S    672  '
lonlat.45 = '  46.38W - 23.13S   1038  '
lonlat.46 = '  54.58W - 25.52S    186  '
lonlat.47 = '  47.05W - 23.72S    845  '
lonlat.48 = '  43.69W - 18.31S   1252  '
lonlat.49 = '  45.77W - 23.54S   1038  '
lonlat.50 = '  39.28W - 14.73S     63  '
lonlat.51 = '  45.23W - 21.23S    845  '
lonlat.52 = '  45.94W - 22.22S    845  '
lonlat.53 = '  46.45W - 23.95S    115  '
lonlat.54 = '  48.80W - 26.42S     63  '
lonlat.55 = '  48.55W - 25.41S      0  '
lonlat.56 = '  52.22W - 31.99S      0  '
lonlat.57 = '  44.69W - 04.22S     20  '
lonlat.58 = '  47.50W - 05.53S    186  '
lonlat.59 = '  43.35W - 04.86S    115  '
lonlat.60 = '  46.03W - 07.53S    277  '
lonlat.61 = '  78.48W - 00.15S   2992  '
lonlat.62 = '  48.80W - 22.59S    672  '
lonlat.63 = '  50.69W - 22.24S    520  '
lonlat.64 = '  51.37W - 29.22S    672  '
lonlat.65 = '  50.93W - 28.51S    845  '
lonlat.66 = '  57.09W - 29.75S     63  '
lonlat.67 = '  55.53W - 30.88S    186  '
lonlat.68 = '  54.10W - 31.33S    277  '
lonlat.69 = '  49.62W - 28.19S   1487  '
lonlat.70 = '  50.61W - 24.32S    672  '
lonlat.71 = '  53.43W - 24.95S    672  '
lonlat.72 = '  51.45W - 25.39S   1038  '
lonlat.73 = '  53.32W - 23.76S    388  '
lonlat.74 = '  51.13W - 23.33S    520  '
lonlat.75 = '  51.38W - 22.12S    520  '
lonlat.76 = '  48.87W - 23.98S    672  '
lonlat.77 = '  54.80W - 22.22S    388  '
lonlat.78 = '  45.59W - 22.74S   1487  '
lonlat.79 = '  46.55W - 21.78S   1252  '
lonlat.80 = '  48.06W - 19.76S    672  '
lonlat.81 = '  48.09W - 19.01S    845  '
lonlat.82 = '  50.32W - 27.81S    845  '
lonlat.83 = '  53.33W - 26.62S    672  '
lonlat.84 = '  49.37W - 20.82S    520  '
lonlat.85 = '  53.91W - 28.39S    277  '
lonlat.86 = '  51.99W - 26.48S   1038  '
lonlat.87 = '  51.15W - 27.00S    845  '
lonlat.88 = '  53.45W - 29.75S    115  '
lonlat.89 = '  68.79W - 32.85S    845  '
lonlat.90 = '  51.25W - 20.97S    388  '
lonlat.91 = '  47.98W - 25.14S      0  '
lonlat.92 = '  47.56W - 22.41S    672  '
lonlat.93 = '  47.85W - 23.35S    672  '
lonlat.94 = '  45.07W - 23.43S      0  '
lonlat.95 = '  49.73W - 20.57S    520  '
lonlat.96 = '  50.13W - 29.98S      0  '
lonlat.97 = '  45.32W - 23.65S      0  '
lonlat.98 = '  44.32W - 23.02S      0  '
lonlat.99 = '  39.96W - 19.55S     20  '
lonlat.100= '  36.32W - 05.10S      0  '
lonlat.101= '  63.17W - 03.93S     20  '
lonlat.102= '  40.02W - 22.42S      0  '
lonlat.103= '  39.00W - 03.00S      0  '
lonlat.104= '  43.35W - 21.77S    845  '
lonlat.105= '  42.83W - 20.74S    672  '
lonlat.106= '  43.65W - 16.67S    672  '
lonlat.107= '  41.94W - 18.85S    277  '
lonlat.108= '  43.50W - 20.29S   1252  '
lonlat.109= '  44.26W - 21.13S    845  '
lonlat.110= '  45.05W - 22.11S   1038  '
lonlat.111= '  46.94W - 19.59S   1038  '
lonlat.112= '  41.11W - 20.84S    277  '
lonlat.113= '  41.33W - 21.75S     20  '
lonlat.114= '  41.80W - 22.35S      0  '
lonlat.115= '  44.56W - 22.49S   1252  '
lonlat.116= '  47.73W - 21.25S    672  '
lonlat.117= '  51.68W - 24.24S    520  '
lonlat.118= '  50.92W - 17.79S    672  '
lonlat.119= '  49.22W - 18.42S    672  '
lonlat.120= '  52.23W - 15.90S    277  '
lonlat.121= '  55.63W - 14.83S    277  '
lonlat.122= '  56.36W - 15.20S    186  '
lonlat.123= '  54.63W - 16.47S    277  '
lonlat.124= '  56.35W - 09.96S    186  '
lonlat.125= '  57.00W - 19.56S     63  '
lonlat.126= '  72.63W - 38.75S     63  '
lonlat.127= '  73.03W - 41.33S     63  '
lonlat.128= '  71.70W - 45.92S    845  '
lonlat.129= '  65.34W - 10.78S    115  '
lonlat.130= '  60.13W - 12.73S    672  '
lonlat.131= '  61.93W - 10.08S    186  '
lonlat.132= '  61.95W - 10.88S    186  '
lonlat.133= '  62.36W - 10.76S    277  '
lonlat.134= '  48.60W - 17.72S    672  '
lonlat.135= '  54.62W - 02.50S     63  '
lonlat.136= '  64.22W - 31.33S    388  '
lonlat.137= '  66.22W - 30.52S    520  '
lonlat.138= '  61.95W - 29.88S     63  '
lonlat.139= '  64.30W - 27.77S    186  '
lonlat.140= '  61.05W - 26.78S    115  '
lonlat.141= '  65.48W - 24.85S   1252  '
lonlat.142= '  60.73W - 22.02S    186  '
lonlat.143= '  59.63W - 18.43S    277  '
lonlat.144= '  64.95W - 14.85S    186  '
lonlat.145= '  64.27W - 22.03S   1252  '
lonlat.146= '  69.20W - 12.63S    186  '
lonlat.147= '  72.67W - 07.63S    186  '
lonlat.148= '  51.51W - 01.71S      0  '
lonlat.149= '  53.58W - 33.63S     20  '
lonlat.150= '  52.35W - 31.75S      0  '
lonlat.151= '  49.07W - 22.32S    520  '
lonlat.152= '  67.53W - 09.59S    186  '
lonlat.153= '  68.65W - 09.06S    186  '
lonlat.154= '  70.49W - 09.43S    186  '
lonlat.155= '  51.29W - 00.11S      0  '
lonlat.156= '  51.83W - 03.84N     63  '
lonlat.157= '  51.17W - 00.86N     20  '
lonlat.158= '  52.00W - 00.89N    115  '
lonlat.159= '  62.92W - 00.97S     20  '
lonlat.160= '  58.44W - 03.14S     20  '
lonlat.161= '  45.24W - 05.50S    115  '
lonlat.162= '  52.35W - 14.67S    388  '
lonlat.163= '  55.50W - 11.86S    388  '
lonlat.164= '  52.20W - 03.20S    115  '
lonlat.165= '  47.18W - 01.19S     20  '
lonlat.166= '  50.10W - 09.29S    186  '
lonlat.167= '  64.22W - 12.44S    115  '
lonlat.168= '  62.46W - 10.44S    186  '
lonlat.169= '  61.90W - 11.17S    186  '
lonlat.170= '  61.12W - 01.81N     63  '
lonlat.171= '  59.62W - 03.88N    115  '
lonlat.172= '  59.83W - 03.36N    115  '
lonlat.173= '  48.20W - 07.19S    277  '
lonlat.174= '  47.29W - 12.22S    672  '
lonlat.175= '  47.72W - 11.71S    388  '
lonlat.176= '  47.41W - 06.32S    186  '
lonlat.177= '  61.14W - 04.43N    845  '
lonlat.178= '  53.18W - 25.30S    672  '
lonlat.179= '  53.30W - 24.80S    672  '
lonlat.180= '  51.97W - 23.42S    388  '
lonlat.181= '  52.68W - 26.12S    672  '
lonlat.182= '  55.72W - 24.78S    277  '
lonlat.183= '  53.07W - 10.78S    186  '
lonlat.184= '  48.49W - 21.84S    609  '
lonlat.185= '  53.69W - 31.49S    204  '
lonlat.186= '  48.78W - 28.48S     47  '
lonlat.187= '  57.57W - 25.26S   87.00  '
lonlat.188= '  55.89W - 26.99S  204.00  '
lonlat.189= '  56.78W - 24.60S   87.00  '
lonlat.190= '  57.08W - 26.69S   87.00  '
lonlat.191= '  55.12W - 26.73S  204.00  '
lonlat.192= '  54.97W - 25.84S  276.00  '
lonlat.193= '  55.27W - 26.15S  354.00  '
lonlat.194= '  55.49W - 26.57S  354.00  '
lonlat.195= '  55.47W - 25.89S  354.00  '
lonlat.196= '  60.64W - 22.96S  140.00  '
lonlat.197= '  60.11W - 22.67S  140.00  '
lonlat.198= '  55.68W - 25.23S  276.00  '
lonlat.199= '  60.47W - 22.09S  140.00  '
lonlat.200= '  59.60W - 23.22S   87.00  '
lonlat.201= '  60.41W - 22.45S  140.00  '
lonlat.202= '  56.83W - 24.26S   87.00  '
lonlat.203= '  59.54W - 22.55S  140.00  '
lonlat.204= '  54.95W - 24.73S  354.00  '
lonlat.205= '  57.68W - 34.33S    47.0  '
lonlat.206= '  56.33W - 34.66S    47.0  '
lonlat.207= '  57.89W - 31.27S    47.0  '
lonlat.208= '  57.13W - 32.00S   140.0  '
lonlat.209= '  54.48W - 33.25S    47.0  '
lonlat.210= '  55.81W - 31.70S   140.0  '
lonlat.211= '  52.09W - 32.03S       0  '

*************************************************
* SIGLAS 
sigla.1=  'XRIO'
sigla.2=  'XARA'
sigla.3=  'XBEH'
sigla.4=  'XBVT'
sigla.5=  'XBRA'
sigla.6=  'XBAS'
sigla.7=  'XBGO'
sigla.8=  'XBEL'
sigla.9=  'XCPA'
sigla.10= 'XCUI'
sigla.11= 'XCAR'
sigla.12= 'XCAY'
sigla.13= 'XCUR'
sigla.14= 'XCPG'
sigla.15= 'XFEN'
sigla.16= 'XFTZ'
sigla.17= 'XFLO'
sigla.18= 'XGEO'
sigla.19= 'XGOI'
sigla.20= 'XJPS'
sigla.21= 'XLPZ'
sigla.22= 'XLIM'
sigla.23= 'XMCO'
sigla.24= 'XMNS'
sigla.25= 'XMAC'
sigla.26= 'XMON'
sigla.27= 'XNAT'
sigla.28= 'XPOA'
sigla.29= 'XPAL'
sigla.30= 'XPAR'
sigla.31= 'XPVH'
sigla.32= 'XPSF'
sigla.33= 'XASN'
sigla.34= 'XREC'
sigla.35= 'XRBR'
sigla.36= 'XSAP'
sigla.37= 'XSLV'
sigla.38= 'XSLU'
sigla.39= 'XSAN'
sigla.40= 'XTER'
sigla.41= 'XVTO'
sigla.42= 'XSJC'
sigla.43= 'XITB'
sigla.44= 'XCPN'
sigla.45= 'XCTR'
sigla.46= 'XFOZ'
sigla.47= 'XCOT'
sigla.48= 'XDIA'
sigla.49= 'XTIE'
sigla.50= 'XILH'
sigla.51= 'XNEP'
sigla.52= 'XPSA'
sigla.53= 'XSNT'
sigla.54= 'XSFS'
sigla.55= 'XPRG'
sigla.56= 'XRGR'
sigla.57= 'XBAC'
sigla.58= 'XIMP'
sigla.59= 'XCAX'
sigla.60= 'XBAL'
sigla.61= 'XQTO'
sigla.62= 'XLEN'
sigla.63= 'XQTA'
sigla.64= 'XG01'
sigla.65= 'XG02'
sigla.66= 'XG03'
sigla.67= 'XG04'
sigla.68= 'XG05'
sigla.69= 'XG06'
sigla.70= 'XG07'
sigla.71= 'XG08'
sigla.72= 'XG09'
sigla.73= 'XG10'
sigla.74= 'XG11'
sigla.75= 'XG12'
sigla.76= 'XG13'
sigla.77= 'XG14'
sigla.78= 'XG15'
sigla.79= 'XG16'
sigla.80= 'XG17'
sigla.81= 'XG18'
sigla.82= 'XG19'
sigla.83= 'XG20'
sigla.84= 'XG21'
sigla.85= 'XG22'
sigla.86= 'XG23'
sigla.87= 'XG24'
sigla.88= 'XG25'
sigla.89= 'XMDZ'
sigla.90= 'XAND'
sigla.91= 'XCAN'
sigla.92= 'XRCL'
sigla.93= 'XTTI'
sigla.94= 'XUBA'
sigla.95= 'XVOT'
sigla.96= 'XTMD'
sigla.97= 'XSSB'
sigla.98= 'XIGR'
sigla.99= 'XREG'
sigla.100='XGUR'
sigla.101='XSOL'
sigla.102='XP01'
sigla.103='XP02'
sigla.104='XJZF'
sigla.105='XVIC'
sigla.106='XMCL'
sigla.107='XGVL'
sigla.108='XOUP'
sigla.109='XSJR'
sigla.110='XSLO'
sigla.111='XARX'
sigla.112='XCIT'
sigla.113='XCMP'
sigla.114='XMCE'
sigla.115='XITT'
sigla.116='XRBP'
sigla.117='XIVA'
sigla.118='XRVD'
sigla.119='XITU'
sigla.120='XARG'
sigla.121='XPRR'
sigla.122='XACO'
sigla.123='XRON'
sigla.124='XCAI'
sigla.125='XPLO'
sigla.126='XTEM'
sigla.127='XPMO'
sigla.128='XBLM'
sigla.129='XGJM'
sigla.130='XVLN'
sigla.131='XRBJ'
sigla.132='XJPR'
sigla.133='XFNS'
sigla.134='XCLN'
sigla.135='XSTR'
sigla.136='XCOR'
sigla.137='XCHA'
sigla.138='XCER'
sigla.139='XSDE'
sigla.140='XPAM'
sigla.141='XSAL'
sigla.142='XMES'
sigla.143='XROB'
sigla.144='XTRI'
sigla.145='XYAQ'
sigla.146='XPMA'
sigla.147='XCRU'
sigla.148='XCXN'
sigla.149='XCHI'
sigla.150='XPEL'
sigla.151='XBAU'
sigla.152='XPAC'
sigla.153='XSMD'
sigla.154='XSRP'
sigla.155='XMZG'
sigla.156='XOPW'
sigla.157='XFRG'
sigla.158='XSRN'
sigla.159='XBCL'
sigla.160='XITC'
sigla.161='XBCD'
sigla.162='XNXV'
sigla.163='XSNP'
sigla.164='XATM'
sigla.165='XCAP'
sigla.166='XSAR'
sigla.167='XCMQ'
sigla.168='XJAR'
sigla.169='XPMC'
sigla.170='XCRI'
sigla.171='XNMD'
sigla.172='XBFM'
sigla.173='XARU'
sigla.174='XCTO'
sigla.175='XNTV'
sigla.176='XTCT'
sigla.177='XPCR'
sigla.178='XCAT'
sigla.179='XCBL'
sigla.180='XMRG'
sigla.181='XPTB'
sigla.182='XTLD'
sigla.183='XFUN'
sigla.184='XGPX'
sigla.185='XOTA'
sigla.186='XLAG'
sigla.187='XASU'
sigla.188='XPAZ'
sigla.189='XFRI'
sigla.190='XSJB'
sigla.191='XJHE'
sigla.192='XUNC'
sigla.193='XRLP'
sigla.194='XCOU'
sigla.195='XPND'
sigla.196='XNEU'
sigla.197='XNE2'
sigla.198='XSOM'
sigla.199='XFER'
sigla.200='XCH1'
sigla.201='XFE2'
sigla.202='XCH2'
sigla.203='XCH3'
sigla.204='XCOO'
sigla.205='XCLE'
sigla.206='XCNB'
sigla.207='XSTG'
sigla.208='XPYG'
sigla.209='XTYT'
sigla.210='XTLM'
sigla.211='XCRG'

'set missconn on'
*'set xlab off'

*************************************************
* Loop das estacoes 
i=1
while (i<=211)
'c'
************************ COM ARQUIVO PROF *******
*'set x 'i
*************************************************

************* COM BIN CTL ***********************
 latlon=lonlat.i
 latit=subwrd(latlon,3)
 sig=substr(latit,6,6)
 lat1=substr(latit,1,5)

 longi=subwrd(latlon,1)
 sig2=substr(longi,6,6)
 lon1=substr(longi,1,5)

 if(sig2=W & sig=S) 
  'set lat -'lat1  
  'set lon -'lon1
   lat2=-lat1
   lon2=-lon1 
 endif

 if(sig2=W & sig=N) 
  'set lat  'lat1 
  'set lon -'lon1
   lat2=lat1
   lon2=-lon1  
  endif

* Escolhendo as estacoes pelo dominio 
 
if (lat2 > latSdominio & lat2 < latNdominio & lon2 > lonOdominio & lon2 < lonLdominio )

Say 'Fazendo: ' lonlat.i '  ' sigla i  

***********************************************
* COLOCA OS TITULOS

  'q dims'
  l = sublin(result,2)
  l2 = subwrd(l,6)
  'set string 8 l 6'
  'set strsiz .13 .14'
  'draw string 4.3 10.7 'cid.i
  'set strsiz .13 .14'
  'draw string 4.0 10.48 'lonlat.i'm'
  'set t 1'
  'q dims'
  t = sublin(result,5)
  t2 = subwrd(t,6)
  t3 = substr(t2,1,12)
  dia = substr(t3,4,2)
  mes = substr(t3,6,3)
  ano = substr(t3,11,2)
  hora = substr(t3,1,3)
  ano=ano+2000
  ano2=ano
  ano=ano', '
  'draw string 1.0 10.70  Hourly from 'dia mes ano hora
  'set string 4 c 6'
  'set strsiz .15 .14'
  'draw string 4.25  10.9 Eta Model - ${_DATAT}UTC'
  'set t 1  '_tf
  'set vpage 0 8.5 10.2 11'
  'set grads off'

* UMIDADE RELATIVA

   'set vpage 0 8.5 5.25 7.20'
   'set parea 0.7 8.22 0.35 1.6'
*   'set lev 850'
   'set grads off'
   'set gxout line'
   'set vrange 0 100'
   'set ylint 20'
   'set cmark 0'
   'set ccolor 4'
   'd ur2m'
   'set string 2 l 6'
   'set strsiz .12 .13'
   'draw string 1 1.7 2-m Relative Humidity (%)'
*   escala()


*************************************************

* PRECIPITACAO

   'set vpage 0 8.5 8.75 10.70'
   'set parea 0.7 8.22 0.35 1.6'
   'set string 2 l 6'
   'set strsiz .12 .13'
   'draw string 1 1.7 Precipitation (mm/h)'
   'set grads off'
   'set gxout bar'
   'set bargap 40'
   'set barbase 0'
   'set baropts filled'
   escalaY(prec)
   'set vrange 0 8'
   'set ylint 1'
   'set grads off'
   'set ccolor 4'
   'd  prec*1000'
*   'set ccolor 3'
*   'd  prcv*1000'
   
*   escala()
''

*************************************************

* OCIS (RADIACAO)

   'set vpage 0 8.5 1.75 3.70'
   'set parea 0.7 8.22 0.35 1.6'
   'set grads off'
   'set gxout line'
   'set string 2 l 6'
   'set strsiz .12 .13'
   'draw string 1 1.7 Downward Shortwave Radiation Flux at the Surface (W/m2)'
   escalaY(ocis)
   'set vrange '_min ' '_max
   'set cmark 0'
   'set ylint '_int
   'set ccolor 4'
   'set grads off'
   'd ocis '
*   escala()


*************************************************

* PRESSAO AO NIVEL DO MAR

*   'set vpage 0 8.5 1.75 3.70'
*   'set parea 0.7 8.22 0.35 1.6'
*   'set grads off'
*   'set gxout line'
*   'set string 2 l 6'
*   'set strsiz .12 .13'
*   'draw string 1 1.7 Mean Sea Level Pressure (hPa)'
*   escalaY(pslm)
*   'set vrange '_min ' '_max
*   'set cmark 0'
*   'set ylint '_int
*   'set ccolor 4'
*   'set grads off'
*   'd pslm '
*   escala()


*************************************************

*  VENTO DO ANEMOMETRO
   'set vpage 0 8.5 3.50 5.45'
   'set parea 0.7 8.22 0.35 1.6'
   'set grads off'
   'set string 2 l 6'
   'set strsiz .12 .13'
   'draw string 1 1.7 10-m Wind (m/s)'
*   'set lev 1002 998'
   'set ylab off'
   'set digsize 0.10'
   'set grads off'
   'd const(u10m,0);u10m;v10m'   
   'set lev 1000'
   'define aux=mag(u10m,v10m)'
   escalaY(aux)
   'set vrange '_min ' '_max
   'set ylint  ' _int
   'set ylab on'
   'set cmark 0'
   'set ccolor 4'
   'd mag(u10m,v10m)'
*   escala()


*************************************************

*  TEMPERATURA DO abrigo
   'set vpage 0 8.5 7.00 8.95'
   'set parea 0.7 8.22 0.35 1.6'
   'set grads off'
   'set gxout line'
   'set ylab on'
   'set string 2 l 6'
   'set strsiz .12 .13'
   'draw string 1 1.7 2-m Temperature (C)'
   'define aux=tp2m-273.16'
   escalaY(aux)
   'set vrange '_min ' '_max
   'set cmark 0'
   'set ccolor 4'
   'set ylint  ' _int
   'd tp2m-273.16'
*   escala()

*************************************************

*  NEBULOSIDADE
   'set vpage 0 8.5 0.0 1.95'
   'set parea 0.7 8.22 0.35 1.6'
   'set gxout bar'
   'set baropts outline'
   'set barbase 0'
   'set grads off'
   'set strsiz .12 .13'
   'set vrange 0 100'
   'set ylint 20'
   'set string 2 l 6'
   'draw string 1 1.7 Cloud Cover (%)'
   'set bargap 50'
   'set ccolor 8'
   'set parea 0.7 8.22 0.35 1.6'
   'set annot 0'
   'd  hinv*100 '
   'set line 8'
   'draw rec  5.84 1.65  5.9 1.8 '
   'set string 8 l 6'
   'draw string 5.96 1.7  high clouds'
   'set parea 0.7 8.22 0.35 1.6'
   'set bargap 30'
   'set ccolor 3'
   'd mdnv*100'
   'set line 3'
   'draw rec  4.17 1.65  4.25 1.8 '
   'set string 3 l 6'
   'draw string 4.31 1.7 middle clouds'
   'set parea 0.7 8.22 0.35 1.6'
   'set bargap 10'
   'set ccolor 4'
   'set annot 1'
   'd lwnv*100'
   'set string 4 l 6'
   'set line 4'
   'draw  rec  2.8 1.65  2.9 1.8 '
   'draw string 2.96 1.7 low clouds'
   'set parea 0.7 8.22 0.35 1.6'

  'printim  ${_dirFIG}/meteogramas/'sigla.i'_${_RunDate}.png png x612 y792 white'
'c'
'set vpage off'

else

Say 'Ponto fora do dominio: ' lonlat.i '  ' sigla i  

endif

i=i+1
endwhile
  
  
'quit'
return
end

function escalaY(var)
   'set t 1'
   'd  min('var',t=1,t='_tf')'
   min=sublin(result,2)
   min=subwrd(min,4)
   _min=math_nint(min)

   if (min!=0) ; min=min-1; endif

   'd  max('var',t=1,t='_tf')'
   max=sublin(result,2)
   max=subwrd(max,4)
   _max=math_nint(max) + 1
   'set t 1 '_tf
   _int=math_nint((_max-_min)/5)
return

function escala()
   k=1
   xi=0.703
   dx=0.682855
   'set string 1 c 5'
   'set strsiz 0.09'
   while(k<=_c)
    'draw string 'xi' 0.21   ' _ddmm.k
    'draw string 'xi' 0.08   ' _hh.k
    'set strsiz 0.03'
    'draw string 'xi' 0.34  |'
    'set line 15 3'
    'draw line  'xi' 0.34 ' xi ' 1.6 '
    'set strsiz 0.09'
    xi=xi+dx
    k=k+2
   endwhile
return
EOF
if [[ ${_gsmeteog} = "ON" ]] ; then
${_gradsr} -bpc "run ${_dirtrab}/gera_meteograma.gs"
fi









