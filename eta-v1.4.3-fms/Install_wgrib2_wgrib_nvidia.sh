#!/bin/bash
#

function clrd(){
    local color=$1;
    local backgnd=$2;
    case $color in
        0) tput setaf 0 ;; # black
        1) tput setaf 1 ;; # red
        2) tput setaf 2 ;; # green
        3) tput setaf 3 ;; # yellow
        4) tput setaf 4 ;; # blue
        5) tput setaf 5 ;; # magenta
        6) tput setaf 6 ;; # cyan
        7) tput setaf 7 ;; # white
        *) tput sgr0;
    esac
    case $backgnd in
        0) tput setb 0 ;; # black
        1) tput setb 1 ;; # red
        2) tput setb 2 ;; # green
        3) tput setb 3 ;; # yellow
        4) tput setb 4 ;; # blue
        5) tput setb 5 ;; # magenta
        6) tput setb 6 ;; # cyan
        7) tput setb 7 ;; # white
    esac
}


#
###########################################################################################################
###########################################################################################################
#
tput clear

divider=================================
divider=$divider$divider$divider$divider$divider$divider$divider
width=`tput cols`
#
Dir_root=`pwd`
ltest="true"
while [ ${ltest} = "true" ] ;do
  printf "%s" "    wgrib   $(clrd 2)Yes$(clrd)/$(clrd 1)No$(clrd): "
  read ans
  ansUP=`echo ${ans}| tr '[:lower:]' '[:upper:]'`
  if [ ${ansUP} = "YES" ] ; then
    Wgrib_Install="Yes"
  else
    Wgrib_Install="No"
  fi
  printf "%s" " "
  printf "%s" "   wgrib2  $(clrd 2)Yes$(clrd)/$(clrd 1)No$(clrd): "
  read ans
  ansUP=`echo ${ans}| tr '[:lower:]' '[:upper:]'`
  if [ ${ansUP} = "YES" ] ; then
    Wgrib2_Install="Yes"
  else
    Wgrib2_Install="No"
  fi
  printf "%s" " "
  printf "%s" "nvidia sdk $(clrd 2)Yes$(clrd)/$(clrd 1)No$(clrd): "
  read ans
  ansUP=`echo ${ans}| tr '[:lower:]' '[:upper:]'`
  if [ ${ansUP} = "YES" ] ; then
    Nvidia_Install="Yes"
  else
    Nvidia_Install="No"
  fi
  printf "%s" " "
  printf "%s\n" "   wgrib : "${Wgrib_Install}
  printf "%s\n" "   wgrib2 : "${Wgrib2_Install}
  printf "%s\n" "Nvidia sdk: "${Nvidia_Install}
  printf "%s" "Install Packages?  [$(clrd 2)Yes$(clrd)/$(clrd 1)No$(clrd)] "
  read ans
  ansUP=`echo ${ans}| tr '[:lower:]' '[:upper:]'`
  if [ ${ansUP} = "YES" ] ; then
    break 
  fi
done

if [ ! -d ~/Softwares ] ; then
  mkdir ~/Softwares
fi
###########################################################################################################
###########################################################################################################
tput clear
if [ ${Wgrib_Install} = "Yes" ] ; then
##################################
text="wgrib INSTALL"
##################################
lng=`expr length "${text}"`
printf "%$width.${width}s\n\v" "$(clrd 2)$divider$(clrd)"
printf "%`expr \( $width / 2 \) + \( $lng / 2 \)`s\n\v" "$(clrd 2)${text}$(clrd)"
printf "%$width.${width}s\n" "$(clrd 2)$divider$(clrd)"
printf "$(clrd)"

# Instalacao wgrib
cd ~/Softwares
mkdir wgrib
cd wgrib/
wget -c -t 10 --no-check-certificate https://ftp.cpc.ncep.noaa.gov/wd51we/wgrib/wgrib.tar
tar -xvf wgrib.tar
make

cp ./wgrib ${Dir_root}/Eta_support_data/util/.
cp ./wgrib ${Dir_root}/datain/util/.
#rm -Rf ~/Softwares/wgrib
fi

###########################################################################################################
###########################################################################################################
tput clear
# Instalacao wgrib2

if [ ${Wgrib2_Install} = "Yes" ] ; then
cd ~/Softwares
##################################
text="wgrib2 INSTALL"
###############################${Dir_root}###
lng=`expr length "${text}"`
printf "%$width.${width}s\n\v" "$(clrd 2)${divider}$(clrd)"
printf "%`expr \( $width / 2 \) + \( $lng / 2 \)`s\n\v" "$(clrd 2)${text}$(clrd)"
printf "%$width.${width}s\n" "$(clrd 2)${divider}$(clrd)"
printf "$(clrd)"
wget -c -t 10 http://ftp1.cptec.inpe.br/pesquisa/grpeta/VII-WorkEta/util/wgrib2.tgz
#Descompactacao
tar -zxvf wgrib2.tgz
#Entrar no diretorio grib2
cd grib2
#executar os comandos abaixo
export CC=gcc
export FC=gfortran
make
cp ~/Softwares/grib2/wgrib2/wgrib2 ${Dir_root}/Eta_support_data/util/.
cp ~/Softwares/grib2/wgrib2/wgrib2 ${Dir_root}/datain/util/.
#rm -f ~/Softwares/wgrib2.tgz
#rm -Rf ~/Softwares/grib2
fi
#

###########################################################################################################
###########################################################################################################
tput clear

#Instalacao nvidia sdk
if [ ${Nvidia_Install} = "Yes" ] ; then

##################################
text="nvidia SDK INSTALL "
text1="The nvidia sdk install will ask for root permission"
##################################
lng=`expr length "${text}"`
lng1=`expr length "${text1}"`
printf "%$width.${width}s\n\v" "$(clrd 2)$divider$(clrd)"
printf "%`expr \( $width / 2 \) + \( $lng / 2 \)`s\n\v" "$(clrd 2)${text}$(clrd)"
printf "%`expr \( $width / 2 \) + \( $lng1 / 2 \)`s\n\v" "$(clrd 2)${text1}$(clrd)"
printf "%$width.${width}s\n" "$(clrd 2)$divider$(clrd)"
printf "$(clrd)"

cd ~/Softwares
wget -c -t 10 https://developer.download.nvidia.com/hpc-sdk/22.7/nvhpc_2022_227_Linux_x86_64_cuda_11.7.tar.gz
#Descompactar e instalar o fortran
tar -zxvf nvhpc_2022_227_Linux_x86_64_cuda_11.7.tar.gz
cd nvhpc_2022_227_Linux_x86_64_cuda_11.7
sudo ./install
rm -f nvhpc_2022_227_Linux_x86_64_cuda_11.7.tar.gz
rm -Rf nvhpc_2022_227_Linux_x86_64_cuda_11.7

printf "%s" "Inclui os caminhos no .bashrc"
printf "%s"  "export MANPATH=\$MANPATH:/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/man"
printf "%s"  "export PATH=.:/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/bin:\$PATH"
printf "%s"  "export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/comm_libs/mpi/bin:\$PATH"
printf "%s"  "export MANPATH=\$MANPATH:/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/comm_libs/mpi/man"
printf "%s\n" " "
fi
if [ ! -s ${Dir_root}/datain/util/wgrib ] ; then
  printf "%s\n" "$(clrd 1)Executar o script novamente e instalar wgrib$(clrd)"
fi
if [ ! -s ${Dir_root}/datain/util/wgrib2 ] ; then
  printf "%s\n" "$(clrd 1)Executar o script novamente e instalar wgrib2$(clrd)"
fi
if [ ! -s /opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/bin/nvfortran ] ; then
  printf "%s\n" "$(clrd 1)Executar o script novamente e instalar nvidia sdk$(clrd)"
fi

