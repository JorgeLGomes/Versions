#!/bin/bash
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
#
tput clear
divider=================================
divider=$divider$divider$divider$divider$divider$divider$divider
width=`tput cols`

##################################
text="DEFINE ENVIRONMENT INSTALL"
##################################
lng=`expr length "${text}"`
printf "%$width.${width}s\n\v" "$divider"
printf "%`expr \( $width / 2 \) + \( $lng / 2 \)`s\n\v" "${text}"
printf "%$width.${width}s\n" "$divider"
text="OPTIONS:"
printf "%s\n" "$(clrd ${clr1})${text} ${msg1}$(clrd)"
Machines=`ls -d ./machines/*`
NumberOfMachines=`echo ${Machines}| wc -w`
for i in `basename -a ${Machines}`
do
 printf "%s\n" "   $(clrd ${clr1})${i} ${msg1}$(clrd)"
done
ltest="true"
while [ ${ltest} = "true" ] ;do
printf "%s" "Choose: "
read Machine
export ${Machine}
Occurrence=`echo ${Machines}|grep ${Machine}|wc -w`
if ((${Occurrence}>0));then
   ltest="false"
else
   echo "ERROR: '${Machine}': unknown input. "
   echo "Note: Case sensitive option!"
fi
done
cd libraries
make_all_libs ${Machine}
pwd
cd ../datain/dprep/install
build_dprep ${Machine}
cd ../../scripts/BAMCLIMT126
configure
cd ../CFS 
configure ${Machine}
cd ../BESMT062
configure
cd ../gfs2gr0.25
configure ${Machine}
cd ../ERA5 
configure ${Machine}

