#!/bin/bash

#####################################################################
# The main root directory of the analysis.
#####################################################################
ROOTANALYSISDIR=`pwd`

#####################################################################
# External shell scripts
#####################################################################
source sh/conf.sh
source sh/global-variable.sh
source sh/init-file-system.sh

source sh/choose-species.sh
source sh/sim-fsc.sh
source sh/sim-imr.sh
source sh/sim1-fsc.sh
source sh/sim2-fsc.sh
source sh/sim1-imr.sh

source sh/convert-fsc-to-imr.sh
# source sh/xxx.sh

#####################################################################
# Read configuration file
#####################################################################
conf

#####################################################################
# Read directories
#####################################################################
SPECIESS=$(ls species|grep -v ^sim|grep -v ^README)
SIMULATIONS=$(ls species|grep sim)

#####################################################################
# Menus
#####################################################################
PS3="Select what you want to do with : "
CHOICES=( init-file-system \
          choose-species \
          ---REALDATA--- \
          prior-file
          ---TEST-SIMULATION--- \
          sim-fsc \
          sim-imr \
          ---EXTRA--- \
          convert-fsc-to-imr \
          xxx
          )
select CHOICE in ${CHOICES[@]}; do 
  if [ "$CHOICE" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  elif [ "$CHOICE" == "init-file-system" ]; then $CHOICE; break
  elif [ "$CHOICE" == "choose-species" ]; then $CHOICE; break
  elif [ "$CHOICE" == "sim-fsc" ]; then $CHOICE; break
  elif [ "$CHOICE" == "sim-imr" ]; then $CHOICE; break
  elif [ "$CHOICE" == "convert-fsc-to-imr" ]; then $CHOICE; break
  elif [ "$CHOICE" == "xxx" ]; then $CHOICE; break
  else
    echo -e "You need to enter something\n"
    continue
  fi
done
