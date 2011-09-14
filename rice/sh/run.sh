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
source sh/mkdir-species.sh

source sh/choose-species.sh
source sh/summarize-sequence-data.sh
source sh/examine-ima2-data.sh
source sh/prior-file.sh
source sh/run-ima2.sh
source sh/ima2-prepare.sh
source sh/ima2-receive.sh

source sh/xyz.sh


#####################################################################
# Read configuration file
#####################################################################
conf

#####################################################################
# Read directories
#####################################################################
SPECIESS=$(ls species|grep -v ^sim)

#####################################################################
# Menus
#####################################################################
PS3="Select what you want to do with : "
CHOICES=( init-file-system \
          choose-species \
          summarize-sequence-data \
          examine-ima2-data \
          prior-file \
          run-ima2 \
          ---IMA---\
          ima2-prepare \
          ima2-receive \
          ---TEMPLATE---\
          xyz
          )
select CHOICE in ${CHOICES[@]}; do 
  if [ "$CHOICE" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  elif [ "$CHOICE" == "init-file-system" ]; then $CHOICE; break
  elif [ "$CHOICE" == "choose-species" ]; then $CHOICE; break
  elif [ "$CHOICE" == "summarize-sequence-data" ]; then $CHOICE; break
  elif [ "$CHOICE" == "examine-ima2-data" ]; then $CHOICE; break
  elif [ "$CHOICE" == "prior-file" ]; then $CHOICE; break
  elif [ "$CHOICE" == "run-ima2" ]; then $CHOICE; break
  elif [ "$CHOICE" == "ima2-prepare" ]; then $CHOICE; break
  elif [ "$CHOICE" == "ima2-receive" ]; then $CHOICE; break
  elif [ "$CHOICE" == "xyz" ]; then $CHOICE; break
  else
    echo -e "You need to enter something\n"
    continue
  fi
done

