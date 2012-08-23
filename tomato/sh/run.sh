#!/bin/bash
###############################################################################
# Copyright (C) 2012 Sang Chul Choi
#
# This file is part of RNAseq Analysis.
# 
# RNAseq Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RNAseq Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RNAseq Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

#####################################################################
# External shell scripts
#####################################################################
source sh/copyright.sh
# source sh/utility.sh
# source sh/conf.sh
# source sh/global-variable.sh
# source sh/read-species.sh
# source sh/init-file-system.sh
source sh/find-position.sh

#####################################################################
# Read configuration file
#####################################################################
# conf

#####################################################################
# Read directories
#####################################################################
# SPECIESS=$(ls $ROOTANALYSISDIR/species|grep -v ^sim)

#####################################################################
# Menus
#####################################################################
PS3="Select the menu : "
CHOICES=( find-position \
          ---BATCH---\
          batch \
          warranty \
          copyright \
          quit )
select CHOICE in ${CHOICES[@]}; do 
  if [ "$CHOICE" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  elif [ "$CHOICE" == "find-position" ]; then $CHOICE; break
  elif [ "$CHOICE" == "batch" ]; then $CHOICE; break
  elif [ "$CHOICE" == "warranty" ]; then $CHOICE; break
  elif [ "$CHOICE" == "copyright" ]; then $CHOICE; break
  elif [ "$CHOICE" == "quit" ]; then break
  else
    echo -e "You need to enter something\n"
    continue
  fi
done
