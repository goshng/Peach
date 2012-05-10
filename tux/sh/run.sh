#!/bin/bash
###############################################################################
# Copyright (C) 2012 Sang Chul Choi
#
# This file is part of tux.
# 
# tux is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tux is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Mauve Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

ROOTDIR=`pwd`

#####################################################################
# External shell scripts
#####################################################################
source sh/noweb-weave.sh
source sh/noweb-tangle.sh

#####################################################################
# Read configuration file
#####################################################################
# conf

#####################################################################
# Read directories
#####################################################################
# SPECIESS=$(ls species|grep -v ^sim)

#####################################################################
# Menus
#####################################################################
CHOICES=( noweb-weave
          noweb-tangle \
          warranty \
          copyright \
          quit )
if [ "$1" != "" ]; then
  NUMCHOICE=$1
  NUMCHOICE=$((NUMCHOICE - 1))
  ${CHOICES[$NUMCHOICE]}
  exit
fi

PS3="Select the menu : "
select CHOICE in ${CHOICES[@]}; do 
  if [ "$CHOICE" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  elif [ "$CHOICE" == "noweb-weave" ]; then $CHOICE; break
  elif [ "$CHOICE" == "noweb-tangle" ]; then $CHOICE; break
  elif [ "$CHOICE" == "warranty" ]; then $CHOICE; break
  elif [ "$CHOICE" == "copyright" ]; then $CHOICE; break
  elif [ "$CHOICE" == "quit" ]; then break
  else
    echo -e "You need to enter something\n"
    continue
  fi
done

