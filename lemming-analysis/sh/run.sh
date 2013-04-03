#!/bin/bash
###############################################################################
# Copyright (C) 2013 Sang Chul Choi
#
# This file is part of Lemming Analysis.
# 
# Lemming Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Lemming Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Lemming Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

#####################################################################
# External shell scripts
#####################################################################
source sh/copyright.sh
source sh/conf.sh
source sh/utility.sh
source sh/global-variable.sh
source sh/read-species.sh
source sh/init-file-system.sh
source sh/project.sh

#####################################################################
# Read configuration file
#####################################################################
conf

#####################################################################
# Read directories
#####################################################################
PROJECTS=$(ls $PROJECTDIR|grep -v ^sim)

#####################################################################
# Menus
#####################################################################
PS3="Select the menu : "
CHOICES=( project \
          batch \
          warranty \
          copyright \
          quit )
if [ $# -eq 1 ]; then
  CHOICEID=$(($1 - 1))
  CHOICE=${CHOICES[$CHOICEID]}
  $CHOICE
else
  select CHOICE in ${CHOICES[@]}; do 
    if [ "$CHOICE" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else
      $CHOICE
      break
    fi
  done
fi
