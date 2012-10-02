#!/bin/bash
###############################################################################
# Copyright (C) 2012 Sang Chul Choi
#
# This file is part of Tomato Analysis.
# 
# Tomato Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Tomato Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Tomato Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

#####################################################################
# External shell scripts
#####################################################################
source sh/copyright.sh
source sh/batch.sh
source sh/summary.sh
source sh/figures.sh
source sh/make-soft-links.sh

#####################################################################
# Menus
#####################################################################
PS3="Select the menu : "
CHOICES=( initialize \
          test \
          batch \
          summary \
          go \
          figures \
          package \
          warranty \
          copyright \
          quit )
if [ $# -eq 1 ]; then
  CHOICEID=$(($1 - 1))
  CHOICE=${CHOICES[$CHOICEID]}
  if [ "$CHOICE" == "test" ]; then 
    RAW=smallraw; NUMBERCHROM=2; batch
  elif [ "$CHOICE" == "batch" ]; then 
    RAW=raw; NUMBERCHROM=12; $CHOICE
  elif [ "$CHOICE" == "summary" ]; then 
    RAW=raw; NUMBERCHROM=12; $CHOICE
  elif [ "$CHOICE" == "go" ]; then 
    RAW=raw; NUMBERCHROM=12; $CHOICE
  elif [ "$CHOICE" == "figures" ]; then 
    RAW=raw; NUMBERCHROM=12; $CHOICE
  else
    $CHOICE
  fi
else
  select CHOICE in ${CHOICES[@]}; do 
    if [ "$CHOICE" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    elif [ "$CHOICE" == "initialize" ]; then $CHOICE; break
    elif [ "$CHOICE" == "test" ]; then 
      RAW=smallraw; NUMBERCHROM=2; batch; break
    elif [ "$CHOICE" == "batch" ]; then 
      RAW=raw; NUMBERCHROM=12; $CHOICE; break
    elif [ "$CHOICE" == "summary" ]; then 
      RAW=raw; NUMBERCHROM=12; $CHOICE; break
    elif [ "$CHOICE" == "go" ]; then 
      RAW=raw; NUMBERCHROM=12; $CHOICE; break
    elif [ "$CHOICE" == "figures" ]; then 
      RAW=raw; NUMBERCHROM=12; $CHOICE; break
    elif [ "$CHOICE" == "package" ]; then $CHOICE; break
    elif [ "$CHOICE" == "warranty" ]; then $CHOICE; break
    elif [ "$CHOICE" == "copyright" ]; then $CHOICE; break
    elif [ "$CHOICE" == "quit" ]; then $CHOICE; break
    else
      echo -e "You need to enter something\n"
      continue
    fi
  done
fi
