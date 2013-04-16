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

function project {
  PS3="Choose your project: "
  select PROJECT in ${PROJECTS[@]}; do 
    if [ "$PROJECT" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      project-repetition
      global-variable
      
      if [ ! -d "$DATADIR" ]; then
        echo -n "$FUNCNAME: creating data, bwa, and analysis directories ... "
        mkdir -p $DATADIR
        mkdir -p $BWADIR
        mkdir -p $RUNANALYSIS
        echo "done!"
      fi

      echo "  Project Data Directory: $PROJECTDATADIR"
      echo "  Data Directory: $DATADIR"

      break
    fi
  done
}
