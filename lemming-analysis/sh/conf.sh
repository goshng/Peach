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

function conf {
  CONFFILE=data/conf/default
  PROJECTNAME=$(eval echo $(grep ^PROJECTNAME\= $CONFFILE | cut -d"=" -f2))
  PROJECTDIR=$(eval echo $(grep ^PROJECTDIR\= $CONFFILE | cut -d"=" -f2))
  OUTPUTDIR=$(eval echo $(grep ^OUTPUTDIR\= $CONFFILE | cut -d"=" -f2))

  ROOTANALYSISDIR=$(eval echo $(grep ^ROOTANALYSISDIR\= $CONFFILE | cut -d"=" -f2))

  # The main base directory contains all the subdirectories.
  ROUTPUTDIR=$(eval echo $(grep ^ROUTPUTDIR\= $CONFFILE | cut -d"=" -f2))
}
