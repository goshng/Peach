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

function global-variable {
  # PROJECT must be set before the call of this bash function. 
  PROJECTFILE=$PROJECTDIR/$PROJECT

  # The subdirectory output contains directories named after the project file.
  # The output project directory would contain all of the results from the
  # analysis.
  BASEDIR=$OUTPUTDIR/$PROJECT
  BASERUNANALYSIS=$BASEDIR/run-analysis
  NUMBERDIR=$BASEDIR/$REPETITION
  DATADIR=$NUMBERDIR/data
  BWADIR=$NUMBERDIR/bwa
  RUNANALYSIS=$NUMBERDIR/run-analysis

  # The cluster has almost the same file system. I used to used Samba client to
  # use the file system of the cluster. This stopped working. I did not know the
  # reason, which I did not want to know. Since then, I use scp command.
  # Note that the cluster base directory does not contain run-analysis. The
  # basic analysis is done in the local machine.
  CAC_BASEDIR=$CAC_OUTPUTDIR/$PROJECT
  CAC_NUMBERDIR=$CAC_OUTPUTDIR/$PROJECT/$REPETITION
  CAC_DATADIR=$CAC_NUMBERDIR/data
  CAC_RUNANALYSIS=$CAC_NUMBERDIR/run-analysis
}

