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

function noweb-tangle {
OUTSRC=output/src
NOWEBFILE=$ROOTDIR/noweb/main.nw
#CPPFILE=( tux.cpp \
#          tuxHelp.cpp \
#          tuxManager.cpp \
#          tuxModel.cpp \
#          tuxMcmc.cpp \
#          tuxXyz.cpp \
#          gsl.cpp \
#          fileman.cpp \
#          yamlEmitter.cpp \
#          monsters.cpp )
#HPPFILE=( tuxHelp.h \
#          tuxManager.h \
#          tuxModel.h \
#          tuxMcmc.h \
#          tuxXyz.h \
#          config.h.in )

CPPFILE=( tux.cpp \
          tuxManager.cpp \
          tuxRandomVariable.cpp \
          tuxParameter.cpp \
          tuxData.cpp \
          tuxSystem.cpp \
          tuxProbability.cpp \
          tuxProbTxBoundaries.cpp \
          tuxProbTxExpression.cpp \
          tuxLikelihood.cpp \
          tuxMoverManager.cpp \
          tuxMover.cpp \
          tuxMoverTxBoundaries.cpp \
          tuxMoverTxExpression.cpp \
          tuxSummarizer.cpp \
          tuxSummary.cpp \
          tuxTxBoundaries.cpp \
          tuxTxExpression.cpp \
          tuxReads.cpp \
          ezlogger.cpp \
          tuxMcmc.cpp \
          tuxChainManager.cpp \
          tuxChain.cpp )
HPPFILE=( tuxManager.h \
          tuxRandomVariable.h \
          tuxParameter.h \
          tuxData.h \
          tuxSystem.h \
          tuxProbability.h \
          tuxProbTxBoundaries.h \
          tuxProbTxExpression.h \
          tuxLikelihood.h \
          tuxMoverManager.h \
          tuxMover.h \
          tuxMoverTxBoundaries.h \
          tuxMoverTxExpression.h \
          tuxSummarizer.h \
          tuxSummary.h \
          tuxTxBoundaries.h \
          tuxTxExpression.h \
          tuxReads.h \
          tuxMcmc.h \
          tuxChainManager.h \
          gtypes.h \
          tuxChain.h )
notangle -RCMakeLists.txt $NOWEBFILE > $OUTSRC/CMakeLists.txt

for f in ${CPPFILE[@]}; do
  notangle -L -R$f $NOWEBFILE > $OUTSRC/$f
done

for f in ${HPPFILE[@]}; do
  notangle -L -R$f $NOWEBFILE > $OUTSRC/$f
  # notangle -R$f $NOWEBFILE | cpif $OUTSRC/$f
done

sed s/IMRBUILDDATETIME/$(date +%y%m%d%H%M%S)/g \
  < $OUTSRC/tuxHelp.cpp > $OUTSRC/tuxHelp.cpp.tmp
mv $OUTSRC/tuxHelp.cpp.tmp $OUTSRC/tuxHelp.cpp
}
