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
CPPFILE=( tux.cpp \
          tux-v1.cpp \
          fileman.cpp \
          monsters.cpp
          yamlEmitter.cpp \
          tux-gsl.cpp \
          tux-fixed.cpp \
          tux-tx.cpp \
          tux-probtx.cpp \
          tux-prob.cpp \
          tux-move.cpp \
          tux-chain.cpp \
          tux-cm.cpp \
          tux-mcmc.cpp \
          tux-single.cpp \
          tux-likelihood.cpp \
          tux-sum.cpp \
          tux-log.cpp \
          stl-transform.cpp \
          tuxSignal.cpp \
          tuxError.cpp \
          tuxDefault.cpp \
          tuxGsl.cpp \
          tuxGslRng.cpp \
          tuxHelp.cpp \
          tuxManager.cpp \
          tuxSystem.cpp \
          tuxProbTxBoundaries.cpp \
          tuxProbTxExpression.cpp \
          tuxLikelihood.cpp \
          tuxTotalLikelihood.cpp \
          tuxMoverManager.cpp \
          tuxMover.cpp \
          tuxMoverTxBoundaries.cpp \
          tuxMoverTxExpression.cpp \
          tuxMoverTxBSingle.cpp \
          tuxSummarizer.cpp \
          tuxSummary.cpp \
          tuxTxBoundaries.cpp \
          tuxTxExpression.cpp \
          tuxTxParameter.cpp \
          tuxReads.cpp \
          tuxFixed.cpp \
          ezlogger.cpp \
          tuxMcmc.cpp \
          tuxMcmcSingle.cpp \
          tuxMcmcMpi.cpp \
          tuxChainManager.cpp \
          tuxChain.cpp \
          MathFunctions/MakeTable.cpp \
          MathFunctions/mysqrt.cpp )
HPPFILE=( tuxManager.h \
          tuxHelp.h \
          tuxSignal.h \
          tuxError.h \
          tuxDefault.h \
          tuxGsl.h \
          tuxGslRng.h \
          benLogDouble.h \
          tuxRandomVariable.h \
          tuxParameter.h \
          tuxData.h \
          tuxSystem.h \
          tuxProbability.h \
          tuxProbTxBoundaries.h \
          tuxProbTxExpression.h \
          tuxLikelihood.h \
          tuxTotalLikelihood.h \
          tuxMoverManager.h \
          tuxMover.h \
          tuxMoverTxBoundaries.h \
          tuxMoverTxExpression.h \
          tuxMoverTxBSingle.h \
          tuxSummarizer.h \
          tuxSummary.h \
          tuxTxBoundaries.h \
          tuxTxExpression.h \
          tuxTxParameter.h \
          tuxReads.h \
          tuxFixed.h \
          tuxMcmc.h \
          tuxMcmcSingle.h \
          tuxMcmcMpi.h \
          tuxChainManager.h \
          tuxChain.h \
          gtypes.h \
          config.h.in )
notangle -Rruntux $NOWEBFILE > $OUTSRC/b/runtux; chmod +x $OUTSRC/b/runtux
notangle -Rtuxbatch.sh $NOWEBFILE > $OUTSRC/b/tuxbatch.sh
notangle -RplotTemperature.sh $NOWEBFILE > $OUTSRC/b/plotTemperature.sh
notangle -RplotTemperature.R $NOWEBFILE > $OUTSRC/b/plotTemperature.R

notangle -RCMakeLists.txt $NOWEBFILE > $OUTSRC/CMakeLists.txt
notangle -RMathFunctions/CMakeLists.txt $NOWEBFILE > $OUTSRC/MathFunctions/CMakeLists.txt

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
