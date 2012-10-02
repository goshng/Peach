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
function figures {
  START_TIME=`date +%s`

  # The circos figures generation takes the most of the time.
  rm -rf output/circos; mkdir output/circos
  rm output/chr-il.txt
  perl pl/prepare-circos-data.pl circos \
    -infile data/$RAW/chr-il.txt \
    -nchr $NUMBERCHROM \
    -outdir output/circos 2>> output/chr-il.txt
  perl pl/prepare-circos-data.pl circos2 \
    -ilpositionfile data/$RAW/ilmap2.txt \
    -nchr $NUMBERCHROM \
    -chril data/$RAW/order.csv \
    -outdir output/circos

  rm -rf output/circosout; mkdir output/circosout
  for g in $(eval echo {1..$NUMBERCHROM}); do
    mkdir output/circosout/$g
    circos -conf output/circos/chr$g.conf -outputdir output/circosout/$g
  done 

  for i in `ls output/position2`; do
    ILNUM=$(echo $i | sed 's/\.csv//')
    CHRNUM=$(echo $i | sed 's/IL\([0-9]*\)-.*/\1/')
    if [[ "$i" =~ IL ]]; then
      mkdir output/circosout/$ILNUM
      circos -conf output/circos/$ILNUM.conf -outputdir output/circosout/$ILNUM
    fi
  done 

  END_TIME=`date +%s`
  ELAPSED=`expr $END_TIME - $START_TIME + 18000`
  echo "FINISHED at "`date`" Elapsed time: "`date -r $ELAPSED +%H:%M:%S`
}

