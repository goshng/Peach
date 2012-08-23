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
function find-position {
  mkdir output/count

  # Count SNPs (output/IL_SNP_gene/pen_SNP) in types in the genomes
  # We use only pen_SNP not others.
  mkdir output/IL_SNP_gene
  for i in `ls input/raw/IL_SNP_gene`; do
    cut -f 1,3,4 input/raw/IL_SNP_gene/$i > output/IL_SNP_gene/$i
  done 
  Rscript R/count-snp-type.R > output/count/pen-snp-type.txt

  rm output/count/snp-uniq.txt 
  CHRNUM=12
  for g in $(eval echo {1..$CHRNUM}); do
    for h in M I D; do
      v=$(cut -f 1,2,3,4 input/raw/IL_SNP_gene/IL$g-* |sort|uniq|cut -f 1|grep "^$h"|wc -l)
      echo $g $h $v >> output/count/snp-uniq.txt 
    done
  done 

  mkdir output/position
  Rscript R/find-position-mrna.R > output/count/mrna-eqtl.txt

  mkdir output/position2
  Rscript R/find-position-mrna2.R > output/count/mrna-eqtl2.txt

  mkdir output/position3
  Rscript R/find-position-mrna3.R > output/count/mrna-eqtl3.txt

  mkdir output/pen
  perl pl/prepare-circos-data.pl pen \
    -infile input/raw/IL_SNP_gene/pen_SNP \
    -nchr 12 \
    -outdir output/pen > output/count/pen.txt

  mkdir output/snp
  mkdir output/othersnp
  rm output/count/snp.txt
  rm output/count/othersnp.txt
  for i in `ls input/raw/IL_SNP_gene`; do
    CHRNUM=$(echo $i | sed 's/IL\([0-9]*\)-.*/\1/')
    perl pl/prepare-circos-data.pl snp -in input/raw/IL_SNP_gene/$i \
         -chr $CHRNUM \
         -out output/snp/$i 2>> output/count/snp.txt

    perl pl/prepare-circos-data.pl othersnp -in input/raw/IL_SNP_gene/$i \
         -chr $CHRNUM \
         -out output/othersnp/$i 2>> output/count/othersnp.txt
  done

  mkdir output/summary
  Rscript R/summarize-count.R

  # We have manually edited the files.
  #mkdir output/ilposition
  #for i in `ls output/snp`; do
    #echo $i
    #Rscript R/find-position-il.R output/snp/$i
  #done

  mkdir output/links
  for i in `ls output/position2`; do
    ILNUM=$(echo $i | sed 's/\.csv//')
    CHRNUM=$(echo $i | sed 's/IL\([0-9]*\)-.*/\1/')
    perl pl/prepare-circos-data.pl links -in output/position2/$i \
              -chr $CHRNUM \
              -il $ILNUM \
              -ilpositiondir output/ilposition \
              -out output/links/$i
  done

  mkdir output/exp
  for i in `ls output/position2`; do
    ILNUM=$(echo $i | sed 's/\.csv//')
    CHRNUM=$(echo $i | sed 's/IL\([0-9]*\)-.*/\1/')
    perl pl/prepare-circos-data.pl exp -in output/position2/$i \
              -chr $CHRNUM \
              -il $ILNUM \
              -ilpositiondir output/ilposition \
              -out output/exp/$i
  done 

  mkdir output/inter
  for i in `ls output/position2`; do
    ILNUM=$(echo $i | sed 's/\.csv//')
    CHRNUM=$(echo $i | sed 's/IL\([0-9]*\)-.*/\1/')
    perl pl/prepare-circos-data.pl inter -in output/position2/$i \
              -chr $CHRNUM \
              -il $ILNUM \
              -ilpositiondir output/ilposition \
              -out output/inter/$i
  done

  mkdir output/exppen
  perl pl/prepare-circos-data.pl exppen -infile output/position2/PENeQTL.csv \
    -nchr 12 \
    -outdir output/exppen

  mkdir output/circos
  rm output/chr-il.txt
  perl pl/prepare-circos-data.pl circos \
    -infile input/chr-il.txt \
    -nchr 12 \
    -outdir output/circos 2>> output/chr-il.txt

  CHRNUM=12
  mkdir output/circosout
  for g in $(eval echo {1..$CHRNUM}); do
    mkdir output/circosout/$g
    circos -conf output/circos/chr$g.conf -outputdir output/circosout/$g
  done 


}
