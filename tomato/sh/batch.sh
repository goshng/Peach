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
function batch {
  date
#  RAW=smallraw
#  NUMBERCHROM=2
  RAW=raw
  NUMBERCHROM=12
  rm -rf output/count; mkdir output/count

  # Number of unique SNPs in IL lines per chromosome. 
  rm output/count/snp-uniq.txt 
  for g in $(eval echo {1..$NUMBERCHROM}); do
    for h in M I D; do
      v=$(cut -f 1,2,3,4 data/$RAW/IL_SNP_gene/IL$g-* |sort|uniq|cut -f 1|grep "^$h"|wc -l)
      echo $g $h $v >> output/count/snp-uniq.txt 
    done
  done 

  # Count SNPs (output/IL_SNP_gene/pen_SNP) in types in the genomes
  # We use only pen_SNP not others.
  rm -rf output/IL_SNP_gene; mkdir output/IL_SNP_gene
  for i in `ls data/$RAW/IL_SNP_gene`; do
    cut -f 1,3,4 data/$RAW/IL_SNP_gene/$i > output/IL_SNP_gene/$i
  done
  # Slow
  Rscript R/count-snp-type.R $RAW $NUMBERCHROM > output/count/pen-snp-type.txt
  

  # Circos data for the outer circle.
  rm -rf output/pen; mkdir output/pen
  perl pl/prepare-circos-data.pl pen \
    -infile data/$RAW/IL_SNP_gene/pen_SNP \
    -nchr $NUMBERCHROM \
    -outdir output/pen > output/count/pen.txt

  rm -rf output/snp; mkdir output/snp
  rm -rf output/othersnp; mkdir output/othersnp
  rm output/count/snp.txt
  rm output/count/othersnp.txt
  for i in `ls data/$RAW/IL_SNP_gene`; do
    ILNUM=$(echo $i | sed 's/\_SNP//')
    CHRNUM=$(echo $i | sed 's/IL\([0-9]*\)-.*/\1/')
    if [[ "$i" =~ IL ]]; then
      perl pl/prepare-circos-data.pl snp -in data/$RAW/IL_SNP_gene/$i \
           -chr $CHRNUM \
           -ilpositiondir data/$RAW/ilposition \
           -il $ILNUM \
           -out output/snp/$i 2>> output/count/snp.txt

      perl pl/prepare-circos-data.pl othersnp -in data/$RAW/IL_SNP_gene/$i \
           -chr $CHRNUM \
           -ilpositiondir data/$RAW/ilposition \
           -il $ILNUM \
           -out output/othersnp/$i 2>> output/count/othersnp.txt
    fi
  done


  # We have manually edited the files.
  #mkdir output/ilposition
  #for i in `ls output/snp`; do
    #echo $i
    #Rscript R/find-position-il.R output/snp/$i
  #done

  # for expression
  # Slow
  rm -rf output/position; mkdir output/position
  Rscript R/find-position-mrna.R $RAW > output/count/mrna-eqtl.txt
  # This depends on R/find-position-mrna.R
  rm -rf output/position2; mkdir output/position2
  Rscript R/find-position-mrna2.R > output/count/mrna-eqtl2.txt
  # This depends on R/find-position-mrna.R
  rm -rf output/position3; mkdir output/position3
  Rscript R/find-position-mrna3.R > output/count/mrna-eqtl3.txt
  # This create product description.
  rm -rf output/note; mkdir output/note
  Rscript R/find-note-mrna.R $RAW
  # This creates Gene Ontology analysis.
  rm -rf output/go; mkdir output/go
  Rscript R/goseq.R $RAW

  rm -rf output/links; mkdir output/links
  rm -rf output/exp; mkdir output/exp
  rm -rf output/inter; mkdir output/inter
  for i in `ls output/position2`; do
    ILNUM=$(echo $i | sed 's/\.csv//')
    CHRNUM=$(echo $i | sed 's/IL\([0-9]*\)-.*/\1/')
    if [[ "$i" =~ IL ]]; then
      perl pl/prepare-circos-data.pl links -in output/position2/$i \
                -chr $CHRNUM \
                -il $ILNUM \
                -ilpositiondir data/$RAW/ilposition \
                -out output/links/$i
      perl pl/prepare-circos-data.pl exp -in output/position2/$i \
                -chr $CHRNUM \
                -il $ILNUM \
                -ilpositiondir data/$RAW/ilposition \
                -out output/exp/$i
      perl pl/prepare-circos-data.pl inter -in output/position2/$i \
                -chr $CHRNUM \
                -il $ILNUM \
                -ilpositiondir data/$RAW/ilposition \
                -out output/inter/$i
    fi
  done

  rm -rf output/exppen; mkdir output/exppen
  perl pl/prepare-circos-data.pl exppen -infile output/position2/PENeQTL.csv \
    -nchr $NUMBERCHROM \
    -outdir output/exppen

  rm -rf output/summary; mkdir output/summary
  Rscript R/summarize-count.R $RAW

  rm -rf output/ilname; mkdir output/ilname
  for i in `ls data/$RAW/ilposition`; do 
    perl pl/s.pl data/$RAW/ilposition/$i > output/ilname/$i; 
  done

  rm -rf output/circos; mkdir output/circos
  rm output/chr-il.txt
  perl pl/prepare-circos-data.pl circos \
    -infile data/$RAW/chr-il.txt \
    -nchr $NUMBERCHROM \
    -outdir output/circos 2>> output/chr-il.txt

  rm -rf output/circosout; mkdir output/circosout
  for g in $(eval echo {1..$NUMBERCHROM}); do
    mkdir output/circosout/$g
    circos -conf output/circos/chr$g.conf -outputdir output/circosout/$g
  done 

  date
}
