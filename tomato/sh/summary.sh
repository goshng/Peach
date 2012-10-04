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
function summary {
  START_TIME=`date +%s`

  # Create an R table of the ilmap file.
  perl pl/prepare-circos-data.pl parse-ilposition \
    -rfriendly \
    -ilpositionfile data/$RAW/ilmap2.txt \
    -chril data/$RAW/order.csv \
    > output/ilmap2.dat

  ##############################################################################
  # BEGIN: for SNPs
  rm -rf output/count; mkdir output/count

  # Numbers of unique SNPs in chromosome. 
  rm -f output/count/snp-uniq.txt 
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

  ##############################################################################
  # BEGIN: for expression
  # Slow
  rm -rf output/position; mkdir output/position
  Rscript R/find-position-mrna.R $RAW > output/count/mrna-eqtl.txt
  # This depends on R/find-position-mrna.R
  rm -rf output/position2; mkdir output/position2
  Rscript R/find-position-mrna2.R > output/count/mrna-eqtl2.txt
  # This depends on R/find-position-mrna.R
  rm -rf output/position3; mkdir output/position3
  Rscript R/find-position-mrna3.R output/ilmap2.dat
  # This create product description.
  rm -rf output/note; mkdir output/note
  Rscript R/find-note-mrna.R $RAW

  rm -rf output/position4; mkdir output/position4
  rm -rf output/cis2; mkdir output/cis2
  rm -rf output/trans2; mkdir output/trans2
  rm -rf output/links2; mkdir output/links2
  rm -rf output/ilname; mkdir output/ilname
  rm -rf output/cissnpinsideil ; mkdir output/cissnpinsideil
  rm -rf output/cissnpoutsideil ; mkdir output/cissnpoutsideil
  rm -rf output/transsnpil ; mkdir output/transsnpil
  for i in `ls output/position3`; do
    ILNUM=$(echo $i | sed 's/\.csv//')
    CHRNUM=$(echo $i | sed 's/IL\([0-9]*\)-.*/\1/')
    if [[ "$i" =~ IL ]]; then
      Rscript R/eqtl-cis-trans.R $NUMBERCHROM \
        data/karyotype/karyotype.sol.txt \
        output/position3/$ILNUM.csv \
        output/IL_SNP_gene/${ILNUM}_SNP \
        output/ilmap2.dat \
        output/cis2 output/trans2 output/links2 output/ilname \
        output/cissnpinsideil output/cissnpoutsideil output/transsnpil
    fi
  done

  rm -rf output/exppen; mkdir output/exppen
  perl pl/prepare-circos-data.pl exppen -infile output/position2/PENeQTL.csv \
    -nchr $NUMBERCHROM \
    -outdir output/exppen

  # Summary
  rm -rf output/summary; mkdir output/summary
  Rscript R/summarize-count.R $RAW > output/summary/log
  Rscript R/find-unique-eqtl.R $RAW
  # See output/summary/ven-snp.txt
  cut -f 1,2,3,4 data/$RAW/IL_SNP_gene/IL* | sort | uniq | grep -v ^type > output/count/ven-il-snp.txt
  cut -f 1,2,3,4 data/$RAW/IL_SNP_gene/pen_SNP | sort | uniq | grep -v ^type > output/count/ven-pen-snp.txt
  Rscript R/ven-snp.R
  # See output/summary/ven-eqtl.csv
  Rscript R/ven-eqtl.R $RAW 

  # IL Map figure
  perl pl/prepare-circos-data.pl ilmap \
    -ilpositionfile data/$RAW/ilmap2.txt \
    -chril data/$RAW/order.csv \
    > output/ilmap.tex
  pdflatex -output-directory output tex/ilmap.tex

  # Correlation of genes and metabolites.
  # eQTL gene-to-gene - done by correlation-il.R
  # total gene-to-gene - done by ...
  # gene-to-metabolite - done by correlation-metabolite.R
  # metabolite-to-metabolite - done by correlation-metabolite.R
  rm -rf output/cor; mkdir output/cor
  rm -rf output/metabolite; mkdir output/metabolite
  rm -rf output/metabolite-cor; mkdir output/metabolite-cor
  rm -rf output/cor-ab; mkdir output/cor-ab
  Rscript R/correlation-il.R $RAW
  for i in `ls data/$RAW/metabolite`; do
    Rscript R/correlation-metabolite.R $RAW $i
  done
  Rscript R/correlation-ab.R $RAW IL-expressionforcorrelation.csv for-apriori-subnetwork-a.csv for-apriori-subnetwork-b.csv ab

  for i in `ls output/cor/*-pval.csv`; do Rscript R/postprocess-il.R $i; done 
  
  END_TIME=`date +%s`
  ELAPSED=`expr $END_TIME - $START_TIME + 18000`
  echo "FINISHED at "`date`" Elapsed time: "`date -r $ELAPSED +%H:%M:%S`
}

function go {
  START_TIME=`date +%s`

  # This creates Gene Ontology analysis.
  rm -rf output/go; mkdir output/go
  Rscript R/goseq.R $RAW
  for i in output/go/*; do echo $i; cat $i; done > output/go-all.csv

  END_TIME=`date +%s`
  ELAPSED=`expr $END_TIME - $START_TIME + 18000`
  echo "FINISHED at "`date`" Elapsed time: "`date -r $ELAPSED +%H:%M:%S`
}

