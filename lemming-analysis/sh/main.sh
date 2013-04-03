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

function main {
  PS3="Choose the project for $FUNCNAME: "
  select PROJECT in ${PROJECTS[@]}; do 
  if [ "$PROJECT" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  else  
    project-repetition
    main-variable
    # mainExtractPartOfPhylip
    mainJModelTest
    break
  fi
  done

}

# Set the shell variables for main script
function main-variable {
  # Local output directories.
  BASEDIR=$PROJECTOUTPUTDIR/$PROJECT
  NUMBERDIR=$BASEDIR/$REPETITION
  DATADIR=$NUMBERDIR/data
  BWADIR=$NUMBERDIR/bwa
  RUNANALYSIS=$NUMBERDIR/run-analysis

  # Remote output directories.
  RBASEDIR=$ROUTPUTDIR/$SPECIES
  RNUMBERDIR=$ROUTPUTDIR/$SPECIES/$REPETITION
  RDATADIR=$RNUMBERDIR/data
  RBWADIR=$RNUMBERDIR/bwa
  RRUNANALYSIS=$RNUMBERDIR/run-analysis

  # Compute node output directories.
  CBASEDIR=output/$SPECIES
  CNUMBERDIR=$CBASEDIR/$REPETITION
  CDATADIR=$CNUMBERDIR/data
  CBWADIR=$CNUMBERDIR/bwa
  CRUNANALYSIS=$CNUMBERDIR/run-analysis
}

function mainExtractPartOfPhylip {
  PHYLIPFILE=$PROJECTDATADIR/Mega_mtGenome_23112011_Fin.phylip
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/ND1.fas -start 2718 -end 3677
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/ND2.fas -start 3882 -end 4916
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/COI.fas -start 5293 -end 6837
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/COII.fas -start 6976 -end 7659
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/ATP8.fas -start 7727 -end 7930
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/ATP6.fas -start 7888 -end 8568
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/COIII.fas -start 8568 -end 9351
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/ND3.fas -start 9421 -end 9768
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/ND4L.fas -start 9839 -end 10135
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/ND4.fas -start 10129 -end 11506
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/ND5.fas -start 11704 -end 13515
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/cytB.fas -start 14111 -end 15253
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/ND6.fas -start 13512 -end 14036
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/12SrRNA.fas -start 70 -end 1016
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/16SrRNA.fas -start 1088 -end 2645
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
    -out $DATADIR/Dloop.fas -start 15386 -end 16348
}

function mainJModelTest {
  GENES=(ND1 ND2 COI COII ATP8 ATP6 COIII ND3 ND4L ND4 ND5 cytB ND6 12SrRNA 16SrRna Dloop)
  for i in "${GENES[@]}"; do
    java -jar $JMODELTEST -d $DATADIR/$i.fas -o $DATADIR/$i.jmt -s 11 -p -AICc
  done 
  tail -n 1 $DATADIR/*.jmt
}


