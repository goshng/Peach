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
    mainExtractPartOfPhylip
    # mainJModelTest
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
  echo -n "  extracting parts of the PHYLIP alignment file ... "
  GENES=(ND1 ND2 COI COII ATP8 ATP6 COIII ND3 ND4L ND4 ND5 cytB ND6 12SrRNA 16SrRna Dloop)
  GENESTART=(2718 3882 5293 6976 7727 7888 8568 9421 9839 10129 11704 14111 13512 70 1088 15386)
  GENEEND=(3677 4916 6837 7659 7930 8568 9351 9768 10135 11506 13515 15253 14036 1016 2645 16348)
  PHYLIPFILE=$PROJECTDATADIR/Mega_mtGenome_23112011_Fin.phylip

  y=$((${#GENES[*]} - 1))
  for i in $(eval echo {0..$y}); do
    perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format nexus \
      -out $DATADIR/${GENES[$i]}.nex -start ${GENESTART[$i]} -end ${GENEEND[$i]}

    perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format fasta \
      -out $DATADIR/${GENES[$i]}.fas -start ${GENESTART[$i]} -end ${GENEEND[$i]}
  done 

  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format nexus \
    -subsetOfSequences Dt4Y,Dt1Y,Dt10Y,Dt5Y,Dt2Y,Dt3Y,Dt6Y,Dt7Y,Dt8Y,Dt9Y \
    -out $DATADIR/region-n1.nex

  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format nexus \
    -subsetOfSequences Dt1T,Dt14T,Dt2T,Dt3T,Dt4T,Dt5T,Dt6T,Dt7T,Dt8T,Dt9T,Dt10T,Dt11T,Dt12T,Dt13T \
    -out $DATADIR/region-n2.nex

  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format nexus \
    -subsetOfSequences Dt1LK,Dt3LK,Dt5LK,Dt2LK,Dt4LK,Dt6LK,Dt7LK,Dt8LK,Dt9LK,Dt10LK \
    -out $DATADIR/region-n3.nex
   
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format nexus \
    -subsetOfSequences Dt1WB,Dt2WB,Dt3WB,Dt4WB,Dt5WB,Dt6WB,Dt7WB,Dt9WB,Dt8WB,Dt10WB,Dt11WB,Dt13WB,Dt12WB,Dt14WB,Dt15WB,Dt16WB \
    -out $DATADIR/region-n4.nex

  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format nexus \
    -subsetOfSequences Dt4Y,Dt1Y,Dt10Y,Dt5Y,Dt2Y,Dt3Y,Dt6Y,Dt7Y,Dt8Y,Dt9Y \
    -subsetOfSequences Dt1T,Dt14T,Dt2T,Dt3T,Dt4T,Dt5T,Dt6T,Dt7T,Dt8T,Dt9T,Dt10T,Dt11T,Dt12T,Dt13T \
    -subsetOfSequences Dt1LK,Dt3LK,Dt5LK,Dt2LK,Dt4LK,Dt6LK,Dt7LK,Dt8LK,Dt9LK,Dt10LK \
    -subsetOfSequences Dt1WB,Dt2WB,Dt3WB,Dt4WB,Dt5WB,Dt6WB,Dt7WB,Dt9WB,Dt8WB,Dt10WB,Dt11WB,Dt13WB,Dt12WB,Dt14WB,Dt15WB,Dt16WB \
    -subsetOfSequences Dg1A,Dg2A,Dg4A,Dg6A,Dg3C,Dg5A \
    -out $DATADIR/regions-5.nex

  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format genepop \
    -subsetOfSequences Dt4Y,Dt1Y,Dt10Y,Dt5Y,Dt2Y,Dt3Y,Dt6Y,Dt7Y,Dt8Y,Dt9Y \
    -subsetOfSequences Dt1T,Dt14T,Dt2T,Dt3T,Dt4T,Dt5T,Dt6T,Dt7T,Dt8T,Dt9T,Dt10T,Dt11T,Dt12T,Dt13T \
    -subsetOfSequences Dt1LK,Dt3LK,Dt5LK,Dt2LK,Dt4LK,Dt6LK,Dt7LK,Dt8LK,Dt9LK,Dt10LK \
    -subsetOfSequences Dt1WB,Dt2WB,Dt3WB,Dt4WB,Dt5WB,Dt6WB,Dt7WB,Dt9WB,Dt8WB,Dt10WB,Dt11WB,Dt13WB,Dt12WB,Dt14WB,Dt15WB,Dt16WB \
    -subsetOfSequences Dg1A,Dg2A,Dg4A,Dg6A,Dg3C,Dg5A \
    -out $DATADIR/regions-5.genepop

  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format genepop \
    -subsetOfSequences Dt4Y,Dt1Y,Dt10Y,Dt5Y,Dt2Y,Dt3Y,Dt6Y,Dt7Y,Dt8Y,Dt9Y \
    -out $DATADIR/region-n1.genepop

  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format genepop \
    -subsetOfSequences Dt1T,Dt14T,Dt2T,Dt3T,Dt4T,Dt5T,Dt6T,Dt7T,Dt8T,Dt9T,Dt10T,Dt11T,Dt12T,Dt13T \
    -out $DATADIR/region-n2.genepop

  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format genepop \
    -subsetOfSequences Dt1LK,Dt3LK,Dt5LK,Dt2LK,Dt4LK,Dt6LK,Dt7LK,Dt8LK,Dt9LK,Dt10LK \
    -out $DATADIR/region-n3.genepop
   
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format genepop \
    -subsetOfSequences Dt1WB,Dt2WB,Dt3WB,Dt4WB,Dt5WB,Dt6WB,Dt7WB,Dt9WB,Dt8WB,Dt10WB,Dt11WB,Dt13WB,Dt12WB,Dt14WB,Dt15WB,Dt16WB \
    -out $DATADIR/region-n4.genepop


  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format migrate-n \
    -subsetOfSequences Dt4Y,Dt1Y,Dt10Y,Dt5Y,Dt2Y,Dt3Y,Dt6Y,Dt7Y,Dt8Y,Dt9Y \
    -out $DATADIR/region-n1.migrate-n

  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format migrate-n \
    -subsetOfSequences Dt1T,Dt14T,Dt2T,Dt3T,Dt4T,Dt5T,Dt6T,Dt7T,Dt8T,Dt9T,Dt10T,Dt11T,Dt12T,Dt13T \
    -out $DATADIR/region-n2.migrate-n

  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format migrate-n \
    -subsetOfSequences Dt1LK,Dt3LK,Dt5LK,Dt2LK,Dt4LK,Dt6LK,Dt7LK,Dt8LK,Dt9LK,Dt10LK \
    -out $DATADIR/region-n3.migrate-n
   
  perl pl/phyilp-part.pl phylip -in $PHYLIPFILE -format migrate-n \
    -subsetOfSequences Dt1WB,Dt2WB,Dt3WB,Dt4WB,Dt5WB,Dt6WB,Dt7WB,Dt9WB,Dt8WB,Dt10WB,Dt11WB,Dt13WB,Dt12WB,Dt14WB,Dt15WB,Dt16WB \
    -out $DATADIR/region-n4.migrate-n

  echo "done!" 
}

function mainJModelTest {
  GENES=(ND1 ND2 COI COII ATP8 ATP6 COIII ND3 ND4L ND4 ND5 cytB ND6 12SrRNA 16SrRna Dloop)
  for i in "${GENES[@]}"; do
    java -jar $JMODELTEST -d $DATADIR/$i.fas -o $DATADIR/$i.jmt -s 11 -p -AICc
  done 
  tail -n 1 $DATADIR/*.jmt
}


