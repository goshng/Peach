function ima2-prepare {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -n "What repetition do you wish to run? (e.g., 1) "
      read REPETITION
      global-variable $SPECIES $REPETITION

      WALLTIME=$(grep REPETITION${REPETITION}-WALLTIME species/$SPECIES | cut -d":" -f2)
      BURNIN=$(grep REPETITION${REPETITION}-BURNIN species/$SPECIES | cut -d":" -f2)
      SAMPLESIZE=$(grep REPETITION${REPETITION}-SAMPLESIZE species/$SPECIES | cut -d":" -f2)
      THIN=$(grep REPETITION${REPETITION}-THIN species/$SPECIES | cut -d":" -f2)
      HEATN=$(grep REPETITION${REPETITION}-HEATN species/$SPECIES | cut -d":" -f2)
      HEATA=$(grep REPETITION${REPETITION}-HEATA species/$SPECIES | cut -d":" -f2)
      HEATB=$(grep REPETITION${REPETITION}-HEATB species/$SPECIES | cut -d":" -f2)
      OTHER=$(grep REPETITION${REPETITION}-OTHER species/$SPECIES | cut -d":" -f2)
      PRIORFILE=$(grep REPETITION${REPETITION}-PRIORFILE species/$SPECIES | cut -d":" -f2)
      INFILE=$(grep REPETITION${REPETITION}-INFILE species/$SPECIES | cut -d":" -f2)
      OUTFILE=$(grep REPETITION${REPETITION}-OUTFILE species/$SPECIES | cut -d":" -f2)
      mkdir-species
      echo "  Creating directories..."
      
      scp -q $DATADIR/* $CAC_USERHOST:$CAC_DATADIR
      echo "  Copying input files..."
      copy-batch-sh-run-ima2 
      echo "Go to CAC's $SPECIES ima2, and execute nsub batch.sh"
      break
    fi
  done
}

function copy-batch-sh-run-ima2 {
cat>$IMA2DIR/batch.sh<<EOF
#!/bin/bash
#PBS -l walltime=${WALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N ${PROJECTNAME}-${SPECIES}
#PBS -q ${QUEUENAME}
#PBS -m e
#PBS -M ${BATCHEMAIL}

WORKDATADIR=\$PBS_O_WORKDIR/../data
IMA2=\$HOME/${BATCHIMA2}

IMA2DIR=\$TMPDIR/output/$SPECIES/$REPETITION/ima2
DATADIR=\$TMPDIR/output/$SPECIES/$REPETITION/data

mkdir -p \$IMA2DIR \$DATADIR
cp \$IMA2 \$TMPDIR/
cp \$WORKDATADIR/* \$DATADIR/

cd \$TMPDIR

# Commands

# First
./IMa2 -b$BURNIN -l$SAMPLESIZE -d$THIN -hn$HEATN -ha$HEATA -hb$HEATB $OTHER \\
  -g\$DATADIR/${PRIORFILE}1 \\
  -i\$DATADIR/${INFILE}1 \\
  -o\$IMA2DIR/$OUTFILE-1-1 &
sleep 5

# Second
./IMa2 -b$BURNIN -l$SAMPLESIZE -d$THIN -hn$HEATN -ha$HEATA -hb$HEATB $OTHER \\
  -g\$DATADIR/${PRIORFILE}1 \\
  -i\$DATADIR/${INFILE}1 \\
  -o\$IMA2DIR/$OUTFILE-1-2 &
sleep 5

# First
./IMa2 -b$BURNIN -l$SAMPLESIZE -d$THIN -hn$HEATN -ha$HEATA -hb$HEATB $OTHER \\
  -g\$DATADIR/${PRIORFILE}2 \\
  -i\$DATADIR/${INFILE}2 \\
  -o\$IMA2DIR/$OUTFILE-2-1 &
sleep 5

# Second
./IMa2 -b$BURNIN -l$SAMPLESIZE -d$THIN -hn$HEATN -ha$HEATA -hb$HEATB $OTHER \\
  -g\$DATADIR/${PRIORFILE}2 \\
  -i\$DATADIR/${INFILE}2 \\
  -o\$IMA2DIR/$OUTFILE-2-2 &
sleep 5

# First
./IMa2 -b$BURNIN -l$SAMPLESIZE -d$THIN -hn$HEATN -ha$HEATA -hb$HEATB $OTHER \\
  -g\$DATADIR/${PRIORFILE}3 \\
  -i\$DATADIR/${INFILE}3 \\
  -o\$IMA2DIR/$OUTFILE-3-1 &
sleep 5

# Second
./IMa2 -b$BURNIN -l$SAMPLESIZE -d$THIN -hn$HEATN -ha$HEATA -hb$HEATB $OTHER \\
  -g\$DATADIR/${PRIORFILE}3 \\
  -i\$DATADIR/${INFILE}3 \\
  -o\$IMA2DIR/$OUTFILE-3-2 &
sleep 5

wait

cp -r \$IMA2DIR/ \$PBS_O_WORKDIR/
cd
rm -rf \$TMPDIR
EOF
  chmod a+x $IMA2DIR/batch.sh
  scp -q $IMA2DIR/batch.sh $CAC_USERHOST:$CAC_IMA2DIR
}
