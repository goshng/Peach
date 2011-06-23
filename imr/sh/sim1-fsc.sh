function sim1-fsc {
      REPETITION=1
      global-variable $SPECIES $REPETITION

      # Create output directories.
      mkdir -p $NUMBERDIR
      mkdir -p $DATADIR
      mkdir -p $IMRDIR
      mkdir -p $FSCDIR

      # Read options from sim1 file to run fsc.
      OUTFILE=$(grep OUTFILE $SPECIESFILE | cut -d":" -f2)
      SEED=$(grep SEED $SPECIESFILE | cut -d":" -f2)
      NUMREP=$(grep NUMREP $SPECIESFILE | cut -d":" -f2)
      NUMPOP=$(grep NUMPOP $SPECIESFILE | cut -d":" -f2)
      POPSIZE=$(grep POPSIZE $SPECIESFILE | cut -d":" -f2)
      NUMSAMPLE=$(grep NUMSAMPLE $SPECIESFILE | cut -d":" -f2)
      GROWTHRATE=$(grep GROWTHRATE $SPECIESFILE | cut -d":" -f2)
      NUMMIG=$(grep NUMMIG $SPECIESFILE | cut -d":" -f2)
      NUMEVENT=$(grep NUMEVENT $SPECIESFILE | cut -d":" -f2)
      NUMLOCI=$(grep NUMLOCI $SPECIESFILE | cut -d":" -f2)
      NUMLINKBLOCK=$(grep NUMLINKBLOCK $SPECIESFILE | cut -d":" -f2)
      PERBLOCK=$(grep PERBLOCK $SPECIESFILE | cut -d":" -f2)

      # Create a par file for FSC.
      FSCPAR=$FSCDIR/$OUTFILE.par
      echo -e "//\n$NUMPOP\n//\n$POPSIZE\n//\n$NUMSAMPLE\n//\n$GROWTHRATE\n//\n$NUMMIG\n//" > $FSCPAR
      echo -e "$NUMEVENT\n//\n$NUMLOCI\n//\n$NUMLINKBLOCK\n//\n$PERBLOCK" >> $FSCPAR

      $FSC -i $FSCPAR -n $NUMREP
}
