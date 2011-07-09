function convert-fsc-to-imr {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SIMULATIONS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      REPETITION=1
      global-variable $SPECIES $REPETITION

      # Create output directories.
      #mkdir -p $NUMBERDIR
      #mkdir -p $DATADIR
      #mkdir -p $IMRDIR
      #mkdir -p $FSCDIR

      # Read options from sim1 file to run fsc.
      NUMREP=$(grep NUMREP $SPECIESFILE | cut -d":" -f2)
      OUTFILE=$(grep OUTFILE $SPECIESFILE | cut -d":" -f2)
      BLOCKLENGTH=$(grep BLOCKLENGTH $SPECIESFILE | cut -d":" -f2)

      # Simulate data using the par file.
      NUMREP=2
      echo perl pl/convert-fsc.pl -in $FSCDIR/${OUTFILE}_1_${NUMREP}.arp \
        -length $BLOCKLENGTH \
        -repetition $NUMREP \
        -inTree $FSCDIR/${OUTFILE}_1_true_trees.trees \
        -out $IMRDIR/${OUTFILE}_1_${NUMREP}.imr
      echo "Check $IMRDIR/${OUTFILE}_1_${NUMREP}.imr"
      break
    fi
  done

}
