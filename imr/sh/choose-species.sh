function choose-species {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -n "What repetition do you wish to run? (e.g., 1) "
      read REPETITION
      global-variable $SPECIES $REPETITION

      mkdir -p $NUMBERDIR
      mkdir -p $DATADIR
      mkdir -p $IMRDIR
      mkdir -p $FSCDIR

      ssh -x $CAC_USERHOST mkdir -p $CAC_NUMBERDIR \
                                    $CAC_DATADIR \
                                    $CAC_IMRDIR \
                                    $CAC_FSCDIR
      break
    fi
  done
}
