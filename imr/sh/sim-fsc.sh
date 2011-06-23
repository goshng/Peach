function sim-fsc {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SIMULATIONS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      $SPECIES-fsc
      break
    fi
  done
}
