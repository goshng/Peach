# Author: Sang Chul Choi
# Date  : Tue Jun 14 10:28:12 EDT 2011

function sim-imaa {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SIMULATIONS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      $SPECIES-imaa
      break
    fi
  done
}
