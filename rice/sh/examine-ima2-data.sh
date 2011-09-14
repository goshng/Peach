# Author: Sang Chul Choi
# Date  : Sat Jun 11 13:08:14 EDT 2011

function examine-ima2-data {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -n "What repetition do you wish to run? (e.g., 1) "
      read REPETITION
      global-variable $SPECIES $REPETITION

      perl pl/$FUNCNAME.pl \
        -ima2data $DATADIR/summarize-sequence-data.txt

      break
    fi
  done
}
