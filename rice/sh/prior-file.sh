# Author: Sang Chul Choi
# Date  : Sat Jun 11 16:10:01 EDT 2011

function prior-file {
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
        -out $DATADIR/$FUNCNAME.txt
      echo "Check $DATADIR/$FUNCNAME.txt"
      break
    fi
  done
}
