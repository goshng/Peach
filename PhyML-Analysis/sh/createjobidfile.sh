for i in `ls input`; do 
  N=$(head -n 1 input/$i | cut -f 2 -d" ") 
  if [ ! "$N" == "2" ]; then
    echo ./phyml -b 500 -m GTR -f e -c 4 -a e -s NNI --quiet -i input/$i
  fi
done 
