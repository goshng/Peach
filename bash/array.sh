declare -a array1=( \
  zero1 \
  one1 \
  two1 )

for (( i = 0 ; i < 3 ; i++ ))
do
    echo "Element [$i]: ${array1[$i]}"
done
