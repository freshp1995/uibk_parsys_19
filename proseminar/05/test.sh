#space
P=10000
#objects
N=1000

for (( i=10; i<$P; i+=10 ))
do
	for (( j=2; j<$N && j < $i; j+=10 ))
	do
        result="$(./n-body $i $j 10000)"
		echo $i";"$j";"$result";">> results.txt
	done
    echo  $i
done
echo 