#!/bin/bash

declare -a arr=("13.txt" "24.txt" "34.txt" "55.txt" "96.txt" "99.txt" "LIN/12.stp" "LIN/22.stp" "LIN/27.stp") 
TIMEFORMAT=%R
line="____________________________________________"
break="---------------------------------------------------"
echo $line `date` $line >> time_oak.log
echo $line `date` $line >> output_oak.log
for i in  "${arr[@]}"
do
	echo "$i"
	echo $break >> time_oak.log
	echo $break >> output_oak.log
	echo -n "$i" >> time_oak.log
	echo -n "$i" >> output_oak.log
	echo -n "    (V,E)--> " >> time_oak.log
	echo -n "    (V,E(--> " >> output_oak.log
	size=$(head -n 1 grouped/$i)
	V=`echo $size java | cut -d' ' -f1`
	echo $V
	(echo $size) >> time_oak.log
	(echo $size) >> output_oak.log
	echo $break >> time_oak.log
	echo $break >> output_oak.log

	for p in {4,8,16,32,64,128,256,384,512}
	do
		echo -n $p": "  >> time_oak.log
		if ((V > p))
		then
			(time mpirun -n $p --hostfile oakonly ./twostar < grouped/"$i" >> output_oak.log) &>> time_oak.log
		else
			echo "# of nodes is more than needed" >> time_oak.log
		fi
	done
	printf "\n\n"  >> time_oak.log
	printf "\n\n"  >> output_oak.log
done
