#!/bin/bash

declare -a arr=("13.txt" "24.txt" "34.txt" "55.txt" "96.txt" "99.txt" "LIN/12.stp" "LIN/22.stp" "LIN/27.stp") 
TIMEFORMAT=%R
line="____________________________________________"
break="---------------------------------------------------"
echo $line `date` $line >> log/time_oak.log
echo $line `date` $line >> log/output_oak.log
for i in  "${arr[@]}"
do
	echo "$i"
	echo $break >> log/time_oak.log
	echo $break >> log/output_oak.log
	echo -n "$i" >> log/time_oak.log
	echo -n "$i" >> log/output_oak.log
	echo -n "    (V,E)--> " >> log/time_oak.log
	echo -n "    (V,E(--> " >> log/output_oak.log
	size=$(head -n 1 grouped/$i)
	V=`echo $size java | cut -d' ' -f1`
	echo $V
	(echo $size) >> log/time_oak.log
	(echo $size) >> log/output_oak.log
	echo $break >> log/time_oak.log
	echo $break >> log/output_oak.log

	for p in {4,8,16,32,64,128}
	do
		echo -n $p": "  >> log/time_oak.log
		if ((V > p))
		then
			(time mpirun -n $p --hostfile hosts/oakonly ./twostar < grouped/"$i" >> log/output_oak.log) &>> log/time_oak.log
		else
			echo "# of nodes is more than needed" >> log/time_oak.log
		fi
	done
	printf "\n\n"  >> log/time_oak.log
	printf "\n\n"  >> log/output_oak.log
done
