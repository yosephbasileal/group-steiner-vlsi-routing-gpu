#!/bin/bash

declare -a arr=("13.txt" "24.txt" "34.txt" "55.txt" "96.txt" "99.txt" "LIN/12.stp" "LIN/22.stp" "LIN/27.stp") 
TIMEFORMAT=%R
line="____________________________________________"
break="---------------------------------------------------"
echo $line `date` $line >> log/time.log
echo $line `date` $line >> log/output.log
for i in  "${arr[@]}"
do
	echo "$i"
	echo $break >> log/time.log
	echo $break >> log/output.log
	echo -n "$i" >> log/time.log
	echo -n "$i" >> log/output.log
	echo -n "    (V,E)--> " >> log/time.log
	echo -n "    (V,E(--> " >> log/output.log
	size=$(head -n 1 grouped/$i)
	V=`echo $size java | cut -d' ' -f1`
	echo $V
	(echo $size) >> log/time.log
	(echo $size) >> log/output.log
	echo $break >> log/time.log
	echo $break >> log/output.log

	for p in {4,8,16,32,64,128,256,384,512,1024}
	do
		echo -n $p": "  >> log/time.log
		if ((V > p))
		then
			(time mpirun -n $p --hostfile hosts/allhosts ./twostar < grouped/"$i"  >> log/output.log) &>> log/time.log
		else
			echo "# of nodes is more than needed" >> log/time.log
		fi
	done
	printf "\n\n"  >> log/time.log
	printf "\n\n"  >> log/output.log
done
