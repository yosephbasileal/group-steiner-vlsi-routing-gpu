#!/bin/bash

declare -a arr=("13.txt" "24.txt" "34.txt" "55.txt" "96.txt" "99.txt" "LIN/12.stp" "LIN/22.stp" "LIN/27.stp") 
TIMEFORMAT=%R
line="____________________________________________"
break="---------------------------------------------------"
echo $line `date` $line >> time.log
echo $line `date` $line >> output.log
for i in  "${arr[@]}"
do
	echo "$i"
	echo $break >> time.log
	echo $break >> output.log
	echo -n "$i" >> time.log
	echo -n "$i" >> output.log
	echo -n "    (V,E)--> " >> time.log
	echo -n "    (V,E(--> " >> output.log
	size=$(head -n 1 grouped/$i)
	V=`echo $size java | cut -d' ' -f1`
	echo $V
	(echo $size) >> time.log
	(echo $size) >> output.log
	echo $break >> time.log
	echo $break >> output.log

	for p in {4,8,16,32,64,128,256,384,512,1024}
	do
		echo -n $p": "  >> time.log
		if ((V > p))
		then
			(time mpirun -n $p --hostfile /home/basileal/hosts ./twostar < grouped/"$i" >> output.log) &>> time.log
		else
			echo "# of nodes is more than needed" >> time.log
		fi
	done
	printf "\n\n"  >> time.log
	printf "\n\n"  >> output.log
done
