#!/bin/zsh
START=1
END=303
STEP=1
RUNDIR="BGS_200826c"

cd /Users/tatsurotanioka/Desktop/Project/mse/run/$RUNDIR/wrappers
for i in $(seq $START $STEP $END) ; do	
	let input=$i
	echo $input
	matlab -nodesktop -nodisplay -r "mse2($input,'$RUNDIR');quit;" 
done

