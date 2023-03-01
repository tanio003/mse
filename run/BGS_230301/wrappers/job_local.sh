#!/bin/zsh
START=1
END=303
STEP=1
RUNDIR="BGS_230301"

cd /Users/tatsurotanioka/Desktop/Project/mse/run/$RUNDIR/wrappers
matlab -nodesktop -nodisplay -r "BatchSetup_metagem;quit;"

for i in $(seq $START $STEP $END) ; do	
	let input=$i
	echo $input
	matlab -nodesktop -nodisplay -r "mse2($input,'$RUNDIR');quit;" 
done

