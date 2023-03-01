# readme.txt
How to run mse on the local machine:
1. Make a new run directory (e.g., BGS_XXXX) in "Path_to_mse/mse/run/"
2. Make wrapper directory in the run directory:
    "mkdir Pathtomse/mse/run/BGS_XXXX/wrappers"
3. Copy wrapper functions in to the wrapper directory: 
    "cp -r Pathtomse/mse/toolbox/wrappers_metagem Pathtomse/mse/run/BGS_XXXX/wrappers"
4. Edit MATLAB file "wrappers/BatchSetup_metagem.m"
5. Edit and run the shell script "wrappers/job_local.sh"
6. All the solutions will be stored in the directory "Pathtomse/mse/run/BGS_XXXX/data/output/Solution"
