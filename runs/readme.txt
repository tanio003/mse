# readme.txt
How to run mse on the local machine:
1. Make a new run directory (e.g., BGS_XXXX) in "Path_to_mse/mse/run/"
2. Copy wrapper functions in to the run directory: 
    "cp -r Pathtomse/mse/toolbox/wrappers_metagem Pathtomse/mse/run/BGS_XXXX/wrappers"
3. Edit the MATLAB file "wrappers/BatchSetup_metagem.m"
4. Edit and run the shell script "wrappers/job_local.sh"
5. All the solutions will be stored in the directory "Pathtomse/mse/run/BGS_XXXX/data/output/Solution"