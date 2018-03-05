
#parallel -j 12 "python ~/bin/sfs_from_slim_update_bootstrap_oneSlim.py -i $1/Rho0.009/runs/Rho0.009.{}.out.gz --gz -w 100 -b 70000  -t 0 -o $1/Rho0.009/sfs/Rho0.009.{}" ::: $(seq 1 1000)
parallel "python ~/bin/sfs_from_slim_update_bootstrap_oneSlim.py -i $1/Rho0.0045/runs/Rho0.0045.{}.out --gz -w 100 -b 70000 -t 0 -o $1/Rho0.0045/sfs/Rho0.0045.{}" ::: $(seq 1 1000)
#parallel -j 12 "python ~/bin/sfs_from_slim_update_bootstrap_oneSlim.py -i $1/Rho0.0009/runs/Rho0.0009.{}.out.gz --gz -w 100 -b 70000 -t 0 -o $1/Rho0.0009/sfs/Rho0.0009.{}" ::: $(seq 1 1000)
