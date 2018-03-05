for i in *DFE*
do
echo $i

python /home/booker/project/15.Estimate_Selection_From_Trough/bin/funcSfsFromSims.py -i $i/Rho0.009/runs -o $i/Rho0.009.sfs

done
