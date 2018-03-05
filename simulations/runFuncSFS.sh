for i in *DFE*
do
echo $i
echo Rho0.009
python /home/booker/project/15.Estimate_Selection_From_Trough/bin/funcSfsFromSimsDirectory.py -i $i/Rho0.009/runs -o $i/Rho0.009.sfs
echo Rho0.0045
python /home/booker/project/15.Estimate_Selection_From_Trough/bin/funcSfsFromSimsDirectory.py -i $i/Rho0.0045/runs -o $i/Rho0.0045.sfs
echo Rho0.0009
python /home/booker/project/15.Estimate_Selection_From_Trough/bin/funcSfsFromSimsDirectory.py -i $i/Rho0.0009/runs -o $i/Rho0.0009.sfs
done
