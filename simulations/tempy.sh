for i in Exp.Nes10_*/
do
echo $i
python ~/bin/analyse_sfs_from_slim_directory.py -i $i/Rho0.009/sfs/ -o $i/full.Rho0.009 -r 1000&
python ~/bin/analyse_sfs_from_slim_directory.py -i $i/Rho0.0045/sfs/ -o $i/full.Rho0.0045 -r 1000&
python ~/bin/analyse_sfs_from_slim_directory.py -i $i/Rho0.0009/sfs/ -o $i/full.Rho0.0009 -r 1000

done
