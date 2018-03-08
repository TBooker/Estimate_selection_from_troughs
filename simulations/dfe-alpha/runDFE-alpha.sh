for i in Nes200_dDFE.Exponential_div10 Nes200_dDFE.Exponential Nes10_dDFE.Exponential_div10 Nes10_dDFE.Exponential 
do
#mkdir $i
#mkdir $i/pointEstimate
#mkdir $i/pointEstimate/1-epoch
#mkdir $i/pointEstimate/2-epoch
#mkdir $i/pointEstimate/3-epoch
#python ../../bin/generateInputFiles.py -i ../$i/$i.func.sfs -b point -o $i/pointEstimate/sfs.txt -s m6
#echo $i
#sh ../../bin/runDFEalphaOnSFS.sh $i
mkdir $i/selected

    for j in $(seq 1 10)
    do
    python ../../bin/generateInputFiles.py -i ../$i/$i.func.sfs -b $j -o sfs.txt -s m6
    # Use the three epoch model  to correct the sfs
    ~/PK_modules/correct-sfs.pl -e $i/pointEstimate/3-epoch/neut_exp_obs_allele_freqs.csv -s sfs.raw.txt -o sfs.txt
    # Run dfe-alpha 5 times on the data to achieve convergence
    python  ~/project/8.SFS_inference/bin/dfe_alpha_random_starts.py -n $i/pointEstimate/3-epoch/est_dfe.out -d 0 -a 1 -r 5 -o TEMP_DFE --gamma --path_to ~/project/8.SFS_inference/dfe-alpha-release-2.16/est_dfe
	# Now extract the best fitting model from the
	python ~/project/8.SFS_inference/bin/get_ML_runs_one_set.py TEMP_DFE.dict.pkl $i/selected/$j.est_dfe.out
    done
done


for i in BimodalDFE BimodalDFE_div10 BimodalDFE_div100
do
mkdir $i/selected
    for j in $(seq 1 10)
    do
    python ../../bin/generateInputFiles.py -i ../$i/$i.func.sfs -b $j -o sfs.txt -s m6
    # Use the three epoch model  to correct the sfs
    ~/PK_modules/correct-sfs.pl -e $i/pointEstimate/3-epoch/neut_exp_obs_allele_freqs.csv -s sfs.raw.txt -o sfs.txt
    # Run dfe-alpha 5 times on the data to achieve convergence
    python  ~/project/8.SFS_inference/bin/dfe_alpha_random_starts.py -n $i/pointEstimate/3-epoch/est_dfe.out -d 0 -a 2 -r 5 -o TEMP_DFE --gamma --path_to ~/project/8.SFS_inference/dfe-alpha-release-2.16/est_dfe
	# Now extract the best fitting model from the
	python ~/project/8.SFS_inference/bin/get_ML_runs_one_set.py TEMP_DFE.dict.pkl $i/selected/$j.est_dfe.out
    done
done
