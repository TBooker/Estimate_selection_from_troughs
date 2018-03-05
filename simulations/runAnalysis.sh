for i in Exp*
#for i in Nes200_dDFE.Exponential Nes200_dDFE.Exponential_div10
do
echo $i
python ../minimisationInPython/BootstrapMouseTroughFit.py -i $i/full.Rho0.009.b.csv.gz $i/full.Rho0.0045.b.csv.gz  $i/full.Rho0.0009.b.csv.gz -b ../BGS/ExpDFE_24.25/Rho0.009.csv ../BGS/ExpDFE_24.25/Rho0.0045.csv ../BGS/ExpDFE_24.25/Rho0.0009.csv -r 1.0 0.5 0.1 --boots 100 --model e -o $i.analysis.txt --label $i --bgs_null  &
done
#full.Rho0.009.b.csv

