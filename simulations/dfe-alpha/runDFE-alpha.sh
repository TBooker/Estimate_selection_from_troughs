for i in Nes200_dDFE.Exponential_div10 Nes200_dDFE.Exponential Nes10_dDFE.Exponential_div10 Nes10_dDFE.Exponential Exp.Nes10_dDFE.Exponential_div10 Exp.Nes10_dDFE.Exponential Exp.Nes100_dDFE.Exponential_div10 Exp.Nes100_dDFE.Exponential
do

mkdir $i
mkdir $i/pointEstimate
mkdir $i/pointEstimate/1-epoch
mkdir $i/pointEstimate/2-epoch
mkdir $i/pointEstimate/3-epoch
python ../../bin/generateInputFiles.py -i ../$i/$i.func.sfs -b point -o $i/pointEstimate/sfs.txt -s m6
echo $i
sh ../../bin/runDFEalphaOnSFS.sh $i

#    for j in $(seq 1 10)
#    do
#    python ../../bin/generateInputFiles.py -i ../$i/$i.func.sfs -b $j -o sfs.txt -s m7

#    done
done
