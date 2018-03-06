## The only input that this script needs is the directory for the data
cp $1/pointEstimate/sfs.txt ./


number1=100
RANGE=1000
number2=$RANDOM
number3=$RANDOM
number4=$RANDOM
number2=123
number3=321
number4=554

let 'number2%=RANGE'
let 'number3%=RANGE'
let 'number4%=RANGE'

# Do the 1-epoch model
# 
# ~/project/8.SFS_inference/dfe-alpha-release-2.16/est_dfe << EOF
# 0
# 0
# 1
# 100
# 0
# EOF
# 
# mv est_dfe.out $1/pointEstimate/1-epoch/
# mv neut_* $1/pointEstimate/1-epoch/
# ~/
# Do the 2-epoch model
# 
# ~/project/8.SFS_inference/dfe-alpha-release-2.16/est_dfe << EOF
# 0
# 0
# 2
# 100
# $number2
# 1
# $number3
# 0
# 1
# 1
# 1000
# EOF
# 
# mv est_dfe.out $1/pointEstimate/2-epoch/
# mv neut_* $1/pointEstimate/2-epoch/

Do the 3-epoch model

~/project/8.SFS_inference/dfe-alpha-release-2.16/est_dfe << EOF
0
0
3
100
$number2
1
$number3
$number4
1
$number5
0
EOF

mv est_dfe.out $1/pointEstimate/3-epoch/
mv neut_* $1/pointEstimate/3-epoch/



