# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cn.mops_cnv

echo ""
echo "cn.mops CNV DETECTION start"
echo ""

rm -Rf cn.mops_output
mkdir cn.mops_output

R -q --vanilla < ~/WCC/cn.mops_cnv.r

echo ""
echo "cn.mops CNV DETECTION job done!"
echo ""