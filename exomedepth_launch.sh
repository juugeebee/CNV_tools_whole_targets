# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh
conda activate exomedepth_env

echo ""
echo "ExomeDepth CNV DETECTION start"
echo ""

rm -Rf exomedepth_output
mkdir exomedepth_output

R -q --vanilla < ~/WCC/exomedepth_cnv.r

echo ""
echo "ExomeDepth CNV DETECTION job done!"
echo ""
