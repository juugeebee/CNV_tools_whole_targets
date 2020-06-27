#!/bin/bash

# Author: Julie BOGOIN

echo ""
echo "pipeline_csv.sh start"
echo ""

source ~/miniconda3/etc/profile.d/conda.sh

# sex determination
conda activate sex_env

python3 ~/WCC/sex_determination.py

conda deactivate

# tools cnv detection
bash ~/WCC/cn.mops_launch.sh
bash ~/WCC/cnvkit_detection.sh
bash ~/WCC/excavator2_detection.sh
bash ~/WCC/exomedepth_launch.sh
bash ~/WCC/gatk4_CNV_detection.sh

# tools results generation
conda activate results_cnv

cd cn.mops_output
python ~/WCC/cn.mops_results.py
cd cnvkit_output
python ~/WCC/cnvkit_results.py
cd excavator2_output
python ~/WCC/excavator2_results.py
cd exomedepth_output
python ~/WCC/exomedepth_results.py
cd gatk_output
python ~/WCC/gatk_results.py
cd ..

# results summary
python ~/WCC/cnv_results.py
python ~/WCC/cnv_interval_object_sample.py
python ~/WCC/cnv_interval_object_run.py

conda deactivate

# annotation
conda activate annot_env

bash ~/WCC/annotation.sh

conda deactivate

echo ""
echo "pipeline_csv.sh job done!"
echo ""

