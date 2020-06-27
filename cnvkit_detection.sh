# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cnvkit_env

REF="/media/Data1/jbogoin/ref/hg38_Mlast/hg38_GenDev.fa"
TARGET="/media/Data1/jbogoin/ref/gencode/v34_hg38/gencode.v34.basic.annotation.CDS.merged.4fields.interval_list"
DATA=$PWD

echo ""
echo "CNVKIT DETECTION start"
echo ""

rm -Rf cnvkit_output
mkdir  cnvkit_output

SAMPLES=$(ls *.dedup.bam)

#create a pooled reference by running batch command specifying only normal samples
~/cnvkit/cnvkit.py batch -n $SAMPLES \
	-m hybrid \
    -f $REF \
	--targets $TARGET \
	--output-reference $DATA/cnvkit_output/pooled-reference.cnn \
	--short-names \
	-d $DATA/cnvkit_output \
	-p 12

# Run WGS batch pipeline using the pooled reference file
~/cnvkit/cnvkit.py batch $SAMPLES \
	-m hybrid \
	-r $DATA/cnvkit_output/pooled-reference.cnn \
	-d $DATA/cnvkit_output \
	-p 12 \
	--scatter \
	--diagram

rm $DATA/cnvkit_output/*.target.bed;

echo ""
echo "CNVKIT DETECTION job done!"
echo ""

