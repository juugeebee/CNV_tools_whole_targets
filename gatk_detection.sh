# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gatk_env

REF="/media/Data1/jbogoin/ref/hg38_Mlast/hg38_GenDev.fa"
DIC="/media/Data1/jbogoin/ref/hg38_Mlast/hg38_GenDev.dict"
TARGETS="/media/Data1/jbogoin/ref/gencode/v34_hg38/gencode.v34.basic.annotation.CDS.merged.4fields.interval_list"
PLOIDY_PRIORS_TABLE="/media/Data1/jbogoin/ref/contig_ploidy_priors/ploidy_priors_table.tsv"
CENTROMETIC="/media/jbogoin/Data1/jbogoin/ref/centromeric_regions.bed"

echo ""
echo "GATK4 CNV DETECTION start"
echo ""

rm -rf gatkcnv_output
mkdir gatkcnv_output

# Define the resolution of the analysis with a genomic intervals list
gatk PreprocessIntervals \
    -R  $REF \
    -L  $TARGETS \
    -XL $CENTROMETIC \
    --bin-length 0 \
    --padding 50 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O gatkcnv_output/targets.preprocessed.interval_list \
    --verbosity ERROR

for sample_id in *.dedup.bam; 
do SAMPLE=${sample_id%%.dedup.bam}; \

# Collect raw integer counts data
gatk CollectReadCounts \
-L gatkcnv_output/targets.preprocessed.interval_list \
    -XL $CENTROMETIC \
    -R $REF \
    --interval-merging-rule OVERLAPPING_ONLY \
    -I $SAMPLE.dedup.bam \
    --format TSV \
    -O gatkcnv_output/$SAMPLE.tsv \
    --verbosity ERROR;

done

# AnnotateIntervals with GC content
gatk AnnotateIntervals \
    -L gatkcnv_output/targets.preprocessed.interval_list \
    -XL $CENTROMETIC \
    -R $REF \
    -imr OVERLAPPING_ONLY \
    -O gatkcnv_output/targets.annotated.tsv \
    --verbosity ERROR


cd gatkcnv_output

COUNTS_LIST=""
for counts_file in *.tsv; do 
        if [ $(grep -v targets.annotated.tsv <<<$counts_file) ]; then
        FILE=${counts_file%%.tsv};
        COUNTS_LIST+="-I $FILE.tsv ";
        fi 
done

# FilterIntervals based on GC-content and cohort extreme counts
gatk FilterIntervals \
        -L targets.preprocessed.interval_list \
        -XL $CENTROMETIC \
        --annotated-intervals targets.annotated.tsv \
        $COUNTS_LIST \
        -imr OVERLAPPING_ONLY \
        -O targets.cohort.gc.filtered.interval_list \
        --verbosity ERROR


# DetermineGermlineContigPloidy in COHORT MODE
gatk DetermineGermlineContigPloidy \
        -L targets.cohort.gc.filtered.interval_list \
        --interval-merging-rule OVERLAPPING_ONLY \
        $COUNTS_LIST \
        --contig-ploidy-priors $PLOIDY_PRIORS_TABLE \
        --output . \
        --output-prefix ploidy \
        --verbosity ERROR


# GermlineCNVCaller in COHORT MODE
gatk GermlineCNVCaller \
        --run-mode COHORT \
        -L targets.cohort.gc.filtered.interval_list\
        $COUNTS_LIST \
        --contig-ploidy-calls ploidy-calls \
        --annotated-intervals targets.annotated.tsv \
        --interval-merging-rule OVERLAPPING_ONLY \
        --output . \
        --output-prefix cohort \
        --verbosity ERROR

cd ..

index=0

for sample_id in *.dedup.bam; 
do SAMPLE=${sample_id%%.dedup.bam};

# PostprocessGermlineCNVCalls COHORT MODE
gatk PostprocessGermlineCNVCalls \
        --model-shard-path gatkcnv_output/cohort-model \
        --calls-shard-path gatkcnv_output/cohort-calls \
        --allosomal-contig chrX --allosomal-contig chrY \
        --contig-ploidy-calls gatkcnv_output/ploidy-calls \
        --sample-index $index \
        --output-genotyped-intervals gatkcnv_output/genotyped-intervals.$SAMPLE.vcf.gz \
        --output-genotyped-segments gatkcnv_output/genotyped-segments.$SAMPLE.vcf.gz \
        --output-denoised-copy-ratios gatkcnv_output/denoised-copy-ratios.$SAMPLE.tsv \
        --sequence-dictionary $DIC \
        --verbosity ERROR;

let "index+=1";

done

cd 

echo ""
echo "GATK4 CNV DETECTION job done!"
echo ""
