# Author: Julien BURATTI

source ~/miniconda3/etc/profile.d/conda.sh
conda activate fastq_bam_env

echo ""
echo "fastq_to_bam.sh start"
echo ""

REF="/media/Data1/jbogoin/ref/hg38_Mlast/hg38_GenDev.fa"

### REMOVE previous dedup.bam files ###
rm -f *.bam*
rm -Rf multiqc*

## MAPPING BWA SEQUENTIAL 
for R1 in *_R1_*.fastq.gz; \
do R2=${R1/_R1_/_R2_}; \
SAMPLE=${R1%%_*}; \
FLOWCELL="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 3)"; \
DEVICE="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 1 | cut -d "@" -f 2)"; \
BARCODE="$(zcat $R1 | head -1 | awk '{print $2}' | cut -d ":" -f 4)"; \
RG=$(echo "\"@RG\tID:${DEVICE}.${FLOWCELL}.${SAMPLE}\tPU:${FLOWCELL}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}-${BARCODE}\"") \
MAPPING_CMD=$(echo "time bwa mem -M -t 18 -R $RG $REF $R1 $R2 | samtools sort -@ 4 -o ${SAMPLE}.bam -"); \
eval $MAPPING_CMD;\
done;

### MARK DUPLICATES
time for i in *.bam; do SAMPLE=${i%.*}; time sambamba markdup -t 12 ${SAMPLE}.bam ${SAMPLE}.dedup.bam; done

echo ""
echo "fastq_to_bam.sh job done!"
echo ""
