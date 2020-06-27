# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh
conda activate fastq_bam_env

echo ""
echo "bam_to_fastq.sh start"
echo ""

for bam_name in *.dedup.bam; \

do SAMPLE=${bam_name%%.dedup.bam} \

samtools bam2fq $bam_name -@12 > $SAMPLE.fastq; 

# split a single .fastq file of paired-end reads into two separated files
# extracting reads ending with '/1' or '/2'
cat $SAMPLE.fastq | grep '^@.*/1$' -A 3 --no-group-separator > ${SAMPLE}_R1_001.fastq;
gzip ${SAMPLE}_R1_001.fastq;
cat $SAMPLE.fastq | grep '^@.*/2$' -A 3 --no-group-separator > ${SAMPLE}_R2_001.fastq;
gzip  ${SAMPLE}_R2_001.fastq;

rm $SAMPLE.fastq;

done

echo ""
echo "bam_to_fastq.sh job done!"
echo ""
