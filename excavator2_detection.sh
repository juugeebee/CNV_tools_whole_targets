# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh
conda activate excavator2_env

TARGETS="/media/Data1/jbogoin/ref/gencode/v34_hg38/gencode.v34.basic.annotation.CDS.merged.bed"
DATA=$PWD

echo ""
echo "EXCAVATOR2 CNV DETECTION start"
echo ""

rm -Rf excavator2_output/
mkdir excavator2_output

# rm -Rf ~/EXCAVATOR2_Package_v1.1.2/data/targets/hg38/MyTarget_w10K/

# BAM list file
ls *.dedup.bam > $DATA/excavator2_output/bam_list.txt

# # Target initialization
# cd ~/EXCAVATOR2_Package_v1.1.2
# perl TargetPerla.pl SourceTarget.txt $TARGETS MyTarget_w10K 10000 hg38

# # Data Prepare Module
cd $DATA

for sample_id in *.dedup.bam; \
do SAMPLE=${sample_id%%.dedup.bam}

echo $DATA/$SAMPLE.dedup.bam $DATA/excavator2_output/$SAMPLE $SAMPLE >> $DATA/excavator2_output/ExperimentalFilePrepare.w10K.txt;

done

# RC calculations
cd ~/EXCAVATOR2_Package_v1.1.2
perl EXCAVATORDataPrepare.pl $DATA/excavator2_output/ExperimentalFilePrepare.w10K.txt --processors 12 --target MyTarget_w10K --assembly hg38

# Experimental analysis file
labels=("C1" "C2" "C3" "C4" "C5" "C6" "C7" "C8" "C9" "C10" "C11")

rm $DATA/excavator2_output/ExperimentalFileAnalysis.w10K.* 

cd $DATA

for sample_id in *.dedup.bam; do 
	SAMPLE=${sample_id%%.dedup.bam};	
	
	echo 'T1' $DATA/excavator2_output/$SAMPLE $SAMPLE > $DATA/excavator2_output/ExperimentalFileAnalysis.w10K.$SAMPLE.txt;
	
	others=$(ls *.dedup.bam | grep -v $SAMPLE);
	
	# TRansformer la liste others en tableau indicable table
	table=( ${others// / } )
	table_clean=${table[@]/.dedup.bam/}
	NORMAL=( ${table_clean// / } )

	for i in `seq 0 10`; do
		
		echo ${labels[$i]} $DATA/excavator2_output/${NORMAL[$i]} ${NORMAL[$i]}\
			>> $DATA/excavator2_output/ExperimentalFileAnalysis.w10K.$SAMPLE.txt;
	
	done

	# Segmentation of the WMRC
	cd ~/EXCAVATOR2_Package_v1.1.2;
	
	perl EXCAVATORDataAnalysis.pl $DATA/excavator2_output/ExperimentalFileAnalysis.w10K.$SAMPLE.txt\
		--processors 12 --target MyTarget_w10K\
	       	--assembly hg38\
		--output $DATA/excavator2_output/w10K_results.$SAMPLE\
	       	--mode pooling;
	
	cd $DATA	

done

echo ""
echo "EXCAVATOR2 CNV DETECTION job done!"
echo ""
