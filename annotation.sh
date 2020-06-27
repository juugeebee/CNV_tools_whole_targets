# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh
conda activate annot_env

echo ""
echo "annotation.sh start"
echo ""

rm -Rf annovar_output
mkdir -p annovar_output

## Download annotation databases from ANNOVAR or UCSC and save to humandb/ directory
#perl ./annotate_variation.pl -downdb avdblist
#perl annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene humandb

## Retirer les noms de colonnes au fichier interval_results.txt > /annovar_output/ex1.avinput
sed '1d' interval_run_results.txt > ex1.avinput

## Remplacer XF par X et XM par X
sed -i -e "s/XF/X/g" ex1.avinput
sed -i -e "s/XM/X/g" ex1.avinput

sudo perl ~/annovar/annotate_variation.pl \
    -out annovar_output/annotation \
    -build hg38 ex1.avinput\
    ~/annovar/humandb/

 sudo mv ex1.avinput annovar_output/    

#first output file: .variant_function
#contains annotation for all variants, by adding two columns to the beginning of each input line

#second output file: .exonic_variant_function
#contains the amino acid changes as a result of the exonic variant. 

echo ""
echo "annotation.sh job done!"
echo ""




