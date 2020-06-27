# Author: Julie BOGOIN

echo ""
echo "REFERENCE INDEXATION start"
echo ""

echo ""
echo "Merci d'indiquer le nom de l'archive .fa.gz à traiter (sans extension):"
read reference 

gunzip $reference.fa.gz

bwa index -a bwtsw $reference.fa

gatk CreateSequenceDictionary -R $reference.fa -O $reference.dict

samtools faidx $reference.fa

echo ""
echo "REFERENCE INDEXATION job done!"
echo ""