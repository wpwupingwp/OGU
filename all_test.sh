# mac
rm -rf rbcL_Poaceae Rosa_ITS Lamiaceae_cp zea_mays Oryza_cp
seq_n=11
python=python3
out=./
$python -m OGU.gb2fasta -gene rbcL -taxon Poaceae -out ${out}rbcL_Poaceae -seq_n $seq_n
$python -m OGU.gb2fasta -query "internal transcribed spacer" -taxon Rosa -out ${out}Rosa_ITS -uniq no -seq_n $seq_n
$python -m OGU -og cp -refseq yes -taxon Lamiaceae -out ${out}Lamiaceae_cp -seq_n $seq_n
$python -m OGU -taxon Zea_mays -min_len 100 -max_len 3000 -out ${out}Zea_mays -primer -seq_n $seq_n
$python -m OGU -taxon Oryza -og mt -min_len 50000 -max_len 200000 -uniq longest -out ${out}Oryza_cp -refseq yes -primer -seq_n $seq_n
