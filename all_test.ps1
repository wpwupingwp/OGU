# mac
$folder="rbcL_Poaceae","Rosa_ITS","Lamiaceae_cp","zea_mays","Oryza_cp"
foreach ($f in $folder)
{
    Remove-Item -Recurse -Force R:\test\$f
}
New-Item -p R:\test
$seq_n="101"
$python="python"
$out="R:\test\"
&$python -m BarcodeFinder.gb2fasta -gene rbcL -taxon Poaceae -out ${out}rbcL_Poaceae -seq_n $seq_n
&$python -m BarcodeFinder.gb2fasta -query "internal transcribed spacer" -taxon Rosa -out ${out}Rosa_ITS -uniq no -seq_n $seq_n
&$python -m BarcodeFinder -og cp -refseq yes -taxon Lamiaceae -out ${out}Lamiaceae_cp -seq_n $seq_n
&$python -m BarcodeFinder -taxon Zea_mays -min_len 100 -max_len 3000 -out ${out}Zea_mays -primer -seq_n $seq_n
&$python -m BarcodeFinder -taxon Oryza -og mt -min_len 50000 -max_len 200000 -uniq longest -out ${out}Oryza_cp -refseq yes -primer -seq_n $seq_n
