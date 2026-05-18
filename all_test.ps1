# mac
$folder="rbcL_Poaceae","Rosa_ITS","zea_mays","Oryza_cp", "Lamiaceae_cp"
foreach ($f in $folder)
{
    Remove-Item -Recurse -Force R:\test\$f
}
New-Item R:\test -ItemType Directory
$count="4"
$python="python"
$out="R:\test\"
&uv run $python -m OGU.gb2fasta -gene rbcL -taxon Poaceae -out ${out}rbcL_Poaceae -count ${count}
&uv run $python -m OGU.gb2fasta -query "internal transcribed spacer" -taxon Rosa -out ${out}Rosa_ITS -uniq no -count $count
&uv run $python -m OGU -taxon Zea_mays -min_len 100 -max_len 3000 -out ${out}Zea_mays -primer -count 10
&uv run $python -m OGU -taxon Oryza -og mt -min_len 50000 -max_len 200000 -uniq longest -out ${out}Oryza_cp -refseq yes -primer -count $count
&uv run $python -m OGU -og cp -refseq yes -taxon Lamiaceae -out ${out}Lamiaceae_cp -count $count
