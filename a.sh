#!/bin/bash
#wpwupingwp@outlook.com
#2014-01-19
mv *.fasta 01.fasta
sed -i 's/ /_/g' 01.fasta
mkdir 02_cluster
./usearch -cluster_fast ./01.fasta -id 0.6 -clusters ./02_cluster/
cd 02_cluster
for i in *
do
		echo $i>>list.txt
		head -1 $i>>list.txt
done
cd ..
mv ./02_cluster/list.txt ./
cp -r 02_cluster 03_rename
cd 03_rename
rename 's/$/\.fasta/' *
for i in *.fasta 
do 
	awk -F '_' '{if($1~/gi/) print $1,$2,$3;else print $0}' $i>$i.1
	awk -F '|' '{if($1~/gi/) printf ">%s_%s_%s\n",$5,$3,$4;else print $0}' $i.1>$i.2
	sed -i 's/\n//' $i.2
done
rm *.fasta
rm *.1
rename 's/\.2//' *.2
cd ..
cp -r 03_rename 04_alignment
cd 04_alignment
for j in *.fasta
do
		mafft --reorder --adjustdirection --auto $j>$j.aln
		sed -i 's/ /_/g' $j
done
rm *.fasta
rename 's/\.aln//' *.aln
cd ..
exit
