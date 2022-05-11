# mac
seq_n=100
python3 -m BarcodeFinder.gb2fasta -gene rbcL -taxon Poaceae -out rbcL_Poaceae -seq_n 2000
python3 -m BarcodeFinder.gb2fasta -query internal transcribed spacer -taxon Rosa -out Rosa_ITS -uniq no -seq_n 2000
python3 -m BarcodeFinder -og cp -refseq -taxon Lamiaceae -out Lamiaceae_cp -seq_n 2000
python3 -m BarcodeFinder -taxon Zea_mays -min_len 100 -max_len 3000 -out Zea_mays -primer -seq_n 2000
python3 -m BarcodeFinder -taxon Oryza -og mt -min_len 50000 -max_len 200000 -uniq longest -out Oryza_cp -refseq yes -primer -seq_n 2000
