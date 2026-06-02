#!/bin/sh

date=$(date +"%d-%m-%Y")

if [ ! -d "test_${date}" ]; then
	dir="test_${date}"
else
	n=1
	while [ -d "test_${date}_${n}" ]; do
		((n++))
	done
	dir="test_${date}_${n}"
fi

mkdir $dir
mv *k.nwk $dir/
mv traverse_* $dir/
cd $dir

touch $dir.txt

python3 traverse_biopython.py > biopython.txt
python3 traverse_cogent3.py > cogent3.txt
python3 traverse_dendropy.py > dendropy.txt
python3 traverse_ete3.py > ete3.txt
python3 traverse_toytree.py > toytree.txt

julia -t 1 traverse_phylo.jl > phylojl.txt

Rscript traverse_ape.r "$(pwd)" > ape.txt
Rscript traverse_castor.r "$(pwd)" > castor.txt

mv *k.nwk ..

cd traverse_phylors
cargo run

mv traverse_* ..


