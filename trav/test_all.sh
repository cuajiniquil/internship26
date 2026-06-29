#!/bin/sh

root_dir=$(pwd)
date=$(date +"%d-%m-%Y")

if [ ! -d "test_${date}" ]; then
	dir="$root_dir/test_${date}"
else
	n=1
	while [ -d "${root_dir}/test_${date}_${n}" ]; do
		n=$((n + 1))
	done
	dir="$root_dir/test_${date}_${n}"
fi

cleanup() {
	for path in "$dir"/*k.nwk; do
		[ -e "$path" ] || continue
		mv "$path" "$root_dir"/
	done

	for path in "$dir"/traverse_*; do
		[ -e "$path" ] || continue
		mv "$path" "$root_dir"/
	done
}

trap cleanup EXIT

mkdir "$dir"
mv "$root_dir"/*k.nwk "$dir"/
mv "$root_dir"/traverse_* "$dir"/
cd "$dir"

resfile="$dir/res_traverse.txt"
: > "$resfile"

write_section() {
	printf '\n#### %s\n\n' "$1" >> "$resfile"
}

write_section "biopython"
python3 traverse_biopython.py >> "$resfile" 2>&1

write_section "cogent3"
python3 traverse_cogent3.py >> "$resfile" 2>&1

write_section "dendropy"
python3 traverse_dendropy.py >> "$resfile" 2>&1

write_section "ete3"
python3 traverse_ete3.py >> "$resfile" 2>&1

write_section "toytree"
python3 traverse_toytree.py >> "$resfile" 2>&1

write_section "Phylo.jl"
julia -t 1 traverse_phylo.jl >> "$resfile" 2>&1

write_section "ape"
Rscript traverse_ape.r "$(pwd)" >> "$resfile" 2>&1

write_section "castor"
Rscript traverse_castor.r "$(pwd)" >> "$resfile" 2>&1

#write_section "traverse_phylors"
#( cd traverse_phylors && cargo run ) >> "$resfile" 2>&1

cd "$root_dir"


