#!/bin/sh

root_dir=$(pwd)
date=$(date +"%d-%m-%Y")

if [ ! -d "test2_${date}" ]; then
	dir="$root_dir/test2_${date}"
else
	n=1
	while [ -d "${root_dir}/test2_${date}_${n}" ]; do
		n=$((n + 1))
	done
	dir="$root_dir/test2_${date}_${n}"
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
mv "$root_dir"/coall_* "$dir"/
cd "$dir"

resfile="$dir/res_coall.txt"
: > "$resfile"

write_section() {
	printf '\n#### %s\n\n' "$1" >> "$resfile"
}

write_section "Phylo.jl"
julia -t 1 coall_phylov2.jl >> "$resfile" 2>&1

write_section "toytree"
python3 coall_toytree.py  >> "$resfile" 2>&1

write_section "cogent3"
python3 coall_cogent.py  >> "$resfile" 2>&1

write section "castor"
Rscript coall_castor.r "$(pwd)" >> "$resfile" 2>&1

cd "$root_dir"
