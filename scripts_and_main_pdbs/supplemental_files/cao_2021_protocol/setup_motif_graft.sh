#!/bin/bash

prepare_run=$2
motifs_file=$1
shift
shift

mkdir motif_runs
rm motif_runs/*.list > /dev/null 2>&1

IFS=$'\n'

for line in $(cat $motifs_file); do
    motif=$(echo $line | awk '{print $1}')
    hotspots=$(echo $line | awk '{print $2}')

    motif_tag=$(basename $motif .pdb.gz)

    $prepare_run -flags="-parser:script_vars motifpdb=$motif hotspots=$hotspots -out:prefix ${motif_tag}_" -destination motif_runs -suffix _$motif_tag "$@"

done

unset IFS

cat motif_runs/*.list > motif_runs_commands.list



