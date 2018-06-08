#!/bin/bash

# usage: 'fastaflatter.sh input.fasta > output_flat.fasta'

less $1 | awk -v RS='>'     -v FS="\n"     -v OFS=""     -v ORS="" '{ if (NR > 1) { printf ">%s\n",$1; $1=""; printf "%s\n",$0 } }'