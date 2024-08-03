#!/bin/bash

input_fasta=$1  # First argument is the input FASTA file

awk '
BEGIN { FS = "\n"; RS = ">" }
NR > 1 {
    header = $1
    sequence = ""
    for (i = 2; i <= NF; i++) {
        sequence = sequence $i
    }
    print header, length(sequence)
}
' "$input_fasta"