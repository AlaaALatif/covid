#!/bin/bash

# Assume path split by "/" and sample name if of format 46idSEARCH0139NBG_S46_L001_R2_001.fastq.gz

echo -e "forward\treverse\tsample\tsample_library\tsequencing_tech"
indir=$1
find $indir -name "*.fastq.gz" | sort | xargs -n 2 | tr ' ' '\t' | awk 'BEGIN{FS="\t";OFS="\t"}{split($1, n, "/");split(n[6], s, "_");print $1,$2,s[1],s[1]"_"s[3],"illumina"}END{}'
