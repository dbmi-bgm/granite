#!/bin/bash

# variables from command line
inputbam=$1
reference=$2
chromosomefile=$3
nthreads=$4

# self variables
directory=mpileups/

# setting up output directory
mkdir -p $directory

# creating basic command
command="mpileup_parser.py -i $inputbam -r $reference --region {} -o ${directory}bam_{}.out"

# running command
cat $chromosomefile | parallel --halt 2 --jobs $nthreads $command || exit 1

# merging the results
array=(${directory}*.out)

IFS=$'\n' sorted=($(sort -V <<<"${array[*]}"))
unset IFS

grep "^#" ${sorted[0]} > mpileup.out

for filename in ${sorted[@]};
  do
    if [[ $filename =~ "M" ]]; then
      chr_M=$filename
    else
      grep -v "^#" $filename >> mpileup.out
      rm -f $filename
    fi
  done

if [[ -v  chr_M  ]]; then
  grep -v "^#" $chr_M >> mpileup.out
  rm -f $chr_M
fi

# compress and index mpileup.out
bgzip mpileup.out
tabix -p vcf mpileup.out.gz
