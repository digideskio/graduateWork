#!/bin/bash

module load bio/samtools

# For each sam file in our directory
for file in /home/wetherc/unmapped/sam/*.unmapped.sam
do
	# Extract the genome name from the full file name;
	# Print it to output
	printf "${file:36:7}\t" >> output.txt

	# Count the number of unmapped reads for that
	# genome; print it to output
	samtools view -c -fox4 $file >> output.txt
done