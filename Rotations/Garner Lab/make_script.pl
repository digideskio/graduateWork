#!/usr/bin/perl

#use strict;
#use warnings "all";
use POSIX;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
use Thread;
use Cwd;

# I could have technically looped through
# each downloaded genome automatically,
# but given constraints on the computing
# resources I can tie up at once,
# there would have needed to be
# a daisy chain from hell of jobs
# waiting to execute pending completion
# of previous jobs.
#
# This seemed a simpler solution.
GetOptions("input=s" => \$input);

# Specify how qsub script will be named
my $qsub = "/home/wetherc/results/$input/qsub.$input.sh";
my $out;

# Open the qsub script
open($out, ">$qsub") || die "Could not create a qsub script. Do I have proper write permissions?\n";

print $out "#!/bin/bash\n";
print $out "#PBS -W group_list=sfx\n";
print $out "#PBS -q sfx_q\n";
print $out "#PBS -N $input\n";
print $out "#PBS -r y\n";
print $out "#PBS -j oe\n";
print $out "#PBS -l walltime=100:00:00\n";
print $out "#PBS -l nodes=1:ppn=5\n";
print $out "#PBS -d /home/wetherc/results/$input\n";

# Define run_cmd function
print $out "run_cmd () {\n";
print $out "    echo; echo `date`; echo \$cmd\n";
print $out "    if !(eval \$cmd); then exit 1; fi;\n";
print $out "}\n\n";

# Load all requisite modules
print $out "module load bio/samtools\n";
print $out "module load bio/bwa\n";
print $out "module load bio/picard\n";
print $out "module load devel/java\n";

# Unzip everything
print $out "cmd=\"gunzip *.gz\"\n";
print $out "run_cmd\n\n";

# Concatenate all of our forward read pairs
print $out "cmd=\"cat *_1.filt.fastq > $input.pairs_1.fastq\"\n";
print $out "run_cmd\n\n";

# Concatenate all of our reverse read pairs
print $out "cmd=\"cat *_2.filt.fastq > $input.pairs_2.fastq\"\n";
print $out "run_cmd\n\n";

# Remove old fastq files
print $out "cmd=\"rm *.filt.fastq\"\n";
print $out "run_cmd\n\n";

# Stitch two fastq files into an output sam file
print $out "cmd=\"bwa mem -aY -t 5 -P -v 3 -U 25 /home/wetherc/reference/hg19.fa /home/wetherc/results/$input/$input.pairs_1.fastq /home/wetherc/results/$input/$input.pairs_2.fastq  > $input.sam\"\n";
print $out "run_cmd\n\n";

# Purge all (old) fastq files
print $out "cmd=\"rm *.fastq\"\n";
print $out "run_cmd\n\n";

# Extract/print all alignments or subalignments in SAM format
# -h : include header information
# -S : input is in SAM format
# -b : output in BAM format
# 
# Sort alignments by leftmost coordinates
print $out "cmd=\"samtools view -hS $input.sam -b | samtools sort - $input\"\n";
print $out "run_cmd\n\n";

# Remove PCR duplicates with samtools tools if it's not PCR free:
print $out "cmd=\"samtools rmdup $input.bam $input.rdup.bam\"\n";
print $out "run_cmd\n\n";

# Index sorted alignment for fast random access
# Creates .bai index file
print $out "cmd=\"samtools index $input.rdup.bam\"\n";
print $out "run_cmd\n\n";

# Replace all read groups in the INPUT file with a new read group and
# assigns all reads to this read group in the OUTPUT
print $out "cmd=\"java -Xmx5g -jar /apps/packages/bio/picard/1.92/bin/AddOrReplaceReadGroups.jar I=$input.rdup.bam O=$input.rdup.GATK.bam SORT_ORDER=coordinate RGID=$input RGLB=bwa-mem RGPL=illumina RGSM=$input CREATE_INDEX=true RGPU=RGPU VALIDATION_STRINGENCY=SILENT\"\n";
print $out "run_cmd\n\n";

# Clean up old files
print $out "cmd=\"rm *.sam\"\n";
print $out "run_cmd\n\n";

# Extract/print all alignments or subalignments in SAM format
# -u : output to uncompressed bam
# -f : only output alignments with all bits in INT present in the FLAG field
#
# Output to sam file containing only unmapped reads from the original genome
print $out "cmd=\"samtools view -u -f 4 $input.rdup.GATK.bam > $input.unmapped.sam\"\n";
print $out "run_cmd\n\n";

# Clean up old files
print $out "cmd=\"rm *.ba*\"\n";
print $out "run_cmd\n\n";

# Move the file to its appropriate location
print $out "cmd=\"mv $input.unmapped.sam /home/wetherc/unmapped/sam/$input.unmapped.sam\"\n" ;
print $out "run_cmd\n\n";

######################################################

print $out "echo end $outputFile\n";
close($out);

# Submit script to job manager
print "$qsub\n";
system("chmod gu+x $qsub\n");
system("qsub -m bae -M wetherc\@vbi.vt.edu /home/wetherc/results/$input/qsub.$input.sh\n");
print "The job was submitted. The final output will be '$input.unmapped.Bam.file'.\n";