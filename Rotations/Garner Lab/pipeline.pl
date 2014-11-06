#!/usr/bin/perl

#use strict;
#use warnings "all";
use POSIX;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
use Thread;
use Cwd;

my ($cmd, $helpFlag, $noGapFlag, $minRate_indelReplace);
my ($fastq1, $fastq2, $referenceGenome, $outputFile, $micFn, $flankingLength, $directory);

$micFn = "regions.strincter.intervals.result";
$flankingLength = 5;

# Define the options for our function
GetOptions(
    "h|?|help"           => \$helpFlag,
    "output=s"           => \$outputFile,
    "reference=s"        => \$referenceGenome,
    "bam=s"              => \$bamFile,
    "genome-directory=s" => \$directory
) || help(1);


############################################################
# DEFINE HELP OUTPUTS
############################################################

# Define help output if no reference genome is supplied
if(!defined $referenceGenome) {
	print STDERR "\nYou forgot to specify a reference sequence!\n\n";
	help(1);
}

# Define help output if no output file is supplied
if(!defined $outputFile) {
	print STDERR "\nYou forgot to specify how you would like your output files to be named!\n\n";
	help(1);
}

# Define help output if no bam file is supplied
if(!defined $bamFile) {
	print STDERR "\nYou forgot to specify a bam file for us to work with!\n\n";
	help(1);
}

# Define help output if no genome directory is supplied
if(!defined $bamFile) {
	print STDERR "\nYou forgot to specify a directory where we can find your data!\n\n";
	help(1);
}

help(0) if defined $helpFlag;

sub help {
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR "Usage: $0 --output myOutputFile --reference myReferenceGenome --bam myBamFile\n\n";
	print STDERR "e.g., $0 --output myOutput --reference hg19.fa --bam HG00504.mapped.bam --genome-directory HG00504 \n\n";
	print STDERR "(The working directory is /home/wetherc/. The reference genome is assumed to be in ./reference/; your novel genome in ./genomes/)\n";

	exit($return);
}

############################################################
# SET ADDITIONAL VARIABLES
############################################################

system("mkdir results/$directory\n");

# Specify the present working directory
my $pwd = "/home/wetherc/";
my $EXE_DIR = "/home/wetherc/";

#Specify names for fastq outputs of SamToFastq.jar
my $fastq1 = "$outputFile.01.fastq";
my $fastq2 = "$outputFile.02.fastq";

# Specify how qsub script will be named
my $qsub = "results/$directory/qsub.$outputFile.sh";
my $out;

# Open the qsub script
open($out, ">$qsub") || die "Could not create a qsub script. Do I have proper write permissions?\n";

############################################################
# CONSTRUCT QSUB SCRIPT!
############################################################

# Spcify qsub script options
print $out "#!/bin/bash\n";
print $out "#PBS -W group_list=sfx\n";
print $out "#PBS -q sfx_q\n";
print $out "#PBS -N $outputFile\n";
print $out "#PBS -j oe\n";
print $out "#PBS -l walltime=100:00:00\n";
print $out "#PBS -l nodes=1:ppn=5\n";
print $out "#PBS -d results/$directory\n";

####################

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
print $out "module load devel/python\n";

# Convert initial bam file to two fastq files
print $out "cmd=\"java -Xmx12g -jar /apps/packages/bio/picard/1.92/bin/SamToFastq.jar VERBOSITY=DEBUG VALIDATION_STRINGENCY=LENIENT INPUT=/home/wetherc/genomes/$directory/$bamFile FASTQ=$fastq1 SECOND_END_FASTQ=$fastq2\"\n";
print $out "run_cmd\n\n";

# Stitch two fastq files into an output sam file
print $out "cmd=\"bwa mem -t 5 -P -v 3 /home/wetherc/reference/$referenceGenome $fastq1 $fastq2 > $outputFile.sam\"\n";
print $out "run_cmd\n\n";

# Clean up old files if we haven't run into errors yet...
print $out "cmd=\"rm -r /home/wetherc/genomes/$genome\"\n";
print $out "run_cmd\n\n";

# Extract/print all alignments or subalignments in SAM format
# -h : include header information
# -S : input is in SAM format
# -b : output in BAM format
# 
# Sort alignments by leftmost coordinates
print $out "cmd=\"samtools view -hS $outputFile.sam -b | samtools sort - $outputFile\"\n";
print $out "run_cmd\n\n";

# Clean up old files
print $out "cmd=\"rm ./*.fastq\"\n";
print $out "run_cmd\n\n";

# Remove PCR duplicates with samtools tools if it's not PCR free:
print $out "cmd=\"samtools rmdup $outputFile.bam $outputFile.rdup.bam\"\n";
print $out "run_cmd\n\n";

# Index sorted alignment for fast random access
# Creates .bai index file
print $out "cmd=\"samtools index $outputFile.rdup.bam\"\n";
print $out "run_cmd\n\n";

# Replace all read groups in the INPUT file with a new read group and
# assigns all reads to this read group in the OUTPUT
print $out "cmd=\"java -Xmx5g -jar /apps/packages/bio/picard/1.92/bin/AddOrReplaceReadGroups.jar I=$outputFile.rdup.bam O=$outputFile.rdup.GATK.bam SORT_ORDER=coordinate RGID=$outputFile RGLB=bwa-mem RGPL=illumina RGSM=$outputFile CREATE_INDEX=true RGPU=RGPU VALIDATION_STRINGENCY=SILENT\"\n";
print $out "run_cmd\n\n";

# Clean up old files
print $out "cmd=\"rm ./*.sam\"\n";
print $out "run_cmd\n\n";

# Extract/print all alignments or subalignments in SAM format
# -u : output to uncompressed bam
# -f : only output alignments with all bits in INT present in the FLAG field
#
# Output to sam file containing only unmapped reads from the original genome
print $out "cmd=\"samtools view -u -f 4 $outputFile.rdup.GATK.bam > $outputFile.unmapped.sam\"\n";
print $out "run_cmd\n\n";

# Clean up old files
print $out "cmd=\"rm ./*.rdup.ba*\"\n";
print $out "run_cmd\n\n";


######################################################


print $out "echo end $outputFile `date`\n";
close($out);

# Submit script to job manager
print " $qsub\n";
system("chmod gu+x $qsub");
system("qsub -m bae -M wetherc\@vbi.vt.edu $qsub");
print "The job was submitted. The final output will be '$outputFile.unmapped.Bam.file'.\n";
