#/usr/bin/perl

#use strict;
#use warnings "all";
use POSIX;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
use Thread;
use Cwd;

my ($cmd, $micFn, $flankingLength);

$micFn = "regions.strincter.intervals.result";
$flankingLength = 5;

# Specify the present working directory
my $pwd = "/home/wetherc/";
my $EXE_DIR = "/home/wetherc/";

@genomes = (
	"HG03673","HG03680","HG03681","HG03684","HG03685",
	"HG03686","HG03687","HG03689","HG03690","HG03691",
	"HG03693","HG03694","HG03695","HG03696","HG03697");

sub qsub{
	system("mkdir /home/wetherc/results/$_[0]\n");

	# Specify how qsub script will be named
	my $qsub = "/home/wetherc/results/$_[0]/qsub.$_[0].sh";
	my $out;

	# Open the qsub script
	open($out, ">$qsub") || die "Could not create a qsub script. Do I have proper write permissions?\n";

	print $out "#!/bin/bash\n";
	print $out "#PBS -W group_list=sfx\n";
	print $out "#PBS -q sfx_q\n";
	print $out "#PBS -N $_[0]\n";
	print $out "#PBS -j oe\n";
	print $out "#PBS -l walltime=100:00:00\n";
	print $out "#PBS -l nodes=1:ppn=5\n";
	print $out "#PBS -d /home/wetherc/results/$_[0]\n";

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

	print $out "cmd=\"wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/$_[0]/sequence_read/*\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"gunzip *.gz\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"cat *_1.filt.fastq > $_[0]_pairs_1.fastq\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"cat *_2.filt.fastq > $_[0]_pairs_2.fastq\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"rm *_1.filt.fastq\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"rm *_2.filt.fastq\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"cat *.filt.fastq > $_[0]_singles.fastq\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"rm *.filt.fastq\"\n";
	print $out "run_cmd\n\n";

	# Stitch two fastq files into an output sam file
	print $out "cmd=\"bwa mem -aY -t 5 -P -v 3 -U 25 /home/wetherc/reference/hg19.fa /home/wetherc/results/$_[0]/$_[0]_pairs_1.fastq /home/wetherc/results/$_[0]/$_[0]_pairs_2.fastq  > $_[0].sam\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"rm *.fastq\"\n";
	print $out "run_cmd\n\n";

	# Extract/print all alignments or subalignments in SAM format
	# -h : include header information
	# -S : input is in SAM format
	# -b : output in BAM format
	# 
	# Sort alignments by leftmost coordinates
	print $out "cmd=\"samtools view -hS $_[0].sam -b | samtools sort - $_[0]\"\n";
	print $out "run_cmd\n\n";

	# Remove PCR duplicates with samtools tools if it's not PCR free:
	print $out "cmd=\"samtools rmdup $_[0].bam $_[0].rdup.bam\"\n";
	print $out "run_cmd\n\n";

	# Index sorted alignment for fast random access
	# Creates .bai index file
	print $out "cmd=\"samtools index $_[0].rdup.bam\"\n";
	print $out "run_cmd\n\n";

	# Replace all read groups in the INPUT file with a new read group and
	# assigns all reads to this read group in the OUTPUT
	print $out "cmd=\"java -Xmx5g -jar /apps/packages/bio/picard/1.92/bin/AddOrReplaceReadGroups.jar I=$_[0].rdup.bam O=$_[0].rdup.GATK.bam SORT_ORDER=coordinate RGID=$_[0] RGLB=bwa-mem RGPL=illumina RGSM=$_[0] CREATE_INDEX=true RGPU=RGPU VALIDATION_STRINGENCY=SILENT\"\n";
	print $out "run_cmd\n\n";

	# Clean up old files
	print $out "cmd=\"rm *.sam\"\n";
	print $out "run_cmd\n\n";

	# Extract/print all alignments or subalignments in SAM format
	# -u : output to uncompressed bam
	# -f : only output alignments with all bits in INT present in the FLAG field
	#
	# Output to sam file containing only unmapped reads from the original genome
	print $out "cmd=\"samtools view -u -f 4 $_[0].rdup.GATK.bam > $_[0].unmapped.sam\"\n";
	print $out "run_cmd\n\n";

	# Clean up old files
	print $out "cmd=\"rm *.ba*\"\n";
	print $out "run_cmd\n\n";


	######################################################


	print $out "echo end $outputFile\n";
	close($out);

	# Submit script to job manager
	print " $qsub\n";
	system("chmod gu+x $qsub");
	system("qsub -m bae -M wetherc\@vbi.vt.edu $qsub");
	print "The job was submitted. The final output will be '$_[0].unmapped.Bam.file'.\n";
}

my $counter = 0;
sleep(7200);

while($counter < @genomes) {
	for (my $i = 0; $i < 5; $i++) {
		qsub($genomes[$counter]);
		$counter++;
	}
	sleep(86400);
}