#!/usr/bin/perl

#use strict;
#use warnings "all";
use POSIX;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
use Thread;
use Cwd;

# We'll run these through velvet in batches
# of 20 to keep from running out of memory
#
# This is an easy parameter to adjust if
# you need to do so
my $velvetLim = 10; 

# Set all of our working directories, etc.
# Not that it matters...
my $pwd = "/home/wetherc/";
my $EXE_DIR = "/home/wetherc/";

# Where are our scripts contained?
my $trf = "/home/wetherc/software/trf407b.linux64";
my $scripts = "/home/wetherc/scripts/detectMicsUnmapped";

# Declare a couple arrays used later
# on to parse through our sam/bams
# and send them through Velvet in
# batches of size N
my @files;
my @bams;
my $dir = '/home/wetherc/unmapped/sam/';


# Specify how qsub script will be named
my $qsub = "/home/wetherc/unmapped/qsub.unmapped.sh";
my $out;

# Open the qsub script
open($out, ">$qsub") || die "Could not create a qsub script. Do I have proper write permissions?\n";

print $out "#!/bin/bash\n";
print $out "#PBS -W group_list=sfx\n";
print $out "#PBS -q sfx_q\n";
print $out "#PBS -N unmapped\n";
print $out "#PBS -r y\n";
print $out "#PBS -j oe\n";
print $out "#PBS -l walltime=100:00:00\n";
print $out "#PBS -l nodes=4:ppn=12\n";
print $out "#PBS -d /home/wetherc/unmapped/\n";

# Define run_cmd function
print $out "run_cmd () {\n";
print $out "    echo; echo `date`; echo \$cmd\n";
print $out "    if !(eval \$cmd); then exit 1; fi;\n";
print $out "}\n\n";

# Load all requisite modules
print $out "module load bio/samtools\n";
print $out "module load bio/picard\n";
print $out "module load bio/bwa\n";
print $out "module load bio/velvet\n";
print $out "module load devel/java\n";
print $out "module load bio/blast\n\n";

# Open our SAM directory containing
# all of our unmapped reads
opendir(DIR, $dir) or die $!;

# While we have files
while (my $file = readdir(DIR)) {

    # We only want files; exclude
    # directories
    next unless (-f "$dir/$file");

    # Use a regular expression to find
    # files ending in .sam for obvious
    # reasons...
    next unless ($file =~ m/\.sam$/);

    # Push the file name onto the
    # @files array
    push (@files, $file); 
}

for(my $i = 0; $i < @files; $i++) {

	# Chomp the full file name for my own sanity
	# This leaves us with only the genome
	# name, sans the .unmapped.sam extension
	my $file = substr($files[$i],0,7);

	# Construct a list of input bams
	# for Picard's lovely syntax
	push (@bams, "INPUT=bam/$file.unmapped.bam");

	# Convert the file from a SAM to BAM
	# Don't redo the operation if the file
	# already exists
	if(! -f "/home/wetherc/unmapped/bam/$file.unmapped.bam") {
		print $out "cmd=\"java -Xmx12g -jar /apps/packages/bio/picard/1.92/bin/SamFormatConverter.jar INPUT=sam/$files[$i] OUTPUT=bam/$file.unmapped.bam\"\n";
		print $out "run_cmd\n\n";
	}
}

# Bin our operations into batches of
# size N, defined above via $velvetLim
for(my $i = 0; $i < @bams / $velvetLim; $i++) {

	if(! -e "/home/wetherc/unmapped/velvet/$i") {
		print $out "cmd=\"mkdir ./velvet/$i/\"\n";
		print $out "run_cmd\n\n";
	}

	# Stitch together all of our BAM files
	# In the range of [i*n, i*n+(n-1)]
	# where n is the group size specified
	# above by $velvetLim
	print $out "cmd=\"java -Xmx12g -jar /home/wetherc/software/picard/picard.jar GatherBamFiles @bams[$i*$velvetLim .. $i*$velvetLim + ($velvetLim-1)] OUTPUT=bam/merged_$i.bam\"\n";
	print $out "run_cmd\n\n";
}

for(my $i = 0; $i < @bams / $velvetLim; $i++) {
	# Run velvet
	#
	# Sacrifice a goat to your choice diety
	# and pray that there's no seg fault
	#
	# Karthik's script specified a contig length
	# of 71 for this. Not sure why---the docs
	# recommend smaller lengths. But why not?
	print $out "cmd=\"velveth ./velvet/$i/ 71 -bam ./bam/merged_$i.bam\"\n";
	print $out "run_cmd\n\n";
}

for(my $i = 0; $i < @bams / $velvetLim; $i++) {
	# More Velvet!
	print $out "cmd=\"velvetg ./velvet/$i/ -scaffolding no -read_trkg yes -amos_file yes -unused_reads yes > ./velvet/$i/velvet_71.log\"\n";
	print $out "run_cmd\n\n";
}

for(my $i = 0; $i < @bams / $velvetLim; $i++) {
	###################################
	# Caveat emptor
	###################################
	#
	# Here on out, I haven't interrogated
	# Karthik's scripts to make sure they
	# do what they should. I'm not
	# too worried about that, but just
	# can't 100% verify their accuracy
	print $out "cmd=\"$trf ./velvet/$i/contigs.fa 2 7 5 80 10 14 6 -h\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"mv ./contigs.fa.2.7.5.80.10.14.6.dat ./velvet/$i/\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"perl $scripts/makeListFromTRF.pl -r ./velvet/$i/contigs.fa -t ./velvet/$i/contigs.fa.2.7.5.80.10.14.6.dat -j 0 -o ./velvet/$i/contigs.trf.noJoin.lst\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"python $scripts/addMotifCount.py ./velvet/$i/ contigs.trf.noJoin.lst\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"mv ./velvet/$i/contigs.trf.noJoin.lst.tmp ./velvet/$i/contigs.trf.noJoin.lst -f\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"perl $scripts/report_motifs_ReadCov.pl ./velvet/$i/contigs.trf.noJoin.lst ./velvet/$i/velvet_asm.afg > ./velvet/$i/motifs.txt\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"perl $scripts/makeListFromTRF.pl -r ./velvet/$i/contigs.fa -t ./velvet/$i/contigs.fa.2.7.5.80.10.14.6.dat -j 10 -o ./velvet/$i/contigs.trf.j10.lst\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"perl $scripts/addCovReadsNum2MicInContigs.pl -c ./velvet/$i/velvet_asm.afg -m ./velvet/$i/contigs.trf.j10.lst -o ./velvet/$i/contigs.trf.j10.readCov.lst\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"perl $scripts/getMicFlankingFromContigs.pl -c ./velvet/$i/contigs.fa -m ./velvet/$i/contigs.trf.j10.readCov.lst -full -o ./velvet/$i/contigs.mic.full.fa\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"perl $scripts/getMicFlankingFromContigs.pl -c ./velvet/$i/contigs.fa -m ./velvet/$i/contigs.trf.j10.readCov.lst -full -o ./velvet/$i/contigs.flanking.fa\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"blastn -db $scripts/nt_db/nt -evalue 0.001 -outfmt 6 -query ./velvet/$i/contigs.flanking.fa -out ./velvet/$i/contigs.flanking.nt.table\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"blastn -db $scripts/human_genome_blast_db/human_genomic -evalue 0.001 -outfmt 6 -query ./velvet/$i/contigs.flanking.fa -out ./velvet/$i/contigs.flanking.hg19.table\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"perl $scripts/getUndetectedMic.pl -f ./velvet/$i/contigs.flanking.fa -m ./velvet/$i/contigs.trf.j10.readCov.lst -b ./velvet/$i/contigs.flanking.hg19.table -nt ./velvet/$i/contigs.flanking.nt.table -o ./velvet/$i/contigs.flanking.txt\"\n";
	print $out "run_cmd\n\n";
}
######################################################


print $out "echo end $outputFile\n";
close($out);

# Submit script to job manager
print " $qsub\n";
system("chmod gu+x $qsub\n");
system("qsub -m bae -M wetherc\@vbi.vt.edu $qsub\n");