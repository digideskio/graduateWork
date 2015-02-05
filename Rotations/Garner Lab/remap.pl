#!/usr/bin/perl

#use strict;
#use warnings "all";
use POSIX;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
use Thread;
use Cwd;

# Declare a couple arrays used later
# on to parse through our sam/bams
# and send them through Velvet in
# batches of size N
my @files;
my @bams;

# Point this to the directory containing
# all of your SAM files
my $dir = '/home/wetherc/unmapped/sam/';

# Specify how qsub script will be named
my $qsub = "/home/wetherc/unmapped/remapped/qsub.remap.sh";
my $out;

open($out, ">$qsub") || die "Could not create a qsub script. Do I have proper write permissions?\n";

print $out "#!/bin/bash\n";
print $out "#PBS -W group_list=sfx\n";
print $out "#PBS -q sfx_q\n";
print $out "#PBS -N unmapped\n";
print $out "#PBS -r y\n";
print $out "#PBS -j oe\n";
print $out "#PBS -l walltime=100:00:00\n";
print $out "#PBS -l nodes=1:ppn=5\n";
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
print $out "module load devel/java\n";

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
	push (@bams, $file);
}

print "What velvetLim would you like to use?\n";
print "(Use the value specified by MST.pl)\n>>> ";
my $velvetLim = <>;

# Construct a new "reference" genome from
# a sampling of our Velvet outputs
#
# If you adjust the upper bound of this,
# make sure to edit $velvetLim so that
# you're not later trying to map a genome
# against itself (or, rather, a reference
# of which it's already part)
print $out "cmd=\"cat velvet/[0-5]/contigs.fa > remapped/reference.fa\"\n";
print $out "run_cmd\n\n";

print $out "cmd=\"bwa index -a bwtsw remapped/reference.fa\"\n";
print $out "run_cmd\n\n";

# Here we use $velvetLim to exclude
# the genomes that we used to construct
# our new reference
for(my $i = $velvetLim * 6; $i < @bams; $i++) {

	# We have to use a more recent version of
	# Picard for this. Otherwise it yells at
	# us about the presence of unpaired reads.
	#
	# (There aren't any, but something's
	# wonky and this at least appeases
	# it without much fussing around)
	print $out "cmd=\"java -Xmx5g -jar /home/wetherc/software/picard/picard.jar SamToFastq INPUT=bam/$bams[$i].unmapped.bam FASTQ=temp/$bams[$i].unmapped.1.fa SECOND_END_FASTQ=temp/$bams[$i].unmapped.2.fa UNPAIRED_FASTQ=temp/$bams[$i].unmapped.singles.fa VALIDATION_STRINGENCY=SILENT\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"bwa mem -aY -t 5 -P -v 3 -U 25 remapped/reference.fa temp/$bams[$i].unmapped.1.fa temp/$bams[$i].unmapped.2.fa > temp/$bams[$i].sam\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"samtools view -hS temp/$bams[$i].sam -b | samtools sort - remapped/$bams[$i]\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"samtools index remapped/$bams[$i].bam\"\n";
	print $out "run_cmd\n\n";

	print $out "cmd=\"rm temp/*.fa remapped/$bams[$i].sam\"\n";
	print $out "run_cmd\n\n";
}

print $out "echo end\n";
close($out);

# Submit script to job manager
print "$qsub\n";
system("chmod gu+x $qsub\n");
system("qsub -m bae -M wetherc\@vbi.vt.edu $qsub\n");
