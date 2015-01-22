#!/usr/bin/perl

#use strict;
#use warnings "all";
use POSIX;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
use Thread;
use Cwd;

# GetOptions("input=s" => \$input)


my $pwd = "/home/wetherc/";
my $EXE_DIR = "/home/wetherc/";
my $trf = "/home/wetherc/software/trf407b.linux64";
my $scripts = "/home/wetherc/scripts/detectMicsUnmapped";

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
print $out "#PBS -l nodes=1:ppn=10\n";
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
print $out "module load bio/blast\n";


my @files;
my @bams;
my $dir = '/home/wetherc/unmapped/sam/';

opendir(DIR, $dir) or die $!;

while (my $file = readdir(DIR)) {

    # We only want files
    next unless (-f "$dir/$file");

    # Use a regular expression to find files ending in .txt
    next unless ($file =~ m/\.sam$/);

    push @files, $file; 
}

for(my $i = 0; $i < @files; $i++) {

	# Chome the full file name for my own sanity
	my $file = substr($files[$i],0,7);

	# Construct a list of input bams
	# for Picard's lovely syntax
	push @bams, "INPUT=concatenated/$files[$i]";

	# Convert the file from a SAM to BAM
	print $out "cmd=\"java -Xmx12g -jar /apps/packages/bio/picard/1.92/bin/SamFormatConverter.jar INPUT=sam/$files[$i] OUTPUT=bam/$file.unmapped.bam\"\n";
	print $out "run_cmd\n\n";

	# Convert our merged sam file into paired-end
	# fastq files. (Velvet is a memory hog
	# and freaks out if we try to pass
	# the original sam files its way)
	#print $out "cmd=\"java -Xmx12g -jar /apps/packages/bio/picard/1.92/bin/SamToFastq.jar INCLUDE_NON_PF_READS=True VALIDATION_STRINGENCY=SILENT INPUT=concatenated/$files[$i] FASTQ=concatenated/$file\_1.fa SECOND_END_FASTQ=concatenated/$file\_2.fa\"\n";
	#print $out "run_cmd\n\n";

	# Knit together the paired-end fasta
	# files to a single output for generally
	# appeasing Velvet
	# print $out "cmd=\"/apps/packages/bio/velvet/current/contrib/shuffleSequences_fasta/shuffleSequences_fasta.pl concatenated/$file\_1.fa concatenated/$file\_2.fa concatenated/$file_unmapped.fa\"\n";
	# print $out "run_cmd\n\n";

	# print $out "cmd=\"rm concatenated/$file\_[1,2].fa\"\n";
	# print $out "run_cmd\n\n";

	# print $out "cmd=\"mkdir velvet/$file/\"\n";
	# print $out "run_cmd\n\n";
}

# Stitch together all of our BAM files
print $out "cmd=\"java -Xmx12g -jar /home/wetherc/software/picard/picard.jar GatherBamFiles @bams OUTPUT=bam/merged.bam\"\n";
print $out "run_cmd\n\n";

# Run velvet
#
# Sacrifice a goat to your choice diety
# and pray that there's no seg fault
#
# Karthik's script specified a contig length
# of 71 for this. Not sure why---the docs
# recommend smaller lengths. But why not?
print $out "cmd=\"velveth ./velvet/ 71 -bam ./bam/merged.bam\"\n";
print $out "run_cmd\n\n";

# More Velvet!
print $out "cmd=\"velvetg ./velvet/ -scaffolding no -read_trkg yes -amos_file yes -unused_reads yes > ./velvet/velvet_71.log\"\n";
print $out "run_cmd\n\n";


print $out "cmd=\"$trf ./velvet/contigs.fa 2 7 5 80 10 14 6 -h\"\n";
print $out "run_cmd\n\n";

print $out "cmd=\"mv ./contigs.fa.2.7.5.80.10.14.6.dat ./velvet/\"\n";
print $out "run_cmd\n\n";

print $out "cmd=\"perl $scripts/makeListFromTRF.pl -r ./velvet/contigs.fa -t ./velvet/contigs.fa.2.7.5.80.10.14.6.dat -j 0 -o ./velvet/contigs.trf.noJoin.lst\"\n";
print $out "run_cmd\n\n";

print $out "cmd=\"python $scripts/addMotifCount.py ./velvet/ contigs.trf.noJoin.lst\"\n";
print $out "run_cmd\n\n";

print $out "mv ./velvet/contigs.trf.noJoin.lst.tmp ./velvet/contigs.trf.noJoin.lst -f\n";
print $out "run_cmd\n\n";

print $out "perl $scripts/report_motifs_ReadCov.pl ./velvet/contigs.trf.noJoin.lst ./velvet/velvet_asm.afg > ./velvet/motifs.txt\n";
print $out "run_cmd\n\n";

print $out "perl $scripts/makeListFromTRF.pl -r ./velvet/contigs.fa -t ./velvet/contigs.fa.2.7.5.80.10.14.6.dat -j 10 -o ./velvet/contigs.trf.j10.lst\n";
print $out "run_cmd\n\n";

print $out "perl $scripts/addCovReadsNum2MicInContigs.pl -c ./velvet/velvet_asm.afg -m ./velvet/contigs.trf.j10.lst -o ./velvet/contigs.trf.j10.readCov.lst\n";
print $out "run_cmd\n\n";

print $out "perl $scripts/getMicFlankingFromContigs.pl -c ./velvet/contigs.fa -m ./velvet/contigs.trf.j10.readCov.lst -full -o ./velvet/contigs.mic.full.fa\n";
print $out "run_cmd\n\n";

print $out "perl $scripts/getMicFlankingFromContigs.pl -c ./velvet/contigs.fa -m ./velvet/contigs.trf.j10.readCov.lst -full -o ./velvet/contigs.flanking.fa\n";
print $out "run_cmd\n\n";

print $out "blastn -db $scripts/nt_db/nt -evalue 0.001 -outfmt 6 -query ./velvet/contigs.flanking.fa -out ./velvet/contigs.flanking.nt.table\n";
print $out "run_cmd\n\n";

print $out "blastn -db $scripts/human_genome_blast_db/human_genomic -evalue 0.001 -outfmt 6 -query ./velvet/contigs.flanking.fa -out ./velvet/contigs.flanking.hg19.table\n";
print $out "run_cmd\n\n";

print $out "perl $scripts/getUndetectedMic.pl -f ./velvet/contigs.flanking.fa -m ./velvet/contigs.trf.j10.readCov.lst -b ./velvet/contigs.flanking.hg19.table -nt ./velvet/contigs.flanking.nt.table -o ./velvet/contigs.flanking.txt\n";
print $out "run_cmd\n\n";

######################################################


print $out "echo end $outputFile\n";
close($out);

# Submit script to job manager
print " $qsub\n";
system("chmod gu+x $qsub\n");
system("qsub -m bae -M wetherc\@vbi.vt.edu $qsub\n");