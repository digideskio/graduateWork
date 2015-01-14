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

# Merge all sam files
#
# These are the remaining reads that were
# unmapped following the first step of
# out pipeline where we used bwa
# to map against hg18
print $out "cmd=\"samtools merge -nr concatenated/unmapped.concatenated.sam  concatenated/HG*.unmapped.sam\"\n";
print $out "run_cmd\n\n";

# Convert our merged sam file into paired-end
# fastq files. (Velvet is a memory hog
# and freaks out if we try to pass
# the original sam files its way)
print $out "cmd=\"java -Xmx12g -jar /apps/packages/bio/picard/1.92/bin/SamToFastq.jar INPUT=concatenated/unmapped.concatenated.sam FASTQ=concatenated/unmapped_1.fa SECOND_END_FASTQ=concatenated/unmapped_2.fa\"\n";
print $out "run_cmd\n\n";

# Knit together the paired-end fasta
# files to a single output for generally
# appeasing Velvet
print $out "cmd=\"/apps/packages/bio/velvet/current/contrib/shuffleSequences_fasta/shuffleSequences_fasta.pl concatenated/unmapped_1.fa concatenated/unmapped_2.fa unmapped.fa\"\n";
print $out "run_cmd\n\n";

# Run velvet
#
# Sacrifice a goat to your choice diety
# and pray that there's no seg fault
#
# Karthik's script specified a contig length
# of 71 for this. Not sure why---the docs
# recommend smaller lengths. But why not?
print $out "cmd=\"velveth ./velvet 71 -fastq unmapped.fa\"\n";
print $out "run_cmd\n\n";

print $out "cmd=\"velvetg ./velvet -scaffolding no -read_trkg yes -amos_file yes -unused_reads yes > ./velvet/velvet_71.log\"\n";
print $out "run_cmd\n\n";

print $out "cmd=\"$trf ./unmapped_velvet/contigs.fa 2 7 5 80 10 14 6 -h\"\n";
print $out "run_cmd\n\n";

print $out "cmd=\"perl $scripts/makeListFromTRF.pl -r ./unmapped_velvet/contigs.fa -t ./unmapped_velvet/contigs.fa.2.7.5.80.10.14.6.dat -j 0 -o contigs.trf.noJoin.lst\"\n";
print $out "run_cmd\n\n";

# print $out "cmd=\"python $scripts/addMotifCount.py ./contigs.trf.noJoin.lst\"\n";
# print $out "run_cmd\n\n";

# print $out "cmd=\"mv ./contigs.trf.noJoin.lst.tmp contigs.trf.noJoin.lst -f\"\n";
# print $out "run_cmd\n\n";

# print $out "cmd=\"perl $scripts/report_motifs_ReadCov.pl ./contigs.trf.noJoin.lst ./unmapped_velvet/velvet_asm.afg\"\n";
# print $out "run_cmd\n\n";

# print $out "cmd=\"perl $scripts/makeListFromTRF.pl -r ./unmapped_velvet/contigs.fa -t ./unmapped_velvet/contigs.fa.2.7.5.80.10.14.6.dat -j 10 -o contigs.trf.j10.lst\"\n";
# print $out "run_cmd\n\n";

# print $out "cmd=\"perl $scripts/addCovReadsNum2MicInContigs.pl -c ./unmapped_velvet/velvet_asm.afg -m contigs.trf.j10.lst -o contigs.trf.j10.readCov.lst\"\n";
# print $out "run_cmd\n\n";

# print $out "cmd=\"perl $scripts/getMicFlankingFromContigs.pl -c ./unmapped_velvet/contigs.fa -m contigs.trf.j10.readCov.lst -full -o contigs.mic.full.fa\"\n";
# print $out "run_cmd\n\n";

# print $out "cmd=\"blastn -db $scripts/nt_db/nt -evalue 0.001 -outfmt 6 -query contigs.flanking.fa -out contigs.flanking.nt.table\"\n";
# print $out "run_cmd\n\n";

# print $out "cmd=\"blasatn -db $scripts/human_genome_blast_db/human_genomic_evalue 0.001 -outfmt 6 -query contigs.flanking.fa -out fontigs.flanking.hg19.table\"\n";
# print $out "run_cmd\n\n";

# print $out "cmd=\"perl $scripts/getUndetectedMic.pl -f contigs.flanking.fa -m contigs.trf.j10.readCov.lst -b contigs.flanking.hg19.table -nt contigs.flanking.nt.table -o contigs.flanking.txt\"\n";
# print $out "run_cmd\n\n";

######################################################


print $out "echo end $outputFile\n";
close($out);

# Submit script to job manager
print " $qsub\n";
system("chmod gu+x $qsub\n");
system("qsub -m bae -M wetherc\@vbi.vt.edu $qsub\n");