#!/usr/bin/perl

#use strict;
#use warnings "all";
use POSIX;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
use Thread;
use Cwd;

# Specify how qsub script will be named
my $qsub = "/home/wetherc/unmapped/remapped/qsub.snp.sh";
my $out;

open($out, ">$qsub") || die "Could not create a qsub script. Do I have proper write permissions?\n";

print $out "#!/bin/bash\n";
print $out "#PBS -W group_list=sfx\n";
print $out "#PBS -q sfx_q\n";
print $out "#PBS -N snps\n";
print $out "#PBS -r y\n";
print $out "#PBS -j oe\n";
print $out "#PBS -l walltime=100:00:00\n";
print $out "#PBS -l nodes=1:ppn=5\n";
print $out "#PBS -d /home/wetherc/unmapped/remapped\n";

# Define run_cmd function
print $out "run_cmd () {\n";
print $out "    echo; echo `date`; echo \$cmd\n";
print $out "    if !(eval \$cmd); then exit 1; fi;\n";
print $out "}\n\n";

# Load all requisite modules
print $out "module load bio/samtools\n";

# Generate BCF or pileup for one or multiple BAM files
print $out "cmd=\"samtools mpileup -uf /home/wetherc/unmapped/remapped/reference.fa /home/wetherc/unmapped/remapped/HG*.bam | bcftools view -vcg > var.raw.bcf\"\n";
print $out "run_cmd\n\n";

# call variant candidates and estimate allele frequencies
print $out "cmd=\"bcftools view var.raw.bcf | /apps/packages/bio/samtools/current/bcfutils/vcfutils.pl varFilter -D100 > var.flt.vcf\"\n";
print $out "run_cmd\n\n";

# print $out "cmd=\"samtools mpileup -uf /home/wetherc/unmapped/remapped/reference.fa /home/wetherc/unmapped/remapped/HG*.bam | java -jar /home/wetherc/software/VarScan.jar copynumber output.basename --mpileup 1\"\n";
# print $out "run_cmd\n\n";

print $out "echo end\n";
close($out);

# Submit script to job manager
print "$qsub\n";
system("chmod gu+x $qsub\n");
system("qsub -m bae -M wetherc\@vbi.vt.edu $qsub\n");