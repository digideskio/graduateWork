#!/usr/bin/perl

#use strict;
#use warnings "all";
use POSIX;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
use Thread;
use Cwd;

my $input;

GetOptions(
    "input=s" => \$input
)

my $pwd = "/home/wetherc/";
my $EXE_DIR = "/home/wetherc/";
my $trf = "/home/wetherc/software/trf407b.linux64";
my $reference = "/home/wetherc/reference/hg19.fa";
my $scripts = "/home/wetherc/scripts/detectMicsUnmapped";
my $fq1 = "unmapped_1.fastq";
my $fq2 = "unmapped_2.fastq";

# Specify how qsub script will be named
my $qsub = "/home/wetherc/unmappedCleaned/qsub.$input.sh";
my $outdir = "/home/wetherc/unmappedCleaned";
my $out;

# Open the qsub script
open($out, ">$qsub") || die "Could not create a qsub script. Do I have proper write permissions?\n";

print $out "#!/bin/bash\n";
print $out "#PBS -W group_list=sfx\n";
print $out "#PBS -q sfx_q\n";
print $out "#PBS -N $input\n";
print $out "#PBS -j oe\n";
print $out "#PBS -l walltime=100:00:00\n";
print $out "#PBS -l nodes=1:ppn=10\n";
print $out "#PBS -d /home/wetherc/unmappedCleaned/\n";

# Load all requisite modules
print $out "module load bio/samtools\n";
print $out "module load bio/picard\n";
print $out "module load bio/bwa\n";
print $out "module load bio/velvet\n";
print $out "module load devel/java\n";

# Convert our unmapped sam file into forward and reverse read fastq
print $out "java -Xmx12g -jar /apps/packages/bio/picard/1.92/bin/SamToFastq.jar VERBOSITY=DEBUG VALIDATION_STRINGENCY=LENIENT INPUT=$input FASTQ=$fq1 SECOND_END_FASTQ=$fq2\n";

print $out "bwa mem -t 10 $reference $f1 $f2 > $outdir/unmappedCleaned.sam\n";

print $out "samtools view -hS $outdir/unmappedCleaned.sam -b | samtools sort - unmappedCleaned\n";

print $out "samtools index $outdir/unmappedCleaned.bam $outdir/unmappedCleaned.bai\n";

print $out "java -Xmx12g -jar /apps/packages/bio/picard/1.92/bin/AddOrReplaceReadGroups.jar I=unmappedCleaned.bam O=unmappedCleaned.rg.bam SORT_ORDER=coordinate RGID=unmappedCleaned RGLB=bwa-mem RGPL=illumina RGSM=unmappedCleaned CREATE_INDEX=true RGPU=RGPU VALIDATION_STRINGENCY=SILENT\n";

print $out "java -Xmx12g -jar $scripts/GenomeAnalysisTK.jar -allowPotentiallyMisencodedQuals -T IndelRealigner -R $reference -targetIntervals $scripts/mic.intervals -I $outdir/unmappedCleaned.rg.bam -o unmappedCleaned.GATK.bam\n";

print $out "perl $scripts/filterUnmappedFromSAM.pl unmappedCleaned.GATK.bam -fa -o unmappedCleaned.fa\n";

print $out "sh $scripts/process.sh $outdir/unmappedCleaned.fa $outdir/unmappedCleaned_N_filtered.fa\n";

print $out "mkdir $outdir/unmapped_velvet/\n";

print $out "sh $scripts/run_velvet_se.sh 71 $outdir/unmapped_velvet/ $outdir/unmappedCleaned_N_filtered.fa\n";

print $out "$trf $outdir/unmapped_velvet/contigs.fa 2 7 5 80 10 14 6 -h\n";

print $out "perl $scripts/makeListFromTRF.pl -r $outdir/unmapped_velvet/contigs.fa -t $outdir/unmapped_velvet/contigs.fa.2.7.5.80.10.14.6.dat -j 0 -o contigs.trf.noJoin.lst\n";

print $out "python $scripts/addMotifCount.py $outdir/contigs.trf.noJoin.lst\n";

print $out "mv $outdir/contigs.trf.noJoin.lst.tmp contigs.trf.noJoin.lst -f\n";

print $out "perl $scripts/report_motifs_ReadCov.pl $outdir/contigs.trf.noJoin.lst $outdir/unmapped_velvet/velvet_asm.afg\n";

print $out "perl $scripts/makeListFromTRF.pl -r $outdir/unmapped_velvet/contigs.fa -t $outdir/unmapped_velvet/contigs.fa.2.7.5.80.10.14.6.dat -j 10 -o contigs.trf.j10.lst\n";

print $out "perl $scripts/addCovReadsNum2MicInContigs.pl -c $outdir/unmapped_velvet/velvet_asm.afg -m contigs.trf.j10.lst -o contigs.trf.j10.readCov.lst\n";

print $out "perl $scripts/getMicFlankingFromContigs.pl -c $outdir/unmapped_velvet/contigs.fa -m contigs.trf.j10.readCov.lst -full -o contigs.mic.full.fa\n";

print $out "blastn -db $scripts/nt_db/nt -evalue 0.001 -outfmt 6 -query contigs.flanking.fa -out contigs.flanking.nt.table\n";

print $out "blasatn -db $scripts/human_genome_blast_db/human_genomic_evalue 0.001 -outfmt 6 -query contigs.flanking.fa -out fontigs.flanking.hg19.table\n";

print $out "perl $scripts/getUndetectedMic.pl -f contigs.flanking.fa -m contigs.trf.j10.readCov.lst -b contigs.flanking.hg19.table -nt contigs.flanking.nt.table -o contigs.flanking.txt\n";

######################################################


print $out "echo end $outputFile\n";
close($out);

# Submit script to job manager
print " $qsub\n";
system("chmod gu+x $qsub\n");
system("qsub -m bae -M wetherc\@vbi.vt.edu $qsub\n");