#!/usr/bin/perl
#PBS -W group_list=sfx
#PBS -q sfx_q
#PBS -N genomes
#PBS -r y
#PBS -j oe
#PBS -l walltime=200:00:00
#PBS -l nodes=1:ppn=1
#PBS -d /home/wetherc/

#use strict;
#use warnings "all";
use POSIX;
use File::Basename;
#use PBS::Client;
use Getopt::Long qw(:config no_ignore_case);
use Thread;
use Cwd;

@genomes = (...);

my $counter = 0;

while($counter < @genomes) {
	system("mkdir /home/wetherc/results/$genomes[$counter]\n");
	system("wget -P /home/wetherc/results/$genomes[$counter]/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/$genomes[$counter]/sequence_read/*\n");
	$counter++;
}