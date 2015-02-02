#!/usr/bin/perl

#use strict;
#use warnings "all";
use POSIX;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
use Thread;
use Cwd;

# Populate with your choice
# array of genome names (HG...)
@genomes = ("HG01048","HG01095","HG01286","HG02700","HG02775",
	"HG02819","HG02820","HG02851","HG02852","HG02854",
	"HG02855","HG02884","HG02885","HG02890","HG02891",
	"HG03024","HG03025","HG03066","HG03069","HG03240",
	"HG03241","HG03246","HG03247","HG03279","HG03294",
	"HG03342","HG03437","HG03458","HG03491","HG03538",
	"HG03539","HG03571","HG03583","HG03615","HG03643",
	"HG03644","HG03652","HG03653","HG03667","HG03668",
	"HG03672","HG03673","HG03679","HG03680","HG03681",
	"HG03684","HG03686","HG03687","HG03689","HG03690",
	"HG03692","HG03693","HG03694","HG03695","HG03696",
	"HG03697","HG03713","HG03714","HG03715","HG03716",
	"HG03717","HG03720","HG03722","HG03727","HG03729",
	"HG03730","HG03731","HG03743","HG03770","HG03771",
	"HG03772","HG03773","HG03775","HG03777","HG03781",
	"HG03784","HG03786","HG03787","HG03955","HG03965",
	"HG03967","HG03969","HG03971","HG03973","HG03974",
	"HG03976","HG03977","HG03985","HG03986","HG03989",
	"HG03990","HG03991","HG03995","HG03998","HG03999",
	"HG04001","HG04002","HG04003","HG04006","HG04017",
	"HG04018","HG04023","HG04029","HG04033","HG04035",
	"HG04038","HG04039","HG04042","HG04047","HG04056",
	"HG03740","HG02814","HG02839","HG02870","HG02888",
	"HG02946","HG02970","HG02981","HG03016","HG03027",
	"HG00105","HG00112","HG00122","HG00132","HG00148",
	"HG00177","HG00235","HG00253","HG00269","HG00346",
	"HG00356","HG00365","HG00382","HG00445","HG00457",
	"HG00566","HG00717","HG00881","HG01046","HG01085",
	"HG01104","HG01125","HG01136","HG01148","HG01200",
	"HG01269","HG01280","HG01305","HG02816","HG02840");

my $counter = 0;
my $lim = @genomes / 5;

for(my $i = 0; $i < $lim; $i++) {
	# Specify how qsub script will be named
	my $qsub = "/home/wetherc/qsub.wget.$i.sh";
	my $out;

	# Open the qsub script
	open($out, ">$qsub") || die "Could not create a qsub script. Do I have proper write permissions?\n";

	print $out "#!/bin/bash\n";
	print $out "#PBS -W group_list=sfx\n";
	print $out "#PBS -q sfx_q\n";
	print $out "#PBS -N $i\n";
	print $out "#PBS -r y\n";
	print $out "#PBS -j oe\n";
	print $out "#PBS -l walltime=100:00:00\n";
	print $out "#PBS -l nodes=1:ppn=1\n";
	print $out "#PBS -d /home/wetherc/results/\n";

	# Define run_cmd function
	print $out "run_cmd () {\n";
	print $out "    echo; echo `date`; echo \$cmd\n";
	print $out "    if !(eval \$cmd); then exit 1; fi;\n";
	print $out "}\n\n";

	while($counter < $lim + $i * $lim) {
		print $out "cmd=\"mkdir /home/wetherc/results/$genomes[$counter]\"\n";
		print $out "run_cmd\n\n";

		print $out "cmd=\"wget -P /home/wetherc/results/$genomes[$counter]/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/$genomes[$counter]/sequence_read/*\n";
		print $out "run_cmd\n\n";

		$counter++;
	}

	close($out);

	# Submit script to job manager
	print "$qsub\n";
	system("chmod gu+x $qsub\n");
	system("qsub -m bae -M wetherc\@vbi.vt.edu /home/wetherc/results/qsub.$i.sh\n");
	print "The job was submitted.\n";
}