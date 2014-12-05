#/usr/bin/perl
#PBS -W group_list=sfx
#PBS -q sfx_q
#PBS -N $_[0]
#PBS -j oe
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=1
#PBS -d /home/wetherc/

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

@genomes = ("HG00577","HG00578","HG00733","HG00982",
	"HG04026","HG04063","HG04093","HG00271","HG00288",
	"HG00337","HG00368","HG00376","HG00383","HG03108",
	"HG00406","HG00422","HG00437","HG00443","HG00478",
	"HG00684","HG00689","HG01412","HG02024","HG02265",
	"HG02266","HG02274","HG02275","HG02425","HG03009",
	"HG03488","HG03660","HG03663","HG01405","HG01486",
	"HG02408","HG02410","HG02879","HG03012","HG03520",
	"HG03585","HG03589","HG03593","HG03812","HG03817",
	"HG03902","HG03905","HG03910","HG03911","HG03916",
	"HG03919","HG03928","HG03931","HG03934","HG03937",
	"HG04146","HG04159","HG04180","HG04131","HG04134",
	"HG04140","HG04141","HG04144","HG04152","HG04153",
	"HG04155","HG04158","HG04161","HG04162","HG04164",
	"HG04171","HG04173","HG04176","HG04177","HG04182",
	"HG04183","HG04185","HG04186","HG04188","HG04189",
	"HG04194","HG04195","HG04198","HG04200","HG04202",
	"HG04206","HG04209","HG04210","HG04211","HG04212",
	"HG04214","HG04216","HG04219","HG04222","HG04225",
	"HG04227","HG04229","HG04235","HG04238","HG04239",
	"HG00621","HG00652","HG00654","HG00694","HG01974",
	"HG01976","HG01980","HG02008","HG02658","HG02660",
	"HG02696","HG02724","HG02733","HG02785","HG02789",
	"HG03760","HG03762","HG03765","HG03767","HG03774",
	"HG03778","HG03779","HG03780","HG03782","HG03785",
	"HG03790","HG03792","HG03802","HG03803","HG03805",
	"HG03808","HG03809","HG03821","HG03823","HG03824",
	"HG03826","HG03829","HG03832","HG03833","HG03836",
	"HG03837","HG03838","HG03844","HG03848","HG03849",
	"HG03850","HG03851","HG03854","HG03856","HG03857",
	"HG03858","HG03861","HG03863","HG03866","HG03868",
	"HG03871","HG03873","HG03875","HG03890","HG03894",
	"HG03895","HG03897","HG03898","HG03899","HG03900",
	"HG03907","HG03913","HG03917","HG03920","HG03922",
	"HG03925","HG03926","HG03940","HG03941","HG03943",
	"HG03944","HG03945","HG03950","HG03951","HG03953",
	"HG03955","HG03965","HG03967","HG03969","HG03971",
	"HG03973","HG03974","HG03976","HG03977","HG03985",
	"HG03986","HG03989","HG03990","HG03991","HG03995",
	"HG03998","HG03999","HG04001","HG04002","HG04003",
	"HG04006","HG04017","HG04018","HG04023","HG04029",
	"HG04033","HG04035","HG04038","HG04039","HG04042",
	"HG04047","HG04054","HG04056","HG04059","HG04060",
	"HG04061","HG04062","HG04070","HG04075","HG04076",
	"HG04080","HG04090","HG04094","HG04096","HG04098",
	"HG04099","HG04100","HG04106","HG04107","HG04118",
	"HG00104","HG00134","HG00135","HG00249","HG00270",
	"HG00359","HG00377","HG00418","HG00427","HG00512",
	"HG00866","HG00983","HG01274","HG01278","HG01347",
	"HG03076","HG03606","HG03633","HG03845","HG04037",
	"HG00477","HG00532","HG00535","HG00538","HG00543",
	"HG00557","HG00609","HG00612","HG00615","HG00628",
	"HG00650","HG00651","HG00653","HG00671","HG01971",
	"HG01979","HG01981","HG01991","HG01992","HG01993",
	"HG01997","HG01998","HG02003","HG02004","HG02089",
	"HG02090","HG02091","HG02104","HG02105","HG02106",
	"HG02146","HG02147","HG02148","HG02259","HG02260",
	"HG02261","HG02271","HG02272","HG02273","HG02277",
	"HG02278","HG02279","HG02285","HG02286","HG02287",
	"HG02301","HG02302","HG02303","HG02600","HG02601",
	"HG02602","HG02603","HG02604","HG02605","HG02654",
	"HG02655","HG02656","HG02657","HG02659","HG02661",
	"HG02662","HG02684","HG02685","HG02686","HG02687",
	"HG02688","HG02689","HG02697","HG02698","HG02725",
	"HG02726","HG02727","HG02728","HG02729","HG02734",
	"HG02735","HG02783","HG02784","HG02786","HG02787",
	"HG02790","HG02791","HG03237","HG03238","HG03239",
	"HG03594","HG03616","HG03625","HG03629","HG03634",
	"HG03636","HG03640","HG03645","HG03646","HG03649",
	"HG03652","HG03653","HG03667","HG03668","HG03672",
	"HG03673","HG03680","HG03681","HG03684","HG03685",
	"HG03686","HG03687","HG03689","HG03690","HG03691",
	"HG03693","HG03694","HG03695","HG03696","HG03697",
	"HG03698","HG03705","HG03706","HG03708","HG03709",
	"HG03711","HG03718","HG03733","HG03736","HG03738",
	"HG03741","HG03744","HG03745","HG03746","HG03750",
	"HG03752","HG03753","HG03754","HG03755","HG03756",
	"HG03757","HG01122","HG01362","HG03702","HG03703",
	"HG04156","HG03054","HG03055","HG03060","HG03061",
	"HG03063","HG03064","HG03072","HG03073","HG03074",
	"HG03077","HG03079","HG03081","HG03082","HG03086",
	"HG03088","HG03091","HG03095","HG03103","HG03105",
	"HG03111","HG03112","HG03117","HG03118","HG03126",
	"HG03127","HG03129","HG03130","HG03139","HG03157",
	"HG03159","HG03160","HG03162","HG03163","HG03166",
	"HG03168","HG03169","HG03172","HG03175","HG03189",
	"HG03190","HG03193","HG03195","HG03196","HG03198",
	"HG03199","HG03202","HG03209","HG03212","HG03224",
	"HG03225","HG03228","HG03229","HG03234","HG03235",
	"HG03265","HG03267","HG03268","HG03270","HG03271",
	"HG03280","HG03291","HG03295","HG03297","HG03298",
	"HG03300","HG03301","HG03303","HG03304","HG03311",
	"HG03313","HG03343","HG03351","HG03352","HG03354",
	"HG03363","HG03366","HG03367","HG03369","HG03370",
	"HG03372","HG03376","HG03378","HG03380","HG03382",
	"HG03385","HG03388","HG03391","HG03394","HG03401",
	"HG03410","HG03419","HG03428","HG03432","HG03433",
	"HG03436","HG03439","HG03442","HG03445","HG03446",
	"HG03449","HG03455","HG03457","HG03461","HG03470",
	"HG03473","HG03478","HG03479","HG03484","HG03485",
	"HG03490","HG03499","HG03511","HG03514","HG03515",
	"HG03517","HG03518","HG03521","HG03547","HG03556",
	"HG03557","HG03558","HG03559","HG03563","HG03565",
	"HG03567","HG03572","HG03575","HG03577","HG03578",
	"HG03595","HG03598","HG03600","HG03603","HG03604",
	"HG03611","HG03619","HG03624","HG03631","HG01970",
	"HG01972","HG01973","HG02291","HG02497","HG02501",
	"HG02502","HG02505","HG02508","HG02511","HG02536",
	"HG02537","HG02541","HG02549","HG02554","HG02555",
	"HG02577","HG02580","HG02597","HG02648","HG02649",
	"HG02651","HG02652","HG02681","HG02682","HG02690",
	"HG02691","HG02694","HG02699","HG02731","HG02736",
	"HG02737","HG02759","HG02760","HG02763","HG02774",
	"HG02778","HG02780","HG02792","HG02793","HG02813",
	"HG02814","HG02816","HG02817","HG02836","HG02837",
	"HG02839","HG02840","HG02851","HG02860","HG02861",
	"HG02870","HG02878","HG02881","HG02882","HG02887",
	"HG02888","HG02895","HG02896","HG02938","HG02941",
	"HG02946","HG02947","HG02952","HG02953","HG02968",
	"HG02970","HG02971","HG02976","HG02977","HG02979",
	"HG02981","HG02982","HG02983","HG03007","HG03015",
	"HG03016","HG03018","HG03019","HG03021","HG03022",
	"HG03027","HG00097","HG00099","HG00101","HG00102",
	"HG00105","HG00106","HG00107","HG00109","HG00110",
	"HG00112","HG00113","HG00115","HG00118","HG00121",
	"HG00122","HG00126","HG00128","HG00129","HG00130",
	"HG00132","HG00139","HG00141","HG00142","HG00143",
	"HG00148","HG00149","HG00150","HG00151","HG00173",
	"HG00177","HG00179","HG00180","HG00231","HG00234",
	"HG00235","HG00237","HG00238","HG00240","HG00252",
	"HG00253","HG00254","HG00255","HG00266","HG00267",
	"HG00269","HG00290","HG00304","HG00332","HG00341",
	"HG00346","HG00349","HG00350","HG00351","HG00355",
	"HG00356","HG00358","HG00360","HG00362","HG00364",
	"HG00365","HG00371","HG00378","HG00379","HG00381",
	"HG00382","HG00384","HG00407","HG00421","HG00442",
	"HG00445","HG00446","HG00448","HG00451","HG00452",
	"HG00457","HG00458","HG00559","HG00560","HG00565",
	"HG00566","HG00592","HG00593","HG00595","HG00596",
	"HG00717","HG00742","HG00743","HG00766","HG00867",
	"HG00881","HG00956","HG00978","HG01028","HG01029",
	"HG01046","HG01058","HG01063","HG01064","HG01077",
	"HG01085","HG01086","HG01088","HG01089","HG01092",
	"HG01104","HG01105","HG01119","HG01121","HG01124",
	"HG01125","HG01130","HG01131","HG01133","HG01134",
	"HG01136","HG01137","HG01139","HG01140","HG01142",
	"HG01148","HG01149","HG01161","HG01162","HG01164",
	"HG01200","HG01256","HG01257","HG01259","HG01260",
	"HG01269","HG01271","HG01272","HG01275","HG01277",
	"HG01280","HG01281","HG01284","HG01302","HG01303",
	"HG01305","HG01308","HG01311","HG01312","HG01323",
	"HG01325","HG01326","HG01344","HG01345","HG01348",
	"HG01356","HG01357","HG01363","HG01369","HG01372",
	"HG01392","HG01393","HG01395","HG01396","HG01398",
	"HG01402","HG01403","HG01413","HG01414","HG01431",
	"HG01432","HG01435","HG01443","HG01444","HG01447",
	"HG01459","HG01468","HG01474","HG01479","HG01485",
	"HG01501","HG01503","HG01504","HG01509","HG01524",
	"HG01525","HG01527","HG01528","HG01530","HG01531",
	"HG01536","HG01537","HG01556","HG01566","HG01571",
	"HG01572","HG01586","HG01589","HG01593","HG01602",
	"HG01603","HG01605","HG01606","HG01607","HG01608",
	"HG01610","HG01612","HG01613","HG01615","HG01617",
	"HG01618","HG01619","HG01620","HG01623","HG01624",
	"HG01625","HG01626","HG01628","HG01630","HG01631",
	"HG01632","HG01668","HG01669","HG01670","HG01672",
	"HG01700","HG01702","HG01765","HG01766","HG01767",
	"HG01768","HG01770","HG01771","HG01773","HG01775",
	"HG01776","HG01777","HG01779","HG01781","HG01783",
	"HG01784","HG01785","HG01786","HG01789","HG01790",
	"HG01860","HG01861","HG01882","HG01883","HG01889",
	"HG01890","HG01894","HG01912","HG01923","HG01924",
	"HG01961","HG01965","HG01988","HG01989","HG01990",
	"HG02002","HG02006","HG02009","HG02010","HG02012",
	"HG02016","HG02017","HG02019","HG02020","HG02031",
	"HG02032","HG02035","HG02050","HG02052","HG02053",
	"HG02054","HG02075","HG02076","HG02078","HG02079",
	"HG02081","HG02082","HG02084","HG02085","HG02086",
	"HG02087","HG02088","HG02102","HG02113","HG02121",
	"HG02122","HG02127","HG02128","HG02138","HG02139",
	"HG02140","HG02141","HG02142","HG02150","HG02152",
	"HG02153","HG02154","HG02155","HG02156","HG02164",
	"HG02165","HG02166","HG02178","HG02179","HG02180",
	"HG02181","HG02182","HG02184","HG02185","HG02186",
	"HG02187","HG02188","HG02190","HG02219","HG02220",
	"HG02221","HG02235","HG02236","HG02252","HG02253",
	"HG02255","HG02256","HG02262","HG02281","HG02282",
	"HG02304","HG02312","HG02317","HG02318","HG02345",
	"HG02348","HG02439","HG02442","HG02455","HG02476",
	"HG02477","HG02479","HG02481","HG02484","HG02485");

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

	close($out);

	# Submit script to job manager
	system("chmod gu+x $qsub");
	system("qsub -m bae -M wetherc\@vbi.vt.edu $qsub");
}

my $counter = 0;

while($counter < @genomes) {
	for (my $i = 0; $i < 5; $i++) {
		qsub($genomes[$counter]);
		$counter++;
	}
	sleep(129600);
	systen("rm -r /home/wetherc/results/HG*/*.o*");
}