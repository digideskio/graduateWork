Pipeline Documentation
======================

General File Hierarchy
----------------------

There are five main directories used in this pipeline:

```
	| - /home
	|   | - wetherc
	|   | - reference
	|   |   | - chromosomal
	|   | - results
	|   | - scripts
	|   |   | - startPipeline
	|   |   | - detectMicsUnmapped
	|   | - software
	|   | - unmapped
	|   |   | - concatenated
	|   |   | - velvet
	|	|	| - unmapped_velvet
```

`./reference` contains the human reference genome (hg19) both as the chromosonal `.fa` files (`reference/chromosomal/`) and as an indexed reference.

`./results` is used to temporarily store the full genome sequences as they are downloaded from 1000 Genomes Project and while they are being mapped against the hg19 human reference. Following completion of this, subdirectories are removed and the files containing only unmapped reads (`HG*.unmapped.sam`) are moved to the `./unmapped` directory.

`./unmapped` contains all unmapped reads against the hg19 human genome (stored in `./unmapped/concatenated`). (The concatenated subdirectory is stupidly named. I had things working a bit differently originally and was too lazy to change the directory name and all references to it in various scripts. So it's stupid, but it stayed...) `./unmapped/velvet` contains the output of both `velveth` and `velvetg`.

Indexing hg19
-------------

To index the reference, we ran:

```
	cd ~/wetherc/reference/
	cat chromosomal/*.fa >> ./hg19.fa

	module load bio/bwa

	bwa index -a bwtsw ./hg19.fa
```

Step 1: Download genomes
------------------------

All genomes were first downloaded using wget via the script `/home/wetherc/scripts/wget.all.qsub`. This serially downloaded all genomes that we wished to analyse from the 1000 Genomes Project. Each genome was downloaded as the original `.fastq.gz` files provided by the 1000 Genomes Project into its own directory at `/home/wetherc/results/HG*` with `HG*` being the sample name for that individual's genome.

Step 2: Mapping against hg19
----------------------------

For each genome downloaded, we ran `/home/wetherc/scripts/startPipeline/make_script.pl` using the syntax:

```
	perl /home/wetherc/scrips/startPipeline/make_script.pl --input HG*
```

where again `HG*` was the complete genome name (e.g., HG04035). This script then created and submitted a qsub job that:

  - unzipped all fastq files in the specified directory;
  - concatenated all of the forward and backward reads, respectively, into `HG*.pairs_1.fastq` and `HG*.pairs_2.fastq`
  - mapped the paired reads against hg19 by `bwa mem -aY -t 5 -P -v 3 -U 25`
  - sorted alignments of the resultant sam file by leftmost coordinates
  - removed PCR duplicates by `samtools rmdup`
  - indexed sorted alignment for fast random access
  - used Picard's AddOrReplaceReadGroups to replace all read groups in the INPUT file with a single new read group and assigns all reads to this read group in the OUTPUT BAM
  - created sam file containing only unmapped reads from the original genome via `samtools view -u -f 4`
  - moved resultant file (HG*.unmapped.sam) to /home/wetherc/unmapped/concatenated

Step 2.5: Count unmapped reads in each genome
---------------------------------------------

To construct a count of the number of unmapped reads against the hg19 human reference for each genome parsed, we run `sh /home/wetherc/unmapped/concatenated/countUnmapped.sh` which executes:

```
	#!/bin/bash

	module load bio/samtools

	for file in /home/wetherc/unmapped/concatenated/*.unmapped.sam
	do
		printf "${file:36:7}\t" >> output.txt
		samtools view -c -fox4 $file >> output.txt
	done
```

This constructs an output.txt file that on each line contains (tab separated) the genome name and the number of unmapped reads against hg19 for that genome.

Step 3: Construct de novo contigs and blast
-------------------------------------------

To construct de novo contigs from our unmapped reads and blast them against known sources, we run `perl /home/wetherc/scripts/detectMicsUnmapped/MST.pl`.

This begins by first mergins all of our HG*.unmapped.sam files (from step 2) into a single sorted sam file using:

```
	samtools merge -nr unmapped.concatenated.sam  /home/wetherc/unmapped/concatenated/HG*.unmapped.sam
```

This merged sam file must then by converted into fasta format. (Technically, this doesn't *have* to be the case, but Velvet is a memory hog and if we try to pass it a sam file, chances are very good that we will then have a seg fault. That's just bad times for everyone.)

To do this, our script runs:

```
	# Convert our sam file into two paired-read .fa files
	java -Xmx12g -jar /apps/packages/bio/picard/1.92/bin/SamToFastq.jar INPUT=unmapped.concatenated.sam FASTQ=unmapped_1.fa SECOND_END_FASTQ=unmapped_2.fa

	# Stitch those two paired-read files into a single
	# output that we can pass to Velvet
	/apps/packages/bio/velvet/current/contrib/shuffleSequences_fasta/shuffleSequences_fasta.pl unmapped_1.fa unmapped_2.fa unmapped.fa
```

- - - - - -
NOTE:

Karthik's original pipeline at this point passed a fasta file through:

```
	perl -e 'my ($name,$seq,$s,$e,$k,$l);while(<>){if(/^>/){$name = $_;}

	else{ chomp;$seq=$_;$l=length($seq);$s=0;$e=$l;$k=0;
	 for($i=0;$i<$l;$i++){ $c=substr($seq, $i, 1); if($c eq "N"){$k=0;$s=$i+1;}else{$k++; if($k==10){last;}}}
	 for($i=$l-1,$k=0;$i>$s;$i--){ $c=substr($seq, $i, 1); if($c eq "N"){$e=$i; $k=0;}else{$k++; if($k==10){last;}}}
	 if($s!=0 || $e != $l){ $seq=substr($seq,$s,$e-$s);}
	 if(length($seq)>50){ print "$name$seq\n";}
	 }
	}' $infile > $outfile
```

That was indecipherable. I elected to not use it. If that's an issue, it should probably be reinserted into the MST.pl script. This function is stored in `~/home/wetherc/scripts/detectMicsUnmapped/process.sh`.
- - - - - -

The script then runs `velveth` on the fasta input (unmapped.fa). Karthik's script specified a hash length of 71. Not sure why; according to the docs anything above 31 will generally be automatically reduced.

We then call `velvetg` for de Bruijn graph construction, error removal and repeat resolution.