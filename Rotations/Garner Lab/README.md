Pipeline Documentation
======================

General File Hierarchy
----------------------

There are five main directories used in this pipeline:

```
| - /home
|   | - wetherc
|   |   | - reference
|   |   |   | - chromosomal
|   |   | - results
|   |   | - scripts
|   |   |   | - startPipeline
|   |   |   | - detectMicsUnmapped
|   |   | - software
|   |   | - unmapped
|   |   |   | - sam
|   |   |   | - bam
|   |   |   | - velvet
|   |   |   | - sample_velvet
|   |   |   | - remapped
```

`./reference` contains the human reference genome (hg19) both as the chromosonal `.fa` files (`reference/chromosomal/`) and as an indexed reference.

`./results` is used to temporarily store the full genome sequences as they are downloaded from 1000 Genomes Project and while they are being mapped against the hg19 human reference. Following completion of this, subdirectories are removed and the files containing only unmapped reads (`HG*.unmapped.sam`) are moved to the `./unmapped` directory.

`./unmapped` contains all unmapped reads against the hg19 human genome (stored in `./unmapped/sam`). The corresponding `./unmapped/bam` folder contains these same sequences in BAM format for appeasing Velvet. `./unmapped/velvet` contains the output of both `velveth` and `velvetg`.

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

(If you want you can pretty easily edit the script to parallelize the downloads across several qsub scripts to make things go a bit faster.)

A list of all genome sample names is contained in `/home/wetherc/sampledGenomes.csv`. The full records for these genomes (as provided by 1000 Genomes Project) are contained in `/home/wetherc/analysis.sequence.index.csv`.

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
  - moved resultant file (HG*.unmapped.sam) to /home/wetherc/unmapped/sam

Step 2.5: Count unmapped reads in each genome
---------------------------------------------

To construct a count of the number of unmapped reads against the hg19 human reference for each genome parsed, we run `sh /home/wetherc/unmapped/sam/countUnmapped.sh` which executes:

```
#!/bin/bash

module load bio/samtools

for file in /home/wetherc/unmapped/sam/*.unmapped.sam
do
	printf "${file:36:7}\t" >> output.txt
	samtools view -c -fox4 $file >> output.txt
done
```

This constructs an output.txt file that on each line contains (tab separated) the genome name and the number of unmapped reads against hg19 for that genome.

Step 3: Construct de novo contigs and blast
-------------------------------------------

To construct de novo contigs from our unmapped reads and blast them against known sources, we run `perl /home/wetherc/scripts/detectMicsUnmapped/MST.pl`.

This begins by first converting all of our HG*.unmapped.sam files (from step 2) into corresponding BAM files:

```
java -Xmx12g -jar /apps/packages/bio/picard/1.92/bin/SamFormatConverter.jar INPUT=sam/$files[$i] OUTPUT=bam/$file.unmapped.bam
```

Technically, we don't *have* to convert these files from SAM to BAM, but Velvet has thrown issues when I've tried feeding it SAM files. Still haven't diagnosed why, but BAM seems to circumvent the issue.

Following this, we merge all of our BAM files into a single input for Velvet:

```
java -Xmx12g -jar /home/wetherc/software/picard/picard.jar GatherBamFiles @bams[$i*$velvetLim .. $i*$velvetLim + ($velvetLim-1)] OUTPUT=bam/merged_$i.bam
```

Strictly speaking, it's several outputs. In our `MST.pl` script, we set `$velvetLim` to an arbitrary integer. This determines how many genomes will be merged together and passed through Velvet at a time. You may have to play around with this until you get a file size that Velvet likes.

Point is, everything from here on out will be run multiple times, once for each velvet input (determined by the value of `$velvetLim` and the number of genomes you have available to use). Our script will create a new velvet subdirectory (`~/velvet/[0..n]`) for each input where all subsequent output will be stored. The `~/unmapped/sample_velvet` directory contains an example of all output that you can expect to see if the process terminates successfully.

<hr />

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

That was indecipherable. I elected to not use it. If that's an issue, it should probably be reinserted into the MST.pl script. This function is stored in `/home/wetherc/scripts/detectMicsUnmapped/process.sh`.

<hr />

The script then runs `velveth` on the BAM input (`merged.bam`). Karthik's script specified a hash length of 71. Not sure why that kmer length, but I kept it...

We then call `velvetg` for de Bruijn graph construction, error removal and repeat resolution.

Step 4: Index and Map Against Contigs
-------------------------------------------

After running Velvet on an arbitrary number of genomes (seriously, no rationale for this. You could probably figure out the inflection point where new information yielded by a single additional genome is insubstantial enough to not matter...) we concatenate all of our `contigs.fa` Velvet outputs to `/home/wetherc/unmapped/remapped/reference.fa`. We then run `bwa index -a bwtsw PATH/reference.fa` on it to index.

Following this, we can map our remaining `unmapped.HG*.sam` files against this new reference, rather than passing them through Velvet. We achieve this by:

```
java -Xmx5g -jar /home/wetherc/software/picard/picard.jar SamToFastq INPUT=/home/wetherc/unmapped/bam/HG*.unmapped.bam FASTQ=/home/wetherc/unmapped/temp/HG*.unmapped.1.fa SECOND_END_FASTQ=/home/wetherc/unmapped/temp/HG*.unmapped.2.fa UNPAIRED_FASTQ=/home/wetherc/unmapped/temp/HG04001.unmapped.singles.fa VALIDATION_STRINGENCY=SILENT

bwa mem -aY -t 5 -P -v 3 -U 25 /home/wetherc/unmapped/remapped/reference.fa /home/wetherc/unmapped/temp/HG*.unmapped.1.fa /home/wetherc/unmapped/temp/HG*.unmapped.2.fa > /home/wetherc/unmapped/remapped/HG*.sam

samtools view -hS /home/wetherc/unmapped/remapped/.sam -b | samtools sort - HG*

samtools rmdup /home/wetherc/unmapped/remapped/HG*.bam /home/wetherc/unmapped/remapped/HG*.rdup.bam

samtools index /home/wetherc/unmapped/remapped/HG*.rmdup.bam

java -Xmx5g -jar /apps/packages/bio/picard/1.92/bin/AddOrReplaceReadGroups.jar I=/home/wetherc/HG*.rmdup.bam O=HG*.rmdup.GATK.bam SORT_ORDER=coordinate RGID=HG* RGLB=bwa-mem RGPL=illumina RGSM=HG* CREATE_INDEX=true RGPU=RGPU VALIDATION_STRINGENCY=SILENT
```

Step 5: Characterize Polymorphisms
----------------------------------------------

To characterize SNPs/INDELs, we use samtool's `mpileup` on our unmapped genomes against the reference that we have compiled in the previous step. To do this, we run (contained in `/home/wetherc/scripts/SNPs.pl`):

```
/home/wetherc/software/samtools/samtools mpileup -uf /home/wetherc/unmapped/remapped/reference.fa /home/wetherc/unmapped/remapped/HG*.bam | /home/wetherc/software/bcftools/bcftools view -vcg > var.raw.bcf

/home/wetherc/software/bcftools/bcftools view var.raw.bcf | /home/wetherc/software/bcftools/vcfutils.pl varFilter -D100 > var.flt.vcf
```

**NOTE:** In this step I compile and use the current version of both samtools and bcftools (v.1.0). Shadowfax only has v.0.1.19 installed; this version contains a bug in `mpileup` regarding BCF/VCF inputs and outputs.