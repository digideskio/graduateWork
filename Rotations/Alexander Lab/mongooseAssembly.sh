#!/bin/bash
#PBS -W group_list=hokiespeed
#PBS -q normal_q
#PBS -N mongooseAssembly
#PBS -r y
#PBS -j oe
#PBS -l walltime=100:00:00
#PBS -l nodes=2:ppn=7
#PBS -d ./

########################################
# Load all requisite modules
########################################
module load bio/abyss

########################################
# Make a project directory
########################################
mkdir mongoose

########################################
# Download FASTA reads
########################################
curl --user hokies:maroon http://mongoose.vbi.vt.edu/blast/db/2012_Illumina_read1.fa > ./mongoose/mongoose_reads_1.fa

curl --user hokies:maroon http://mongoose.vbi.vt.edu/blast/db/2012_Illumina_read2.fa > ./mongoose/mongoose_reads_2.fa

########################################
# Optimize kmer length
########################################

for k in {20..40}; do
	# make new directory for each kmer
	mkdir mongoose/k$k

	# de novo assemble paired-end reads
	# for each kmer length, parallelized
	# across 14 processors
	abyss-pe np=14 -d -C ./mongoose/k$k in='mongoose_reads_1.fa mongoose_reads_2.fa' name=mongoose
done

# Calculate assembly contiguity statistics
# for our various kmer coverages
abyss-fac ./mongoose/k*/mongoose-contigs.fa

########################################
# ftp all those juicy contigs
########################################

# Compress our contigs
# tar -czvf mongoose.tar.gz mongoose

# # host URL
# host='chriswetherill.me'

# # login information
# username='vbi'
# passwd='mongoose4life'

# # file for upload
# PUTFILE='./mongoose.tar.gz'

# # ftp that stuff
# ftp -n -v $host <<SCRIPT
# quote USER $username
# quote PASS $passwd
# mput $PUTFILE 
# quit
# SCRIPT

# rm mongoose.tar.gz
# 
# rm -r mongoose
# 
exit 0