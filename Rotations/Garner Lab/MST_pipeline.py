#09/21/2014
#Written by Karthik


import sys
import getopt
import os


def main(argv):                          
  
	try:                                
		opts, args = getopt.getopt(argv, "h", ["help","fq1=","fq2=","ref=","out-header=","scriptsDir=","java-path="])
	except getopt.GetoptError:                                 
		sys.exit('Program terminated\n') 
			
	for opt, arg in opts:  
  
		if opt in ("-h", "--help"):

			print '\n********************* PROGRAM USAGE *********************\n'
			print 'python MST_pipeline.py --fq1 --fq2 --ref --out-header --scriptsDir --java-path\n'
			print 'Perl programs directory: /home/kvel/natalie/detect_micS_unmapped/\n'
			print 'Please do not forget to use the slash symbol at the end of the perl-directory\n'  
			print 'Also, do not forget to give the absolute path for all input files\n' 
			print 'The scripts directory path must end with a "/" forward slash\n'                
			sys.exit()    
		  
		elif opt == '--fq1':		              
			global f1
			f1 = arg  

		elif opt == '--fq2':                
			global f2
			f2 = arg

		elif opt == '--out-header':                
			global header 
			global pbsfilename
			header = arg	
			pbsfilename = arg


		elif opt == '--ref':
			global reference
			reference = arg

		elif opt == '--scriptsDir':
			global scriptsDir
			scriptsDir = arg
			
		elif opt == '--java-path':
			global javaPath
			javaPath = arg		
	
	print "\nCommand line arguments:\n"
	print 'Path to fastq_read1 file: '+str(f1)
	print 'Path to fastq_read2 file: '+str(f2)
	print 'Head: '+str(header)
	print 'Reference fastq file: '+str(reference)
	print 'Script directory: '+str(scriptsDir)
	print 'Path to Java executable: '+str(javaPath)+'\n'
	 
def usage():
	print 'Error in input\n'


if __name__ == "__main__":
	main(sys.argv[1:])
    

cwd = os.getcwd()

pbsfile = open(cwd+'/'+pbsfilename+'_detectMicS.pbs','w')

pbsfile.write('#!/bin/bash\n')
pbsfile.write('#PBS -l nodes=1:ppn=10\n#PBS -j oe\n')
pbsfile.write('#PBS -o $PBS_JOBID.output\n')
pbsfile.write('#PBS -l walltime=80:00:00\n')
pbsfile.write('#PBS -q sfx_q\n#PBS -W group_list=sfx\n\n\n')

pbsfile.write('module load bio/velvet\n')
pbsfile.write('module load bio/samtools\n')
pbsfile.write('module load bio/bwa\n\n')
pbsfile.write('module load bio/blast\n\n')

#trf location needs to be changed according to where the user's trf program location

trf = '/home/kvel/softwares/trf407b.linux64'

infile_array = f1.split('/')
del infile_array[0]
del infile_array[-1]
outFolder = '/'.join(infile_array)
outFolder = '/'+outFolder+'/'

outpath = outFolder+header

pbsfile.write('mkdir '+outpath+'\n\n')

pbsfile.write('cd '+outpath+'/\n\n')

pbsfile.write('bwa mem -t 10 '+reference+' '+f1+' '+f2+' > '+outpath+'/'+header+'.sam\n\n')

pbsfile.write('samtools view -hS '+outpath+'/'+header+'.sam -b | samtools sort - '+header+'\n\n')

pbsfile.write('samtools index '+outpath+'/'+header+'.bam '+outpath+'/'+header+'.bai\n\n')

pbsfile.write(javaPath+' -Xmx12g -jar /apps/packages/bio/picard/1.92/bin/AddOrReplaceReadGroups.jar I='+header+'.bam O='+header+'.rg.bam SORT_ORDER=coordinate RGID='+header+' RGLB=bwa-mem RGPL=illumina RGSM='+header+' CREATE_INDEX=true RGPU=RGPU VALIDATION_STRINGENCY=SILENT\n\n')

pbsfile.write(javaPath+' -Xmx12g -jar '+scriptsDir+'GenomeAnalysisTK.jar -allowPotentiallyMisencodedQuals -T IndelRealigner -R '+reference+' -targetIntervals '+scriptsDir+'mic.intervals -I '+outpath+'/'+header+'.rg.bam -o '+header+'.GATK.bam\n\n')


pbsfile.write('perl '+scriptsDir+'filterUnmappedFromSAM.pl '+header+'.GATK.bam -fa -o '+header+'_unmapped.fa\n\n')

pbsfile.write('cd '+scriptsDir+'\n\n')

pbsfile.write('./process.sh '+ outpath+'/'+header+'_unmapped.fa '+ outpath+'/'+header+'_unmapped_N_filtered.fa\n\n')

pbsfile.write('mkdir '+outpath+'/'+header+'_velvet/\n\n')

pbsfile.write('./run_velvet_se.sh 71 '+outpath+'/'+header+'_velvet/ '+outpath+'/'+header+'_unmapped_N_filtered.fa\n\n')

pbsfile.write('./trf407b.linux64 '+outpath+'/'+header+'_velvet/contigs.fa 2 7 5 80 10 14 6 -h\n\n')

pbsfile.write('mv contigs.fa.2.7.5.80.10.14.6.dat '+outpath+'/contigs.fa.2.7.5.80.10.14.6.dat\n\n')

pbsfile.write('cd '+outpath+'/\n\n')

pbsfile.write('perl '+scriptsDir+'makeListFromTRF.pl -r '+outpath+'/'+header+'_velvet/contigs.fa -t contigs.fa.2.7.5.80.10.14.6.dat -j 0 -o contigs.trf.noJoin.lst\n\n')


#ADDING A SIXTH COLUMN TO THE contigs.trf.noJoin.lst FILE i.e. NUMBER OF TIMES THE micS MOTIF IS REPEATED
#pbsfile.write("perl -e 'while(<>){$_ =~ /cov_(\d+)/;chomp;print '$_\t$1\n';}'"+' contigs.trf.noJoin.lst > contigs.trf.noJoin.lst.tmp\n\n') --> HONGSOEK'S CODE

pbsfile.write('python '+scriptsDir+'addMotifCount.py '+outpath+'/ contigs.trf.noJoin.lst\n\n')
	
pbsfile.write('mv contigs.trf.noJoin.lst.tmp contigs.trf.noJoin.lst -f\n\n')

pbsfile.write('perl '+scriptsDir+'report_motifs_ReadCov.pl contigs.trf.noJoin.lst '+header+'_velvet/velvet_asm.afg\n\n')

pbsfile.write('perl '+scriptsDir+'makeListFromTRF.pl -r '+header+'_velvet/contigs.fa -t contigs.fa.2.7.5.80.10.14.6.dat -j 10 -o contigs.trf.j10.lst\n\n')

pbsfile.write('perl '+scriptsDir+'addCovReadsNum2MicInContigs.pl -c '+header+'_velvet/velvet_asm.afg -m contigs.trf.j10.lst -o contigs.trf.j10.readCov.lst\n\n')

pbsfile.write('perl '+scriptsDir+'getMicFlankingFromContigs.pl -c '+header+'_velvet/contigs.fa -m contigs.trf.j10.readCov.lst  -full -o contigs.mic.full.fa\n\n')

pbsfile.write('perl '+scriptsDir+'getMicFlankingFromContigs.pl -c '+header+'_velvet/contigs.fa -m contigs.trf.j10.readCov.lst  -full -o contigs.flanking.fa\n\n')

pbsfile.write('blastn -db '+scriptsDir+'nt_db/nt -evalue 0.001 -outfmt 6 -query contigs.flanking.fa -out contigs.flanking.nt.table\n\n')

pbsfile.write('blastn -db '+scriptsDir+'human_genome_blast_db/human_genomic -evalue 0.001 -outfmt 6 -query contigs.flanking.fa -out contigs.flanking.hg19.table\n\n')

pbsfile.write('perl '+scriptsDir+'getUndetectedMic.pl -f contigs.flanking.fa -m contigs.trf.j10.readCov.lst -b contigs.flanking.hg19.table -nt contigs.flanking.nt.table  -o contigs.flanking.txt\n\n')


pbsfile.close()

print 'Find your pbs file at:\t'+str(cwd)+'\n'


