__author__ = 'Jamie Wangen'

"""
Workflow for generating figures for SRI37240 paper: *** DOI PENDING ***

"""



import sys
import os
import importlib
import argparse
import subprocess
from ftplib import FTP

rootDir = os.path.dirname(os.path.realpath(__file__)) ## add current directory as rootDir
threadNumb = 40



##### 1) Generate Genome Files:
	
	### download hg38 from genecode:

class generateGenomes(object):

	### building all files need for downstream analysis

	def __init__(self, rootDir, threadNumb):
		self.rootDir = rootDir
		self.threadNumb = threadNumb

	def download_hg38(self):
		"""

		"""
		outfile = "%s/genomes/GRCh38.primary_assembly.genome.fa.gz" % self.rootDir
		hg38file = "GRCh38.primary_assembly.genome.fa.gz"

		if not os.path.exists("%s/genomes" % rootDir):	os.makedirs("%s/genomes" % self.rootDir)

		ftp = FTP('ftp.ebi.ac.uk')
		ftp.login()
		ftp.cwd('pub/databases/gencode/Gencode_human/release_30/')

		with open(outfile, 'wb') as fp:
			ftp.retrbinary('RETR %s' % hg38file, fp.write)

		ftp.quit()


	def clean_hg38(self):

		GRCh38 = "%s/genomes/GRCh38.primary_assembly.genome.fa.gz" % self.rootDir

		# checkCmnd = "zcat %s | grep -in '>' " % GRCh38
		# print checkCmnd
		# subprocess.check_output(checkCmnd, shell=True)

		gunzipCmnd = "gunzip %s" % GRCh38
		subprocess.Popen(gunzipCmnd, shell= True).wait()

		### 51471479 last line of valid chromosomes
		cleanCmnd = "sed -n '1,51471479p' %s/genomes/GRCh38.primary_assembly.genome.fa > %s/genomes/hg38.fa" % (self.rootDir, self.rootDir)
		subprocess.Popen(cleanCmnd, shell= True).wait()

		indexCmnd = "samtools faidx %s/genomes/hg38.fa" % (self.rootDir)
		subprocess.Popen(indexCmnd, shell= True).wait()

		TwoBitCmnd = "faToTwoBit %s/genomes/hg38.fa %s/genomes/hg38.2bit" % (self.rootDir, self.rootDir)
		subprocess.Popen(TwoBitCmnd, shell= True).wait()

		removeGRCh38 = "rm %s/genomes/GRCh38.primary_assembly.genome.fa" % (self.rootDir)
		subprocess.Popen(removeGRCh38, shell= True).wait()


	def download_Gencode_annotation(self):

		outfile = "%s/genomes/gencode.v30.annotation.gtf.gz" % (self.rootDir)
		annotationFile = "gencode.v30.annotation.gtf.gz"

		outfile_tRNA = "%s/genomes/gencode.v30.tRNAs.gtf.gz" % (self.rootDir)
		annotationFiletRNA = "gencode.v30.tRNAs.gtf.gz"

		ftp = FTP('ftp.ebi.ac.uk')
		ftp.login()
		ftp.cwd('pub/databases/gencode/Gencode_human/release_30/')

		with open(outfile, 'wb') as fp:
			ftp.retrbinary('RETR %s' % annotationFile, fp.write)

		with open(outfile_tRNA, 'wb') as fp:
			ftp.retrbinary('RETR %s' % annotationFiletRNA, fp.write)		

		ftp.quit()

	def parse_GTF_file(self):

		gtfInFile = "%s/genomes/gencode.v30.annotation.gtf.gz" % (self.rootDir)

		parseCmnd = "python2 %s/utils/GTF_hg38_protCode_termStop_validUTR.py --gtfInFile %s --rootDir %s" % (
			self.rootDir, gtfInFile, self.rootDir)

		subprocess.Popen(parseCmnd, shell= True).wait()

	def build_annotation_files(self):

		gtfInFilePrefix = "%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs" % (self.rootDir)
		twoBitGenome = "%s/genomes/hg38.2bit" % (self.rootDir)

		mRNA_cmnd = "python2 %s/utils/mRNA_sequences_from_gtf.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s" % (
			self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome)
		subprocess.Popen(mRNA_cmnd, shell= True).wait()

		prot_cmnd = "python2 %s/utils/protein_sequences_from_gtf.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s" % (
			self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome)
		subprocess.Popen(prot_cmnd, shell= True).wait()

		utrTable_cmnd = "python2 %s/utils/makeUTRtable.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s" % (
			self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome)
		subprocess.Popen(utrTable_cmnd, shell= True).wait()

		stopFinder_cmnd = "python2 %s/utils/stopcodon_finder.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s" % (
			self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome)
		subprocess.Popen(stopFinder_cmnd, shell= True).wait()

		stopPositions_cmnd = "python2 %s/utils/utr3_stop_positions.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s" % (
			self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome)
		subprocess.Popen(stopPositions_cmnd, shell= True).wait()

		# uORFs_cmnd = "python2 %s/utils/uORF_finder.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s" % (
		# 	self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome)
		# subprocess.Popen(uORFs_cmnd, shell= True).wait()

		# codons_cmnd = "python2 %s/utils/codonFinder.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s" % (
		# 	self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome)
		# subprocess.Popen(uORFs_cmnd, shell= True).wait()

	# def parse_GTF_allTr(self):

	# 	gtfInFile = "%s/genomes/gencode.v30.annotation.gtf.gz" % (self.rootDir)
	# 	parseCmnd = "python2 %s/utils/GTF_hg38_all_tr.py --gtfInFile %s --rootDir %s" % (
	# 		self.rootDir, gtfInFile, self.rootDir)

	# 	subprocess.Popen(parseCmnd, shell= True).wait()

	# def build_annotation_allTR(self):

	# 	gtfInFilePrefix = "%s/genomes/gencodeV30_all_tr" % (self.rootDir)
	# 	twoBitGenome = "%s/genomes/hg38.2bit" % (self.rootDir)

	# 	utrTable_cmnd = "python2 %s/utils/makeUTRtable.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s" % (
	# 		self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome)
	# 	subprocess.Popen(utrTable_cmnd, shell= True).wait()



	def build_ncRNA_depletion(self):

		gtfInFile = "%s/genomes/gencode.v30.annotation.gtf.gz" % (self.rootDir)
		gtfInFileTrna = "%s/genomes/gencode.v30.tRNAs.gtf.gz" % (self.rootDir)
		rRNAncbi = "%s/genomes/refseq_human_rRNA_all.fa" % (self.rootDir)
		twoBitGenome = "%s/genomes/hg38.2bit" % (self.rootDir)

		ncRNA_depletion_cmnd = "python2 %s/utils/rRNA_depletion_hg38.py --gtfInFile %s --gtfInFileTrna %s --rRNAncbi %s --rootDir %s --twoBitGenome %s" % (
			self.rootDir, gtfInFile, gtfInFileTrna, rRNAncbi, self.rootDir, twoBitGenome)
		subprocess.Popen(ncRNA_depletion_cmnd, shell= True).wait()


	def build_STAR_indexes(self):

		### ncRNA
		ncSparsity = 1
		ncGenomeDir = "%s/genomes/star_hg38_ncRNA" % self.rootDir
		ncGenomeFasta = "%s/genomes/gencodeV30_ncRNA_all.fa" % self.rootDir
		ncSjdbGTF = "0"
		ncSAindexNbases = 9

		ncRNA_build_STAR_index = "python2 %s/utils/buildStarIndex.py --rootDir %s --threadNumb %s --STARsparsity %s --genomeDir %s --genomeFastaFiles %s --sjdbGTF %s --SAindexNbases %s" % (
			self.rootDir, self.rootDir, self.threadNumb, ncSparsity, ncGenomeDir, ncGenomeFasta, ncSjdbGTF, ncSAindexNbases)
		subprocess.Popen(ncRNA_build_STAR_index, shell=True).wait()

		### hg38
		## need to unzip GTF file to work with STAR
		unzip_cmnd = "gunzip %s/genomes/gencode.v30.annotation.gtf.gz" % self.rootDir
		subprocess.Popen(unzip_cmnd, shell=True).wait()


		hgSparsity = 1
		hgGenomeDir = "%s/genomes/star_gtf_gencodeV30annotation" % self.rootDir
		hgGenomeFasta = "%s/genomes/hg38.fa" % self.rootDir
		hgSjdbGTF = "%s/genomes/gencode.v30.annotation.gtf" % self.rootDir
		hgSAindexNbases = 14

		hg38_build_STAR_index = "python2 %s/utils/buildStarIndex.py --rootDir %s --threadNumb %s --STARsparsity %s --genomeDir %s --genomeFastaFiles %s --sjdbGTF %s --SAindexNbases %s" % (
			self.rootDir, self.rootDir, self.threadNumb, hgSparsity, hgGenomeDir, hgGenomeFasta, hgSjdbGTF, hgSAindexNbases)
		subprocess.Popen(hg38_build_STAR_index, shell=True).wait()

		gzip_cmnd = "gzip %s/genomes/gencode.v30.annotation.gtf" % self.rootDir
		subprocess.Popen(gzip_cmnd, shell=True).wait()



##### 2) Retrieve RawData

class RawData(object):

	### building all files need for downstream analysis

	def __init__(self, rootDir, threadNumb):
		self.rootDir = rootDir
		self.threadNumb = threadNumb


	# def FASTQ_dump_sequences(self):
	# 	"""
	# 	Download all FASTQ_files in appropriate directories for downstream analysis
	#	*** Pending release of GEO datasets ***
	# 	"""

	# 	fq_dump_cmnd = "python2 %s/utils/fastqDump_sri37240.py --rootDir %s --threadNumb %s" % (
	# 		self.rootDir, self.rootDir, self.threadNumb)
	# 	subprocess.Popen(fq_dump_cmnd, shell=True).wait()


	def merge_allAG_experiment(self):
		"""
		merge fastq files from allAG profiling experiment
		"""

		inputDir = "%s/Data/RPF/FASTQ" % self.rootDir
		outputDir = "%s/Data/RPF/FASTQ/merge" % self.rootDir

		if not os.path.exists(outputDir):	os.makedirs(outputDir)

		mergeFQ_cmnd = "python2 %s/utils/merge_fastq.py --inputDir %s --outputDir %s" % (
			self.rootDir, inputDir, outputDir)
		subprocess.Popen(mergeFQ_cmnd, shell=True).wait()

##### 3) Run main analysis pipeline on ribosome profiling data

class RibosomeProfiling_workflow(object):

	def __init__(self, rootDir, threadNumb, libSetFile):
		self.rootDir = rootDir
		self.threadNumb = threadNumb
		self.libSetFile = libSetFile


	def RPexp(self):
		"""
		run the main ribosome profiling workflow for all samples here
		"""

		riboseq_cmnd = "python2 %s/riboseq/riboseq_main.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			self.rootDir, self.rootDir, self.libSetFile, self.threadNumb)

		subprocess.Popen(riboseq_cmnd, shell=True).wait()

	def RP_raw_countTables(self):
		"""
		build raw count tables for RNA seq files
		"""
		raw_countTables_cmnd = "python2 %s/riboseq/riboseq_build_exp_RAW_countTables.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			self.rootDir, self.rootDir, self.libSetFile, self.threadNumb)
		subprocess.Popen(raw_countTables_cmnd, shell=True).wait()

	def RP_avgene_cdsNorm_start(self):

		avgene_cmnd = "python2 %s/riboseq/riboseq_avggene_cdsNorm_start.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			self.rootDir, self.rootDir, self.libSetFile, self.threadNumb)
		subprocess.Popen(avgene_cmnd, shell=True).wait()

	def RP_avgene_cdsNorm_stop(self):

		avgene_cmnd = "python2 %s/riboseq/riboseq_avggene_cdsNorm_stop.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			self.rootDir, self.rootDir, self.libSetFile, self.threadNumb)
		subprocess.Popen(avgene_cmnd, shell=True).wait()

	def RP_stopCodon_pauseScore(self):

		gtfInFilePrefix = "%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs" % (self.rootDir)

		stop_pause_cmnd = "python2 %s/riboseq/riboseq_stopCodon_pauseScore.py --rootDir %s --libSetFile %s --threadNumb %s --gtfInFilePrefix %s" % (
			self.rootDir, self.rootDir, self.libSetFile, self.threadNumb, gtfInFilePrefix)
		subprocess.Popen(avgene_cmnd, shell=True).wait()

	# def RP_codon_occ(self):

	# 	gtfInFilePrefix = "%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs" % (self.rootDir)

	# 	codon_occ_cmnd = "python2 %s/riboseq/riboseq_codon_occ_workflow.py --gtfInFilePrefix %s --rootDir %s --libSetFile %s --threadNumb %s" % (
	# 		self.rootDir, gtfInFilePrefix, self.rootDir, self.libSetFile, self.threadNumb)
	# 	subprocess.Popen(codon_occ_cmnd, shell=True).wait()

	# def densebuild_allTr(self):

	# 	dballTr_cmnd = "python2 %s/riboseq/riboseq_allTr_wf.py --rootDir %s --libSetFile %s --threadNumb %s" % (
	# 		self.rootDir, self.rootDir, self.libSetFileAllTr, self.threadNumb) 
	# 	subprocess.Popen(dballTr_cmnd, shell=True).wait()



# class RNAseq_workflow(object):

# 	def __init__(self, rootDir, threadNumb, libSetFile):
# 		self.rootDir = rootDir
# 		self.threadNumb = threadNumb
# 		self.libSetFile = libSetFile

# 	def RNAexp(self):

# 		rnaseq_cmnd = "python2 %s/RNAseq/RNAseq_main.py --rootDir %s --libSetFile %s --threadNumb %s" % (
# 			self.rootDir, self.rootDir, self.libSetFile, self.threadNumb)
# 		subprocess.Popen(rnaseq_cmnd, shell=True).wait()

# 	def RNA_raw_countTables(self):
# 		raw_countTables_cmnd = "python2 %s/RNAseq/RNAseq_build_exp_RAW_countTables.py --rootDir %s --libSetFile %s --threadNumb %s" % (
# 			self.rootDir, self.rootDir, self.libSetFile, self.threadNumb)
# 		subprocess.Popen(raw_countTables_cmnd, shell=True).wait()


### plot figures:

class plot_figures(object):

	def __init__(self, rootDir, threadNumb, 
			libSet_RP_all, libSet_RP_merge):
		self.rootDir = rootDir
		self.threadNumb = threadNumb
		self.libSet_RP_all = libSet_RP_all
		self.libSet_RP_merge = libSet_RP_merge


	def plot_figure_6(self):

		fig_cmnd_6A = "python2 %s/figures/figscripts/plot_figure6A.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			self.rootDir, self.rootDir, self.libSet_RP_merge, self.threadNumb)
		subprocess.Popen(fig_cmnd_6A, shell=True).wait()

		fig_cmnd_6B = "python2 %s/figures/figscripts/plot_figure6B.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			self.rootDir, self.rootDir, self.libSet_RP_merge, self.threadNumb)
		subprocess.Popen(fig_cmnd_6B, shell=True).wait()

		fig_cmnd_6C = "python2 %s/figures/figscripts/plot_figure6C.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			self.rootDir, self.rootDir, self.libSet_RP_merge, self.threadNumb)
		subprocess.Popen(fig_cmnd_6C, shell=True).wait()



##### MAIN #####

def main():

	### 1) Generate Genome Files
	genomeGen = generateGenomes(rootDir, threadNumb)
	genomeGen.download_hg38()
	genomeGen.clean_hg38()
	genomeGen.download_Gencode_annotation()
	genomeGen.parse_GTF_file() ## this takes awhile
	genomeGen.build_annotation_files()
	genomeGen.parse_GTF_allTr() ## this takes a very long time
	genomeGen.build_annotation_allTR()
	genomeGen.build_ncRNA_depletion()
	genomeGen.build_STAR_indexes()

	### 2) Process Raw Data
	rawData = RawData(rootDir, threadNumb)
	# rawData.FASTQ_dump_sequences()
	rawData.merge_allAG_experiment()

	### 3) Run Ribosome Profiling analysis pipeline
	RP = RibosomeProfiling_workflow(rootDir, threadNumb, 
		libSetFile="riboseq_libsettings_all")
	RP.RPexp()
	RP.RP_raw_countTables()
	RP.RP_avgene_cdsNorm_stop()

	RP2 = RibosomeProfiling_workflow(rootDir, threadNumb, 
		libSetFile="riboseq_libsettings_merge")
	RP2.RPexp()
	RP2.RP_raw_countTables()
	RP2.RP_avgene_cdsNorm_stop()
	RP2.RP_stopCodon_pauseScore()

	### 4) Plot Figures
	pltFig = plot_figures(rootDir, threadNumb,
		libSet_RP_all = 'riboseq_libsettings_allG418',
		libSet_RP_merge = 'riboseq_libsettings_allAGmerge')
	pltFig.plot_figure_6()


if __name__ == '__main__':
	main()





