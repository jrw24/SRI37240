"""
191210

This script calculates pause scores for the stop codons of all transcripts in an annotation
	The output is a csv file containing transcript names, and stop codon pause scores
"""



import sys
import os
import csv
from Bio import SeqIO
import twobitreader
import struct
import pandas as pd
import importlib
import argparse


### import libsettings file and add to global namespace
parser= argparse.ArgumentParser()
parser.add_argument('--rootDir', help= 'the root directory containing data and scripts')
parser.add_argument('--libSetFile', help= 'riboseq libsettings file to run riboseq_main')
parser.add_argument('--threadNumb', help= 'number of threads')
parser.add_argument('--gtfInFilePrefix', help='prefix used for genome files')
args = parser.parse_args()
sys.path.append("%s/riboseq" % args.rootDir)
sys.path.append("%s/riboseq/libsettings" % args.rootDir)
import rphelper as rph

rootDir = args.rootDir
libsetName = args.libSetFile
libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)
threadNumb = str(args.threadNumb)



### inputs
mRNAdf = pd.read_csv('%s/genomes/%s_mRNAseqs.csv' % (rootDir, args.gtfInFilePrefix), sep=',', index_col=0)

print samplelist

threshold = 10
customSize = 30 ### not used generally
pop = "fl"


assignment = "5"
ribosome_shift = "A"
norm_type = 'rpm'
densitystring = 'Density_rpm'
inset_choice = 'default'
norm = 'eq'


if pop == "fl":
	minlen = str(flmin)
	maxlen = str(flmax)
elif pop == "eA":
	minlen = str(eAmin)
	maxlen = str(eAmax)
elif pop == "eE":
	minlen = str(eEmin)
	maxlen = str(eEmax)
elif pop == "custom":
	minlen = str(customSize)
	maxlen = str(customSize)
else:
	print "read lengths not set"

def load_genomes(UTRfilestring, twobitfile):
	"""
	make this a separate function so that these only need to be loaded a single time
	"""
	UTRdict= rph.readindict(open(UTRfilestring, "rU"))
	genome= twobitreader.TwoBitFile(twobitfile) # do we actually need to load this in here?
	return UTRdict, genome

def build_scDense(UTRdict, threshold, samp, inset_choice):

	print "*** running build_scDense ***"

	fp_assign_path = '%s/FPassignment/%s/%s/%s' % (rootpath, genome_name, experiment, samp)
	trspfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
		fp_assign_path, densitystring, assignment, ribosome_shift, 
		minlen, maxlen, samp, minlen, maxlen, samp, minlen, maxlen)
	totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, samp)
	totreadcountf = open(totreads_countfile, "r")
	totreads = int(totreadcountf.read())
	totreadcountf.close()

	### This is where all the count files are loaded into a dictionary
	trspdict= rph.readcountsf(trspfilestring) ### this takes a minute

	## create a list of 0's that is the length of the region of interest
		## add 3 to account for the stop codon
		## this will be added to for every transcript
	# averagegene= [0 for num in range(0, 3+upstreamNTs+ downstreamNTs)] # add 3 for start or stop codon

	# print averagegene

	## add counters and set to zero
	noUTRentry = 0 ### discard transcripts not in UTRdict
	zeroCdsdense = 0 ### discard transcripts with zero reads in CDS
	lowCdsdense = 0 ## Optional CDS density thresholding, in RPKM
	totalCountedTranscripts = 0 ## number included in final output


	### calculate mRNA-region densities 

	defaultInsets = { 'utr5Inset3' : 6, 'cdsInset5' : 18, 'cdsInset3' : 15, 'utr3Inset5' : 6 }
	zeroInsets    = { 'utr5Inset3' : 0, 'cdsInset5' : 0, 'cdsInset3' : 0, 'utr3Inset5' : 0 }
	customInsets  = { 'utr5Inset3' : 15, 'cdsInset5' : 24, 'cdsInset3' : 15, 'utr3Inset5' : 15 }

	if inset_choice == "default":
		insets = defaultInsets
	elif inset_choice == "zero":
		insets = zeroInsets
	elif inset_choice == "custom":
		insets = customInsets
	else:
		print "Insets were not set"
		sys.exit()

	scDenseDict = {}


	### Iterated through transcripts one at a time, retrieving counts in region of interest:
	for trsp in trspdict:

		# print trsp

		# if trsp != 'ENST00000319248.13': ## testing with PRDX1
		# 	continue

		### Load in count file for the transcript here
		exonsplicedcounts = trspdict[trsp]

		if UTRdict.has_key(trsp)!=True: # check to make sure density file has an annotation in the UTR csv
			noUTRentry +=1
			continue

		mrnalen = int(UTRdict[trsp][3])
		cdslen = int(UTRdict[trsp][4])
		utr5len = int(UTRdict[trsp][5])
		utr3len = int(UTRdict[trsp][6])
		assert mrnalen == cdslen + utr5len + utr3len ## check that this is true

		### define Coding sequence here
		cdsstart = utr5len
		cdsend = len(exonsplicedcounts) - utr3len # cdsend == first position of utr3
		if cdsstart == cdsend:
			print "Error, gene length is 0 for transcript "+ trsp
			sys.exit()

		### Calculate Region Densities ###
		utr5lenMod = utr5len-insets['utr5Inset3']
		cdslenMod = cdslen-insets['cdsInset5']-insets['cdsInset3']
		utr3lenMod = utr3len-insets['utr3Inset5']
		mrnalenMod = utr5lenMod+cdslenMod+utr3lenMod

		utr5Counts = sum(exonsplicedcounts[:cdsstart-insets['utr5Inset3']])
		cdsCounts = sum(exonsplicedcounts[cdsstart+insets['cdsInset5']:cdsend-insets['cdsInset3']])
		utr3Counts = sum(exonsplicedcounts[cdsend+insets['utr3Inset5']:])
		mrnaCounts = utr5Counts+cdsCounts+utr3Counts


		### RAW counts
		RAWutr5Counts = int(utr5Counts*(totreads/1E6))
		RAWutr3Counts = int(utr3Counts*(totreads/1E6))
		RAWcdsCounts = int(cdsCounts*(totreads/1E6))
		RAWmrnaCounts = int(mrnaCounts*(totreads/1E6))

		### denisites
		# mrnaDensity = (mrnaCounts/mrnalenMod) 
		cdsDensity = (cdsCounts/cdslenMod) 
		# utr5Density = (utr5Counts/utr5lenMod)
		# utr3Density = (utr3Counts/utr3lenMod)

		#### RPKM densities
		# mrnaDensity_rpkm = (mrnaCounts/mrnalenMod) * 1000
		cdsDensity_rpkm = (cdsCounts/cdslenMod) * 1000
		# utr5Density_rpkm = (utr5Counts/utr5lenMod) * 1000
		# utr3Density_rpkm = (utr3Counts/utr3lenMod) * 1000

		### throw out zero's
		if cdsDensity == 0:
			zeroCdsdense += 1
			continue

		if cdsDensity*float(1000)< int(threshold):	# Threshold on cds density: (thresholding on "rpkm")
			lowCdsdense += 1
			continue

		### define vector in valid CDS region, normalize by cdsDensity
	#     cdsSplicedCounts = exonsplicedcounts[cdsstart+insets['cdsInset5']:cdsend-insets['cdsInset3']]
		### include FULL CDS for normalization calculation
		cdsSplicedCounts = exonsplicedcounts[cdsstart:cdsend]

		cdsNormCounts = [rpf/cdsDensity for rpf in cdsSplicedCounts] ## just region of cds within insets for counts, sum/len == 1


		assert len(cdsNormCounts) == cdslen
		assert len(exonsplicedcounts) == mrnalen

		# print len(cdsNormCounts)
		# print cdslen
		
		### finally retrieve stop codon normalized counts and calculate the density here
		scNorm = cdsNormCounts[-3:] ## last 3 codons, the stop codon
		scDense = sum(scNorm)/3.0

		### add this to density dictionary
		scDenseDict[trsp] = scDense
		totalCountedTranscripts +=1
		
	### build an output dataframe to write to a csv file
	dfout = pd.DataFrame.from_dict(scDenseDict, orient='index')
	cols = ['scDense']
	dfout.columns = cols
	outdir = "%s/codon" % (fp_assign_path)
	if not os.path.exists(outdir):	os.makedirs(outdir)
	print outdir
	dfout.to_csv("%s/%s_stopCodonDense_cdsNorm.csv" % (outdir, samp), index=True)

	print "*****"
	print "run complete for %s" % samp
	print "no utr entery == %s" % noUTRentry
	print "zero CDS density == %s" % zeroCdsdense
	print "below threshold of %s CDS density == %s" % (threshold, lowCdsdense)
	print "total remaining transcripts == %s" % totalCountedTranscripts
	print "*****"
	print ""


def main():
	print UTRfilestring
	print twobitfile 

	UTRdict, genomes = load_genomes(UTRfilestring, twobitfile)

	for samp in samplelist:
		print samp
		build_scDense(UTRdict, threshold, samp, inset_choice)

if __name__ == '__main__':
	main()

