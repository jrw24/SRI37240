### python script for downloading fastq files

import os, sys
import subprocess
from collections import OrderedDict
from pathos.multiprocessing import ProcessingPool as Pool
import argparse

### import libsettings file and add to global namespace
parser= argparse.ArgumentParser()
parser.add_argument('--rootDir', help= 'the root directory containing data and scripts')
parser.add_argument('--threadNumb', help= 'number of threads')
args = parser.parse_args()

rootDir = args.rootDir
threadNumb = str(args.threadNumb)


### 9 ribosome profiling samples to download 


outdir_RP = '%s/Data/RPF/FASTQ' % rootDir

if not os.path.exists(outdir_RP):	os.makedirs(outdir_RP)


srrList_RP = [
	'SRR10957150',
	'SRR10957151',
	'SRR10957152',
	'SRR10957153',
	'SRR10957154',
	'SRR10957155',
	'SRR10957156',
	'SRR10957157',
	'SRR10957158'
]

samplelist_RP = [
	"1_dmso_A",
	"2_dmso_B",
	"3_dmso_C",
	"4_g418_A",
	"5_g418_B",
	"6_g418_C",
	"7_sri37240_A",
	"8_sri37240_B",
	"9_sri37240_C"
]

RP_dict = OrderedDict(zip(srrList_RP, samplelist_RP))


def get_FASTQ_RP(srr):

	command_to_run = "fastq-dump --readids --split-files --skip-technical --gzip %s --outdir %s" % (srr, outdir_RP)
	print command_to_run
	os.system(command_to_run)


def rename_RP_fastq():
	for srr in RP_dict:
		rename_cmnd = "mv %s/%s_1.fastq.gz %s/%s_RP.fastq.gz" % (
			outdir_RP, srr, outdir_RP, RP_dict[srr])
		os.system(rename_cmnd)

def main():

	print "starting fastq-dump"

	p = Pool(nodes=threadNumb)
	p.map(get_FASTQ_RP, srrList_RP)
	rename_RP_fastq()


if __name__ == '__main__':
	main()






###