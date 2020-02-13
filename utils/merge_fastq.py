### Script for merging fastq files from seperate experiments

import sys
import os
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--inputDir', help= 'directory with fastq files')
parser.add_argument('--outputDir', help = 'directory to send output')

args = parser.parse_args()

inpath = args.inputDir
outpath = args.outputDir

fq1 = [
	"1_dmso_A",
	"4_g418_A",
	"7_sri37240_A"
]

fq2 = [
	"2_dmso_B",
	"5_g418_B",
	"8_sri37240_B"
]

fq3 = [
	"3_dmso_C",
	"6_g418_C",
	"9_sri37240_C"
]


fq_merged = [

	"1_dmso",
	"2_g418",
	"3_sri372340"
]


# if not os.path.exists(FASTQpath): os.makedirs(FASTQpath)

def mergeFastQ(fq1Input, fq2Input, fq3Input, fqOutput):
	fq1 = '%s/%s*.fastq.gz' % (inpath, fq1Input)
	fq2 = '%s/%s*.fastq.gz' % (inpath, fq2Input)
	fg3 = '%s/%s*.fastq.gz' % (inpath, fq3Input)

	fqOut = '%s/%s.fastq.gz' % (outpath, fqOutput)

	merge_command = 'cat %s %s %s > %s' % (fq1, fq2, fq3, fqOut)
	print merge_command
	os.system(merge_command)

for sample in range(len(fq_merged)):
	mergeFastQ(fq1[sample], fq2[sample], fq3[sample], fq_merged[sample])
