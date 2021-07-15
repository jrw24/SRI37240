"""
Plotting stop codon pause scores
"""
## plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

## import dependencies
import sys
import os
import math
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.distributions.empirical_distribution import ECDF
pd.set_option('display.max_columns', 50)
import seaborn as sns
from pylab import *
import argparse
import importlib

### import libsettings file and add to global namespace
parser= argparse.ArgumentParser()
parser.add_argument('--rootDir', help= 'the root directory containing data and scripts')
parser.add_argument('--libSetFile', help= 'riboseq libsettings file to run riboseq_main')
parser.add_argument('--threadNumb', help= 'number of threads')
args = parser.parse_args()
sys.path.append("%s/riboseq" % args.rootDir)
sys.path.append("%s/riboseq/libsettings" % args.rootDir)


rootDir = args.rootDir
libsetName = args.libSetFile
libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)
threadNumb = str(args.threadNumb)
sys.path.append('%s/riboseq' % rootDir)


### set colors
black = '#000000'
orange = '#ffb000'
cyan = '#63cfff'
red = '#eb4300'
green = '#00c48f'
pink = '#eb68c0'
yellow = '#fff71c'
blue = '#006eb9'

colorList = [black, orange, cyan]

def log_trans_b10(x):
    try:
        return math.log(x, 10)
    except:
        return float(-6.00)
#         return float("NaN")
def log_trans_b2(x):
    try:
        return math.log(x, 2)
    except:
        # return float("NaN")
        return float(-15.00) # set arbitrarily low value

def corrfunc(x, y, **kws):
    r, _ = stats.pearsonr(x, y)
    rho, _ = stats.spearmanr(x, y)
    ax = plt.gca()
    ax.annotate("r = {:.3f}".format(r),
                xy=(.1, .9), xycoords=ax.transAxes)
    ax.annotate(u"p = {:.3f}".format(rho),
                xy=(.1, .85), xycoords=ax.transAxes)

def load_countTables():
	cdsThresh = 10 ## threshold used to calc sc Pause score, not currently in use

	FPassignpath = "%s/FPassignment/%s/%s" % (rootpath, genome_name, experiment)
	namelist = []
	dflist= []

	for samp in samplelist:

	    ### scDense Counts
	    dftemp = pd.read_csv('%s/%s/codon/%s_stopCodonDense_cdsNorm.csv' % (FPassignpath, samp, samp),
	                        index_col = 0)
	    dftemp['sampname'] = samp
	    dftemp['scDenseLog2'] = dftemp['scDense'].apply(log_trans_b2)
	    dftemp.rename(columns={"scDense": "scDense_%s" % samp, "scDenseLog2":"scDenseLog2_%s" % samp}, inplace=True)
	    dflist.append(dftemp)
	    namelist.append(samp)
	    
	### combine into one master dataframe
	df = pd.concat(dflist, axis=0, ignore_index = True)
	### merge all dataframes into single data frame, that has only transcripts present in all other dataframes
	counter = 0
	dfm = dflist[0].copy()

	for i in range(len(dflist))[1:]:
	    print len(dfm)
	    dfm =dfm.merge(dflist[i], left_index=True, right_index=True)
	    
	return dfm

def plot_stop_ecdf(dfm):

	figout = "%s/figures/Fig3B.pdf" % rootDir
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))
	counter = 0
	for i in samplelist:
	    ecdf = ECDF(dfm['scDenseLog2_%s' % i])
	    ax.plot(ecdf.x, ecdf.y, color=colorList[counter])
	    counter +=1

	# x0, x1 = ax.get_xlim()
	# y0, y1 = ax.get_ylim()
	# lims = [max(x0, y0), min(x1, y1)]

	ax.plot((-4,6), (1.0,1.0), ':k')
	ax.plot((-4,6), (0.5,0.5), ':k')
	ax.plot((-4,6), (0.0,0.0), ':k')
	ax.set_xlim(-2,6)

	plt.savefig(figout, format='pdf', bbox_inches = "tight")


def main():
	dfm = load_countTables()
	plot_stop_ecdf(dfm)

if __name__ == '__main__':
	main()

