#! /usr/local/Cellar/python/2.7.12/bin/python

import numpy as np
import seaborn as sns
import pandas as pd 
import sys
import os
import re


def main(path):
	df = pd.read_csv(path,header=None)
	sns.distplot(df)
	pname = re.split('.out',path)[0]
	pname = pname.title()
	title = "NLSQ Estimates for " + pname + " Parameter"
	sns.plt.title(title)

	data = np.array(df[0]);
	print "{0} >> MEAN: {1}, MEDIAN: {2}, STD: {3}".format(os.path.splitext(os.path.basename(path))[0], np.mean(data), np.median(data), np.std(data))

	sns.plt.savefig(os.path.splitext(path)[0]+'_plots.png')
	sns.plt.clf()



if __name__ == '__main__':
	for n in range(1,len(sys.argv)):
		main(path=sys.argv[n])