#! /usr/local/Cellar/python/2.7.12/bin/python

import numpy as np
import pandas as pd
import seaborn as sns
import sys



def main():
    if sys.argv.__len__()>2:
        sns.plt.figure(1)
        sns.set_style("whitegrid")

        # sns.plt.plot([0,40],[598.877302,598.877302],'r--',label='Serial Matlab N = 10000')
        sns.plt.plot([0,40],[237.478829,237.478829],'r-',label='Reference: Parallel Matlab (4 Threads)')

        ORIG_DF = pd.read_csv(sys.argv[1])

        df = ORIG_DF.loc[ORIG_DF["N"]==10000]
        print df
        n = np.array(df["NT"])
        t = np.array(df["T"])
        nn = np.linspace(1, np.max(n), 1000)
        p = np.polyfit(np.log10(n),np.log10(t),float(sys.argv[2]))
        fit = np.poly1d(p)
        sns.plt.plot(n, t, c='b', marker='o',linestyle='')
        sns.plt.plot(nn, 10 ** fit(np.log10(nn)), 'b--', label='Parallel OpenMP (1 node)')

        
        ORIG_DF = pd.read_csv('mpi_timing.csv')

        df = ORIG_DF.loc[ORIG_DF["N"]==10000]
        print df
        n = np.array(df["NT"])
        t = np.array(df["T"])
        nn = np.linspace(1, np.max(n), 1000)
        p = np.polyfit(np.log10(n),np.log10(t),float(sys.argv[2]))
        fit = np.poly1d(p)
        sns.plt.plot(n, t, c='k', marker='o',linestyle='')
        sns.plt.plot(nn, 10 ** fit(np.log10(nn)), 'k--', label='Parallel OpenMP & MPI (2 nodes)')

        


        sns.plt.title('Paralleized Monte Carlo Simulations\nN=10,000')
        sns.plt.ylabel('Time (seconds)')
        sns.plt.xlabel('Number of OpenMP Threads per Node')
        sns.plt.legend(loc=5)



        sns.plt.xlim([1,40])
        sns.plt.savefig('./res.pdf')
        sns.plt.show()
    else:
        print "Usage: ./fitting.py path_to_csv fit_order"

if __name__ == '__main__':
    main()

