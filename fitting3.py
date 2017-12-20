#! /usr/local/Cellar/python/2.7.12/bin/python

import numpy as np
import pandas as pd
import seaborn as sns
import sys



def main():
    if sys.argv.__len__()>2:
        sns.plt.figure(1)
        sns.set_style("whitegrid")

        df = pd.read_csv(sys.argv[1])
        n = np.array(df["N"])
        t = np.array(df["CPU"])
        nn = np.linspace(1, np.max(n), 1000)
        p = np.polyfit(np.log10(n),np.log10(t),float(sys.argv[2]))
        fit = np.poly1d(p)
        sns.plt.loglog(n, t, c='k', marker='o',linestyle='',basex=2)
        sns.plt.loglog(nn, 10 ** fit(np.log10(nn)), 'k--', label='CPU',basex=2)

        t = np.array(df["GPUin"])
        nn = np.linspace(1, np.max(n), 1000)
        p = np.polyfit(np.log10(n),np.log10(t),float(sys.argv[2]))
        fit = np.poly1d(p)
        sns.plt.loglog(n, t, c='b', marker='o',linestyle='')
        sns.plt.loglog(nn, 10 ** fit(np.log10(nn)), 'b--', label='HW10GPU Inclusive',basex=2)

        t = np.array(df["GPUex"])
        nn = np.linspace(1, np.max(n), 1000)
        p = np.polyfit(np.log10(n),np.log10(t),float(sys.argv[2]))
        fit = np.poly1d(p)
        sns.plt.loglog(n, t, c='r', marker='o',linestyle='')
        sns.plt.loglog(nn, 10 ** fit(np.log10(nn)), 'r--', label='HW10GPU Exclusive',basex=2)

        t = np.array(df["Thrustin"])
        nn = np.linspace(1, np.max(n), 1000)
        p = np.polyfit(np.log10(n),np.log10(t),float(sys.argv[2]))
        fit = np.poly1d(p)
        sns.plt.loglog(n, t, c='b', marker='*',linestyle='')
        sns.plt.loglog(nn, 10 ** fit(np.log10(nn)), 'b-', label='ThrustGPU Inclusive',basex=2)

        t = np.array(df["Thrustex"])
        nn = np.linspace(1, np.max(n), 1000)
        p = np.polyfit(np.log10(n),np.log10(t),float(sys.argv[2]))
        fit = np.poly1d(p)
        sns.plt.loglog(n, t, c='r', marker='*',linestyle='')
        sns.plt.loglog(nn, 10 ** fit(np.log10(nn)), 'r-', label='ThrustGPU Exclusive',basex=2)


        sns.plt.title('Problem 1 Scaling Analysis \n(Euler99)')
        sns.plt.ylabel('Time (ms)')
        sns.plt.xlabel('N')
        sns.plt.legend(loc=4)



        sns.plt.xlim([2**0,2**24])
        # sns.plt.ylim([0,3000])
        sns.plt.savefig('./res.pdf')
        sns.plt.show()
    else:
        print "Usage: ./fitting.py path_to_csv fit_order"

if __name__ == '__main__':
    main()

