import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as pat

import csv

#  loga, logl, skew, kurto, mean, SDNN, RMSSD;

nr_alpha = []
nr_lambda = []
nr_sdnn = []
nr_rmssd = []
ar_alpha = []
ar_lambda = []
ar_sdnn = []
ar_rmssd = []
with open('normal.csv',newline='') as ncsv:
    spamreader = csv.reader(ncsv, delimiter=' ', quotechar='|')
    for row in spamreader:
        NR = row[0].split(',')
        nr_alpha.append(float(NR[0]))
        nr_lambda.append(float(NR[1]))
        nr_sdnn.append(float(NR[5]) * 1000.)
        nr_rmssd.append(float(NR[6]) * 1000.)

nr = len(nr_alpha)

with open('abnormal.csv') as abcsv:
    spamreader = csv.reader(abcsv, delimiter=' ', quotechar='|')
    for row in spamreader:
        AR= row[0].split(',')
        ar_alpha.append(float(AR[0]))
        ar_lambda.append(float(AR[1]))
        ar_sdnn.append(float(AR[5]) * 1000.)
        ar_rmssd.append(float(AR[6]) * 1000.)

ar = len(ar_alpha)

bins = np.linspace(0,450,100)




# Figure of SDNN
fig1 = plt.figure(1)
nr_patch = pat.Patch(color='red', label='NR')
ar_patch = pat.Patch(color='blue', label='AR')
plt.hist(nr_sdnn, bins,alpha=0.7,color='red')
plt.hist(ar_sdnn, bins,alpha=0.7,color='blue')
plt.legend(handles=[nr_patch, ar_patch])
plt.xlabel("SDNN [ms]")
plt.ylabel("Frequency")

#  Figure of gamma parms.

fig2 , ax = plt.subplots(2,1)
ax[0].scatter(nr_alpha, nr_lambda, color='red', alpha=0.3, marker='.', label='NR',s=20)
ax[0].scatter(ar_alpha, ar_lambda, color='blue',alpha=0.3, marker='x', label='AN',s=20)
ax[0].legend()
ax[0].set_ylabel('log λ')
ax[0].set_xlabel('log α')
ax[1].scatter(nr_alpha, nr_lambda, color='red', alpha=0.7, marker='.', label='NR',s=20)
ax[1].scatter(ar_alpha, ar_lambda, color='blue',alpha=0.7, marker='x', label='AN',s=20)
ax[1].set_ylabel('log λ')
ax[1].set_xlabel('log α')


# Figure of RMSSD

fig3 = plt.figure(3)
nr_patch = pat.Patch(color='red', label='NR')
ar_patch = pat.Patch(color='blue', label='AN')
plt.hist(ar_rmssd, bins,alpha=0.7,color='blue')
plt.hist(nr_rmssd, bins,alpha=0.7,color='red')

plt.legend(handles=[nr_patch, ar_patch])
plt.xlabel("RMSSD [ms]")
plt.ylabel("Frequency")


# Figure of gamma params with RMSSD

fig4 = plt.figure(4)
ax = fig4.add_subplot(111, projection='3d')
ax.scatter(nr_alpha, nr_lambda, nr_rmssd, color='red', alpha=0.3, marker='.', label='NR',s=20)
ax.scatter(ar_alpha, ar_lambda, ar_rmssd, color='blue',alpha=0.3, marker='x', label='AN',s=20)
ax.legend()
ax.set_ylabel('log λ')
ax.set_xlabel('log α')
ax.set_zlabel('RMSSD [ms]')


#  Figure of SDNN and RMSSD

fig2 , ax = plt.subplots(2,1)
ax[0].scatter(nr_sdnn, nr_rmssd, color='red', alpha=0.3, marker='.', label='NR',s=20)
ax[0].scatter(ar_sdnn, ar_rmssd, color='blue',alpha=0.3, marker='x', label='AN',s=20)
ax[0].legend()
ax[0].set_ylabel('log λ')
ax[0].set_xlabel('log α')
ax[1].scatter(nr_sdnn, nr_rmssd, color='red', alpha=0.3, marker='.', label='NR',s=20)
ax[1].scatter(ar_sdnn, ar_rmssd, color='blue',alpha=0.3, marker='x', label='AN',s=20)
ax[1].set_ylabel('log λ')
ax[1].set_xlabel('log α')

plt.show()
