#!/usr/bin/env python3

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

# name of the parameter file 
fsetup='3bdf_bins.dat'
fid=open(fsetup);
nsample=int(fid.readline());
ntype=int(fid.readline());
rsetup=[float(i) for i in fid.readline().split()];
rmin=rsetup[0];
rmax=rsetup[1];
rstep=rsetup[2];
#print(rmin, rmax, rstep);
fid.close();

# name of input files
# default output files of "3bdf2s3"
f3B='S3B.dat'
f3F='S3Fluct.dat'
f3I='S3Info.dat'

S3B=np.loadtxt(f3B);
S3Fluct=np.loadtxt(f3F);
S3Info=np.loadtxt(f3I);

n=1000;
epsilon=2*rstep;
x=np.linspace(0+epsilon, rmax-epsilon, n);
idxes=x>rmin;
spl = interpolate.interp1d(S3B[:,0], S3B[:,1]);
S3B_spl=spl(x);
spl = interpolate.interp1d(S3Fluct[:,0], S3Fluct[:,1]);
S3Fluct_spl = np.zeros(n);
S3Fluct_spl[idxes]=spl(x[idxes]);
spl = interpolate.interp1d(S3Info[:,0], S3Info[:,1]);
S3Info_spl = np.zeros(n);
S3Info_spl[idxes]=spl(x[idxes]);
# plt.figure();
# plt.plot(x, S3B_spl)
# plt.plot(S3B[:,0], S3B[:,1], 'o')
# plt.show()
S3Fluct_spl=S3Fluct_spl+S3B_spl+1.0/6;
S3tot=S3Fluct_spl + S3Info_spl;

[print(x[i], S3tot[i], S3Fluct_spl[i], S3Info_spl[i]) for i in range(n)]
