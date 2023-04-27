#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

fid = open('3bdf_bins.ini', 'r');
line = fid.readline();
rho = float(line.split()[0]);
natom = int(line.split()[1]);
typelist = np.array([int(elem) for elem in fid.readline().split()]);
#print(rho, natom, typelist);
fid.close();

fid = open("3bdf_bins.dat_sum", "r");
nstep = int(fid.readline());
ntype = int(fid.readline());
line = fid.readline().split();
rmin = float(line[0]);
rc = float(line[1]);
nr1,nr2,ncost = [int(elem) for elem in fid.readline().split()]
#print(nstep, ntype, rc);
#print(nr1, nr2, ncost);
fid.close();

bin_hist = np.loadtxt('3bdf_bins.dat_sum', skiprows=4);
g3fluct = np.loadtxt('3bdf_g3_fluct.dat');
g3info = np.loadtxt('3bdf_g3_info.dat');
s3mapfluct = np.loadtxt('3bdf_s3map_fluct.dat');
s3mapinfo = np.loadtxt('3bdf_s3map_info.dat');
#l,w = np.shape(bin_hist);
#print (l, w);

#bin_x3 = np.reshape(bin_x, (nr1,nr2,ncost));
bin_hist3 = np.reshape(bin_hist, (nr1,nr2,ncost));
g3info = np.reshape(g3info, (nr1,nr2,ncost));
g3fluct = np.reshape(g3fluct, (nr1,nr2,ncost));
s3mapinfo = np.reshape(s3mapinfo, (nr1,nr2,ncost));
s3mapfluct = np.reshape(s3mapfluct, (nr1,nr2,ncost));

r1 = np.linspace(rmin,rc,nr1);
r2 = np.linspace(rmin,rc,nr2);
cost = np.linspace(rmin,rc,ncost);
dr = (rc-rmin)/nr1;

nv0 =natom * (rho * 4*np.pi/3*rc**3)**2;

xx = np.linspace(rmin+0.1,rc-0.1,10);
for x in xx:
    plt.figure();
    nrx = int((x-rmin)/(rc-rmin)*nr1);
    #h = plt.contourf(cost, r2, s3mapfluct[nrx,:,:], norm=mcolors.CenteredNorm(0.0),cmap='seismic');
    h = plt.contourf(cost, r2, s3mapfluct[nrx,:,:],cmap='seismic');
    plt.xlabel('R3 [Ang]');
    plt.ylabel('R2 [Ang]');
    plt.title("s3map_fluct_r"+str(x));
    plt.colorbar()
    plt.savefig("s3map_fluct_r"+str(int(x*100))+".png");

for x in xx:
    plt.figure();
    nrx = int((x-rmin)/(rc-rmin)*nr1);
    #h = plt.contourf(cost, r2, s3mapinfo[nrx,:,:], norm=mcolors.CenteredNorm(0.0),cmap='seismic');
    h = plt.contourf(cost, r2, s3mapinfo[nrx,:,:],cmap='seismic');
    plt.xlabel('R3 [Ang]');
    plt.ylabel('R2 [Ang]');
    plt.title("s3map_info_r"+str(x));
    plt.colorbar()
    plt.savefig("s3map_info_r"+str(int(x*100))+".png");
    
for x in xx:
    plt.figure();
    nrx = int((x-rmin)/(rc-rmin)*nr1);
    #h = plt.contourf(cost, r2, g3fluct[nrx,:,:], norm=mcolors.CenteredNorm(0.0),cmap='seismic');
    h = plt.contourf(cost, r2, g3fluct[nrx,:,:],cmap='seismic');
    plt.xlabel('R3 [Ang]');
    plt.ylabel('R2 [Ang]');
    plt.title("g3_fluct_r"+str(x));
    plt.colorbar()
    plt.savefig("g3_fluct_r"+str(int(x*100))+".png");

for x in xx:
    plt.figure();
    nrx = int((x-rmin)/(rc-rmin)*nr1);
    #h = plt.contourf(cost, r2, g3info[nrx,:,:], norm=mcolors.CenteredNorm(0.0),cmap='seismic');
    h = plt.contourf(cost, r2, g3info[nrx,:,:],cmap='seismic');
    plt.xlabel('R3 [Ang]');
    plt.ylabel('R2 [Ang]');
    plt.title("g3_info_r"+str(x));
    plt.colorbar()
    plt.savefig("g3_info_r"+str(int(x*100))+".png");

plt.figure();
ainfo=np.cumsum(np.sum(s3mapinfo[:,:,:],axis=(1,2)));
afluct=np.cumsum(np.sum(s3mapfluct[:,:,:],axis=(1,2)));
#plt.plot(r1, afluct, label="fluct");
#plt.plot(r1, ainfo, label="info");
plt.plot(r1, ainfo+afluct, label="fluct+info");
plt.legend();
plt.savefig("s3_vs_r1.png");

