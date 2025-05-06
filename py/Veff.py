import sys
import os

import numpy as np

from cosmoprimo.fiducial import DESI
distance = DESI(engine='class').comoving_radial_distance

#DESI P(k=0.14) pre recon values
P0_used4veff = {'bgs':9200,'lrg1':8900,'lrg2':8900,'lrg3':8400,'elg1':2600,'elg2':2900}
P0_dr1 = {'bgs':9600,'lrg1':9550,'lrg2':9700,'lrg3':9400,'elg1':2780,'elg2':3100} #directly looking at pk files taking ~mean of k=0.1375 and 0.1425
P0_dr2 = {'bgs':7900,'lrg1':9600,'lrg2':9750,'lrg3':9480,'elg1':2850,'elg2':3100} #directly looking at pk files taking ~mean of k=0.1375 and 0.1425

def calc_veff_int(Ntot,area,zbin_edges,zdist,P0,areafac=1,h=.676):
    #Ntot should be the total number of galaxies; this information is often the easiest to find
    #area should be the sky coverage in deg2
    #zbin_edges contains the bin edges for zdist (should have length 1 greater than zdist)
    #zdist contains fraction of galaxies such that frac*Ntot gives the number in the bin
    fullsky = 360*360/np.pi
    veff = 0
    disl = distance(zbin_edges)
    for zbin in range(0,len(zdist)):
        ngal = zdist[zbin]*Ntot
        fsky = area/fullsky
        
        vshell = 4/3.*np.pi*(disl[zbin+1]**3.-disl[zbin]**3.)

        vtot = vshell*fsky
        nbar = ngal/vtot
        print('np is '+str(round(nbar*P0,2)))
        veff += ((nbar*P0)/(1+nbar*P0))**2.*vtot
    return veff/h**3./1e9

