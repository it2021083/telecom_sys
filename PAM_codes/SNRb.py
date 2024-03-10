#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 21:20:58 2024

@author: pp
"""

from scipy import special as sp
import numpy as np
import matplotlib.pyplot as plt

def SNRb_dB(Pb, M):
    return ((1/3) * (sp.erfcinv((Pb * M * np.log2(M)) / (M - 1))**2) * (((M**2)-1) / (np.log2(M))))

MM=np.array([4, 8, 16, 32, 64]).astype(int)
n=np.arange(-14,-1,0.1)
Pb=np.array( 10.0** n )

SNRb=np.zeros([Pb.size, MM.size])

plt.close('all')
plt.figure()    
for i, M in enumerate(MM):
    print(i,M)
    SNRb[:,i] = 10.0*np.log10(SNRb_dB( Pb, M ))
    plt.semilogx(Pb,SNRb[:,i],label='$M='+str(M)+'$')
    
plt.legend()
plt.title('$P_b$ with respect to $\mathrm{SNR}_b^{min}$ (PAM)')
plt.xlabel('$P_b$')
plt.ylabel('$\mathrm{SNR}_b^{min}$ [dB]')
plt.grid(True)

# Save the plot as a PNG file with M in the filename
plt.savefig('SNRb_M_Pb.png')