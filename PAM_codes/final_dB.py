#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 20:30:24 2024

@author: pp
"""

from scipy import special as sp
import numpy as np
import matplotlib.pyplot as plt

def SNRb_dB(Pb, M):
    return 10 * np.log10((1/3) * (sp.erfcinv((Pb * M * np.log2(M)) / (M - 1))**2) * (((M**2)-1) / (np.log2(M))))

# Generate M values
MM = np.array([2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]).astype(int)

# Store SNRb values for each M
SNRb_values_dB = [SNRb_dB(10 ** -5, M) for M in MM]



plt.close('all')
plt.figure()   

# Plot the SNRb values in decibels
plt.plot(MM, SNRb_values_dB, marker='o', linestyle='-', label='Pb=1e-5')

plt.legend()
plt.title('SNRb with respect to M (PAM)')
plt.xlabel('M')
plt.ylabel('SNRb [dB]')
plt.grid(True)

# Save the plot as a PNG file with M in the filename
plt.savefig('SNRb_M.png')

