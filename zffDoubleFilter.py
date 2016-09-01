# -*- coding: utf-8 -*-
"""
Created on Sat May 30 12:25:32 2015
Description: Finding refined resonant contour for filtering the ZFFS signal
Inputs: Zffs signal, fs, window, FFT points, Hop size
Outputs: Resonant frequency contour for filering the ZFFS signla second time
@author: Gurunath Reddy M

"""

import numpy as np
import dftStft as stft
import matplotlib.pyplot as plt
import plotSpect as pltSpec
import getMusicZff as getZFFS


def filterZffs(zfss, fs, w, N, H, sampleOnsets):

    mx, px = stft.stftAnal(zfss, fs, w, N, H) # mx and px are the magnitude and phase spectrum of ZFFS
    tempMx = np.zeros([mx.shape[0], mx.shape[1]])
    numFrames = np.shape(mx)[0]
    maxFreq = 8000.0
    pltSpec.plotSectrogram(mx, maxFreq, numFrames, fs, N, H)
    
    thresh = 20.0
    for i in range(numFrames):
        maxIndi = np.where(mx[i, :] > thresh)[0]    # Get the bins of mag. spectrum greater than threshold 
        tempMx[i, maxIndi] = mx[i, maxIndi]         # Place mag.'s above threshold into a tempMx
    
    zfsF0 = np.argmax(tempMx, axis=1)               # Bin corresponding to max mag. of each frame
    zfsF0 = zfsF0 * fs/N                            # Bin to frequency in Hz
    
    zffsDoub = getZFFS.getZFFS(zfss, fs, H, sampleOnsets, zfsF0) # ZFFS for each voiced region
    print('Reached double filtering')    
    
    return zffsDoub
