# -*- coding: utf-8 -*-
"""
Created on Sat May 30 17:47:04 2015
Description:
Inputs:
Outputs:
@author: Gurunath Reddy M

"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as matio

def plotSectrogram(mx, maxFreq, numFrames, fs, N, H):
    #maxFreq = 2000 # Hz
    binMax = np.floor(maxFreq * N/fs)
    numFrames = int(mx[:,0].size)
    frmTime = H*np.arange(numFrames)/float(fs)                             
    binFreq = np.arange(binMax)*float(fs)/N                         
    plt.figure()
    plt.pcolormesh(frmTime, binFreq, np.transpose(mx[:, :binMax]))
    #plt.title('mx, M=1001, N=1024, H=256')
    plt.autoscale(tight=True)
    plt.tight_layout()
    #plt.savefig('spectrogram.png')
    plt.show() 
    matio.savemat('/home/gurunath/MS_work/papers/Indicon-2015/happy_synth/code_to _generate_figures/orig/synthSpectTime.mat', mdict={'synthSpectTime':frmTime})
    matio.savemat('/home/gurunath/MS_work/papers/Indicon-2015/happy_synth/code_to _generate_figures/orig/synthSpectBinfreq.mat', mdict={'synthSpectBinfreq':binFreq})
    matio.savemat('/home/gurunath/MS_work/papers/Indicon-2015/happy_synth/code_to _generate_figures/orig/synMagSpect.mat', mdict={'synMagSpect':mx})

