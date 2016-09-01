# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 10:45:24 2015

@author: gurunath
"""
import numpy as np
#import waveio as io
#import matplotlib.pyplot as plt


def remTrend(sig, winSize):
    
    # Added on 2-May-2015
    if np.isnan(winSize):
        winSize = 250 # choose by default 4ms window length
#        print "The default window size is choosen: " + '\t' +  str(winSize) + " samples"
             
    window = np.ones(winSize)
    
    rm = np.convolve(sig, window)
    rm = rm[winSize/2:np.shape(rm)[0] - (winSize/2)+1] 
    
    norm = np.convolve(np.ones(np.shape(sig)[0]), window)
    norm = norm[winSize/2:np.shape(norm)[0]-(winSize/2)+1]
    
    rm = rm/norm # Check whether this is doing proper division
    out = sig - rm
    return out


def zeroFreqFilt(wav, fs, winLen):
    dwav = np.diff(wav)
    dwav = np.append(dwav, dwav[-1])
    dwav = dwav/np.max(np.abs(dwav))
    N = np.shape(dwav)[0]
    
    zfSig = np.cumsum(np.cumsum(np.cumsum(np.cumsum(dwav))))
        
    winLen = np.round((winLen*fs)/1000.0)
    zfSig = remTrend(zfSig, winLen)
    zfSig = remTrend(zfSig, winLen)
    zfSig = remTrend(zfSig, winLen)
    zfSig[N-winLen*2:N] = 0
    zfSig[0:winLen*2] = 0
    return zfSig

def zff(wav, fs, winLen):
#    print "Window length is not calculated using autocor in zff.py"
    zf = zeroFreqFilt(wav, fs, winLen)
#    plt.figure(10)
#    plt.plot(zf)
    zf1 = np.copy(zf)
    zf[zf>0] = 1
    zf[zf<0] = -1
    gci = np.where(np.diff(zf) == 2)[0]
    es = np.abs(zf1[gci+1]-zf1[gci-1])
    T0 = np.diff(gci)
    T0 = T0/float(fs)
    f0 = 1.0/T0
    f0 = np.append(f0, f0[-1])
    return zf1, gci, es, f0

#inputFile = 'out_female.wav'
#(fs, x) = io.wavread(inputFile)
#x = np.array(x, np.float64)
#[zf, gci, es, f0] = zff(x, fs, 5)


    