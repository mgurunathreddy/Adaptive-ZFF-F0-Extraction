# -*- coding: utf-8 -*-
"""
Created on Fri May 22 14:32:14 2015
Description:
Inputs:
Outputs:
@author: Gurunath Reddy M

"""

import numpy as np
from scipy.signal import get_window
import dftStft as stft
from scipy.fftpack import fft
import peakDetectCorrect as pdc

def getSinusoids(x, fs):

    # FFT parameter setting
    M = (20*fs)/1000.0 # 20ms window 
    N = 1024
    H = M/4
    window = 'hamming'
    w = get_window(window, M)
    
    # Computing STFT magnitude and phase of the signal
    mx, px = stft.stftAnal(x, fs, w, N, H) # mx and px are the magnitude and phase spectrum
    mxLinear = 10 ** (mx/20.0) # Linear spectrum may not be needed
    
    # -------------------------------------------------------------------------------- #
    # Hamming window spectrum
    hN = N/2     
    hM = M/2
    fftbuffer = np.zeros(N)
    mX1 = np.zeros(N)
    hamWin = np.hamming(M)
    hamWin = hamWin/sum(hamWin)
    fftbuffer[hN-hM:hN+hM] = hamWin
    
    X = fft(fftbuffer)
    mXDb = 20 * np.log10(abs(X)) 
    mX1[:hN] = mXDb[hN:]
    mX1[N-hN:] = mXDb[:hN]
    mXLinear = 10 ** (mX1/20.0)
    print("---------------------------------------------------")
    print("Hamming window main lobe selection needs to be changed for change in FFT parameters")
    print("---------------------------------------------------")
    hamMainLobe = mXLinear[507:518]
    # -------------------------------------------------------------------------------- #
    # Process the spectrum upto only 5kHz
    maxBinFreq = np.round((5000.0 * N)/fs)
    widthHalfHamMainLobe = hamMainLobe.size/2
    
    sinThresh = 0.6
    t = 0.001
    
    detSinusoids = np.zeros([np.shape(mxLinear)[0], maxBinFreq])
    
    #for i in range(np.shape(mxLinear)[0]):
    for i in range(np.shape(mxLinear)[0]):
        tempMag = mxLinear[i, :maxBinFreq]
        ploc = pdc.peakDetection(tempMag, t)                           # detect peak locations   
    #    collAM = np.array([])
    #    collEm = np.array([])
        collS = np.array([])
        # Peak locations are available, now compute similarity between the measured and the ideal sinusoid(spectrum of hamming window)
        index = np.where(ploc <= widthHalfHamMainLobe)[0] # Remove the index which is less than length of half main lobe for correct indexing
        ploc = np.delete(ploc, index)
        index = np.where(ploc > (maxBinFreq - (widthHalfHamMainLobe+1)))[0] # Remove plocation which is close to maxfrequency bin considered
        ploc = np.delete(ploc, index)    
    
        for j in range(np.size(ploc)):
            begSilce = ploc[j]-widthHalfHamMainLobe
            endSilce = ploc[j]+widthHalfHamMainLobe+1
            measuredSin = tempMag[begSilce:endSilce]
            Am = np.sum(hamMainLobe * measuredSin)
            mainLobeEng = np.sum(np.power(hamMainLobe, 2))
            Am = Am/mainLobeEng
            Em = np.sum(np.power((measuredSin - (Am*hamMainLobe)), 2))
            S = 1 - (Em/np.sum(np.power(measuredSin, 2)))
            #collAM = np.append(collAM, Am)
            #collEm = np.append(collEm, Em)
            collS = np.append(collS, S)
        sinPeaksIndx = np.where(collS > sinThresh)[0]
        sinPeaksBin = ploc[sinPeaksIndx]
        sinPeaksMag = tempMag[sinPeaksBin]        
        #detSinusoids[i, sinPeaksBin] = sinPeaksMag # T0 get only the center frequency of each main lobe
        # Place the entired lobe of the detected sinusoid. O.W SSH creates probelm
        for k in range(sinPeaksBin.size):
            detSinusoids[i, sinPeaksBin[k]-5:sinPeaksBin[k]+5] = tempMag[sinPeaksBin[k]-5:sinPeaksBin[k]+5]
        
        
    return detSinusoids, mxLinear, sinPeaksMag, sinPeaksBin    
