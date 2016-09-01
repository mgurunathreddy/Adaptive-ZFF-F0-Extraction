# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:14:35 2015
Description: Determine the sinusoidal components present in the spectrogram
Inputs: Magnitude spectrum of the music signal
Outputs: The sinusoidal components present in each frame
@author: Gurunath Reddy M

"""

import sys, os
import numpy as np
import waveio as io
from scipy.signal import get_window
import dftStft as stft
from scipy.fftpack import fft
import matplotlib.pyplot as plt
import peakDetectCorrect as pdc
sys.path.append(os.path.join(os.path.dirname(os.path.realpath('__file__')), '/home/gurunath/coursera/audio_signal_peocessing/sms-tools-master/software/models/'))
import utilFunctions as UF
from scipy.signal import medfilt
import f0Detect as f0Detect

#fileName = raw_input('Enter the music file name: ')

fileName = '3'
(fs, x) = io.wavread('segwav/' + fileName +'.wav')
x = np.array(x, np.float64) # Convert the samples to matlab double data type
x = x/(1.01*np.max(np.abs(x)));  # Normalize sample values

## To process only a part of music signal
#dur = np.array([13.377, 13.450]);
#dur = np.round(dur * fs)
#x = x[dur[0]:dur[1]]

# FFT parameter setting
M = (40*fs)/1000.0 # 40ms window 
N = 2048
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
hamMainLobe = mXLinear[1021:1028]
# -------------------------------------------------------------------------------- #
# Process the spectrum upto only 5kHz
maxBinFreq = np.round((5000.0 * N)/fs)
minBinFreq = np.round((100.0 * N)/fs)
widthHalfHamMainLobe = hamMainLobe.size/2

sinThresh = 0.6
t = 0.001
f0et = 5.0
minf0 = 100
maxf0 = 1000
f0t = 0
f0stable = 0
f0 = []

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
    ipmag = 20 * np.log10(sinPeaksMag)
    
    detSinusoids[i, sinPeaksBin] = sinPeaksMag
    
#    if (i <= 2): 
#        plt.figure()
#        #tempMag = tempMag / np.max(tempMag) 
#        freqAxis = (np.arange(np.size(tempMag)) * fs)/N
#        sinPeaksMag = tempMag[sinPeaksBin]
#        sinFreqAxis = sinPeaksBin * (fs/N)
#        plt.plot(freqAxis, tempMag)
#        #plt.stem(ploc, collS, 'c')
#        plt.stem(sinFreqAxis, sinPeaksMag, 'r')
#        #plt.stem(sinFreqAxis, collS[sinPeaksIndx], 'c')
#plt.stem(ploc, tempMag[ploc], 'r')
#plt.stem(ploc, collAM, 'g')
#plt.stem(ploc, collEm, 'y')
#plt.stem(ploc, collS, 'c')
