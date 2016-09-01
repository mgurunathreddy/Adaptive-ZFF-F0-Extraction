# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 17:10:10 2015

@author: gurunath
"""
import sys, os
import numpy as np
import peakDetectCorrect as pdc
sys.path.append(os.path.join(os.path.dirname(os.path.realpath('__file__')), '/home/gurunath/coursera/audio_signal_peocessing/sms-tools-master/software/models/'))
import utilFunctions as UF
def f0Detection(stftMx, stftPx, fs, N, t, minf0, maxf0, f0et):
    """
    Fundamental frequency detection of a sound using twm algorithm
    x: input sound; fs: sampling rate; w: analysis window; 
    N: FFT size; t: threshold in negative dB, 
    minf0: minimum f0 frequency in Hz, maxf0: maximim f0 frequency in Hz, 
    f0et: error threshold in the f0 detection (ex: 5),
    returns f0: fundamental frequency
    """
    if (minf0 < 0):                                            # raise exception if minf0 is smaller than 0
        raise ValueError("Minumum fundamental frequency (minf0) smaller than 0")
    if (maxf0 >= 10000):                                       # raise exception if maxf0 is bigger than fs/2
        raise ValueError("Maximum fundamental frequency (maxf0) bigger than 10000Hz")
    f0 = []                                                    # initialize f0 output
    f0t = 0                                                    # initialize f0 track
    f0stable = 0                                               # initialize f0 stable
    for i in range(np.shape(stftMx)[0]):
        mX = stftMx[i, :]
        pX = stftPx[i, :]
        ploc = pdc.peakDetection(mX, t)                           # detect peak locations   
        iploc, ipmag, ipphase = pdc.peakInterp(mX, pX, ploc)      # refine peak values
        ipfreq = fs * iploc/N                                    # convert locations to Hez
        f0t = UF.f0Twm(ipfreq, ipmag, f0et, minf0, maxf0, f0stable)  # find f0
        if ((f0stable==0)&(f0t>0)) \
    				or ((f0stable>0)&(np.abs(f0stable-f0t)<f0stable/5.0)):
            f0stable = f0t                                         # consider a stable f0 if it is close to the previous one
        else:
            f0stable = 0
        f0 = np.append(f0, f0t)                                  # add f0 to output array
    return f0
