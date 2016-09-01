# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 11:57:20 2015

Description:
        Finds the zffs of each melody segment and returns the concatinated melody 
    
Inputs: 
        x - music signal
        fs - sampling frequency of music signal
        sampOnsets - vector of onsets followed by offsets of voiced note segments
        f0MedFilt - median filtered TWM output for signal x
        
Outputs: 
        addZsp - the combined zffs of all voiced notes
        
@author: gurunath
"""
import os, sys
import numpy as np
import zff as zff
from scipy.signal import get_window
import dftStft as stft
sys.path.append(os.path.join(os.path.dirname(os.path.realpath('__file__')), '/home/gurunath/coursera/audio_signal_processing/sms-tools-master/software/models/utilFunctions_C'))
import utilFunctions_C as UF_C
import peakDetectCorrect as pdc


def getZFFS(x, fs, H, sampOnsets, f0MedFilt):
    
    #f0MedFilt = resntFreq
    addZsp = np.zeros(np.shape(x)[0]) # Add ZFFS of all notes into a single vector
    
    # Process first segment separately
    begSeg = sampOnsets[0]
    begF0 = 1
    
    collectMean = np.array([]) # Collect mean F0 values of all notes
    collectBegF0 = np.array([]) # Collect all beg note f0 frame of TWM
    collectEndF0 = np.array([]) # collect all end note f0 frame of TWM
    collectBegF0 = np.append(collectBegF0, begF0) # Place first frame no. of TWM f0
    
    
    for i in range(np.shape(sampOnsets)[0]-1):
        # To get the timing information, place the note exactly at the instants of original signal of a vector of zeros of length = input music sinal
        # Now it is difficult to store the ZFFS for each segment, hence it is note stored
        wavSeg = np.zeros(np.shape(x)[0])
    
        endSeg = sampOnsets[i+1] # Get sample index of each note
        endF0 = np.floor(endSeg/H) # Find the end frame index of TWM f0 
        collectEndF0 = np.append(collectEndF0, endF0) # Collect all end note frame indicies of TWM f0 for debugging
        # TODO: In future try to take at the stable portions of note 
        meanF0 = np.mean(f0MedFilt[begF0:endF0]) # Find the mean f0 of each note
        collectMean = np.append(collectMean, meanF0) # Again collect the mean f0 for cross validation of results
        
        # Check if meanF0 is 0 to avoid divide by zero error
        if (meanF0 <= 0):
            meanF0 = 250
            
        winLength = (1.0/meanF0)*1000 # Convert window length to ms
        print winLength
        #appZeros = np.zeros(np.floor((winLength*fs)/1000)) # Sudden truncations of signal may create problems, hence append zeros
        #Get the segment of note in samples, appended with zeros on either sides
        #wavSeg = np.concatenate((appZeros, x[begSeg:endSeg])) 
        #wavSeg = np.concatenate((wavSeg, appZeros))
        wavSeg[begSeg:endSeg] = x[begSeg:endSeg]
        zsp, gclocssp, epssp, f0sp = zff.zff(wavSeg, fs, winLength) # Find the ZFF of current note
        print 'Double filtering is not performed'
    
        # Here we need to add instead of concatenation     
        addZsp = addZsp+zsp        
        begSeg = endSeg+1 # Beg note index of next note is the end index + 1 of current note
        begF0  = np.floor(begSeg/H) # Find the f0 frame index for TWM f0
        collectBegF0 = np.append(collectBegF0, endF0) # Collect the frame number
        
    return addZsp
        