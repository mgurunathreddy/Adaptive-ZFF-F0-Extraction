# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 19:11:06 2015
Find the voiced regions in the wave file and then modify pitch of each voiced region
@author: gurunath
"""

import numpy as np
from scipy.signal import _savitzky_golay as savgol_filter
import matplotlib.pyplot as plt
import zff as zff
import waveio as io

def speechVUV(x, epssp1, gclocssp1, fs, lenSig):

    # Read the wave file for processing
#    [fs, x] = io.wavread('../wav/6_syn.wav')
#    x = np.array(x, np.float64) # Datatype conversion is must for python to compatable with matlab
#    x = x/(1.01*np.max(np.abs(x)));
#    lenSig = x.size
#    winLen = 6.0 # Time being, let window lenth to be 4ms. In future compute it by automatic way
#    zf1, gclocssp1, epssp1, f0sp1 = zff.zff(x, fs, winLen)  # Perform ZFF of the signal
    
    epssp1 = epssp1/np.max(epssp1) # Normalize the epoch strengh of ZFF
    #sgfiltEpssp1 = savgol_filter.savgol_filter(epssp1, 21, 2) # Low pass filter the signal to get the envelope of epoch strengh signal
    sgfiltEpssp1 = np.copy(epssp1)
    PlotsgfiltEpssp1 = np.copy(sgfiltEpssp1) # copied for visualization
    
    # TODO: Threshold may not work all the time. Come up with some different solution
    #thresh = np.mean(sgfiltEpssp1) - np.std(sgfiltEpssp1) # Decide thresholding for V/UV classification
    print('V/UV threshold is hardcoded in zffMusicVUV change in future')
    thresh = 0.04
    epstrsp1 = np.zeros(lenSig)
    epstrsp1[gclocssp1] = epssp1 # Place epoch strength at glottal closure instants
    
    sgfiltEpssp1 = sgfiltEpssp1 > thresh # Keep zeros at uv and ones at v regions
    
    #    vgclocssp1 = gclocssp1[sgfiltEpssp1] # Only take the voiced segments indicies
    
    #    timeSeg = gclocssp1/float(fs) # Get the time scale from sample values
    sgfiltEpssp1 = np.array(sgfiltEpssp1, dtype=int) # converting bool to logical to take diff
    
    diffSgfiltEpssp1 = np.diff(sgfiltEpssp1) # Take the difference of signal
    diffSgfiltEpssp1 = np.append(diffSgfiltEpssp1, diffSgfiltEpssp1[-1]) # Adjust the signal length
    
    begVoic = np.where(diffSgfiltEpssp1 == 1)[0] # Values = 1 in diffSgfiltEpssp1 corresponds to the beg. instant of voiced segment 
    
    endVoic = np.where(diffSgfiltEpssp1 == -1)[0] # Values = -1 corresponds to UV regions
    
    # Get the time instants of beg and end of voiced segments
    #    timeBegVoic = timeSeg[begVoic] 
    #    timeEndVoic = timeSeg[endVoic]
    
    # Get the sample index of the begin and end of each voiced segments
    timeBegVoicSamp = gclocssp1[begVoic]
    timeEndVoicSamp = gclocssp1[endVoic]
    
    # Eliminate all false voiced segments based on minimum duration criterian
    
    durVoicSeg = timeEndVoicSamp - timeBegVoicSamp
    
    
    minVoicDur = 40*16 # Minimum voice duration must be around 40ms 
    indFalsVoic = np.where(durVoicSeg <  minVoicDur)[0]   
    timeBegVoicSamp = np.delete(timeBegVoicSamp, indFalsVoic)
    timeEndVoicSamp = np.delete(timeEndVoicSamp, indFalsVoic) 
    
    # Find the unvoiced indicies by using voiced markers
    
    timeBegUVSamp = np.array([])
    timeEndUVSamp = np.array([])
    for i in range(timeEndVoicSamp.size-1):
        timeBegUVSamp = np.append(timeBegUVSamp, timeEndVoicSamp[i])
        timeEndUVSamp = np.append(timeEndUVSamp, timeBegVoicSamp[i+1])
    
    timeBegUVSamp = timeBegUVSamp + 1
    timeEndUVSamp = timeEndUVSamp - 1    
    # --------------------------------------------------------------------------------- #
#    
#    voicMark = np.ones(timeBegVoicSamp.size)
#    
#    #plt.figure()
#    #plt.plot(gclocssp1, PlotsgfiltEpssp1)
#    #threshold = thresh * np.ones(gclocssp1.size)
#    #plt.plot(gclocssp1, threshold)
#    #plt.stem(timeBegVoicSamp, voicMark, 'g')
#    #plt.stem(timeEndVoicSamp, voicMark, 'r')
#    
#    
#    plt.figure()
#    plt.subplot(211)
#    plt.plot(x)
#    #plt.plot(gclocssp1, sgfiltEpssp1)
#    plt.stem(timeBegVoicSamp, voicMark, 'g')
#    plt.stem(timeEndVoicSamp, voicMark, 'r')
#    
#    threshold = thresh * np.ones(gclocssp1.size)
#    plt.subplot(212)
#    plt.plot(gclocssp1, epssp1)
#    plt.stem(timeBegVoicSamp, voicMark, 'g')
#    plt.stem(timeEndVoicSamp, voicMark, 'r')
#    plt.plot(gclocssp1, threshold)    
#        
#    #    uvMark = np.ones(timeBegUVSamp.size)
#    #    #plt.figure()
#    #    #plt.plot(x)
#    #    plt.stem(timeBegUVSamp, uvMark, 'y')
#    #    plt.stem(timeEndUVSamp, uvMark, 'c')
#    #    
        # -------------------------------------------------------------------------------- #
        
    return timeBegVoicSamp, timeEndVoicSamp, timeBegUVSamp, timeEndUVSamp, PlotsgfiltEpssp1
