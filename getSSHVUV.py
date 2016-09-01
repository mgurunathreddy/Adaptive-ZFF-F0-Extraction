# -*- coding: utf-8 -*-
"""
Created on Thu May 28 11:22:29 2015
Description:
Inputs:
Outputs:
@author: Gurunath Reddy M

"""

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 19:11:06 2015
Find the voiced regions in the wave file and then modify pitch of each voiced region
@author: gurunath
"""

import numpy as np
import matplotlib.pyplot as plt

def sshVUV(sshVal, H, fs, N):

    # TODO: Threshold may not work all the time. Come up with some different solution
    print('V/UV threshold is hardcoded in getSHHVUV. Change in future')
    thresh = 0.01
    
    sshVUV = sshVal > thresh # Keep zeros at uv and ones at v regions
    
    sshVUV = np.array(sshVUV, dtype=int) # converting bool to logical to take diff
    
    diffShhVal = np.diff(sshVUV) # Take the difference of signal
    diffShhVal = np.append(diffShhVal, diffShhVal[-1]) # Adjust the signal length
    
    begVoic = np.where(diffShhVal == 1)[0] # Values = 1 in diffSgfiltEpssp1 corresponds to the beg. instant of voiced segment 
    
    endVoic = np.where(diffShhVal == -1)[0] # Values = -1 corresponds to UV regions
    
    begVoicIndx = np.ones(begVoic.size)
    
    
#    begVoic = (begVoic*H)/fs
#    endVoic = (endVoic*H)/fs
#    
#    plt.plot(timeAxis, x)
#    plt.plot(frame2Time, sshVal, 'r')
#    plt.stem(begVoic, begVoicIndx, 'g')
#    plt.stem(endVoic, begVoicIndx, 'r') 
    
    ## Get the time instants of beg and end of voiced segments
    ##    timeBegVoic = timeSeg[begVoic] 
    ##    timeEndVoic = timeSeg[endVoic]
    #
    ## Get the sample index of the begin and end of each voiced segments
    #timeBegVoicSamp = gclocssp1[begVoic]
    #timeEndVoicSamp = gclocssp1[endVoic]
    #
    ## Eliminate all false voiced segments based on minimum duration criterian
    #
    #durVoicSeg = timeEndVoicSamp - timeBegVoicSamp
    #
    #minVoicDur = 40*16 # Minimum voice duration must be around 40ms 
    #indFalsVoic = np.where(durVoicSeg <  minVoicDur)[0]   
    #timeBegVoicSamp = np.delete(timeBegVoicSamp, indFalsVoic)
    #timeEndVoicSamp = np.delete(timeEndVoicSamp, indFalsVoic) 
    #
    ## Find the unvoiced indicies by using voiced markers
    #timeBegUVSamp = np.array([])
    #timeEndUVSamp = np.array([])
    #for i in range(timeEndVoicSamp.size-1):
    #    timeBegUVSamp = np.append(timeBegUVSamp, timeEndVoicSamp[i])
    #    timeEndUVSamp = np.append(timeEndUVSamp, timeBegVoicSamp[i+1])
    #
    #timeBegUVSamp = timeBegUVSamp + 1
    #timeEndUVSamp = timeEndUVSamp - 1    
        
    return begVoic, endVoic