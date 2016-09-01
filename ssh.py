# -*- coding: utf-8 -*-
"""
Created on Fri May 22 17:01:23 2015
Description: Finds the summmation of spectral partials which are detected by similarity measure
Inputs: detected sinusoids, sampling frequency, frame hop in samples, FFT samples 
Outputs: Approximate pitch contour equal to the number of frames
@author: Gurunath Reddy M

"""

import numpy as np
import matplotlib.pyplot as plt

def sumSpectHarm(detSinusoids, fs, H, N):

    minF0 = 50  # Min. F0 considered for SSH
    print 'Adjust the min. and max. F0 values'
    maxF0 = 500 # Max. F0 considered for SSH
    
    minF0Bin = np.floor((minF0 * N) /fs) # Convert frequency to bin
    maxF0Bin = np.floor((maxF0 * N) /fs)
    minF0Bin = int(minF0Bin)
    maxF0Bin = int(maxF0Bin)
        
    numFrames = np.shape(detSinusoids)[0] 
    
    ssh = np.zeros(maxF0Bin) # Remove first few bins which are less than minF0. Lower bins are kept for simplicity
    sshVal = np.zeros(numFrames)
    sshTt = np.zeros([numFrames, maxF0Bin])
    F0s = np.zeros(numFrames)
    
    #    nHar = 10
    
    for frm in range(numFrames):
        for freq in range(minF0Bin, maxF0Bin):
            #ssh[freq] = detSinusoids[frm, freq] + np.sum( detSinusoids[frm, freq*np.arange(2, nHar)] - detSinusoids[frm, np.floor(freq*(np.arange(2, nHar)-0.5))])                    
            ssh[freq] = (detSinusoids[frm, freq] + detSinusoids[frm, 2*freq] + detSinusoids[frm, 3*freq] + detSinusoids[frm, 4*freq] + detSinusoids[frm, 5*freq]) - (detSinusoids[frm, np.int(1.5*freq)] + detSinusoids[frm, np.int(2.5*freq)] + detSinusoids[frm, np.int(3.5*freq)] + detSinusoids[frm, np.int(4.5*freq)] + detSinusoids[frm, np.int(5.5*freq)])                                                      
            #ssh[freq] = (detSinusoids[frm, freq] + detSinusoids[frm, 2*freq] + detSinusoids[frm, 3*freq] + detSinusoids[frm, 4*freq] + detSinusoids[frm, 5*freq])# - (detSinusoids[frm, np.int(1.5*freq)] + detSinusoids[frm, np.int(2.5*freq)] + detSinusoids[frm, np.int(3.5*freq)] + detSinusoids[frm, np.int(4.5*freq)] + detSinusoids[frm, np.int(5.5*freq)])                                                      
        sshTt[frm, :] = ssh
        maxi = np.max(ssh)
        posi = np.where(ssh == maxi)[0]
        posi = posi[0]
        F0s[frm] = posi
        sshVal[frm] = ssh[posi]
    
    pitch = F0s * (fs/N)

    #timeAxis = (np.arange(pitch.size)*H)/fs        
        
    #melodyPitch = np.loadtxt('../wav/6_syn.txt')
    #plt.plot(melodyPitch[:, 0], melodyPitch[:, 1])
    
    #plt.plot(timeAxis, pitch)
    
    return pitch, sshVal

#frame = np.int((1.521*fs)/H)
#for i in range(frame, frame+3):
#    plt.figure()
#    plt.plot(detSinusoids[i, :])
#    
#maxPeakAmp = np.max(detSinusoids, 0)    
#maxPeakBin = np.where(detSinusoids == maxPeakAmp)[0]
    
    