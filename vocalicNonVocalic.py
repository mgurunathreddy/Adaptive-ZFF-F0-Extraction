# -*- coding: utf-8 -*-
"""
Created on Fri May 22 12:26:59 2015
Description: Used to detect the vocalic and non-vocalic regions of the percussion supressed music signal
Inputs: 
Outputs: 
@author: Gurunath Reddy M

"""

import numpy as np
from scipy.signal import _savitzky_golay as savgol_filter
import matplotlib.pyplot as plt
import zff as zff
import waveio as io
import voiceRegMod as speechVuv
import detSinusouds as detSin
import ssh as ssh
import getMusicZff as getZFFS
import smooth as smooth
import getSSHVUV as sshVUV
import dftStft as stft
import zffDoubleFilter as doubleFilt
import plotSpect as pltSpec
import scipy.io as matio

#Read the wave file for processing
fileName = '2harmonicComp'
[fs, x] = io.wavread('../wav/'+fileName+'.wav')
x = np.array(x, np.float64) # Datatype conversion is must for python to compatable with matlab
x = x/(1.01*np.max(np.abs(x)));
lenSig = x.size
timeAxis = np.arange(lenSig)/float(fs)

#TODO: Find the representative F0 for each voiced region and then perform ZFF filteri
detSinusoids, mxLinear, sinPeaksMag, sinPeaksBin  = detSin.getSinusoids(x, fs)

H = 80.0 # 5ms = 80 samples hop size for 16kHz
N = 1024 # FFT size

# Find summation of spectral harmonics and get pitch contour which is used as resonant frequency in ZFF filtering 
resntFreq, sshVal = ssh.sumSpectHarm(detSinusoids, fs, H, N)
sshVal = np.power(sshVal, 2)
sshVal = sshVal/np.max(sshVal)

# Finding V/UV regions based on the power of the sshVal
begVoic, endVoic = sshVUV.sshVUV(sshVal, H, fs, N)

begVoicIndx = np.ones(begVoic.size)
frame2Time = (np.arange(np.size(sshVal)) * H)/float(fs)
begVoicTime = (begVoic*H)/float(fs)
endVoicTime = (endVoic*H)/float(fs)

tempTime = np.arange(np.size(sshVal))
tTime = (tempTime * H)/fs

plt.figure()
plt.plot(timeAxis, x)
plt.plot(tTime, sshVal, 'r')
plt.stem(begVoicTime, begVoicIndx, 'g')
plt.stem(endVoicTime, begVoicIndx, 'r') 
plt.title('Salience based voiced and unvoiced classification')