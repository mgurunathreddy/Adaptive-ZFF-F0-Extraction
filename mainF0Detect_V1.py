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
fileName = '1.7.happy-15'
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
frame2Time = (np.arange(np.size(sshVal)) * float(H))/fs
begVoicTime = (begVoic*H)/fs
endVoicTime = (endVoic*H)/fs

tempTime = np.arange(np.size(sshVal))
tTime = (tempTime * H)/fs

plt.figure()
plt.plot(timeAxis, x)
plt.plot(tTime, sshVal, 'r')
plt.stem(begVoicTime, begVoicIndx, 'g')
plt.stem(endVoicTime, begVoicIndx, 'r') 
plt.title('Salience based voiced and unvoiced classification')

#matio.savemat('/home/gurunath/MS_work/papers/Indicon-2015/happy_synth/code_to _generate_figures/orig/synthtime.mat', mdict={'synthtime':timeAxis})
#matio.savemat('/home/gurunath/MS_work/papers/Indicon-2015/happy_synth/code_to _generate_figures/orig/synx.mat', mdict={'synx':x})
#matio.savemat('/home/gurunath/MS_work/papers/Indicon-2015/happy_synth/code_to _generate_figures/orig/synframe2time.mat', mdict={'synframe2time':frame2Time})
#matio.savemat('/home/gurunath/MS_work/papers/Indicon-2015/happy_synth/code_to _generate_figures/orig/synsshval.mat', mdict={'synsshval':sshVal})
#matio.savemat('/home/gurunath/MS_work/papers/Indicon-2015/happy_synth/code_to _generate_figures/orig/synbegvoictime.mat', mdict={'synbegvoictime':begVoicTime})
#matio.savemat('/home/gurunath/MS_work/papers/Indicon-2015/happy_synth/code_to _generate_figures/orig/synendvoictime.mat', mdict={'synendvoictime':endVoicTime})
    
# ------------------------------------------------------------------------------
# frame to sample number conversion
timeBegVoicSamp = np.floor(begVoic * H)
timeEndVoicSamp = np.floor(endVoic * H)
# ------------------------------------------------------------------------------
sampOnsets = np.sort(np.concatenate((timeBegVoicSamp, timeEndVoicSamp))) # smapleOnsets are in sample number
#matio.savemat('/home/gurunath/MS_work/papers/Indicon-2015/happy_synth/synthAngOnsets.mat', mdict={'synthAngOnsets':sampOnsets})

zffs = getZFFS.getZFFS(x, fs, H, sampOnsets, resntFreq) # ZFFS for each voiced region

# Once the ZFFS is found, get the instants of zero crossings, strength of excitation and f0 at GCI's
zf1 = np.copy(zffs)
zffs[zffs>0] = 1 # To find zero crossings, place value 1 at all locations zffs>0
zffs[zffs<0] = -1 # Place value -1 at all locations zffs<0
gci = np.where(np.diff(zffs) == 2)[0] # Take difference and look for the positions of value 2
es = np.abs(zf1[gci+1]-zf1[gci-1]) # Positions of value 2 are the instants of zero crossings
T0 = np.diff(gci) # Finding period interms of sample number 
T0 = T0/float(fs) # Period in seconds
f0 = 1.0/T0 # Frequency in Hz
f0 = np.append(f0, f0[-1], f0[-1]) # Filling holes created by two difference operation

# Smoothing the melody to remove high frequency contents
f0Smooth = smooth.smooth(f0, 10, 'flat')

# ------------------------------------------------------------------------------- #
# Remove unvoiced portions from F0 contour
tempF0 = np.zeros(lenSig)
tempF0[gci] = f0Smooth

# Get F0 at only voiced regions
voicPitch = np.zeros(np.shape(x)[0])

if(np.shape(timeBegVoicSamp)[0] == np.shape(timeEndVoicSamp)[0]): # Check if number of instants are same
    for i in range(np.shape(timeBegVoicSamp)[0]):
        begInst = timeBegVoicSamp[i] # Get beg instant
        endInst = timeEndVoicSamp[i] # Get end instant
        voicPitch[begInst:endInst] = tempF0[begInst:endInst] # Place voiced portion f0

plt.figure()
plt.subplot(211)
plt.plot(timeAxis, x)
plt.xlim([np.min(timeAxis), np.max(timeAxis)])
plt.subplot(212)
plt.plot(timeAxis, voicPitch, 'r*')
plt.xlim([np.min(timeAxis), np.max(timeAxis)])
plt.ylim([0, 400])
plt.title('F0 extraction by adaptive ZFF')

#matio.savemat('/home/gurunath/MS_work/papers/Indicon-2015/happy_synth/code_to _generate_figures/orig/synpitchbefdouble.mat', mdict={'synpitchbefdouble':voicPitch})


#plt.subplot(212)
#plt.plot(timeAxis, zf1)

#plt.figure()
#plt.plot(timeAxis, x)
#plt.plot(frame2Time, sshVal/np.max(sshVal), 'r')
 
# ----------------------------------------------------------------------------- # 
# Performing ZFF double filtering part
# ----------------------------------------------------------------------------- #
# FFT parameter setting
M = (20*fs)/1000.0 # 20ms window 
N = 512
H = M/4
window = 'hamming'
w = np.hamming(M)

zffsDoub = doubleFilt.filterZffs(zf1, fs, w, N, H, sampOnsets) # Get the F0 of ZFFS signal     
#plt.plot(zfsF0)

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Filter the zero frequnecy filtered signal one more time with the help of resonant frequency contour obtained from ZFFS signal
#zffsDoub = getZFFS.getZFFS(zf1, fs, H, sampOnsets, resntFreq) # ZFFS for each voiced region

# Once the ZFFS is found, get the instants of zero crossings, strength of excitation and f0 at GCI's
zf1Doub = np.copy(zffs)
zffsDoub[zffsDoub>0] = 1 # To find zero crossings, place value 1 at all locations zffs>0
zffsDoub[zffsDoub<0] = -1 # Place value -1 at all locations zffs<0
gci = np.where(np.diff(zffsDoub) == 2)[0] # Take difference and look for the positions of value 2
es = np.abs(zf1Doub[gci+1]-zf1Doub[gci-1]) # Positions of value 2 are the instants of zero crossings
T0 = np.diff(gci) # Finding period interms of sample number 
T0 = T0/float(fs) # Period in seconds
f0 = 1.0/T0 # Frequency in Hz
f0 = np.append(f0, f0[-1], f0[-1]) # Filling holes created by two difference operation

# Smoothing the melody to remove high frequency contents
f0Smooth = smooth.smooth(f0, 10, 'flat')

# ------------------------------------------------------------------------------- #
# Remove unvoiced portions from F0 contour
tempF0 = np.zeros(lenSig)
tempF0[gci] = f0Smooth

# Get F0 at only voiced regions
voicPitch = np.zeros(np.shape(x)[0])

if(np.shape(timeBegVoicSamp)[0] == np.shape(timeEndVoicSamp)[0]): # Check if number of instants are same
    for i in range(np.shape(timeBegVoicSamp)[0]):
        begInst = timeBegVoicSamp[i] # Get beg instant
        endInst = timeEndVoicSamp[i] # Get end instant
        voicPitch[begInst:endInst] = tempF0[begInst:endInst] # Place voiced portion f0

plt.figure()
plt.subplot(211)
plt.plot(timeAxis, x)
plt.xlim([np.min(timeAxis), np.max(timeAxis)])
plt.subplot(212)
plt.plot(timeAxis, voicPitch, 'r*')
plt.xlim([np.min(timeAxis), np.max(timeAxis)])
plt.ylim([0, 400])
plt.title('F0 by double filtering')

#matio.savemat('/home/gurunath/MS_work/papers/Indicon-2015/happy_synth/code_to _generate_figures/orig/synvopitchaftdoub.mat', mdict={'synvoicpitchaftdoub':voicPitch})


#-------------------------------------------------------  Original ZFF ---------------------------------- #

windSize = 3 # Choose ZFF mean subtraction window emperically as 4ms for v/uv classification
[zffs, gci, vuvEpssp1, vuvF0sp1] = zff.zff(x, fs, windSize) # Compute the zff of x
zf1 = np.copy(zffs)
zffs[zffs>0] = 1 # To find zero crossings, place value 1 at all locations zffs>0
zffs[zffs<0] = -1 # Place value -1 at all locations zffs<0
gci = np.where(np.diff(zffs) == 2)[0] # Take difference and look for the positions of value 2
es = np.abs(zf1[gci+1]-zf1[gci-1]) # Positions of value 2 are the instants of zero crossings
T0 = np.diff(gci) # Finding period interms of sample number 
T0 = T0/float(fs) # Period in seconds
f0 = 1.0/T0 # Frequency in Hz
f0 = np.append(f0, f0[-1], f0[-1]) # Filling holes created by two difference operation
f0Smooth = smooth.smooth(f0, 10, 'flat')

# Remove unvoiced portions from F0 contour
tempF0 = np.zeros(lenSig)
tempF0[gci] = f0Smooth

# Get F0 at only voiced regions
voicPitch = np.zeros(np.shape(x)[0])

if(np.shape(timeBegVoicSamp)[0] == np.shape(timeEndVoicSamp)[0]): # Check if number of instants are same
    for i in range(np.shape(timeBegVoicSamp)[0]):
        begInst = timeBegVoicSamp[i] # Get beg instant
        endInst = timeEndVoicSamp[i] # Get end instant
        voicPitch[begInst:endInst] = tempF0[begInst:endInst] # Place voiced portion f0

plt.figure()
plt.subplot(211)
plt.plot(timeAxis, x)
plt.xlim([np.min(timeAxis), np.max(timeAxis)])
plt.subplot(212)
plt.plot(timeAxis, voicPitch, 'r*')
plt.xlim([np.min(timeAxis), np.max(timeAxis)])
plt.ylim([0, 400])
plt.title('Original ZFF')
#matio.savemat('/home/gurunath/MS_work/papers/Indicon-2015/happy_synth/code_to _generate_figures/orig/synorigzffpit.mat', mdict={'synorigzffpit':voicPitch})

# -------------------------------------------------------------------------------------------------------- #

#zfsF0Doub = doubleFilt.filterZffs(zf1Doub, fs, w, N, H) # Get the F0 of ZFFS signal     


#plt.figure()
#plt.plot(x)
#plt.plot(zf1Doub)

#t = 0.001
#import peakDetectCorrect as pdc
#maxBinFreq = np.round((5000.0 * N)/fs)
#tempMag = mxLinear[512, :maxBinFreq]
#ploc = pdc.peakDetection(tempMag, t)          
#plt.figure()
#plt.plot(tempMag)
#plt.stem(ploc, tempMag[ploc], 'r')
#
#plt.figure()
#tempMag = detSinusoids[512, :maxBinFreq]
#ploc = pdc.peakDetection(tempMag, t)          
#plt.plot(tempMag)
#plt.stem(ploc, tempMag[ploc], 'r')




