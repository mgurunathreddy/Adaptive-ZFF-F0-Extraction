# Adaptive-ZFF-F0-Extraction

Authors: Gurunath Reddy M, and K. Sreenivasa Rao

Emotion like joyous happy and laughter sequence is characterised by dynamically varying source and system features from one voiced region to other compared to neutral speech. Since, ZFF [1] depends on mean subtraction filter (MSF) for accurate source feature (SF: F0 and strength of excitation) extraction, a single MSF is not sufficient to obtain the source parameters accurately. In order to obtain the accurate SF's, the trend in the output of the Zeor Frequency Resonator (ZFR) of each voiced segment in case of happy emotive speech and voiced laughter call (voiced segments in the laughter sequence) in the laughter sequence should be removed adaptively with the MSF length corresponding to the average pitch period of the voiced segment under consideration. Hence, a spectro-temporal based method is proposed to extract the SF’s from the happy emotive speech and laughter signal.

The main program is "mainF0Detect_V1.py"

From the Ubuntu termianl, run the source code as 

python mainF0Detect_V1.py

Note: In line number 28 of "mainF0Detect_V1.py", specify the file name of the audio file for which source features needs to be extracted without .wav extension as shown below 

fileName = 'wavefile name'

Dependencies: pyhthon 2.7, numpy, scipy, matplotlib


[1] K. S. R. Murty and B. Yegnanarayana, “Epoch extraction from speech signals,” IEEE Transactions on Audio, Speech, and
Language Processing, vol. 16, no. 8, pp. 1602–1613, 2008.




