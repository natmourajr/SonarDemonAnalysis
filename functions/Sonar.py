""" 
  This file contents implementation of some sonar functions
"""
import numpy as np
import scipy.signal

def DemonAnalysis(data, fs, decimate_ratio1=25, decimate_ratio2=25, n_fft_pts=1024, overlap=0.5):
    """
        Demon Analysis: def
        decimation process in 2 steps
        decimate_ratio1 = first step decimation ratio
        decimate_ratio2 = second step decimation ratio
        n_fft_pts = fft number of points
        ovrlap = time to compute overlap
    """
    
    # demodulation
    data_abs = np.abs(data)
    
    # first step of decimation
    data_decimate1 = scipy.signal.decimate(data_abs,decimate_ratio1,ftype='fir',n=10)
    fs_decimate1 = float(fs)/float(decimate_ratio1)
    
    # second step of decimation
    data_decimate2 = scipy.signal.decimate(data_decimate1,decimate_ratio2,ftype='fir',n=10)
    fs_decimate2 = float(decimate_ratio1)/float(decimate_ratio2)
    
    n_pts_ovrlap = np.floor(n_fft_pts-fs_decimate2*overlap)
    
    f, t, Sxx = scipy.signal.spectrogram((data_decimate2.T
                                          -np.mean(data_decimate2)),
                                         fs=fs_decimate2,
                                         window=scipy.signal.hanning(n_fft_pts),
                                         nperseg=n_fft_pts,
                                         noverlap=n_pts_ovrlap, nfft=n_fft_pts)
    data_demon =  Sxx
    data_demon_abs = np.abs(data_demon)
    
    
    return [data_demon_abs,f,t]