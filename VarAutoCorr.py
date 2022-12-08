# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 10:17:07 2021

@author: marie
"""
import numpy as np
from scipy import signal


def variance(time_series, tau):
    variance_ = []
    for i in range(int(len(time_series) - tau + 1)):
        Y = time_series[i: i+tau-1]
        Y = signal.detrend(Y)
        variance = 0
        Y_mean = np.mean(Y)
        for j in range(0,len(Y)):
            variance += (Y[j]-Y_mean)**2
        variance_.append(variance/len(Y))
        
    return variance_

def autocorrelation(time_series, tau):
    autocorr_ = []
    for i in range(int(len(time_series) - tau + 1)):
        Y = time_series[i: i+tau-1]    
        Y = signal.detrend(Y)
        Y_mean = np.mean(Y)
        autocorr = 0
        for j in range(len(Y)-1):
            autocorr += (Y[j]-Y_mean)*(Y[j+1]-Y_mean)
        vari = 0
        vari1 = 0
        for j in range(len(Y)):
            vari += (Y[j]-Y_mean)**2
            vari1 += (Y[j]-Y_mean)**2
        var = np.sqrt(vari*vari1)
        
        autocorr = autocorr/var
        
        autocorr_.append(autocorr)
        
    return autocorr_

