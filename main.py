# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 08:31:03 2022

@author: Marie
"""
#import stability indicator

import StabInd as SIP

#import data. Note that the stability indicator requires time series of the form [X,t], where t is time and x is univariate 

import numpy as np
time_series = np.genfromtxt('example_timeseries.csv', delimiter=',')

#fix simulation parameters
tau = 350 

#set names for figure if Plot==True, as well as axis name
x = 'x_label'
y = 'y_label'

Figurename = 'example_timeseries_plot'

if __name__ == '__main__':   
    ARIMA_Result = SIP.stability_indicator_parallel(time_series, tau, p_max = 10, q_max =10, wind_shift = 1, Plot_Y = True, xlabel = x, ylabel = y, fig_name = Figurename)
 
#Note that StabInd does not autimatically save the generated data, so that has to be implemented manually    