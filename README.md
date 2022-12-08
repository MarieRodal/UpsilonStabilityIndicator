# The Upsilon Stability Indicator

## Installation

1. Download the Stability Indicator folder from Github
3. Install dependencies. The setup.py file lists the dependencies. Pay particular 
   attention to the installation of the rpy2 library which is used to import the R function auto.arima into python. 
4. Run main.py to see that the installation has worked

For any questions contact:

marie.rodal@uantwerpen.be 


## Description

Provided a window length, tau, the function stability_indicator_parallel in the StabInd module finds the best ARIMA model according to the Baysian Information Criterion (BIC), for each set of overlapping time series intervals. The degree of overlap can be set by the user; by default the time window is moved one datapoint each time. 
From there it computes the Upsilon indicator or the modified Upsilon indicator (see references), depending on the choice of base model. The function found in the module StabInd_AR1 is the original indicator used by Kaiser et al (2020), while the indicator found in the module StabInd_Modified is the indicator used in Rodal et al (2022). StabInd_AR00 uses ARMA(0,0) as the base model when computing the difference in BIC value, and is included for completeness. 

In all cases stationary = False, meaning that auto.arima first determines the differencing order needed to make the time series stationary. 

The function outputs a list object of the form [index, sub sequence , Y, BIC(1,0), BIC(0,0), BIC(p,q), order, persistence, AR1coeff] or [index, sub sequence , Y, BIC(1,0), BIC(p,q), order, persistence, AR1coeff]. Here index is a list of the indicices for all the time intervals, sub sequence is a list of all the considered time intervals, Y is a list of the stability indicator for each time interval, BIC(1,0) , BIC(0,0) and BIC(p,q) are the BIC's of the ARMA(1,0) fit , ARMA(0,0) and the ARMA(p,q) fit (here p and q denote the order of the best fit) and order is a list of tuples of the form (p,d,q), showing the order of the best fitted ARIMA model; d denotes the differencing order. Persistence is the sum of the absolute value of the model coefficients, while AR1coeff is the absolute value of the AR1 coefficient for the ARMA(1,0) fit. See Hyndman et al (2008) for a description of auto.arima. 
If plot_Y is set to True a plot of the stability indicator values overlaid onto the time series is shown. 

The function is written to run in paralell, and a bar should appear on screen to indicate the progression. 

To call the function write

ARIMA_Result = stability_indicator_paralell(time_series, tau, p_max = 10, q_max =10, wind_shift = 1, Plot_Y = True, xlabel = x, ylabel = y, fig_name = Figurename)

Here, time_series is a time series of N data points the form [ [x1 , ... , xN],[t1, ... , tN] ], max_p and max_q are the maximum order of the ARIMA model, and wind_shift describes the degree of overlap between the time intervals; the default of 1 means that the time window is moved one data point each time.

Also included in this folder are the functions used to simulate R-, B- and N-tipping for the 3-box model of the Atlantic Meridional Overturning Circulation (AMOC) in Rodal et al (2022), previusly described in Alkhayoun et al (2019), in addition to the functions used to generate the colored-noise time series. colorplots.py contain the plotting routine used to generate the figures in Rodal et al (2022), while VarAutoCor.py contains functions for the variance and autocorrelation of a time series. 


## Example
```

#import stability indicator

import StabInd_Modified as SIP

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
```


## References 

[1] KAISER A., FARANDA D., KRUMSCHEID S., BELUSIC D., 2020, VERCATUEREN N., "Detecting Regime Transitions of the Nocturnal and Polar Near-SurfaceTemperature Inversion",   American Meteorological Society,  Journal of the Atmospheric Sciences, Vol. 77, p.2921-2940 

[2] HYNDMAN R. J. and KHANDAKAR, Y., 2008, "Automatic Time Series Forecasting: The forecast Package for R", Journal of Statistical Software, Vol. 27, Issue 3. 

[3] RODAL, M., KRUMSCHEID, S., MADAN, G., LaCasce, J and VERCAUTEREN N.,2022, "Dynamical Stability Indicator based on Autoregressive Moving-Average Models: Critical Transitions and the Atlantic Meridional Overturning Circulation," Chaos, Vol.32, Issue 11 
 
[4]  ALKHAYUON H, ASHWIN P, JACKSON LC, QUINN C, WOOD RA. 2019,  "Basin bifurcations, oscillatory instability and rate-induced thresholds for Atlantic meridional overturning circulation in a global oceanic box model", .Proc. R. Soc. A 475:20190051.
