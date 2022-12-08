
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 13:50:06 2020

@author: Marie Rodal
"""

#Uses the modified definition of the Upsilon indicator 


#import libraries
import numpy as np 
import matplotlib.pyplot as plt  #plotting library
import multiprocessing as mp #library for paralell processing
import time #library to check how fast the paralell processing is
import warnings
from tqdm import tqdm #progressbar for parallel processing


        
#create a class to store variables 
class multiResultClass():
    """
    Class to store relevant properties of a time seris to which an arima model is applied. 
    """
    
    def __init__(self, index, x_tau, Y , BIC10, BIC00, BICpq, order, persistence, AR1coeff):
        self.index = index
        self.x_tau = x_tau
        self.Y = Y
        self.BIC10 = BIC10
        self.BIC00 = BIC00
        self.BICpq = BICpq
        self.order = order 
        self.persistence = persistence
        self.AR1coeff = AR1coeff
   
def Y_ARIMA_i(time_series, tau, p_max, q_max, i,wind_shift):
    """
    Function that computes the stability indicator Y for a Gaussian ARMA(p,q) process, for one 
    time window
    time_series: 2D pandas dataframe
    tau: length of time_window
    max_p: maximum value for p. Default should be 5
    max_q: maximum value for q. Default should be 5
    Returns a multiresultClass object
    """
    import rpy2.robjects as robjects
    import rpy2.robjects.packages as rpackages
    from rpy2.robjects import numpy2ri
    from rpy2.rinterface import RRuntimeWarning
    from rpy2.robjects.packages import importr
    from rpy2.robjects.vectors import StrVector
    
    numpy2ri.activate()

    required_packages = ['base', 'forecast'] # list of required R packages 
    if all(rpackages.isinstalled(x) for x in required_packages):
        check_packages = True # True if packages are already installed 
    else:
       check_packages = False # False if packages are not installed 
    if check_packages == False: # Not installed? Then install.
        utils = rpackages.importr('utils')
        utils.chooseCRANmirror(ind=1)
        packages_to_install = [x for x in required_packages if not rpackages.isinstalled(x)]
        if len(packages_to_install) > 0:
            utils.install_packages(StrVector(packages_to_install))
        check_packages = True 
      
    forecast = importr('forecast')

    
    results = multiResultClass(None, None, None, None, None, None, None, None, None)
    #save i
    results.index = i
    
    #slice time series by intervals of length tau, each slice is one column in x_tau 
    results.x_tau = time_series[0][i*wind_shift: i*wind_shift+tau-1]
   
            
    # Estimate ARMA models with all possible combinations of chosen order ranges.
    # Select the best according to BIC
    # do not allow differencing, unit root test is augumented dickey fuller test -> we assume stationarity
    # Auto.Arima requires univariate time series and gaussian white noise (Note: normally distributed errors is equivalent to having normally 
    # distributed observations for any linear time series model) 
    
    warnings.filterwarnings("ignore", category= UserWarning)
    warnings.filterwarnings("ignore", category= RRuntimeWarning)
    
    
    resultsr = robjects.FloatVector(results.x_tau)
    
  
    try:
        model = forecast.auto_arima(resultsr, D = 0, max_p = p_max, max_q = q_max, ic = "bic" ,stepwise=False, allowmean = True, stationary = False, approximation = False, test="adf", allowdrift=False)
    except:
        None

    try:
        results.BICpq = np.array(model.rx('bic'))[0][0]
       
    except:
        None
    
    try:
        results.order = forecast.arimaorder(model)
    except:
        None
      
    #calculate BIC for AR1
    
    try:
        d = results.order[1]  
        order10 = robjects.FloatVector([1,d,0])
        AR1 = forecast.Arima(resultsr, order=order10)
    except:
        None
    
    #calculate BIC for AR0
    
    try:
        d = results.order[1]  
        order00 = robjects.FloatVector([0,d,0])
        AR0 = forecast.Arima(resultsr, order=order00)
    except:
        None

    try:
        results.BIC10 = np.array(AR1.rx('bic'))[0][0]
       
    except: 
        None
   
    #calculate AR1 coeff
    
    try:
        AR1coeff = AR1.rx('coef')[0][0]
        results.AR1coeff = AR1coeff
    except:
        None

    try:
        results.BIC00 = np.array(AR0.rx('bic'))[0][0]
       
    except: 
        None
        
    #calculate stability indicator
    try:
        DeltaBIC10 = results.BICpq-results.BIC10
        DeltaBIC00 = results.BICpq-results.BIC00
        if  DeltaBIC10 >= 0:
                DeltaBIC = np.amin([DeltaBIC10, DeltaBIC00])
        elif DeltaBIC10 < 0:
                DeltaBIC = DeltaBIC00 
        results.Y = 1 - np.exp(-abs(DeltaBIC)/tau)
        
    except:
        None   
        
    #calculate persistence    
    
    try:
        p = results.order[0]
        q = results.order[2]
        coefficients = model.rx('coef')[0]
        coeff = coefficients[:len(coefficients)-1] #remove intercept
        if coeff.size == 0:
            results.persistence = 0
        elif (q==0 and p==0):
            results.persistence = 0 
        else:         
            persistence = np.sum(np.abs(coeff))
            results.persistence = persistence       
    except:
        None
  
    numpy2ri.deactivate()

    return results



def stability_indicator_parallel(time_series, tau, p_max = 5, q_max =5, wind_shift =1, Plot_Y = True, ylabel = 'x', xlabel = 'time', fig_size = (15,5), fig_name='StabilityIndicator'):
    '''
    Computes the stability indicator for all time windows given a certain window lenght. Written to run in parallel
    time_series: time sereies data, numpy array
    max_p: integer, maximum value for p. Default should be 5
    max_q: integer, maximum value for q. Default should be 5
    wind_shift: integer, degree of overlap between time interval. Default is 1
    Plot_Y: Boolean,  If set to true, the function will plot the time series colorcoded accoring to Y
    xlabel: string, label for x-axis
    ylabel: string, label for y-axis
    Returns a list of object attributes
    '''
    from rpy2.rinterface import RRuntimeWarning
    

    start = time.perf_counter()

    pool = mp.Pool(mp.cpu_count())
    
    warnings.filterwarnings("ignore", category= UserWarning)
    warnings.filterwarnings("ignore", category= RRuntimeWarning)

    list_of_result_classes = pool.starmap(Y_ARIMA_i, tqdm([(time_series, tau, p_max, q_max, i, wind_shift) for i in range(int((len(time_series[0]) - tau + 1) / wind_shift))], desc ='Determining Stability Indicator'))
 
    pool.close()

    pool.join()

    end = time.perf_counter()

    
    print(f' Parallel Finished in {round(end - start, 2)} second(s)')

    index = [list_of_result_classes[i].index for i in range(len(list_of_result_classes))]
    subSeq = [list_of_result_classes[i].x_tau for i in range(len(list_of_result_classes))]
    Y = [list_of_result_classes[i].Y for i in range(len(list_of_result_classes))]
    BIC10 = [list_of_result_classes[i].BIC10 for i in range(len(list_of_result_classes))]
    BICpq = [list_of_result_classes[i].BICpq for i in range(len(list_of_result_classes))]
    order = [list_of_result_classes[i].order for i in range(len(list_of_result_classes))]
    persistence = [list_of_result_classes[i].persistence for i in range(len(list_of_result_classes))]
    AR1coeff = [list_of_result_classes[i].AR1coeff for i in range(len(list_of_result_classes))]
    
    Y_data = [time_series,index, subSeq, Y, BIC10, BICpq, order, persistence, AR1coeff]
    
    if Plot_Y == True:
            plot_Y(time_series, Y_data, wind_shift, tau, ylabel, xlabel, fig_size, fig_name)


    return Y_data

def plot_Y(time_series, Y_data, wind_shift, tau, ylabel, xlabel, fig_size, fig_name):
    
    
        subSeq = Y_data[2]
        Y = np.array(Y_data[3])
        points = [subSeq[i][-1] for i in range(len(subSeq))]
       
        time_subSeq = [time_series[1][i*wind_shift: i*wind_shift+tau-1] for i in range(int((len(time_series[0]) - tau + 1) / wind_shift))]
        time_points = [time_subSeq[i][-1] for i in range(len(time_subSeq))]
        
        plt.clf()
        fig = plt.figure(figsize=fig_size)
        
        ax1 = fig.add_subplot(1,1,1)
        ax1.scatter(time_series[1],time_series[0], s= 5, c = 'grey')
        im1 = ax1.scatter(time_points, points, s=6, c = Y, cmap = 'viridis')
        fig.colorbar(im1, ax=ax1)
        
        from matplotlib.colors import LinearSegmentedColormap
        import matplotlib as mpl
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        
        my_colormap = LinearSegmentedColormap.from_list('my colormap' ,  ['#22577a','#38a3a5','#57cc99','#80ed99','#c7f9cc'], N =100)

        plt.clf()
        fig = plt.figure(figsize=fig_size)
        
        points0 = [points[i] if Y[i] is None else None for i in range(0,len(points))]
        
        points1 = [points[i] if (Y[i] is not None and Y[i] <0.2) else None for i in range(0,len(points))]
        
        points2 = [points[i] if (Y[i] is not None and Y[i] > 0.2 and Y[i]<0.4 ) else None for i in range(0,len(points))]
    
        points3 = [points[i] if (Y[i] is not None and Y[i] > 0.4 and Y[i]<0.6 ) else None for i in range(0,len(points))]
   
        points4 = [points[i] if (Y[i] is not None and Y[i] > 0.6 and Y[i]<0.8 ) else None for i in range(0,len(points))]
    
        points5 = [points[i] if (Y[i] is not None and Y[i] > 0.8 ) else None for i in range(0,len(points))]
        
        ax1 = fig.add_subplot(1,1,1)
        ax1.scatter(time_points, points0, s= 5, c = 'grey')
        ax1.scatter(time_points, points1, s= 5, c = '#22577a')
        ax1.scatter(time_points, points2, s= 5, c = '#38a3a5')
        ax1.scatter(time_points, points3, s= 7, c = '#57cc99')        
        ax1.scatter(time_points, points4, s= 7, c = '#80ed99')        
        ax1.scatter(time_points, points5, s= 7, c = '#c7f9cc')
        
        divider = make_axes_locatable(plt.gca())
        ax_cb = divider.new_horizontal(size="2%", pad=0.3)    
        cb1 = mpl.colorbar.ColorbarBase(ax_cb, cmap=my_colormap, orientation='vertical')
        plt.gcf().add_axes(ax_cb)
        cb1.set_label(label = r'$\Upsilon$', rotation = 0, fontsize = 16)
        
        delta = (np.max(points)-np.min(points))/5
        ymax = np.max(points) + delta
        ymin = np.min(points) - delta
        ax1.set_xlabel(f"{xlabel}", fontsize = 16)
        ax1.set_ylabel(f"{ylabel}", fontsize = 16)
        ax1.set_ylim(ymin,ymax)
        ax1.grid()
        
        plt.savefig(fig_name)
