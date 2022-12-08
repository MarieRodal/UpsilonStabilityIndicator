# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 12:35:42 2021

@author: Marie
"""

import matplotlib.pyplot as plt
import numpy as np

#import data

time_series = np.genfromtxt('example_timeseries.csv', delimiter=',')
Y = np.genfromtxt('example_Y.csv', delimiter=',')
BIC10 = np.genfromtxt('example_BIC10.csv', delimiter=',')
BICpq = np.genfromtxt('example_BICpq.csv', delimiter=',')
order = np.genfromtxt('example_order.csv', delimiter=',')
subSeq = np.genfromtxt('example_subSeq.csv', delimiter=',')
Persistence = np.genfromtxt('example_persistence.csv', delimiter=',')


#set parameters from simulation run 

wind_shift = 1
tau = 350

#%%

#create colorplot of q, including q plotted as function of time t
#Use ciscrete color scheme, i.e. N=5 in my_colormap

from matplotlib.colors import LinearSegmentedColormap

my_colormap = LinearSegmentedColormap.from_list('my colormap' ,  ['#22577a','#38a3a5','#57cc99','#80ed99','#c7f9cc'], N =5)


points = [subSeq[i][-1] for i in range(len(subSeq))]

p = [order[i][0] if order[i] is not None else None for i in range(len(order))]  #extract p-values
q = [order[i][2] if order[i] is not None else None for i in range(len(order))]  #extract q-values

      
time_subSeq = [time_series[1][i*wind_shift: i*wind_shift+tau-1] for i in range(int((len(time_series[0]) - tau + 1) / wind_shift))]
time_points = [time_subSeq[i][-1] for i in range(len(time_subSeq))]

fig = plt.figure(figsize=(17,5))

ax1 = fig.add_subplot(1,1,1)
ax1.scatter(time_series[1],time_series[0], s= 4, c = 'grey')
im1 = ax1.scatter(time_points, points, s=0, c = q, cmap = my_colormap)
points1 = [points[i] if (q[i] is not None and q[i]<1.5 ) else None for i in range(0,len(points))]
points2 = [points[i] if (q[i] is not None and q[i] >= 1.5 and q[i]< 2.5 ) else None for i in range(0,len(points))]
points3 = [points[i] if (q[i] is not None and q[i] >= 2.5 and q[i]< 3.5 ) else None for i in range(0,len(points))]   
points4 = [points[i] if (q[i] is not None and q[i] >= 3.5 and q[i]<4.5 ) else None for i in range(0,len(points))]
points5 = [points[i] if (q[i] is not None and q[i] >= 4.5 ) else None for i in range(0,len(points))]
ax1.scatter(time_points, points1, s= 4, c='#22577a')
ax1.scatter(time_points, points2, s=4, c='#38a3a5' )
ax1.scatter(time_points, points3, s=4,c = '#57cc99' )
ax1.scatter(time_points, points4, s =10 ,c ='#80ed99')
ax1.scatter(time_points, points5, s=10,c = '#c7f9cc')
cbar = fig.colorbar(im1, ax=ax1)
cbar.set_label('q', rotation = 0, fontsize = 18)

ax1.set_xlabel('t', fontsize = 18)
ax1.set_ylabel(r'$S_N$', fontsize = 18)
ax1.grid()

ax12=ax1.twinx()
ax12.plot(time_points, q, label = r'$q$', alpha = 0.5)
ax12.spines['right'].set_color('lightskyblue')
ax12.legend(fontsize = 16)


#%%

#create colorplot of p, including q plotted as function of time t
#Use discrete color scheme, i.e. N=5 in my_colormap

fig = plt.figure(figsize=(17,5))

ax2 = fig.add_subplot(1,1,1)
ax2.scatter(time_series[1],time_series[0], s= 4, c = 'grey')
im2 = ax2.scatter(time_points, points, s=0, c = p, cmap = my_colormap)              
points1 = [points[i] if (p[i] is not None and p[i]<1.5 ) else None for i in range(0,len(points))]
points2 = [points[i] if (p[i] is not None and p[i] >= 1.5 and p[i]< 2.5 ) else None for i in range(0,len(points))]
points3 = [points[i] if (p[i] is not None and p[i] >= 2.5 and p[i]< 3.5 ) else None for i in range(0,len(points))]   
points4 = [points[i] if (p[i] is not None and p[i] >= 3.5 and p[i]<4.5 ) else None for i in range(0,len(points))]
points5 = [points[i] if (p[i] is not None and p[i] >= 4.5 ) else None for i in range(0,len(points))]
ax2.scatter(time_points, points1, s= 4, c='#22577a')
ax2.scatter(time_points, points2, s=4, c='#38a3a5' )
ax2.scatter(time_points, points3, s=4,c = '#57cc99' )
ax2.scatter(time_points, points4, s =10 ,c ='#80ed99')
ax2.scatter(time_points, points5, s=10,c = '#c7f9cc')
cbar = fig.colorbar(im2, ax=ax2)
cbar.set_label('p', rotation = 0, fontsize = 18)

ax2.set_xlabel('t', fontsize = 18)
ax2.set_ylabel(r'$S_N$', fontsize = 18)
ax2.grid()

ax22=ax2.twinx()
ax22.plot(time_points, p, label = r'$p$', alpha = 0.5)
ax22.spines['right'].set_color('lightskyblue')
ax22.legend(fontsize = 16)

#%%

fig = plt.figure(figsize=(15,5))

ax1 = fig.add_subplot(1,1,1)
ax1.plot(time_points, Y)
ax1.set_xlabel('t', fontsize = 16)
ax1.set_ylabel('Y', fontsize = 16)
ax1.grid()

#%%

#Create colorplot of deltaBIC values for positive deltaBIC, and include a time series plot of DeltaBIC including negative values
#Use continous color scheme, i.e. N=100 in my_colormap

my_colormap = LinearSegmentedColormap.from_list('my colormap' ,  ['#22577a','#38a3a5','#57cc99','#80ed99','#c7f9cc'], N =100)

points = [subSeq[i][-1] for i in range(len(subSeq))]

#Compute DeltaBIC
DeltaBIC1 = [BIC10[i] - BICpq[i] if (BIC10[i] and BICpq[i]) is not None else None for i in range(0,len(BIC10))]      
time_subSeq = [time_series[1][i*wind_shift: i*wind_shift+tau-1] for i in range(int((len(time_series[0]) - tau + 1) / wind_shift))]
time_points = [time_subSeq[i][-1] for i in range(len(time_subSeq))]

#Find points for which DeltaBIC > 0 
DeltaBIC = [DeltaBIC1[i] if (DeltaBIC1[i] is not None and DeltaBIC1[i]>=0) else None for i in range(len(DeltaBIC1))]

fig = plt.figure(figsize=(17,5))

ax1 = fig.add_subplot(1,1,1)
ax1.scatter(time_series[1],time_series[0], s= 4, c = 'grey')
im1 = ax1.scatter(time_points, points, s=0, c = DeltaBIC, cmap = my_colormap)
points1 = [points[i] if (DeltaBIC[i] is not None  and DeltaBIC[i]<2 ) else None for i in range(0,len(points))] 
points2 = [points[i] if (DeltaBIC[i] is not None and DeltaBIC[i] >= 2 and Y[i]<4 ) else None for i in range(0,len(points))] 
points3 = [points[i] if (DeltaBIC[i] is not None and DeltaBIC[i] >= 4 and Y[i]<6 ) else None for i in range(0,len(points))]   
points4 = [points[i] if (DeltaBIC[i] is not None and DeltaBIC[i] >= 6 and Y[i]<8 ) else None for i in range(0,len(points))]
points5 = [points[i] if (DeltaBIC[i] is not None and DeltaBIC[i] >= 10 ) else None for i in range(0,len(points))]
ax1.scatter(time_points, points1, s= 4, c='#22577a')
ax1.scatter(time_points, points2, s=4, c='#38a3a5' )
ax1.scatter(time_points, points3, s=4,c = '#57cc99' )
ax1.scatter(time_points, points4, s =10 ,c ='#80ed99')
ax1.scatter(time_points, points5, s=10,c = '#c7f9cc')
cbar = fig.colorbar(im1, ax=ax1)
cbar.set_label(r'$\Delta$BIC$>0$', fontsize = 20)

ax1.set_xlabel('t', fontsize = 20)
ax1.set_ylabel(r'$S_N$', fontsize = 20)
ax1.grid()

ax2=ax1.twinx()
ax2.plot(time_points, DeltaBIC1, label = r'$\Delta BIC$', alpha = 0.5)
#ax2.set_ylabel(r'$\Delta BIC$', fontsize= 16)
ax2.spines['right'].set_color('lightskyblue')
ax2.legend(fontsize = 16)



#%%

#Create color plot of Upsilon (Y) and include plot of equilibirum branches in the case of R-, N- and B-tipping in the 
# 3-box model 
#Use continous color scheme, i.e. N=100 in my_colormap

fig = plt.figure(figsize = (17,5))
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel("$t$", fontsize = 18)
ax1.set_ylabel("$S_N$", fontsize = 18)
ax1.grid()


time = time_series[1] #extract time part of time series
X = time_series[1] #extract spatial part of time series

points0 = [points[i] if Y[i] is None else None for i in range(0,len(points))]
points1 = [points[i] if (Y[i] is not None and Y[i] < 0.2) else None for i in range(0,len(points))]
points2 = [points[i] if (Y[i] is not None and Y[i] > 0.2 and Y[i] < 0.4 ) else None for i in range(0,len(points))]
points3 = [points[i] if (Y[i] is not None and Y[i] > 0.4 and Y[i] < 0.6 ) else None for i in range(0,len(points))]   
points4 = [points[i] if (Y[i] is not None and Y[i] > 0.6 and Y[i] < 0.8 ) else None for i in range(0,len(points))]
points5 = [points[i] if (Y[i] is not None and Y[i] > 0.8 ) else None for i in range(0,len(points))]

            
ax1.scatter(time, X, s=5, c = 'grey')
im2 = ax1.scatter( time_points,points, s=0, c = Y, cmap = my_colormap)  
ax1.scatter(time_points, points0, s= 5, c='grey')
ax1.scatter(time_points, points1, s= 5, c='#22577a')
ax1.scatter(time_points, points2, s=5, c='#38a3a5' )
ax1.scatter(time_points, points3, s=5,c = '#57cc99' )
ax1.scatter(time_points, points4, s =10 ,c ='#80ed99')
ax1.scatter(time_points, points5, s=10,c = '#c7f9cc') 


#For R-, B- and N-tipping in the 3-box model, include equilibrium branches

'''
#import equilibirum branch data 
eq_upper = np.genfromtxt('example_equil_upper.csv', delimiter=',')
eq_unstable = np.genfromtxt('example_equil_unstable.csv', delimiter=',')
eq_lower = np.genfromtxt('example_equil_lower.csv', delimiter=',')

ax1.plot(time_points, eq_upper[tau-1:], lw = 3, color ='dimgray')
ax1.plot(time_points, eq_lower[tau-1:], linestyle = '--', lw = 3, color ='dimgray')
ax1.plot(time_points, eq_unstable[tau-1:], lw=3,  color ='dimgray')
'''

cbar = fig.colorbar(im2, ax=ax1)
cbar.set_label(label = r'$\Upsilon$', rotation = 0, fontsize = 18)


#%%

#plot persistence as a function of time

time = time_series[1]  #extract time part of time series

fig = plt.figure(figsize=(15,5))

ax1 = fig.add_subplot(1,1,1)
ax1.plot(time[tau-1:], Persistence)
ax1.set_xlabel('t', fontsize = 18)
ax1.set_ylabel('Persistence', fontsize = 18)
ax1.grid()
