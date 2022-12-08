# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 10:34:04 2021

@author: Marie
"""

#generate data for bifrucation analysis

import numpy as np

import sdeint

import matplotlib.pyplot as plt


#range H: -0.375 to 0.375 
i= 21
H = -0.4  + (i-1)*0.0375  

def global_oceanic_3box_model(X,t):
    #fixed parameters (2 X C02)    
    
    KN = 1.762 *10**6    #gyro strength
    KS = 1.872*10**6
    KIP = 99.977 *10**6
    
    SN_eq = 0.034912
    ST_eq = 0.035435
    SS_eq = 0.034427
    SIP_eq = 0.034668
    SB_eq = 0.034538
    
    Lambda = 1.62*10**7 
    alpha = 0.12
    beta = 790
    gamma = 0.36
    mu = 22*10**(-8)
    
    Y = 100*3.15*10**7
    
    TS = 7.919  #temperatures
    T0 = 3.87  
    
    
    FN = 0.486*10**6  + H*0.1311*10**6 #fluxes
    FT = -0.997 *10**6 + H*0.6961*10**6
    FIP = -0.754*10**6 - H*0.5646*10**6
    
    VN = 0.3683*10**(17) #volumes
    VT = 0.5418*10**(17) 
    VS = 0.6097*10**(17)
    VIP = 1.4860*10**(17)
    VB = 9.9250*10**(17)
    
    
    C = VN*SN_eq+VT*ST_eq+VS*SS_eq+VIP*SIP_eq+VB*SB_eq
    
    
    S0 = 0.035 # reference salinity
    
    SS = (0.034427 - S0)*100.0 # fixed salinities
    SB = (0.034538 - S0)*100.0

    
    SN = X[0];
    ST = X[1]
    
    
    SIP = 100*(C-(VN*SN+VT*ST+VS*SS+VB*SB)/100-S0*(VB+VN+VT+VIP+VS))/VIP
    
    q = Lambda*(alpha*(TS-T0)+beta*(X[0]/100-SS/100))/(1+Lambda*alpha*mu)
    aq = abs(q)
    
    
    z1p = (Y/VN)*(q*(ST/100-SN/100)+KN*(ST/100-SN/100)-FN*S0)
       
    z2p = (Y/VT)*(q*(gamma*SS/100+(1-gamma)*SIP/100-ST/100)+KS*(SS/100-ST/100)+KN*(SN/100-ST/100)-FT*S0)
        
    
    z1n = (Y/VN)*(aq*(SB/100-SN/100)+KN*(ST/100-SN/100)-FN*S0)
        
    z2n = (Y/VT)*(aq*(SN/100-ST/100)+KS*(SS/100-ST/100)+KN*(SN/100-ST/100)-FT*S0)
    
    if q >= 0:
        return np.array([z1p,z2p])
    else:
        return np.array([z1n,z2n])
    
def G(X, t):
    B = np.diag([ 0.004 , 0.004])
    return B


tspan = np.linspace(0,2000,10000) 

#inital salinity values
S_T0 = 0.15
S_N0_upper = 0.035
S_N0_lower = -0.2

if H < 0:
    X0 = [S_N0_lower,S_T0]
elif H > 0:
    X0 = [S_N0_upper,S_T0]

result_SDE = sdeint.itoint(global_oceanic_3box_model, G, X0, tspan)

title = f'2xCO2_H={H}_timeseries_Ntipping.pkl'

time_series =  [result_SDE[:,0],tspan]

data = np.array([tspan,result_SDE[:,0], result_SDE[:,1]])

fig = plt.figure(figsize=(15,5))
#title = 'Atlantic meridional overturning circulation in a global oceanic box model'
#plt.suptitle(title, fontsize=16)

ax1 = fig.add_subplot(1,1,1)
ax1.plot(tspan, result_SDE[:,0], label = '$S_N$', color = 'seagreen')
#ax1.plot(tspan, result_SDE[:,1], label = '$S_T$', color = 'darkcyan')
ax1.set_xlabel("time")
ax1.set_ylabel("Salinity")
ax1.legend()
ax1.grid()

plt.show()
#plt.savefig(title +'.jpg')

