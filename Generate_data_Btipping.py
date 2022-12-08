# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 12:18:49 2021

@author: Marie
"""


#generate data for phase plane analysis, time dependent H non-zero Tris and infinite T_pert

import numpy as np

import sdeint

import matplotlib.pyplot as plt


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
    
    def H(t, H0 = 0, Hpert = 0.5, Tris = 1000):
        T0 = 0
        
        rris    = (Hpert - H0)/Tris;
        
        if t <= T0:
            output = H0
        elif (t>T0 and t<=(Tris+T0)):
            output = (rris * t + (H0 - (rris * T0))) # linear 
        elif (t>(Tris+T0)):
            output = Hpert;
                 
        return output
    
    FN = 0.486*10**6  + H(t)*0.1311*10**6 #fluxes
    FT = -0.997 *10**6 + H(t)*0.6961*10**6
    FIP = -0.754*10**6 - H(t)*0.5646*10**6
    
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
S_N0 = 0.035

X0 = [S_N0,S_T0]

result_SDE = sdeint.itoint(global_oceanic_3box_model, G, X0, tspan)

time_series =  [result_SDE[:,0],tspan]

data = np.array([tspan,result_SDE[:,0], result_SDE[:,1]])


fig = plt.figure(figsize=(15,5))
title = 'Time series'
plt.suptitle(title, fontsize=16)

ax1 = fig.add_subplot(1,1,1)
#ax1.scatter(S_N0, S_T0)
ax1.plot(tspan,result_SDE[:,0], color = 'seagreen')
ax1.set_xlabel("$t$")
ax1.set_ylabel("$S_N$")
ax1.grid()

def H(t, H0 = 0, Hpert = 0.37, Tris = 700):
    T0 = 0
    
    rris    = (Hpert - H0)/Tris;
    
    if t <= T0:
        output = H0
    elif (t>T0 and t<=(Tris+T0)):
        output = (rris * t + (H0 - (rris * T0))) # linear 
    elif (t>(Tris+T0)):
        output = Hpert;
             
    return output

Hosing = [H(t) for t in tspan]


ax2 = ax1.twinx() 
ax2.plot(tspan,Hosing)
ax2.set_ylabel("$H$")

plt.show()
