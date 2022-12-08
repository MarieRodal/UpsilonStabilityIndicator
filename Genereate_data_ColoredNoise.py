# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 15:35:58 2021

@author: marie
"""


import numpy as np
import matplotlib.pyplot as plt


def dW(t,tmax):
    """Sample a random number at each call."""
    b = (10-0.2)/tmax
    sigma = 0.2 + b*t
    randn = np.random.normal(loc=0.0, scale=sigma)
    return randn 

def phi_1(t,tmax):
    """AR(1) coefficient of red noise"""
    b = (0.95-0.1)/tmax
    phi = 0.1 + b*t
    return phi

tmin=0

tmax = 1000

ts =np.linspace(tmin,tmax, 10000)

dt = (tmax-tmin)/len(ts)

X0 = -0.0

Xs = np.zeros(len(ts))
eta = np.zeros(len(ts))

eta[0] = 0

Xs[0] = X0


for i in range(1, ts.size):
    t = (i - 1) * dt
    X = Xs[i - 1]
    eta[i] = phi_1(t,tmax)*eta[i-1] + dW(t,tmax)
    Xs[i] = X - 5*X*dt + eta[i]
    
fig = plt.figure(figsize=(15,5))

ax1 = fig.add_subplot(1,1,1)
ax1.plot(ts, Xs, color = 'seagreen')
ax1.set_xlabel("time")
ax1.set_ylabel("x")
ax1.grid()

plt.show()
