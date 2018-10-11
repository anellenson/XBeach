import numpy as np
import pandas as pd
import datetime as dt
import netCDF4 as nc
import os
import matplotlib.pyplot as pl
from scipy.optimize import fsolve
from math import log10, floor

#First pull in the bathymetry you want to use
os.chdir('Research/XBeach1/python/')
import bathytools as bt

bathyset = bt.retrieveData(2016,3,1)
time,x,transect = bt.transect(bathyset,960)


##beta is the asymptotic offshore slope, assume an x number to be seaward of the active
##bar 

#Shift to shore base coordinate system so h=0 at x=0. Cut the transect at h=0, and then shift 
#x coordinate system so it's 0 where h = 0.
shore_ind = np.where([np.abs(transect)-0] == np.min([np.abs(transect) - 0]))[1][0]

#Find slope seaward of bar movement. 350 is relatively arbitrary, found qualitatively.
offshore = np.where(x>350)[0]
beta = np.divide(np.diff(transect[offshore]),np.diff(x[offshore]))
beta_0 = np.mean(beta)

##Choose a point to find depth h_, randomly chose 30 offshore of the offshore point
h_ = transect[offshore[10]]
x_ = x[offshore[10]]

#Find slope at the shoreline
beta_s = np.divide(np.diff(transect[shore_ind:shore_ind+2]),np.diff(x[shore_ind:shore_ind+2]))[0]

#Now, find gamma and k by solving the two functions (equations 5 and 7) simultaneously
def equations5and7(vars):
    gamma,k = vars
    fun1 = h_ - gamma*(np.exp(-k*x_)-1) - beta_0*x_ 
    fun2 = beta_s + gamma*k - beta_0 
    return fun1,fun2
guess = [-10,10]
sol = fsolve(equations5and7,guess)
gamma = np.float(sol[0])
k = np.float(sol[1])

#Find background transect shape, h0
h0 = [gamma*(np.exp(-k*xx) -1) + beta_0*xx for xx in x]

#add the shoreline back in
h0 = np.append(transect[0:shore_ind],h0)

###########################################
## Find hbar
a = 0.53
b = 0.57
delta = 0.3
c = 0.09
#Cut the transect so it's just between the bar activity
offshore = np.where(x>350)[0][1]
shore_ind = np.where([np.abs(h0)-0] == np.min([np.abs(h0) - 0]))[1][0]
h0_shore = -1*h0[0]
h0_sea = -1*h0[offshore]
x_sea = x[offshore]
x_shore = x[shore_ind]
x_bar = x[shore_ind:offshore]
Smax = 0.2*h0_sea

bar_zone_depth = [(hh - h0_shore)/(h0_sea - h0_shore) for hh in -1*h0]
expfun = [np.exp(-1*((1 - bb)**a - b)**2/c) for bb in bar_zone_depth]
xmax = x[np.where(expfun == np.nanmax(expfun))[0]][0]

S = []
for xi,xx in enumerate(x):
    if isnan(expfun[xi]):
        S.append(delta*(xx/x_sea))
    else:
        S.append(delta*(xx/x_sea) + (Smax - delta*(xmax/x_sea))*expfun[xi])

aL = 100
bL = 0.27
Smax = 0.2*hsea


#Find theta, theta is the integral of 2pi/L from the offshore location to the shoreline
L = [aL*np.exp(bL*hh) for hh in -1*h0]
theta = np.zeros((x.shape))
for ii in range(len(x[:offshore])):
    int_L = L[ii:offshore]
    int_x = x[ii:offshore]
    int_L.reverse()
    int_x = np.flip(int_x,axis = 0)
    integrand = [(2*np.pi)/ll for ll in int_L]
    integral = np.trapz(integrand,int_x)
    theta[ii] = integral 

cos = np.cos(theta)
hbar = -1*np.multiply(S,cos)
transect_simulated = np.nansum(np.dstack((h0,hbar)),2).squeeze()