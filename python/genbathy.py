import numpy as np
import pandas as pd
import datetime as dt
import netCDF4 as nc
import os
import matplotlib.pyplot as pl
from matplotlib import interactive
interactive(True)
from scipy.optimize import fsolve
import bathytools as bt

bathyset = bt.retrieveData(2016,6,1)
pNumber = 914
time,x,transect = bt.transect(bathyset,pNumber)


##beta is the asymptotic offshore slope, assume an x number to be seaward of the active
##bar
#Find slope at the shoreline
shore_range = np.where((transect < 0.5) & (transect > -0.5))[0]
beta_s = np.divide(np.diff(transect[shore_range]),np.diff(x[shore_range]))
beta_s = -1*np.mean(beta_s)


#Shift to shore base coordinate system so h=0 at x=0. Cut the transect at h=0, and then shift
#x coordinate system so it's 0 where h = 0. Also use soundings instead of depths. 
shore_ind = np.where([np.abs(transect)-0] == np.min([np.abs(transect) - 0]))[1][0]
x = x[shore_ind:]-x[shore_ind]
transect = -1*transect[shore_ind:]

#Find slope seaward of bar movement as where h=4.5m. (Holman)
offshore = np.where(np.abs(transect-4.5) == np.min(np.abs(transect-4.5)))[0][0]
beta = np.divide(np.diff(transect[offshore:]),np.diff(x[offshore:]))
beta_0 = np.mean(beta)

##Choose a point to find depth h_, randomly chose 60 offshore of the offshore point
#Criteria is that betaOffshore * x Offshore < hOffshore
h_ = transect[offshore+60]
x_ = x[offshore+60]


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
print(gamma,k)

#Find background transect shape, h0
h0 = [gamma*(np.exp(-k*xx) -1) + beta_0*xx for xx in x]

###########################################
## Find hbar
#Constants for S
a = 0.53
b = 0.57
delta = 0.3
c = 0.09
#Find sea depth and shore depth
h0 = np.array(h0)
offshore = np.where(np.abs(h0-4.5) == np.min(np.abs(h0-4.5)))[0][0]
h0_shore = h0[shore_ind]
h0_sea = h0[offshore]
x_sea = x[offshore]
Smax = 0.2*h0_sea
bar_zone_depth = [(hh - h0_shore)/(h0_sea - h0_shore) for hh in h0]
expfun = [np.exp(-1*((1 - bb)**a - b)**2/c) for bb in bar_zone_depth]
maxInd = np.where(expfun == np.nanmax(expfun))[0][0]
xmax = x[maxInd)[0]][0]

#S is the envelope function constrained to a maximum amplitude of Smax that's offset in order
#to force that S = 0 at x = 0 and S = delta at h0_sea 
S = []
for xi,xx in enumerate(x[:offshore]):
        S.append(xx/x_sea*delta + (Smax - delta*(xmax/x_sea))*expfun[xi])


#Exponential taper offshore for S with the same slope as the offshore point
betaoff = S[offshore-1]/h0[offshore-1]
k = betaoff/S[-1]
for hh in h0[offshore:]:
    S.append(S[-1]*np.exp(-k*(hh-h0[offshore])))

#Bar phase function constants 
aL = 100
bL = 0.27
#Find theta, theta is the integral of 2pi/L from the offshore location to the shoreline
L = [aL*np.exp(bL*hh) for hh in h0]
theta = np.zeros((x.shape))
for ii in range(len(x[:offshore])):
    int_L = L[ii:offshore]
    int_x = x[ii:offshore]
    int_L.reverse()
    int_x = np.flip(int_x,axis = 0)
    integrand = [(2*np.pi)/ll for ll in int_L]
    integral = np.trapz(integrand,int_x)
    theta[ii] = integral

#Cosine function
cos = np.cos(theta - theta[100])
hbar = -1*np.multiply(S,cos)

transect_simulated = h0 + hbar
