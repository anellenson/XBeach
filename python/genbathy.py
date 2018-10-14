import numpy as np
import pandas as pd
import datetime as dt
import netCDF4 as nc
import os
import matplotlib.pyplot as pl
from matplotlib import interactive
interactive(True)
from scipy.optimize import fsolve


#First pull in the bathymetry you want to use
#os.chdir('/Users/ashelyellenson/Research/XBeach/python')
import bathytools as bt

bathyset = bt.retrieveData(2016,3,1)
pNumber = 1006
time,x,transect = bt.transect(bathyset,pNumber)


##beta is the asymptotic offshore slope, assume an x number to be seaward of the active
##bar

#Shift to shore base coordinate system so h=0 at x=0. Cut the transect at h=0, and then shift
#x coordinate system so it's 0 where h = 0.
shore_ind = np.where([np.abs(transect)-0] == np.min([np.abs(transect) - 0]))[1][0]
x = x-x[shore_ind]
#transect_shoreline = transect[shore_ind:]
#Find slope seaward of bar movement. 350 is relatively arbitrary, found qualitatively.
offshore = np.where([-1*transect-4.5] == np.min(np.abs(-1*transect-4.5)))[1][0]
beta = np.divide(np.diff(transect[offshore:]),np.diff(x[offshore:]))
beta_0 = np.mean(beta)

##Choose a point to find depth h_, randomly chose 30 offshore of the offshore point
h_ = transect[offshore+30]
x_ = x[offshore+30]

#Find slope at the shoreline
shore_range = np.where((transect < 0.5) & (transect > -0.5))[0]
beta_s = np.divide(np.diff(transect[shore_range]),np.diff(x[shore_range]))
beta_s = np.mean(beta_s)

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

#add the shoreline back in
#h0 = np.append(transect[0:shore_ind],h0)
pl.plot(x,h0)
###########################################
## Find hbar
a = 0.53
b = 0.57
delta = 0.3
c = 0.09
#Find sea depth and shore depth
d0 = np.multiply(-1,h0)

offshore = np.where(np.abs([d0-4.5]) == np.min(np.abs(d0-4.5)))[1][0]
d0_shore = d0[shore_ind]
d0_sea = d0[offshore]
x_sea = x[offshore]
Smax = 0.2*d0_sea
bar_zone_depth = [(dd - d0_shore)/(d0_sea - d0_shore) for dd in d0]
bar_zone_depth = [bb if bb < 1 else 1 for bb in bar_zone_depth]
expfun = [np.exp(-1*((1 - bb)**a - b)**2/c) for bb in bar_zone_depth]
xmax = x[np.where(expfun == np.nanmax(expfun))[0]][0]


S = []
for xi,xx in enumerate(x):
        S.append((Smax - delta*(xmax/x_sea))*expfun[xi])

aL = 100
bL = 0.27
#Find theta, theta is the integral of 2pi/L from the offshore location to the shoreline
L = [aL*np.exp(bL*dd) for dd in d0]
theta = np.zeros((x.shape))
for ii in range(len(x[:offshore])):
    int_L = L[ii:offshore]
    int_x = x[ii:offshore]
    int_L.reverse()
    int_x = np.flip(int_x,axis = 0)
    integrand = [(2*np.pi)/ll for ll in int_L]
    integral = np.trapz(integrand,int_x)
    theta[ii] = integral


cos = np.cos(theta - theta[220])
hbar = -1*np.multiply(S,cos)
transect_simulated = np.nansum(np.dstack((h0,hbar)),2).squeeze()
pl.plot(x,transect_simulated)
pl.plot(x,transect,label = str(pNumber))
