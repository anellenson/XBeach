import numpy as np 
import matplotlib.pyplot as pl
import netCDF4 as nc
import datetime as dt
import BuoySpec


def JSwapSpec(gamma,beta,alpha,peak_wf,wavefreq):
    sigma_lo = 0.07
    sigma_hi = 0.09
    wavefreq_lo = [x for x in wavefreq if x <= peak_wf]
    wavefreq_hi = [x for x in wavefreq if x > peak_wf]
    
    a_lo = [np.exp(-(x - peak_wf)**2/(2*peak_wf**2*sigma_lo**2)) for x in wavefreq_lo]
    a_hi = [np.exp(-(x - peak_wf)**2/(2*peak_wf**2*sigma_hi**2)) for x in wavefreq_hi]
 
    
    Jswap = []
    for ai, a_val in enumerate(a_lo):
        Jswapval = (alpha * 9.8**2)/wavefreq_lo[ai]**5 * np.exp(-beta * (peak_wf**4)/wavefreq_lo[ai]**4)*gamma**a_val
        Jswap.append(Jswapval) 
    for ai, a_val in enumerate(a_hi):
        Jswapval = (alpha * 9.8**2)/wavefreq_hi[ai]**5 * np.exp(-beta * (peak_wf**4)/wavefreq_hi[ai]**4)*gamma**a_val
        Jswap.append(Jswapval) 

    
    return Jswap

#Find frequency and directional bins
dt_start = dt.datetime(2011,9,30)

filename = 'FRF-ocean_waves_8m-array_' + dt_start.strftime('%Y%m') + '.nc'
ncfile = nc.Dataset('https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/oceanography/waves/8m-array/' + dt_start.strftime('%Y') +'/' + filename)
freq = ncfile['waveFrequency'][:]
dir = ncfile['waveDirectionBins'][:]

dir = np.arange(-180,180,5)
alpha = 6.37E-6
beta  = 1
T = 14.2 #s
peak_wf = 1/T
gamma = 5
incident_dir_nautical = 250
incident_dir = 90 - incident_dir_nautical + 180
H_desired = 1.2
jswap_spec = JSwapSpec(gamma, beta, alpha, peak_wf,freq)


#Check to make sure it's the right wave height

m0 = np.trapz(jswap_spec,freq)
Hm0 = 4*np.sqrt(m0)
beta = np.abs(Hm0-H_desired)

while beta>0.01:
    if Hm0 > H_desired:
        alpha = alpha-0.01E-6
        jswap_spec = JSwapSpec(gamma, beta, alpha, peak_wf,freq)
        m0 = np.trapz(jswap_spec,freq)
        Hm0 = 4*np.sqrt(m0)
        beta = np.abs(Hm0-H_desired)
        print(str(Hm0))
    if Hm0 < H_desired:
        alpha = alpha + 0.01E-6
        jswap_spec = JSwapSpec(gamma, beta, alpha, peak_wf,freq)
        m0 = np.trapz(jswap_spec,freq)
        Hm0 = 4*np.sqrt(m0)
        beta = np.abs(Hm0-H_desired)
        print(str(Hm0))

#Create a matrix
dir_spec = np.zeros((len(freq),len(dir)))

diff_dir = np.unique(np.diff(dir))
jswap_spec = np.divide(jswap_spec,diff_dir)
#Chose an angle of 270
dir_ind = np.where(dir == incident_dir)[0]
dir_spec[:,dir_ind[0]] = jswap_spec



pl.figure()
ax = pl.subplot(111,projection = 'polar')
BuoySpec.polarSpecPlot_NDBC(freq,dir,dir_spec,dt.datetime(2009,11,30),ax)

freq_num = len(freq)
wavedir_num = len(dir)

file = open('/home/server/pi/homes/aellenso/Research/XBeach/runs/FRF_TestCase/jonswap/bathy/extendedforeshore/vardens.txt','w')
file.write('%d\n' %freq_num)
for ff in freq:
    file.write('%7.3f\n' %ff)
file.write('%d\n' %wavedir_num)
for ww in dir:
    file.write('%7.3f\n' %ww)
for line in dir_spec:
    for ll in line:
        file.write('%  11.7e'%ll)
    file.write('\n')
file.close()

