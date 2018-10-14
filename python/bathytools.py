import numpy as np
import pandas as pd
import datetime as dt
import netCDF4 as nc

def toDateTime(d):
    return (dt.datetime(1970,1,1) + dt.timedelta(seconds = d))

def toTimestamp(d):
    return (d-dt.datetime(1970,1,1)).total_seconds()

def retrieveData(Yr,Mo,Da):
#First pull in the bathymetry you want to use
    surveyinfo = pd.read_pickle('../data/surveyinfo.pickle')
    sind = np.where(surveyinfo['dates'] > dt.datetime(Yr,Mo,Da))[0][-1]
    bathyset = nc.Dataset('https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/geomorphology/elevationTransects/survey/' + surveyinfo['bathy_fnames'].iloc[sind] + '.nc')
    return bathyset

def transect(bathyset, pNumber):
    profileNumber = bathyset['profileNumber'][:]
    pinds = np.where(profileNumber== pNumber)[0]
    transect = bathyset['elevation'][pinds]
    x = bathyset['xFRF'][pinds]
    time = bathyset['time'][pinds]
    time = [toDateTime(tt) for tt in time]
    return time,x,transect
