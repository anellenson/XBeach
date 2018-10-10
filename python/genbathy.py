import numpy as np
import pandas as pd
import datetime as dt

#First pull in the bathymetry you want to use
base_dir = '/home/server/pi/homes/aellenso/Research/XBeach1/'
surveyinfo = pd.read_pickle(base_dir + '/data/surveyinfo.pickle')
sdates_inds = np.where(surveyinfo['dates'] > dt.datetime(2016,4,1))[0]

