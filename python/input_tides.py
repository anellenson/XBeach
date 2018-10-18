import pandas as pd
import numpy as np
import datetime as dt
import csv
import requests

datum = 'NAVD'
timeZone = 'lst'
application ='web_services'
Format ='csv'
products = 'hourly_height'
no_days = np.arange(0,2)
start_date = dt.datetime(2017,3,26,15) #Start one hour before the waves
dt_vec = [start_date + dt.timedelta(days = int(x)) for x in no_days]
seconds = []
waterlevel = []

fname = '/Users/ashleyellenson/Research/XBeach/model/run/Hs0.5/tides.txt'

for dd in dt_vec:
    tideurl = 'https://tidesandcurrents.noaa.gov/api/datagetter?begin_date='  + dd.strftime('%Y%m%d') + \
              '&end_date=' + dd.strftime('%Y%m%d') + '&station=8651370&product=' + products + \
              '&datum=' + datum + '&units=metric&time_zone=gmt&application=web_services&format=csv'

    with requests.Session() as s:
        download = s.get(tideurl)

        decoded_content = download.content.decode('utf-8')

        cr = csv.reader(decoded_content.splitlines(), delimiter=',')
        my_list = list(cr)
        for row in my_list[1:]:#First entry is the header, and have to start on the second entry
            dt_row = dt.datetime.strptime(row[0],'%Y-%m-%d %H:%M')
            if dt_row > start_date:
                seconds.append((dt.datetime.strptime(row[0],'%Y-%m-%d %H:%M') - start_date).total_seconds())
                waterlevel.append(np.float(row[1]))

with open(fname, "w") as text_file:
    for ii in range(len(seconds[1:])):
        text_file.write("{0:1d} {1:3f} \n".format(np.int(seconds[ii]), waterlevel[ii]))

text_file.close()
