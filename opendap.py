#!/usr/bin/env python
# coding: utf-8

# ### Importing

# -*- coding: utf-8 -*-
import os
os.environ['PROJ_LIB'] = r'/home/luiz/anaconda3/pkgs/proj-7.1.0-h966b41f_1/share/proj'

from math import trunc
import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import datetime as dtt


# ### Loading data

day = 9
month = 2
year = 2014

# [0:1:0]
time_min=0
time_timestep=1
time_max=0
# [0:1:15999],
lat_min=7000
lat_timestep=1
lat_max=9000
# [0:1:35999]
lon_min=17000
lon_timestep=1
lon_max=19000

date =  dtt.date(year, month, day)
formatted_date = str(date).replace('-', '')

url = "https://podaac-opendap.jpl.nasa.gov/opendap/hyrax/allData/ghrsst/data/L4/GLOB/JPL_OUROCEAN/G1SST/{0}/{1:03d}/{2}-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST.nc.bz2?time[{3}:{4}:{5}],lat[{6}:{7}:{8}],lon[{9}:{10}:{11}],analysed_sst[{3}:{4}:{5}][{6:d}:{7}:{8}][{9}:{10}:{11}]".format(date.year,  date.timetuple().tm_yday, formatted_date, time_min, time_timestep, time_max, lat_min, lat_timestep, lat_max, lon_min, lon_timestep, lon_max)

print('Reading via OpenDAP, {} '.format(url))
dataset = nc.Dataset(url, "r", format="NETCDF4" )

print(dataset)


print(dataset.variables["lat"], dataset.variables["lon"], dataset.variables["analysed_sst"], sep="\n")

lats = dataset.variables["lat"][:]
lons = dataset.variables["lon"][:]


# ### Plot


plt.figure(figsize=(10, 10))
mp = Basemap(projection="mill", llcrnrlon=lons.min(), urcrnrlon=lons.max(), llcrnrlat=lats.min(), urcrnrlat=lats.max(), resolution='i')
# mp = Basemap(projection='robin', lon_0=0, resolution='c')

lon, lat = np.meshgrid(lons, lats)
x, y = mp(lon, lat)

color_scheme = mp.pcolormesh(x, y, np.squeeze(dataset.variables["analysed_sst"][:]), cmap = 'jet', shading='auto')

mp.drawcoastlines()
mp.drawstates()
mp.drawcountries()

mp.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
mp.drawmeridians(np.arange(mp.lonmin,mp.lonmax+30,60),labels=[0,0,0,1])

mp.fillcontinents(color='coral',lake_color='aqua')

color_bar = mp.colorbar(color_scheme, location='right', pad="5%")

plt.title("Ocean temperature on {}".format(date))
plt.show()

dataset.close()

