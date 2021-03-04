import os
import math
import netCDF4 as nc
import numpy as np
import pandas as pd
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import datetime as dtt
import cmocean as cmo


def kelvinToCelsius(k):
    return k - 273.15

# Loading data


def handle_map_region_plot(day, month, year, lat_max, lat_min, lon_min, lon_max):
    # [0:1:0]
    time_min = 0
    time_timestep = 1
    time_max = 0
    # [0:1:15999],
    lat_min_v = int((8000+lat_min*100)+(lat_min/abs(lat_min))
                    ) if lat_min != 0 else 8000
    lat_timestep = 1
    lat_max_v = int((8000+lat_max*100)+(lat_max/abs(lat_max))
                    ) if lat_max != 0 else 8000
    # [0:1:35999]
    lon_min_v = int((18000+lon_min*100)+(lon_min/abs(lon_min))
                    ) if lon_min != 0 else 18000
    lon_timestep = 1
    lon_max_v = int((18000+lon_max*100)+(lon_max/abs(lon_max))
                    ) if lon_max != 0 else 18000

    date = dtt.date(year, month, day)
    formatted_date = str(date).replace('-', '')

    url = "https://podaac-opendap.jpl.nasa.gov/opendap/allData/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/{0}/{1:03d}/{2}090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc?time[{3}:{4}:{5}],lat[{6}:{7}:{8}],lon[{9}:{10}:{11}],analysed_sst[{3}:{4}:{5}][{6:d}:{7}:{8}][{9}:{10}:{11}]".format(date.year,  date.timetuple().tm_yday, formatted_date, time_min, time_timestep, time_max, lat_min, lat_timestep, lat_max, lon_min, lon_timestep, lon_max)

    print('Reading via OpenDAP, {} '.format(url))
    dataset = nc.Dataset(url, "r", format="NETCDF4")

    print(dataset)

    print(dataset.variables["lat"], dataset.variables["lon"],
          dataset.variables["analysed_sst"], sep="\n")

    lats = dataset.variables["lat"][:]
    lons = dataset.variables["lon"][:]

    # Plot
    print(lats.max(), lats.min())
    print(lons.max(), lons.min())

    plt.figure(figsize=(10, 10))
    mp = Basemap(projection="mill", llcrnrlon=lon_min, urcrnrlon=lon_max,
                 llcrnrlat=lat_min, urcrnrlat=lat_max, resolution='i')

    lon, lat = np.meshgrid(lons, lats)
    x, y = mp(lon, lat)

    vfunc = np.vectorize(kelvinToCelsius)
    data = vfunc(np.squeeze(dataset.variables["analysed_sst"][:]))

    average = np.average(data)
    standard_deviation = np.std(data)

    print('average: {}'.format(average), 'std: {}'.format(standard_deviation))

    upper_limit = average + 2*standard_deviation
    lower_limit = average - 2*standard_deviation

    cor = cmo.cm.thermal
    color_scheme = mp.pcolormesh(x, y, data, cmap=cor, shading='auto', vmin=math.floor(
        lower_limit), vmax=math.floor(upper_limit))

    mp.drawcoastlines()
    mp.drawstates()
    mp.drawcountries()

    mp.drawparallels(np.arange(lat_min, lat_max+2, 2), labels=[1, 0, 0, 0])
    mp.drawmeridians(np.arange(lon_min, lon_max+2, 2), labels=[0, 0, 0, 1])

    mp.fillcontinents(color='coral', lake_color='aqua')

    color_bar = mp.colorbar(color_scheme, location='right', extend='both')
    color_bar.ax.set_ylabel('Celsius (°C)')

    plt.title("Ocean temperature on {}".format(date))

    cur_dir = os.path.abspath('')
    plt.savefig(cur_dir + '/maps_region_GHRSST/{1:03d}-{2}-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST.png'.format(
        date.year,  date.timetuple().tm_yday, formatted_date))
    plt.show()
    dataset.close()


start_date = '02/09/2014'
end_date = '02/12/2014'

for timestamp in pd.date_range(start=start_date, end=end_date):
    date = timestamp.date()
    handle_map_region_plot(date.day, date.month, date.year, 0, -20, -40, -20)
