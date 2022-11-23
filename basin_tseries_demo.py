"""
Small script to demonstrate how to create and visualize a time series of 
basin-averaged temperature using pyshedtools.

In this example, we are looking to compare the basin-average temp and
the temp at a Water Survey of Canada station of interest.
"""

# load requirements
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# load pyshed utilities
from pyshed.nhs_polys import WSCStation
from pyshed.grid import Grid, VoronoiMask


#------------------------
# Create a station object
#------------------------

# for this example, we're going to look at WSC station, 01BE001
station_num = '01BE001'

# instantiate a station object
stn = WSCStation(station_num)

# check the station's major drainage area
basin_code, basin = stn.drainage_basin()

# grab the lat/lon coordinates of the station
tst_lat, tst_lon = stn.location()


#------------------
# Read gridded data
#------------------

var = 'RDRS_v2.1_P_TT_1.5m'

# create path to sample files
gdata_path = os.path.join(os.getcwd(), 'sample_data', 'gridded_data')
files = os.listdir(gdata_path)

arrys = []
for file in files:
    tmp_path = os.path.join(os.getcwd(), 'sample_data', 'gridded_data', file)
    
    # load data from sample file into DataArray
    ds_tmp = xr.open_dataset(tmp_path)
    arrys.append(ds_tmp)
    
ds = xr.concat([ds_tmp for ds_tmp in arrys], dim='time')

# instantiate Grid object
grid = Grid(ds.lon.values, ds.lat.values)


#-------------------------------------------------------------------
# Mask the gridded data with basin shapefile for station of interest
#-------------------------------------------------------------------

mask = VoronoiMask(grid=grid, stn='01BE001')

# calculate the spatial avg for the basin containing the sample station
avg = np.array([mask.spatial_avg(ds[var][i]) for i in range(len(ds[var]))])


#-------------------------
# Creating the time series
#-------------------------

# tst_lon = -50.712
# tst_lat = 47.5615

# find grid cell that contains test lat/lon from above
cell_x, cell_y = grid.find_nearest_cell(tst_lon, tst_lat)

# extract the data at the grid point corresponding to the sample station
stn_data = np.array([ds[var][i][cell_x, cell_y].data for i in range(len(ds[var]))])

# create new dataset with time series of station data and basin-avg data
ds_ts = xr.Dataset(
    {
     "stn_data" : (["time"], stn_data),
     "basin_data" : (["time"], avg),
     },
    coords={
        "time" : ds.time.data,
        },
    
    )


#----------------------------
# Visualizing the time series
#----------------------------

fig, ax = plt.subplots(figsize=(9.4, 4.8), dpi=200)

ds_ts.plot.scatter(x="time", 
                   y="stn_data", 
                   label="At WSC station: "+station_num, 
                   ax=ax)
ds_ts.plot.scatter(x="time", 
                   y="basin_data",
                   label="Basin average",
                   ax=ax)

ax.set_ylabel(var)
plt.legend(loc='lower left')
plt.show()




