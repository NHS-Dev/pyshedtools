"""
Small script to demonstrate the utilities available in pyshedtools.
"""

# load requirements
import os
import xarray as xr
# load pyshed utilities
from pyshed.nhs_polys import WSCStation
from pyshed.grid import Grid, VoronoiMask


#---------------
# nhs_poly tools
#---------------

# for this example, we're going to look at WSC station, 01BE001
station_num = '01BE001'

# instantiate a station object
stn = WSCStation(station_num)

# check the station's major drainage area
basin_code, basin = stn.drainage_basin()

# grab the lat/lon coordinates of the station
tst_lat, tst_lon = stn.location()


#-----------
# Grid tools
#-----------

# create path to sample file
test_file = os.path.join(os.getcwd(), 'sample_data', 'gridded_data', '1980010112.nc')

# load data from sample file into DataArray
ds = xr.open_dataset(test_file)

# instantiate Grid object
test_grid = Grid(ds.lon.values, ds.lat.values)

# check the grid shape
grd_shape = test_grid.grid_shape

# find grid cell that contains test lat/lon from above
cell_x, cell_y = test_grid.find_nearest_cell(tst_lon, tst_lat)

# check the lat/lon at cell_x, cell_y to see how close they are to the actual lat/lon
cell_lat = ds.lat[cell_x, cell_y].data
cell_lon = ds.lon[cell_x, cell_y].data


#------------------
# VoronoiMask tools
#------------------

var = 'RDRS_v2.1_P_TT_1.5m'

# choose whether the modules will use the sample data in ../polygons/samples/
# or fetch the necessary data from the internet
sample_data = False

# let's mask the grid with a shapefile for our basin of interest
if sample_data:
    # first, create the shapefile path
    test_shape = os.path.join(os.getcwd(), 'sample_data', 'polygons', station_num,
                              f'{station_num}_DrainageBasin_BassinDeDrainage.shp')
    # mask the grid with the shapefile
    mask = VoronoiMask(grid=test_grid, path_to_shp=test_shape)
else:
    # mask the grid by specifying the station (shapefile fetched from internet)
    mask = VoronoiMask(grid=test_grid, stn=station_num)

# find the spatial average of air temperature for the basin corresponding ...
# to our station of interest
avg = mask.spatial_avg(ds[var][0])

# view voronoi grid with shapefile mask
mask.view_mask()

# view mask intersection with grid
mask.view_itx()

# view data with mask applied
mask.view_masked_data(ds[var][0])





