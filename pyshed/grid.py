import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import logging as log

from numba import jit
from shapely.geometry import Polygon

from pyshed.nhs_polys import WSCStation


def convert_lon_to_360(x):
    """ Converts longitude values from [-180, 180] to [0, 360] """
    return np.mod(x, 360)


def convert_lon_to_180(x):
    """ Converts longitude values from [0, 360] to [-180, 180] """
    return (x + 180) % 360 - 180


class Grid(object):
    
    
    def __init__(self, x, y, crs=None):
        
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        
        self.ny, self.nx = self.x.shape
        
        if crs is None:
            self.crs = 'EPSG:4326'
            log.info("No CRS given. Using default crs: EPSG:4326")
        else:
            self.crs = crs
            
        
    @property
    def grid_shape(self):
        """ Shape of the grid, taken from the x array """
        return self.x.shape


    def find_nearest_cell(self, X, Y):
        """ Finds index position of grid cell closest to point (X, Y) """
        abs_x = np.abs(self.x - X)
        abs_y = np.abs(self.y - Y)
        c = np.maximum(abs_x, abs_y)
        ([xloc], [yloc]) = np.where(c == np.min(c))
        return xloc, yloc
    
    
    # def view_stn(self):
   
    
class VoronoiMask:
    
    
    def __init__(self, grid=None, path_to_shp=None, stn=None, 
                 path_to_vor=None, path_to_clp=None):
        
        if path_to_vor is not None and path_to_clp is not None:
            self._gdf_vor = gpd.read_file(path_to_vor)
            self._gdf_clp = gpd.read_file(path_to_clp)
        elif path_to_vor or path_to_clp is not None:
            msg = ("Only one of path_to_vor or path_to_clp provided. "
                   "To create a VoronoiMask from file, both path_to_vor and "
                   "path_to_clp must be provided."
                   )
            raise TypeError(msg)
        else:
            if grid is None:
                raise TypeError("Missing grid object to create Voronoi mask.")
            else:
                self.grid = grid
            
            self._gdf_vor = self._transform_grid()
            self._gdf_clp = self._clip_grid(path_to_shp, stn)
        
        self._gdf_clp = self._gdf_clp.to_crs(self.grid.crs)


    @classmethod
    def from_grid_params(cls, x, y, crs_grd=None, path_to_shp=None, 
                         stn=None, darea=None):
        grid = Grid(x, y, crs_grd)
        return cls(grid, path_to_shp, stn, darea)


    @classmethod
    def from_files(cls, vor, clp):
        return cls(path_to_vor=vor, path_to_clp=clp)

    
    @property
    def clip(self):
        """ GeoDataFrame of masked Voronoi grid"""
        return self._gdf_clp
    
    
    @property
    def voronoi_gpd(self):
        """ GeoDataFrame of input grid as Voronoi polygons """
        return self._gdf_vor
    
    
    def _clip_grid(self, path_to_shp=None, stn=None):
        if path_to_shp is None and stn is None:
            raise TypeError("Missing argument path_to_shp or stn to clip Voronoi grid.")
        elif path_to_shp and stn is None:
            gdf_clp = self._from_shape(path_to_shp)
        elif stn and path_to_shp is None:
            gdf_clp = self._from_station(stn)
        else:
            raise ValueError("Too many arguments to clip Voronoi grid. Include only one of path_to_shp or stn.")
        return gdf_clp
    
    
    def _make_canvas(self):
        plt.figure(figsize=(6.4, 4.8), dpi=200)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.OCEAN)
        ax.coastlines(resolution='10m', color='black', linewidth=1)
        ax.gridlines(draw_labels=True, dms=True, alpha=0)
        
        minlon = self._gdf_clp.x.min()
        maxlon = self._gdf_clp.x.max()
        minlat = self._gdf_clp.y.min()
        maxlat = self._gdf_clp.y.max()
        
        ax.set_extent([minlon-0.5, maxlon+0.5, minlat-0.5, maxlat+0.5], 
                      crs=ccrs.PlateCarree())
        
        return plt, ax
    
    
    def view_mask(self):
        
        plt, ax = self._make_canvas()
        
        ax.add_geometries(self._gdf_clp.geometry, 
                          crs=ccrs.PlateCarree(), 
                          facecolor='blue',
                          color='blue',
                          alpha=0.3,
                          zorder=1)
        plt.show()
        
        
    def view_itx(self):
        
        plt, ax = self._make_canvas()
        
        gdf_clp_plt = self._gdf_clp
        gdf_vor_plt = self._gdf_vor
        
        ax.add_geometries(gdf_vor_plt.geometry, 
                          crs=ccrs.PlateCarree(), 
                          facecolor='black',
                          color='black',
                          alpha=0.1,
                          linestyle='--',
                          zorder=0)
        
        ax.add_geometries(gdf_clp_plt.geometry, 
                          crs=ccrs.PlateCarree(), 
                          facecolor='blue',
                          color='blue',
                          alpha=0.3,
                          zorder=1)
        
        for i in range(len(self._gdf_clp)):
        
            ax.scatter(gdf_clp_plt.iloc[i].x, gdf_clp_plt.iloc[i].y, 
                       transform=ccrs.PlateCarree(),
                       marker='o', 
                       color='gold', 
                       edgecolors='black',
                       linewidths=0.25,
                       zorder=2)
        
        ax.add_feature(cfeature.OCEAN, zorder=3)
        ax.coastlines(resolution='10m', color='black', linewidth=1, zorder=4)
        
        plt.show()
        
    
    def view_masked_data(self, ds):
        
        plt, ax = self._make_canvas()
        
        ax.contourf(ds.lon, ds.lat, ds.data, transform=ccrs.PlateCarree(),
                    zorder=0)
        
        symdiff = self._gdf_vor.overlay(self._gdf_clp, how='symmetric_difference',
                                      keep_geom_type=False)
        
        ax.add_geometries(symdiff.geometry, 
                          crs=ccrs.PlateCarree(), 
                          facecolor='white',
                          alpha=1.0,
                          zorder=1)
        
        ax.add_feature(cfeature.OCEAN, zorder=2)
        ax.coastlines(resolution='10m', color='black', linewidth=1, zorder=3)
        
        plt.show()
        
        
    def spatial_avg(self, ds):
        # !!!
        # transform gdf_clp and ds to crs with units of meters?
        gdf_clp = self._gdf_clp.to_crs('EPSG:26910')
        gdf_wt = gdf_clp.area / gdf_clp.area.sum()
        sa = 0
        for i in gdf_clp.index:
            sa += ds.data[gdf_clp.ind_y[i], gdf_clp.ind_x[i]] * gdf_wt[i]
        return sa
    
    
    def _from_shape(self, path_to_shp):
        """ Returns GeoDataFrame of Voronoi masked with input shapefile """
        gdf_shp = gpd.read_file(path_to_shp)
        gdf_shp = gdf_shp.to_crs(self.grid.crs)
        gdf_clp = gpd.clip(self._gdf_vor, gdf_shp)
        return gdf_clp
    
    
    def _from_station(self, stn):
        """ 
        Returns GeoDataFrame of Voronoi masked with watershed 
        shapefile corresponding to WCS station.
        """
        gdf_shp = WSCStation(stn)._read_shps_from_tmp(shp_type='d_basin')
        gdf_shp = gdf_shp.to_crs(self.grid.crs)
        gdf_clp = gpd.clip(self._gdf_vor, gdf_shp)
        return gdf_clp

    
    def _create_vor_gpd(self, ind_x, ind_y, x_out, y_out, poly):
        
        polygon = [Polygon(coords) for coords in poly]

        df_coords = pd.DataFrame({'ind_y': ind_y, # lat
                                  'ind_x': ind_x, # lon 
                                  'y': y_out, # lat
                                  'x': x_out}) # lon

        gdf_coords = gpd.GeoDataFrame(df_coords, geometry=polygon)
        gdf_coords = gdf_coords.set_crs(self.grid.crs)
        
        return gdf_coords
    
    
    def _validate_lon(self):
        """ Converts longitude values from [0, 360] to [-180, 180] if necessary """
        if np.any(self.grid.x > 180):
            self.grid.x = convert_lon_to_180(self.grid.x)
        else:
            pass

    
    def _transform_grid(self):
        self._validate_lon()
        ind_x, ind_y, x_out, y_out, poly = _create_vor_polys(self.grid.ny,
                                                             self.grid.nx,
                                                             self.grid.y,
                                                             self.grid.x)
        gdf_vor = self._create_vor_gpd(ind_x, ind_y, x_out, y_out, poly)
        return gdf_vor
 

@jit(nopython=True)
def _create_vor_polys(ny, nx, y, x):
    ind_x = []
    ind_y = []
    x_out = []
    y_out = []
    poly = []
    for i in range(1, ny-1):
        
        for j in range(1, nx-1):
            
            ind_x.append(i)
            ind_y.append(j)
            x_out.append(x[i][j])
            y_out.append(y[i][j])
            
            p0 = (((x[i-1][j-1]) + (x[i][j]))/2, ((y[i-1][j-1]) + (y[i][j]))/2)
            p1 = (((x[i-1][j]) + (x[i][j]))/2, ((y[i-1][j]) + (y[i][j]))/2)
            p2 = (((x[i-1][j+1]) + (x[i][j]))/2, ((y[i-1][j+1]) + (y[i][j]))/2)
            p3 = (((x[i][j+1]) + (x[i][j]))/2, ((y[i][j+1]) + (y[i][j]))/2)
            p4 = (((x[i+1][j+1]) + (x[i][j]))/2, ((y[i+1][j+1]) + (y[i][j]))/2)
            p5 = (((x[i+1][j]) + (x[i][j]))/2, ((y[i+1][j]) + (y[i][j]))/2)
            p6 = (((x[i+1][j-1]) + (x[i][j]))/2, ((y[i+1][j-1]) + (y[i][j]))/2)
            p7 = (((x[i][j-1]) + (x[i][j]))/2, ((y[i][j-1]) + (y[i][j]))/2)

            poly.append([p0, p1, p2, p3, p4, p5, p6, p7, p0])
    
    return ind_x, ind_y, x_out, y_out, poly