import os

import geopandas as gpd

from io import BytesIO
from tempfile import TemporaryDirectory
from zipfile import ZipFile
from urllib.request import urlopen


wsc_dareas = {1 : 'Maritime Provinces',
              2 : 'St. Lawrence',
              3 : 'Northern Quebec and Labrador',
              4 : 'Southwestern Hudson Bay',
              5 : 'Nelson River',
              6 : 'Western and Northern Hudson Bay',
              7 : 'Great Slave Lake',
              8 : 'Pacific',
              9 : 'Yukon River',
              10 : 'Arctic',
              11 : 'Mississippi River'}


class WSCStation:
    
    def __init__(self, stn_num):
        
        if isinstance(stn_num, str):
            pass
        else:
            stn_num = str(stn_num)
        self.stn_num = stn_num
        
        
    # @classmethod
    # def from_name(cls, stn_name):
    #     return cls(stn_num)
    
    def location(self):
        gdf_shp = self._read_shps_from_tmp(shp_type='station')
        lon = gdf_shp.geometry.x[0]
        lat = gdf_shp.geometry.y[0]
        return lat, lon
        
    
    def drainage_basin(self):
        code = int(self.stn_num[0:2])
        return code, wsc_dareas[code]
        
    
    def download_shps(self, path, shp_type):
        
        os_path = os.path.join(path)
        stn = self.stn_num
        zipfile = self._grab_zip()
        shp_typ_names = self._nhs_net_fnames()
        
        if shp_type not in ('d_basin', 'pour_pnt', 'station', 'all'):
            raise ValueError(shp_type)
        else:
            pass
        
        prefix = []
        if shp_type == 'all':
            for key in shp_typ_names.keys():
                prefix.append(f"{stn}/{stn}_{shp_typ_names[shp_type]}")
        else:
            prefix.append(f"{stn}/{stn}_{shp_typ_names[shp_type]}")
        
        all_files = zipfile.namelist()
        stn_files = []
        for pre in prefix:
            for f in all_files:
                if pre in f:
                    stn_files.append(f)

        for f in stn_files:
            zipfile.extract(f, path=os_path)
        
    
    def _read_shps_from_tmp(self, shp_type):
        
        stn = self.stn_num
        zipfile = self._grab_zip()
        shp_typ_names = self._nhs_net_fnames()
        
        if shp_type not in ('d_basin', 'pour_pnt', 'station'):
            raise ValueError(shp_type)
        else:
            prefix = f"{stn}/{stn}_{shp_typ_names[shp_type]}"
        
        all_files = zipfile.namelist()
        stn_files = [f for f in all_files if prefix in f]
        
        with TemporaryDirectory() as tmpdirname:
            for f in stn_files:
                zipfile.extract(f, path=tmpdirname)
            
            tfile = os.path.join(tmpdirname, prefix+'.shp')
            gdf_shp = gpd.read_file(tfile)
        
        return gdf_shp
    
    
    def _nhs_net_fnames(self):
        names = {'d_basin' : "DrainageBasin_BassinDeDrainage",
                 'pour_pnt' : "PourPoint_PointExutoire",
                 'station' : "Station"}
        return names
    
    
    def _grab_zip(self):
        parent = "https://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/HydrometricNetworkBasinPolygons/"
        zipfolder = parent + self.stn_num[0:2] + '.zip'
        resp = urlopen(zipfolder)
        zipfile = ZipFile(BytesIO(resp.read()))
        return zipfile
    
    