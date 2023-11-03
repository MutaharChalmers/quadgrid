#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Class to generate a quadcell Uniform Resolution Grid (URG) at arbitrary 
resolutions and to perform simple operations.
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely as shp
import shapely.vectorized as shpv
from .convert import lls2qids
from .distance import dmat


# Authalic radius of Earth in kilometres
R = 6371.007


class QuadGrid():
    def __init__(self, res, lon_from=-180, lon_to=180, lat_from=-90, lat_to=90):
        self.res = res
        self.lon_from = lon_from
        self.lon_to = lon_to
        self.lat_from = lat_from
        self.lat_to = lat_to
        self.lons_1d = np.arange(lon_from+res/2, lon_to+res/2, res)
        self.lats_1d = np.arange(lat_from+res/2, lat_to+res/2, res)
        lons_2d, lats_2d = np.meshgrid(self.lons_1d, self.lats_1d, indexing='xy')
        self.lons, self.lats = lons_2d.ravel(), lats_2d.ravel()
        self.qids = lls2qids(self.lons, self.lats, res)
        self.mask = np.full(self.qids.shape, True, dtype=bool)
        self.dmat = None
        self.mix = pd.MultiIndex.from_arrays([self.lats, self.lons], 
                                             names=['lat','lon'])
        
        # Approximate quadcell areas in km2 assuming spherical Earth.
        # See: https://gis.stackexchange.com/questions/29734/
        res_rad = np.deg2rad(self.res)
        areas = (np.sin(np.deg2rad(self.lats_1d+self.res/2)) - 
                 np.sin(np.deg2rad(self.lats_1d-self.res/2))) * res_rad * R**2
        self.areas = np.repeat(areas, self.lons_1d.size)
        
    def __repr__(self):
        return f'QuadGrid ({self.res}) ' \
               f'{self.lon_from}<=lon<={self.lon_to} | ' \
               f'{self.lat_from}<=lat<={self.lat_to}'
        
    def apply_mask(self, geom, buff=None):
        """Generate boolean mask for all quadcells within geom."""
        if buff is None:
            buff = np.sqrt(2*self.res**2)/2
        self.mask = shpv.contains(geom.buffer(buff), self.lons, self.lats)
        
    def distance(self, lon, lat):
        """Distance matrix in km between a single point and all quadcells."""
        return pd.Series(dmat(self.lons, self.lats, lon, lat, R).ravel(), 
                            index=self.mix, name=f'distance_km')

    def to_pandas(self):
        """Convert to pandas DataFrame."""
        return pd.DataFrame({'lat': self.lats, 'lon': self.lons, 
                             'qid': self.qids, 'area': self.areas,
                             'mask': self.mask}).set_index(['lat','lon']
                                                          ).sort_index()
        
    def to_geopandas(self):
        """Convert to geopandas GeoDataFrame."""
        geoms = [shp.geometry.Polygon([(lon+self.res/2, lat+self.res/2),
                                       (lon+self.res/2, lat-self.res/2),
                                       (lon-self.res/2, lat-self.res/2),
                                       (lon-self.res/2, lat+self.res/2)]) 
                 for lon, lat in zip(self.lons, self.lats)]
        return gpd.GeoDataFrame({'lat': self.lats, 'lon': self.lons, 
                                 'qid': self.qids, 'area': self.areas,
                                 'mask': self.mask, 'geometry': geoms}, 
                                crs='epsg:4326')
    
    def to_xarray(self):
        """Convert to xarray Dataset."""
        attrs = {'Type': 'quadgrid', 'Resolution': self.res, 'Area units': 'km2'}
        ds = self.to_pandas().to_xarray().reindex(lon=self.lons_1d, lat=self.lats_1d)
        return ds.assign_attrs(**attrs)
