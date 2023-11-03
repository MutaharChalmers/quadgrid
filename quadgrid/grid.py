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
from .qtree import lls2qids
from .distance import dmat


# Authalic radius of Earth in kilometres
R = 6371.007


class QuadGrid():
    def __init__(self, res, lon_bounds=(-180,180), lat_bounds=(-90,90)):
        """Converts a single (lon, lat) quadcell centroid to qid.

        Parameters
        ----------
            res : float
                Quadgrid resolution in decimal degrees.
            lon_bounds : (float, float)
                Longitude bounds to be covered by quadgrid.
            lat_bounds : (float, float)
                Latitude bounds to be covered by quadgrid.
        """
        self.res = res
        self.lon_bounds = lon_bounds
        self.lat_bounds = lat_bounds

        # Global centroid lons and lats
        lons_1d = np.arange(-180+res/2, 180+res/2, res)
        lats_1d = np.arange(-90+res/2, 90+res/2, res)

        # Apply bounds
        self.lons_1d = lons_1d[(lons_1d>(lon_bounds[0]-res/2)) &
                               (lons_1d<(lon_bounds[1]+res/2))]
        self.lats_1d = lats_1d[(lats_1d>(lat_bounds[0]-res/2)) &
                               (lats_1d<(lat_bounds[1]+res/2))]
        lons_2d, lats_2d = np.meshgrid(self.lons_1d, self.lats_1d, indexing='xy')
        self.lons, self.lats = lons_2d.ravel(), lats_2d.ravel()

        # Generate qids, initial mask and MultiIndex
        self.qids = lls2qids(self.lons, self.lats, res)
        self.mask = np.full(self.qids.shape, True, dtype=bool)
        self.mix = pd.MultiIndex.from_arrays([self.lats, self.lons], 
                                             names=['lat','lon'])
        
        # Approximate quadcell areas in km2 assuming spherical Earth.
        # See: https://gis.stackexchange.com/questions/29734/
        res_rad = np.deg2rad(self.res)
        areas = (np.sin(np.deg2rad(self.lats_1d+self.res/2)) - 
                 np.sin(np.deg2rad(self.lats_1d-self.res/2))) * res_rad * R**2
        self.areas = np.repeat(areas, self.lons_1d.size)
        
    def __repr__(self):
        return f'QuadGrid {self.res} deg | ' \
               f'{self.lon_bounds[0]}<=lon<={self.lon_bounds[1]} | ' \
               f'{self.lat_bounds[0]}<=lat<={self.lat_bounds[1]}'
        
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
        attrs = {'Type': 'quadgrid', 'Resolution': f'{self.res} deg', 'Area units': 'km2'}
        ds = self.to_pandas().to_xarray().reindex(lon=self.lons_1d, lat=self.lats_1d)
        return ds.assign_attrs(**attrs)
