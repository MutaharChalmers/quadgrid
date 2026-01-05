#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Class to generate a quadtree-based Uniform Resolution Grid (URG) at arbitrary
resolutions and to perform simple operations.
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely as shp
from .qtree import QTree
from .distance import dmat


# Authalic radius of Earth in kilometres
R = 6371.007


class QuadGrid():
    def __init__(self, res, lon_bounds=(-180,180), lat_bounds=(-90,90)):
        """Class constructor for QuadGrid.
        Coordinates are defined internally in terms of milliarcseconds.

        Parameters
        ----------
        res : float
            Quadgrid resolution in decimal degrees.
        lon_bounds : (float, float)
            Longitude bounds to be covered by quadgrid.
        lat_bounds : (float, float)
            Latitude bounds to be covered by quadgrid.
        """

        self.mas_per_degree = 3_600_000
        self.res = res
        self.res_mas = int(self.res * self.mas_per_degree)
        self.lon_bounds = lon_bounds
        self.lat_bounds = lat_bounds

        # Global centroid lons and lats in milliarcseconds
        dlons = np.arange(360*self.mas_per_degree//self.res_mas)*self.res_mas
        dlats = np.arange(180*self.mas_per_degree//self.res_mas)*self.res_mas
        _lons = -180*self.mas_per_degree + self.res_mas//2 + dlons
        _lats = -90*self.mas_per_degree + self.res_mas//2 + dlats

        # Apply bounds and create object attributes
        self._lons = _lons[(_lons>=int(lon_bounds[0]*self.mas_per_degree)) &
                           (_lons<=int(lon_bounds[1]*self.mas_per_degree))]
        self._lats = _lats[(_lats>=int(lat_bounds[0]*self.mas_per_degree)) &
                           (_lats<=int(lat_bounds[1]*self.mas_per_degree))]

        # Make mesh, ravel and create object attributes
        _lons_2d, _lats_2d = np.meshgrid(self._lons, self._lats, indexing='xy')
        self._lons_2d, self._lats_2d = _lons_2d.ravel(), _lats_2d.ravel()

        # Create object attributes in decimal degrees for convenience
        self.lons = self._lons / self.mas_per_degree
        self.lats = self._lats / self.mas_per_degree
        self.lons_2d = self._lons_2d / self.mas_per_degree
        self.lats_2d = self._lats_2d / self.mas_per_degree

        # Generate qids, initial mask and MultiIndex
        self.qt = QTree(res, mas=True)
        self.qids = self.qt.lls2qids(self.lons_2d, self.lats_2d)
        self.base_mask = np.full(self.qids.shape, True, dtype=bool)
        self.mix = pd.MultiIndex.from_arrays([self.lats_2d, self.lons_2d],
                                             names=['lat','lon'])

        # Approximate quadcell areas in km2 assuming spherical Earth.
        # See: https://gis.stackexchange.com/questions/29734/
        res_rad = np.deg2rad(self.res)
        self.areas = (np.sin(np.deg2rad(self.lats_2d+self.res/2)) -
                      np.sin(np.deg2rad(self.lats_2d-self.res/2)))*res_rad*R**2

        # Create geopandas GeoDataFrame as the reference version of the grid
        lons_1 = self.lons_2d + self.res/2
        lats_1 = self.lats_2d + self.res/2
        lons_0 = self.lons_2d - self.res/2
        lats_0 = self.lats_2d - self.res/2
        geoms = np.empty((self.lons_2d.size*4, 2), dtype=np.float64)
        geoms[::4,0], geoms[::4,1] = lons_1, lats_1
        geoms[1::4,0], geoms[1::4,1] = lons_0, lats_1
        geoms[2::4,0], geoms[2::4,1] = lons_0, lats_0
        geoms[3::4,0], geoms[3::4,1] = lons_1, lats_0
        geoms = shp.polygons(geoms.reshape((self.lons_2d.size, 4, 2)))

        self.grid = gpd.GeoDataFrame({'lat': self.lats_2d, 'lon': self.lons_2d,
                                      'qid': self.qids, 'area': self.areas,
                                      'mask': self.base_mask, 'res': self.res,
                                      'geometry': geoms}, crs='epsg:4326')

    def __repr__(self):
        return f'QuadGrid({self.res}°) | ' \
               f'{self.lon_bounds[0]}°<=lon<={self.lon_bounds[1]}° | ' \
               f'{self.lat_bounds[0]}°<=lat<={self.lat_bounds[1]}°'

    def apply_mask(self, gdf, from_base=True):
        """Create an intersection boolean mask from another GeoDataFrame."""

        # Check CRS is EPSG:4326
        if gdf.crs != 'epsg:4326':
            gdf = gdf.to_crs('epsg:4326')

        # Do the spatial join
        grid_mask = gpd.sjoin(self.grid, gdf[['geometry']], how='left'
                              ).drop_duplicates('qid')

        # Determine the reference mask - if from_base, apply new mask to base,
        # otherwise, apply the new mask to the current mask
        if from_base:
            ref_mask = grid_mask['mask']
        else:
            ref_mask = grid_mask
        grid_mask['mask'] = ref_mask & ~grid_mask['index_right'].isna()
        self.grid = grid_mask.drop('index_right', axis=1)

    def query(self, lons, lats):
        """QuadcellID lookup."""
        return self.qt.lls2qids(np.atleast_1d(lons), np.atleast_1d(lats))

    def distance(self, lon, lat):
        """Distance matrix in km between a single point and all quadcells."""
        return pd.Series(dmat(self.lons_2d, self.lats_2d, lon, lat, R).ravel(),
                         index=self.mix, name=f'distance_km')

    def to_xarray(self, masked=True):
        # Create xarray Dataset
        attrs = {'Resolution': f'{self.res}°', 'Area units': 'km2'}
        df = self.grid.drop('geometry', axis=1
                            ).set_index(['lat','lon']).sort_index()
        return df.to_xarray().reindex(lon=self.lons, lat=self.lats,
                                      method='nearest').assign_attrs(**attrs)

    def to_geojson(self):
        """Convert to GeoJSON."""
        return self.grid.to_json(show_bbox=True, drop_id=True, to_wgs84=True)
