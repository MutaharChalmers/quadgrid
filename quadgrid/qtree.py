#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class with associated functions to convert between quadcell centroid (lon, lat)
coordinates and quadtreeIDs (qids). Coordinates must be in decimal degrees, and
longitudes must be defined from -180 to 180 degrees. Objects have four methods:

    ll2qid - converts a single (lon, lat) pair to qid
    qid2ll - converts a single qid to a (lon, lat) pair
    lls2qids - converts numpy arrays of lons and lats to an array of qids
    qids2lls - converts a numpy array of qids to arrays of lons and lats

From v0.2.0, robust to edge-case floating point round-off errors by working
in milliarcseconds. Pre-v0.2.0 functions are still available in the package
namespace, but it is recommended to use the QTree class in future.
"""

import numpy as np


class QTree():
    def __init__(self, res, mas=False):
        """Create QTree object for geospatial analysis.

        Parameters
        ----------
        res : float
            Target resolution in decimal degrees.
        """

        self.res = res
        self.mas = mas
        self.mas_per_degree = 3_600_000
        self.res_mas = int(res * self.mas_per_degree)

        # Determine number of levels given target resolution
        self.i_max = int(np.ceil(np.log2(180/self.res)))

    def __repr__(self):
        mas_str = '' if not self.mas else '[mas]'
        return f'QuadTree({self.res}°{mas_str})'

    def ll2qid(self, lon, lat, verbose=False):
        """Converts a lon and lat pair to a qid.

        Parameters
        ----------
        lon : float
            Longitude in decimal degrees.
        lat : float
            Latitude in decimal degrees.
        verbose : bool, optional
            Output extra information about the cell centroid.

        Returns
        -------
        qid : int
            Quadcell id.
        """

        # Discretise lon and lat to integer number of cell widths
        if not self.mas:
            lon_int = int(np.ceil(lon/self.res))
            lat_int = int(np.ceil(lat/self.res))
            shift = 1 << self.i_max
        # Convert lon and lat to milliarcseconds
        else:
            lon_int = int(np.round(lon*self.mas_per_degree))
            lat_int = int(np.round(lat*self.mas_per_degree))
            shift = self.res_mas << self.i_max

        # Initialise origin and qid
        origin_lon_int, origin_lat_int, qid = 0, 0, 0

        for i in range(self.i_max+1):
            delta = 4**(self.i_max-i)
            shift >>= 1
            right = lon_int > origin_lon_int
            top = lat_int > origin_lat_int
            if top and right:
                # Do nothing to the qid
                origin_lon_int += shift
                origin_lat_int += shift
            elif top and not right:
                qid += delta
                origin_lon_int -= shift
                origin_lat_int += shift
            elif not top and not right:
                qid += 2*delta
                origin_lon_int -= shift
                origin_lat_int -= shift
            else:
                qid += 3*delta
                origin_lon_int += shift
                origin_lat_int -= shift

        if verbose:
            if self.mas:
                print(f'({lon}, {lat}) -> {qid}')
                centroid_lon = origin_lon_int/self.mas_per_degree
                centroid_lat = origin_lat_int/self.mas_per_degree
            else:
                print(f'({lon}, {lat}) -> {qid}')
                dlon = self.res/2 if lon_int > origin_lon_int else -self.res/2
                dlat = self.res/2 if lat_int > origin_lat_int else -self.res/2
                centroid_lon = origin_lon_int * self.res + dlon
                centroid_lat = origin_lat_int * self.res + dlat
            print(f'Centroid=({centroid_lon}, {centroid_lat})')
        return qid

    def qid2ll(self, qid):
        """Converts a single qid to quadcell centroid lon and lat.

        Parameters
        ----------
        qid : int
            Quadcell qid.

        Returns
        -------
        lon : float
            Quadcell centroid longitude in decimal degrees.
        lat : float
            Quadcell centroid latitude in decimal degrees.
        """

        # Initialise integer lon, lat and shift variables
        lon_int, lat_int, shift = 0, 0, 1

        for i in range(self.i_max+1):
            qid, mod = divmod(qid, 4)
            if mod == 0:
                lon_int += shift
                lat_int += shift
            elif mod == 1:
                lon_int -= shift
                lat_int += shift
            elif mod == 2:
                lon_int -= shift
                lat_int -= shift
            else:
                lon_int += shift
                lat_int -= shift
            shift <<= 1

        lon = lon_int * self.res/2
        lat = lat_int * self.res/2
        return lon, lat

    def lls2qids(self, lons, lats):
        """Converts arrays of lons and lats to an array of qids.

        Parameters
        ----------
        lons : ndarray
            Array of longitudes in decimal degrees.
        lats : ndarray
            Array of latitudes in decimal degrees.

        Returns
        -------
        qids : ndarray
            Array of qids.
        """

        # Discretise lons and lats to integer numbers of cell widths
        if not self.mas:
            lons_int = np.ceil(lons/self.res).astype(np.int64)
            lats_int = np.ceil(lats/self.res).astype(np.int64)
            shift = 1 << self.i_max
        # Convert lons and lats to milliarcseconds
        else:
            lons_int = np.round(lons*self.mas_per_degree).astype(np.int64)
            lats_int = np.round(lats*self.mas_per_degree).astype(np.int64)
            shift = self.res_mas << self.i_max

        # Initialise qid array, origins and shift
        qids = np.zeros_like(lons, dtype=np.int64)
        origin_lons_int = np.zeros_like(lons, dtype=np.int64)
        origin_lats_int = np.zeros_like(lats, dtype=np.int64)

        for i in range(self.i_max+1):
            delta = 4**(self.i_max-i)
            shift >>= 1
            right = lons_int > origin_lons_int
            top = lats_int > origin_lats_int

            # Define disjoint masks by quadrant for all points in the arrays
            mask0 = top & right
            mask1 = top & ~right
            mask2 = ~top & ~right
            mask3 = ~top & right

            # Shift origins by quadrant
            # Do nothing to the qids for mask0
            origin_lons_int[mask0] += shift
            origin_lats_int[mask0] += shift
            qids[mask1] = qids[mask1] + delta
            origin_lons_int[mask1] -= shift
            origin_lats_int[mask1] += shift
            qids[mask2] = qids[mask2] + 2*delta
            origin_lons_int[mask2] -= shift
            origin_lats_int[mask2] -= shift
            qids[mask3] = qids[mask3] + 3*delta
            origin_lons_int[mask3] += shift
            origin_lats_int[mask3] -= shift
        return qids

    def qids2lls(self, qids):
        """Converts array of qids to arrays of quadcell centroids lons and lats.

        Parameters
        ----------
        qids : ndarray
            Array of qids.

        Returns
        -------
        lons : ndarray
            Array of quadcell centroid longitudes in decimal degrees.
        lats : ndarray
            Array of quadcell centroid latitudes in decimal degrees.
        """

        # Initialise integer lon and lat arrays
        lons_int = np.zeros_like(qids, dtype=np.int64)
        lats_int = np.zeros_like(qids, dtype=np.int64)
        shift = 1

        for i in range(self.i_max+1):
            qids, mods = np.divmod(qids, 4)
            mask0 = mods == 0
            mask1 = mods == 1
            mask2 = mods == 2
            mask3 = mods == 3
            lons_int[mask0] += shift
            lats_int[mask0] += shift
            lons_int[mask1] -= shift
            lats_int[mask1] += shift
            lons_int[mask2] -= shift
            lats_int[mask2] -= shift
            lons_int[mask3] += shift
            lats_int[mask3] -= shift
            shift <<= 1

        lons = lons_int * self.res/2
        lats = lats_int * self.res/2
        return lons, lats
