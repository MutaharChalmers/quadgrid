#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functions to calculate distances relative to the grid.
"""

import numpy as np


def dmat(lons1, lats1, lons2, lats2, R=6371.007):
    """Great circle (haversine) distance matrix in km from two different
    pairs of lon, lat arrays. Assumes both arrays in decimal degrees.

    Sources:
        http://en.wikipedia.org/wiki/Haversine_formula
        http://stackoverflow.com/questions/34502254

    Parameters
    ----------
        lons1 : numpy array
            Longitudes in decimal degrees.
        lats1 : numpy array
            Latitudes in decimal degrees.
        lons2 : numpy array
            Longitudes in decimal degrees.
        lats2 : numpy array
            Latitudes in decimal degrees.
        R : float
            Earth's radius in km.

    Returns
    -------
        dmat : numpy array
            Distance matrix in km.
    """

    lons1_rad, lats1_rad = np.deg2rad(lons1), np.deg2rad(lats1)
    lons2_rad, lats2_rad = np.deg2rad(lons2), np.deg2rad(lats2)
    
    dlon = lons1_rad[:, None] - lons2_rad
    dlat = lats1_rad[:, None] - lats2_rad

    d = (np.sin(dlat/2)**2 + np.cos(lats1_rad[:, None]) * 
         np.cos(lats2_rad) * np.sin(dlon/2)**2)
    return 2*R*np.arcsin(np.sqrt(d))
