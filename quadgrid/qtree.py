#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to convert between quadcell centroid (lon, lat) coordinates and
quadtreeIDs (qids). Coordinates must be in decimal degrees, and longitudes
must be defined from -180 to 180 degrees.

    ll2qid - converts a single (lon, lat) pair to qid
    qid2ll - converts a single qid to a (lon, lat) pair
    lls2qids - converts numpy arrays of lons and lats to an array of qids
    qids2lls - converts a numpy array of qids to arrays of lons and lats
    
The functions all start by classifying a (lon, lat) point into one of the 
main quadrants on the Earth's surface. The origin (initially {0,0}) is then 
moved to the centre of that quadrant ({res/2, res/2} in quadrant coordinates),
new child quadrants are defined, and the point is classified into one of the
child quadrants. The process continues recursively until the target resolution
is reached. The initial resolution is chosen such that it is:
  (a) larger by a factor of a power of 2 than the grid resolution
  (b) >=180 degrees
"""

import numpy as np


def ll2qid(lon, lat, res_target, verbose=False):
    """Converts a single (lon, lat) quadcell centroid to qid.

    Parameters
    ----------
        lon : float
            Quadcell centroid longitude in decimal degrees.
        lat : float
            Quadcell centroid latitude in decimal degrees.
        res_target : float
            Target resolution.

    Returns
    -------
        qid : int
            Quadcell qid.
    """

    # Determine resolution of top-level quadrants, given target resolution
    i_max = int(np.ceil(np.log2(180/res_target)))
    res = res_target * 2**i_max

    i, qid, origin_lon, origin_lat = 0, 0, 0, 0
    while res >= res_target:
        delta = 4**(i_max-i)
        if lon >= origin_lon and lat >= origin_lat:
            # Do nothing to the qid
            origin_lon, origin_lat = (origin_lon+res/2, origin_lat+res/2)
        elif lon < origin_lon and lat >= origin_lat:
            qid = qid + delta
            origin_lon, origin_lat = (origin_lon-res/2, origin_lat+res/2)
        elif lon < origin_lon and lat < origin_lat:
            qid = qid + 2*delta
            origin_lon, origin_lat = (origin_lon-res/2, origin_lat-res/2)
        elif lon >= origin_lon and lat < origin_lat:
            qid = qid + 3*delta
            origin_lon, origin_lat = (origin_lon+res/2, origin_lat-res/2)
        else:
            print(f'Error:\n lon: {lon}\n lat: {lat}')
            return None

        # Halve the resolution and update counter
        res /= 2
        i += 1

    if verbose:
        print(f'({lon}, {lat}) -> {qid}')
        print(f'Resolution={2*res} | centroid=({origin_lon}, {origin_lat})')
    return qid


def qid2ll(qid, res_target, verbose=False):
    """Converts a single qid to the (lon, lat) of the centroid of the quadcell.

    Parameters
    ----------
        qid : int
            Quadcell qid.
        res_target : float
            Target resolution.

    Returns
    -------
        lon : float
            Quadcell centroid longitude in decimal degrees.
        lat : float
            Quadcell centroid latitude in decimal degrees.

    """

    qid0, lon, lat, delta = qid*1, 0, 0, res_target/2

    while delta <= 180:
        qid, mod = divmod(qid, 4)
        if mod == 0:
            lon += delta
            lat += delta
        elif mod == 1:
            lon -= delta
            lat += delta
        elif mod == 2:
            lon -= delta
            lat -= delta
        elif mod == 3:
            lon += delta
            lat -= delta
        delta *= 2

    if verbose:
        print(f'{qid0} -> ({lon}, {lat})')
        print(f'Resolution={res_target}')
    return lon, lat


def lls2qids(lons, lats, res_target):
    """Converts arrays of (lon, lat) quadcell centroids to qids.

    Parameters
    ----------
        lons : ndarray
            1d ndarray of quadcell centroid longitudes in decimal degrees.
        lats : ndarray
            1d ndarray of quadcell centroid latitudes in decimal degrees.
        res_target : float or int
            Target resolution.

    Returns
    -------
        qids : ndarray
            1d ndarray of qids.
    """

    # Determine resolution of top-level quadrants, given target resolution
    i_max = int(np.ceil(np.log2(180/res_target)))
    res = res_target * 2**i_max

    i = 0
    qids = np.zeros_like(lons, dtype=np.int64)
    origin_lons = np.zeros_like(lons)
    origin_lats = np.zeros_like(lats)

    while res >= res_target:
        delta = 4**(i_max-i)

        # Define disjoint masks by quadrant for all points in the arrays
        mask0 = (lons >= origin_lons) & (lats >= origin_lats)
        mask1 = (lons < origin_lons) & (lats >= origin_lats)
        mask2 = (lons < origin_lons) & (lats < origin_lats)
        mask3 = (lons >= origin_lons) & (lats < origin_lats)

        # Shift origins by quadrant
        # Do nothing to the qids for mask0
        origin_lons[mask0] = origin_lons[mask0] + res/2
        origin_lats[mask0] = origin_lats[mask0] + res/2

        qids[mask1] = qids[mask1] + delta
        origin_lons[mask1] = origin_lons[mask1] - res/2
        origin_lats[mask1] = origin_lats[mask1] + res/2

        qids[mask2] = qids[mask2] + 2*delta
        origin_lons[mask2] = origin_lons[mask2] - res/2
        origin_lats[mask2] = origin_lats[mask2] - res/2

        qids[mask3] = qids[mask3] + 3*delta
        origin_lons[mask3] = origin_lons[mask3] + res/2
        origin_lats[mask3] = origin_lats[mask3] - res/2

        # Halve the resolution and update counter
        res /= 2
        i += 1

    return qids
 

def qids2lls(qids, res_target):
    """Converts arrays of qids to (lon, lat) arrays of quadcell centroids.

    Parameters
    ----------
        qids : ndarray
            1d ndarray of qids.
        res_target : float or int
            Target resolution.

    Returns
    -------
        lons : ndarray
            1d ndarray of quadcell centroid longitudes in decimal degrees.
        lats : ndarray
            1d ndarray of quadcell centroid latitudes in decimal degrees.
    """

    qids0 = qids*1
    lons = np.zeros_like(qids, dtype=np.float64)
    lats = np.zeros_like(qids, dtype=np.float64)
    delta = res_target/2

    while delta <= 180:
        qids, mods = np.divmod(qids, 4)

        mask0 = mods == 0
        mask1 = mods == 1
        mask2 = mods == 2
        mask3 = mods == 3

        lons[mask0] += delta
        lats[mask0] += delta
        lons[mask1] -= delta
        lats[mask1] += delta
        lons[mask2] -= delta
        lats[mask2] -= delta
        lons[mask3] += delta
        lats[mask3] -= delta

        delta *= 2

    return lons, lats
