#!/usr/bin/env python
#
# cartesian_map.py
# Implements a class for querying dust maps that are stored in an 
# Equirectangular projection.
#
# Copyright (C) 2019  Gregory M. Green
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from __future__ import print_function, division

import numpy as np
import astropy.coordinates as coordinates
import astropy.units as units
from scipy.spatial import cKDTree as KDTree
from astropy.coordinates import Longitude

from .map_base import DustMap, ensure_flat_coords


class EquirectangularDustMap(DustMap):
    """
    A class for querying dust maps stored in an Equirectangular
    projection. The maps may optionally include distances as well.
    """

    def __init__(self, pix_values, lon0, lon1, lat0, lat1,
                       dist0=None, dist1=None,
                       axis_order=('lon','lat','dist'),
                       frame='galactic',
                       dist_interp='linear'):
        """
        Args:
            pix_values (:obj:`np.ndarray`): An array containing the pixel
                values. For a 3D dust map (with distance), the array should be
                3D, with the order of the axes corresponding to
                :obj:`axis_order`. For a 2D dust map, the array should be 2D.
            lon0 (float): The lower limiting longitude of the map, in deg.
            lon1 (float): The upper limiting longitude of the map, in deg.
            lat0 (float): The lower limiting latitude of the map, in deg.
            lat1 (float): The upper limiting latitude of the map, in deg.
            dist0 (Optional[:obj:`astropy.units.Quantity`]): The lower
                limiting distance of the map. If :obj:`dist0` has units of
                distance of the map. If :obj:`dist0` has units of
                :obj:`'mag'`, then it is assumed to be a distance modulus,
                and the distance bins are assumed to be spaced linearly
                in distance modulus, instead of distance.
            dist1 (Optional[:obj:`astropy.units.Quantity`]): The upper
                limiting distance of the map. If :obj:`dist1` has units of
                :obj:`'mag'`, then it is assumed to be a distance modulus,
                and the distance bins are assumed to be spaced linearly
                in distance modulus, instead of distance.
            axis_order (tuple of str): The order of the axes in 
                :obj:`pix_values`. Defaults to :obj:`('lon','lat','dist')`.
                For 2D maps, do not include :obj:`'dist'`.
            frame (Optional[:obj:`str`]): The coordinate frame to which the
                longitudes and latitudes correspond. Must be a frame
                understood by :obj:`astropy.coordinates.SkyCoord`.
                Defaults to :obj:`'galactic'`.
            dist_interp (:obj:`str`): How to interpolate between distance
                slices in the map. Valid choices are :obj:`'step'` and
                :obj:`'linear'`. Defaults to :obj:`'linear'`.
        """
        self._frame = frame

        # Read (lon, lat, dist) bounds
        self._wrap_angle = lon0 * units.deg
        self._lon_lim = (
            Longitude(lon0, unit='deg', wrap_angle=self._wrap_angle),
            Longitude(lon1, unit='deg', wrap_angle=self._wrap_angle),
        )

        self._lat_lim = (lat0, lat1)
        
        if dist0 is not None:
            if dist0.unit == 'mag':
                self._dm = True
                self._dist_lim = (dist0.value, dist1.value)
            else:
                self._dm = False
                self._dist_lim = (
                    dist0.to('kpc').value,
                    dist1.to('kpc').value
                )

        # Read mapping from axis -> (lon, lat, dist)
        self._axis_dist = None

        if len(axis_order) < 2:
            raise ValueError("axis_order must have at least " +
                             "two entries ('lon' and 'lat').")

        for a in ('lon', 'lat'):
            if a not in axis_order:
                raise ValueError("'{}' must be in axis_order.".format(a))

        for i,a in enumerate(axis_order):
            if a == 'lon':
                self._axis_lon = i
            elif a == 'lat':
                self._axis_lat = i
            elif a == 'dist':
                self._axis_dist = i
            else:
                raise ValueError("Unknown entry in axis_order: " +
                                 "'{}'".format(a))
        
        # Get (lon,lat,dist) grid properties
        self._n_lon = pix_values.shape[self._axis_lon]
        self._n_lat = pix_values.shape[self._axis_lat]

        if self._lon_lim[1] == self._lon_lim[0]:
            self._dlon = 360. / self._n_lon
        else:
            self._dlon = (self._lon_lim[1] - self._lon_lim[0]) / self._n_lon
            self._dlon = self._dlon.to('deg').value
        self._dlat = (self._lat_lim[1] - self._lat_lim[0]) / self._n_lat

        if self._axis_dist is not None:
            self._n_dist = pix_values.shape[self._axis_dist]
            self._ddist = (self._dist_lim[1] - self._dist_lim[0]) / self._n_dist
            
            if self._dm:
                # Physical distance of each slice (in kpc)
                self._ddist_phys = 10.**(
                    0.2 * np.linspace(
                        self._dist_lim[0],
                        self._dist_lim[1],
                        self._n_dist
                    ) - 2.
                )

        # Get pixel values
        self._pix_values = pix_values
        self._n_axes = len(pix_values.shape)

        self._dist_interp = dist_interp

    def _coords2idx(self, coords, diff=False):
        c = coords.transform_to(self._frame).represent_as('spherical')

        idx = np.empty(coords.shape + (self._n_axes,), dtype='i4')

        lon = Longitude(c.lon, wrap_angle=self._wrap_angle)
        lon_idx = (lon - self._lon_lim[0]).to('deg').value / self._dlon
        lat_idx = (c.lat.deg - self._lat_lim[0]) / self._dlat

        lon_idx = np.floor(lon_idx).astype('i4')
        lat_idx = np.floor(lat_idx).astype('i4')

        mask = (
            (lon_idx < 0) | (lon_idx >= self._n_lon) |
            (lat_idx < 0) | (lat_idx >= self._n_lat)
        )

        if self._axis_dist is not None:
            dist = c.distance.to('kpc').value

            if self._dm:
                dist = 5. * (np.log10(dist) + 2.)

            dist_idx = (dist - self._dist_lim[0]) / self._ddist
            dist_idx_close = np.floor(dist_idx).astype('i4')
            
            # Differential extinction
            if diff:
                dist_idx_far = (dist_idx_close + 1).clip(0, self._n_dist-1)
                close_mask = (dist_idx_close < 0)
                far_mask = (dist_idx_close >= self._n_dist)
            elif self._dist_interp == 'step':
                dist_idx_close = dist_idx_close.clip(-1, self._n_dist-1)
                close_mask = (dist_idx_close == -1)
            elif self._dist_interp == 'linear':
                far_weight = dist_idx - dist_idx_close

                # Beyond farthest distance slice
                far_mask = (dist_idx_close >= self._n_dist-1)
                if np.any(far_mask):
                    dist_idx_close[far_mask] = self._n_dist - 2
                    far_weight[far_mask] = 1.0

                # Closer than closest distance slice
                close_mask = (dist_idx_close < 0)
                dist_idx_close[close_mask] = -1
                if np.any(close_mask):
                    if self._dm:
                        far_weight[close_mask] = (
                            10.**(0.2 * (dist[close_mask]-self._dist_lim[0]))
                        )
                    else:
                        far_weight[close_mask] = (
                            dist[close_mask] / self._dist_lim[0]
                        )

            #mask |= (dist_idx < 0) | (dist_idx >= self._n_dist)

        
        if np.any(mask):
            lon_idx[mask] = -1
            lat_idx[mask] = -1

            if self._axis_dist is not None:
                dist_idx_close[mask] = -1

        idx[:,self._axis_lon] = lon_idx
        idx[:,self._axis_lat] = lat_idx

        # Include distance indices
        if self._axis_dist is not None:
            if diff:
                # For differential reddening, need:
                #   1. Index of distance bin just beyond d
                #   2. Mask of out-of-bounds (lon, lat)
                #   3. Mask of which pixels are closer than closest slice
                #   3. Mask of which pixels are farther than farthest slice
                idx[:,self._axis_dist] = dist_idx_far
                return idx, mask, (close_mask, far_mask)
            else:
                idx[:,self._axis_dist] = dist_idx_close

            if self._dist_interp == 'step':
                # For step-function interpolation, need:
                #   1. Index of distance bin just inside d
                #   2. Mask of out-of-bounds (lon, lat)
                #   3. Mask of which pixels are closer than closest slice
                return idx, mask, close_mask
            elif self._dist_interp == 'linear':
                # For piecewise-linear interpolation, need:
                #   1. Index of distance bin just inside d
                #   2. Mask of out-of-bounds (lon, lat)
                #   3. For each pixel, weight to apply to next distance bin
                return idx, mask, far_weight
        
        return idx, mask, None

    @ensure_flat_coords
    def query(self, coords, diff=False):
        """
        Returns the cumulative reddening or reddening density (in mag/kpc)
        at the given coordinates.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): Coordinates at which
                to query the reddening. If the map is 3D, these coordinates
                must include distance. If the map is 2D, then the distance
                of the coordinates will be ignored.
            diff (bool): If :obj:`False` (the default), then the cumulative
                reddening is returned. If :obj:`True`, then the reddening
                density (in mag/kpc) is returned. This parameter is ignored
                for 2D maps.

        Returns:
            Either the cumulative reddening or reddening density, as a
            numpy array or float, with the same shape as the input
            :obj:`coords`.
        """
        idx,mask,dist_info = self._coords2idx(coords, diff=diff)
        idx_tuple = tuple([idx[...,i] for i in range(idx.shape[-1])])

        v = self._pix_values[idx_tuple]

        if self._axis_dist is not None:
            if diff:
                if self._dm:
                    # Convert distance modulus of closest slice to
                    # physical distance
                    d0 = 10.**(0.2 * self._dist_lim[0] - 2.)
                else:
                    d0 = self._dist_lim[0]
                
                close_mask, far_mask = dist_info

                if np.any(close_mask):
                    # Inside closest distance slice, interpolate linearly
                    # between reddening = 0 at distance = 0 and the reddening
                    # of the first distance slice.
                    v[close_mask] /= d0

                if np.any(far_mask):
                    # Beyond farthest distance bin, no differential reddening.
                    v[far_mask] = 0.

                middle_mask = ~close_mask & ~far_mask

                if np.any(middle_mask):
                    # Subtract reddening in previous bin, and divide
                    # by width of distance bins.
                    if self._dm:
                        ddist = (
                            self._ddist_phys[idx_tuple[self._axis_dist]] -
                            self._ddist_phys[idx_tuple[self._axis_dist]-1]
                        )[middle_mask]
                    else:
                        ddist = self._ddist

                    idx_tuple[self._axis_dist][:] -= 1
                    v[middle_mask] -= self._pix_values[idx_tuple][middle_mask]
                    v[middle_mask] /= ddist

                v *= units.mag / units.kpc
            elif self._dist_interp == 'step':
                close_mask = dist_info
                if np.any(close_mask):
                    # Zero reddening nearer than closest distance slice.
                    v[close_mask] = 0.
            elif self._dist_interp == 'linear':
                close_mask = (idx_tuple[self._axis_dist] == -1)
                if np.any(close_mask):
                    # Insert reddening = 0 slice at distance = 0
                    v[close_mask] = 0.
                
                # Weight near slice and add in weighted far slice
                v *= (1. - dist_info)
                idx_tuple[self._axis_dist][:] += 1
                v += dist_info * self._pix_values[idx_tuple]

        if np.any(mask):
            # Set reddening to NaN for out-of-bounds (lon, lat)
            v[mask] = np.nan

        return v
