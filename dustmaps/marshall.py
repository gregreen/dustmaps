#!/usr/bin/env python
#
# iphas.py
# Reads the Marshall et al. (2006) 2MASS-based dust extinction map.
#
# Copyright (C) 2016  Gregory M. Green
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
from scipy.spatial import cKDTree as KDTree
import h5py
import os

import astropy.coordinates as coordinates
import astropy.units as units

from .std_paths import *
from .map_base import DustMap, ensure_flat_galactic
from . import fetch_utils


class MarshallQuery(DustMap):
    """
    Galactic-plane 3D dust map of Marshall et al. (2006), based on 2MASS
    photometry.
    """

    def __init__(self, map_fname=None):
        """
        Args:
            map_fname (Optional[str]): Filename at which the map is stored.
                Defaults to `None`, meaning that the default filename is used.
        """
        if map_fname is None:
            map_fname = os.path.join(data_dir(), 'marshall', 'marshall.h5')

        with h5py.File(map_fname, 'r') as f:
            self._l = f['l'][:]
            self._b = f['b'][:]
            self._A = f['A'][:]
            self._sigma_A = f['sigma_A'][:]
            self._dist = f['dist'][:]
            self._sigma_dist = f['sigma_dist'][:]

        # self._l.shape = (self._l.size,)
        # self._b.shape = (self._b.size,)
        # self._A.shape = (self._A.shape[0], self._A.shape[1]*self._A.shape[2])

        # Shape of the (l,b)-grid
        self._shape = self._l.shape

        # Number of distance bins in each sightline
        self._n_dists = np.sum(np.isfinite(self._dist), axis=0)

        # idx = ~np.isfinite(self._dist)
        # if np.any(idx):
        #     self._dist[idx] = np.inf

        self._l_bounds = (-100., -100.) # min,max Galactic longitude, in deg
        self._b_bounds = (-10., 10.)    # min,max Galactic latitude, in deg
        self._inv_pix_scale = 4.        # 1 / (pixel scale, in deg)

    def _gal2idx(self, gal):
        """
        Converts from Galactic coordinates to pixel indices.

        Args:
            gal (``astropy.coordinates.SkyCoord``): Galactic coordinates. Must
                store an array of coordinates (i.e., not be scalar).

        Returns:
            ``j, k, mask`` - Pixel indices of the coordinates, as well as a mask
            of in-bounds coordinates. Outputs have the same shape as the input
            coordinates.
        """

        # Make sure that l is in domain [-180 deg, 180 deg)
        l = coordinates.Longitude(gal.l, wrap_angle=180.*units.deg)

        j = (self._inv_pix_scale * (l.deg - self._l_bounds[0])).astype('i4')
        k = (self._inv_pix_scale * (gal.b.deg - self._b_bounds[0])).astype('i4')

        idx = (j < 0) | (j >= self._shape[0]) | (k < 0) | (k >= self._shape[1])

        if np.any(idx):
            j[idx] = -1
            k[idx] = -1

        return j, k, ~idx

    @ensure_flat_galactic
    def query(self, coords, return_sigma=False):
        """
        Returns 2MASS Ks-band extinction at the given coordinates.

        Args:
            coords (`astropy.coordinates.SkyCoord`): The coordinates to query.
                Must contain distances.
            return_sigma (Optional[bool]): If True, return the uncertainty in
                extinction as well. Defaults to False.

        Returns:
            Extinction at the specified coordinates, in mags of 2MASS Ks-band
            extinction. If ``return_sigma`` is ``True``, then the uncertainty
            in reddening is also returned, so that the output is
            ``(A, sigma_A)``, where both ``A`` and ``sigma_A`` have the same
            shape as the input coordinates.
        """

        # Ensure that distance has been requested
        has_dist = hasattr(coords.distance, 'kpc')
        if not has_dist:
            raise ValueError('Input `coords` must specify distance.')

        # Convert coordinates to pixel indices
        j, k, mask_idx = self._gal2idx(coords)

        # Which distances to extract
        d = coords.distance.kpc
        dist_idx_ceil = np.sum(d > self._dist[:, j, k], axis=0)

        # Initialize return arrays
        A_ret = np.full(coords.shape, np.nan, dtype='f4')
        if return_sigma:
            sigma_ret = np.full(coords.shape, np.nan, dtype='f4')

        # d < d(nearest distance slice)
        idx_near = (dist_idx_ceil == 0) & mask_idx

        if np.any(idx_near):
            a = d[idx_near] / self._dist[0, j[idx_near], k[idx_near]]
            A_ret[idx_near] = a * self._A[0, j[idx_near], k[idx_near]]

            if return_sigma:
                sigma_ret[idx_near] = self._sigma_A[0, j[idx_near], k[idx_near]]

        # d > d(farthest distance slice)
        idx_far = (dist_idx_ceil == self._n_dists[j,k]) & mask_idx

        if np.any(idx_far):
            A_ret[idx_far] = (
                self._A[self._n_dists[j[idx_far],k[idx_far]]-1, j[idx_far], k[idx_far]])

            if return_sigma:
                sigma_ret[idx_far] = (
                    self._sigma_A[self._n_dists[j[idx_far],k[idx_far]]-1, j[idx_far], k[idx_far]])

        # d(nearest distance slice) < d < d(farthest distance slice)
        idx_btw = (~idx_near & ~idx_far) & mask_idx

        if np.any(idx_btw):
            d_ceil = self._dist[dist_idx_ceil[idx_btw], j[idx_btw], k[idx_btw]]
            d_floor = (
                self._dist[dist_idx_ceil[idx_btw]-1, j[idx_btw], k[idx_btw]])
            a = (d_ceil - d[idx_btw]) / (d_ceil - d_floor)

            A_ret[idx_btw] = (
                (1.-a) * self._A[dist_idx_ceil[idx_btw], j[idx_btw], k[idx_btw]]
                + a * self._A[dist_idx_ceil[idx_btw]-1, j[idx_btw], k[idx_btw]]
            )

            if return_sigma:
                w0 = (1.-a)**2
                w1 = a**2
                norm = 1. / (w0 + w1)
                w0 *= norm
                w1 *= norm
                sigma_ret[idx_btw] = np.sqrt(
                    w0 * self._sigma_A[dist_idx_ceil[idx_btw], j[idx_btw], k[idx_btw]]**2
                    + w1 * self._sigma_A[dist_idx_ceil[idx_btw]-1, j[idx_btw], k[idx_btw]]**2
                )

        if return_sigma:
            return A_ret, sigma_ret

        return A_ret


def fetch():
    """
    Downloads the Marshall et al. (2006) dust map, which is based on 2MASS
    stellar photometry.
    """

    raise NotImplementedError('Automatic downloading of the Marshal et al. '
                              '(2006) map not yet implemented.')
