#!/usr/bin/env python
#
# drimmel.py
# Reads the Drimmel dust reddening maps, described in
# http://adsabs.harvard.edu/abs/2003A%26A...409..205D.
# But queeried via the https://github.com/jobovy/mwdust interface
#
# Copyright (C) 2016-2019  Gregory M. Green
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

import os
import h5py
import numpy as np

import astropy.coordinates as coordinates
import astropy.units as units
import h5py
import healpy as hp

from .std_paths import *
from .map_base import DustMap, WebDustMap, ensure_flat_galactic
from . import fetch_utils

from time import time


def lb2pix(nside, l, b, nest=True):
    """
    Converts Galactic (l, b) to HEALPix pixel index.

    Args:
        nside (:obj:`int`): The HEALPix :obj:`nside` parameter.
        l (:obj:`float`, or array of :obj:`float`): Galactic longitude, in degrees.
        b (:obj:`float`, or array of :obj:`float`): Galactic latitude, in degrees.
        nest (Optional[:obj:`bool`]): If :obj:`True` (the default), nested pixel ordering
            will be used. If :obj:`False`, ring ordering will be used.

    Returns:
        The HEALPix pixel index or indices. Has the same shape as the input :obj:`l`
        and :obj:`b`.
    """

    theta = np.radians(90. - b)
    phi = np.radians(l)

    if not hasattr(l, '__len__'):
        if (b < -90.) or (b > 90.):
            return -1

        pix_idx = hp.pixelfunc.ang2pix(nside, theta, phi, nest=nest)

        return pix_idx

    idx = (b >= -90.) & (b <= 90.)

    pix_idx = np.empty(l.shape, dtype='i8')
    pix_idx[idx] = hp.pixelfunc.ang2pix(nside, theta[idx], phi[idx], nest=nest)
    pix_idx[~idx] = -1

    return pix_idx


class DrimmelQuery(DustMap):
    """
    Queries the Drimmel 3D dust maps Drimmel+ 2003. The map is all-sky.
    """

    def __init__(self, map_fname=None):
        """
        Args:
            map_fname (Optional[:obj:`str`]): Filename of the Bayestar map. Defaults to
                :obj:`None`, meaning that the default location is used.
        """

        if map_fname is None:
            map_fname = os.path.join(data_dir(), 'drimmel.h5')

        t_start = time()
        
        with h5py.File(map_fname, 'r') as f:
            # Load pixel information
            print('Loading pixel_info ...')
            self._pixel_info = f['/pixel_info'][:]
            self._DM_bin_edges = np.linspace(4,19,31)
            self._n_distances = len(self._DM_bin_edges)
            self._n_pix = self._pixel_info.size
            
            t_pix_info = time()

            # Load reddening
            print('Loading best_fit ...')
            self._best_fit = f['/best_fit'][:]
            
            t_best = time()

        # Reshape best fit
        s = self._best_fit.shape
        self._best_fit.shape = (s[0], 1, s[1])  # (pixels, samples=1, distances)

        # Get healpix indices at each nside level
        print('Sorting pixel_info ...')
        sort_idx = np.argsort(self._pixel_info, order=['nside', 'healpix_index'])
        
        t_sort = time()

        self._nside_levels = np.unique(self._pixel_info['nside'])
        self._hp_idx_sorted = []
        self._data_idx = []

        start_idx = 0

        print('Extracting hp_idx_sorted and data_idx at each nside ...')
        for nside in self._nside_levels:
            print('  nside = {}'.format(nside))
            end_idx = np.searchsorted(self._pixel_info['nside'], nside,
                                      side='right', sorter=sort_idx)

            idx = sort_idx[start_idx:end_idx]

            self._hp_idx_sorted.append(self._pixel_info['healpix_index'][idx])
            self._data_idx.append(idx)

            start_idx = end_idx
        
        t_finish = time()
        
        print('t = {:.3f} s'.format(t_finish - t_start))
        print('  pix_info: {: >7.3f} s'.format(t_pix_info-t_start))
        print('      best: {: >7.3f} s'.format(t_best-t_pix_info))
        print('      sort: {: >7.3f} s'.format(t_sort-t_best))
        print('       idx: {: >7.3f} s'.format(t_finish-t_sort))

    def _find_data_idx(self, l, b):
        pix_idx = np.empty(l.shape, dtype='i8')
        pix_idx[:] = -1

        # Search at each nside
        for k,nside in enumerate(self._nside_levels):
            ipix = lb2pix(nside, l, b, nest=True)

            # Find the insertion points of the query pixels in the large, ordered pixel list
            idx = np.searchsorted(self._hp_idx_sorted[k], ipix, side='left')

            # Determine which insertion points are beyond the edge of the pixel list
            in_bounds = (idx < self._hp_idx_sorted[k].size)

            if not np.any(in_bounds):
                continue

            # Determine which query pixels are correctly placed
            idx[~in_bounds] = -1
            match_idx = (self._hp_idx_sorted[k][idx] == ipix)
            match_idx[~in_bounds] = False
            idx = idx[match_idx]

            if np.any(match_idx):
                pix_idx[match_idx] = self._data_idx[k][idx]

        return pix_idx

    def get_query_size(self, coords, return_flags=False, pct=None):
        # Check that the query mode is supported

        # Validate percentile specification

        n_coords = np.prod(coords.shape, dtype=int)

        n_samples = self._n_samples

        if hasattr(coords.distance, 'kpc'):
            n_dists = 1
        else:
            n_dists = self._n_distances

        return n_coords * n_samples * n_dists

    @ensure_flat_galactic
    def query(self, coords):
        """
        Returns reddening at the requested coordinates. There are several
        different query modes, which handle the probabilistic nature of the map
        differently.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.
        Returns:
            Reddening at the specified coordinates, in magnitudes of reddening.
        """

        # Get number of coordinates requested
        n_coords_ret = coords.shape[0]

        # Determine if distance has been requested
        has_dist = hasattr(coords.distance, 'kpc')
        d = coords.distance.kpc if has_dist else None

        # Extract the correct angular pixel(s)
        # t0 = time.time()
        pix_idx = self._find_data_idx(coords.l.deg, coords.b.deg)
        in_bounds_idx = (pix_idx != -1)

        # t1 = time.time()

        samp_idx = slice(None)
        n_samp_ret = 1
        
        # t2 = time.time()

        val = self._best_fit


        # Extract the correct distance bin (possibly using linear interpolation)
        if has_dist: # Distance has been provided
            # Determine ceiling bin index for each coordinate
            dm = 5. * (np.log10(d) + 2.)
            bin_idx_ceil = np.searchsorted(self._DM_bin_edges, dm)

            # Create NaN-filled return arrays
            if isinstance(samp_idx, slice):
                ret = np.full((n_coords_ret, n_samp_ret), np.nan, dtype='f4')
            else:
                ret = np.full((n_coords_ret,), np.nan, dtype='f4')

            # d < d(nearest distance slice)
            idx_near = (bin_idx_ceil == 0) & in_bounds_idx
            if np.any(idx_near):
                a = 10.**(0.2 * (dm[idx_near] - self._DM_bin_edges[0]))
                if isinstance(samp_idx, slice):
                    ret[idx_near] = (
                        a[:,None]
                        * val[pix_idx[idx_near], samp_idx, 0])
                else:
                    # print('idx_near: {} true'.format(np.sum(idx_near)))
                    # print('ret[idx_near].shape = {}'.format(ret[idx_near].shape))
                    # print('val.shape = {}'.format(val.shape))
                    # print('pix_idx[idx_near].shape = {}'.format(pix_idx[idx_near].shape))

                    ret[idx_near] = (
                        a * val[pix_idx[idx_near], samp_idx[idx_near], 0])

            # d > d(farthest distance slice)
            idx_far = (bin_idx_ceil == self._n_distances) & in_bounds_idx
            if np.any(idx_far):
                # print('idx_far: {} true'.format(np.sum(idx_far)))
                # print('pix_idx[idx_far].shape = {}'.format(pix_idx[idx_far].shape))
                # print('ret[idx_far].shape = {}'.format(ret[idx_far].shape))
                # print('val.shape = {}'.format(val.shape))
                if isinstance(samp_idx, slice):
                    ret[idx_far] = val[pix_idx[idx_far], samp_idx, -1]
                else:
                    ret[idx_far] = val[pix_idx[idx_far], samp_idx[idx_far], -1]

            # d(nearest distance slice) < d < d(farthest distance slice)
            idx_btw = ~idx_near & ~idx_far & in_bounds_idx
            if np.any(idx_btw):
                DM_ceil = self._DM_bin_edges[bin_idx_ceil[idx_btw]]
                DM_floor = self._DM_bin_edges[bin_idx_ceil[idx_btw]-1]
                a = (DM_ceil - dm[idx_btw]) / (DM_ceil - DM_floor)
                if isinstance(samp_idx, slice):
                    ret[idx_btw] = (
                        (1.-a[:,None])
                        * val[pix_idx[idx_btw], samp_idx, bin_idx_ceil[idx_btw]]
                        + a[:,None]
                        * val[pix_idx[idx_btw], samp_idx, bin_idx_ceil[idx_btw]-1]
                    )
                else:
                    ret[idx_btw] = (
                        (1.-a) * val[pix_idx[idx_btw], samp_idx[idx_btw], bin_idx_ceil[idx_btw]]
                        +    a * val[pix_idx[idx_btw], samp_idx[idx_btw], bin_idx_ceil[idx_btw]-1]
                    )

        else:   # No distances provided
            ret = val[pix_idx, samp_idx, :]   # Return all distances
            ret[~in_bounds_idx] = np.nan

        # t4 = time.time()

        # Reduce the samples in the requested manner
        s = ret.shape
        ret.shape = s[:1] + s[2:]
        return ret

    @property
    def distances(self):
        """
        Returns the distance bin edges that the map uses. The return type is
        :obj:`astropy.units.Quantity`, which stores unit-full quantities.
        """
        d = 10.**(0.2*self._DM_bin_edges - 2.)
        return d * units.kpc

    @property
    def distmods(self):
        """
        Returns the distance modulus bin edges that the map uses. The return
        type is :obj:`astropy.units.Quantity`, with units of mags.
        """
        return self._DM_bin_edges * units.mag


def fetch(clobber = False):
    """
    Downloads the Drimmel+ 2003 dust map.

    Args:
        clobber (Optional[:obj:`bool`]): If ``True``, any existing file will be
            overwritten, even if it appears to match. If ``False`` (the
            default), :obj:`fetch()` will attempt to determine if the dataset
            already exists. This determination is not 100\% robust against data
            corruption.
    """

    dest_dir = fname_pattern = os.path.join(data_dir(), 'drimmel')
    url = 'https://keeper.mpdl.mpg.de/f/90fd101e3e1844d4aa20/?dl=1'
    h5_fname = os.path.join(dest_dir, 'drimmel.h5')
    md5 = 'b0ae9437f1aad695bcc689cd78b76db5'

    # Check if file already exists
    if not clobber:
        if os.path.isfile(h5_fname):
            print('File appears to exist already. Call `fetch(clobber=True)` '
                  'to force overwriting of existing file.')
            return

    # Download the table
    print('Downloading {}'.format(url))
    fetch_utils.download_and_verify(url, md5, fname=h5_fname)

