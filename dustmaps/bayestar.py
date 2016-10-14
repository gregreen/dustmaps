#!/usr/bin/env python
#
# bayestar.py
# Reads the Bayestar dust reddening map, described in
# Green, Schlafly, Finkbeiner et al. (2015).
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

import os
import h5py
import numpy as np

import astropy.coordinates as coordinates
import h5py
import healpy as hp

from std_paths import *
from map_base import DustMap, ensure_flat_galactic
import fetch_utils


def lb2pix(nside, l, b, nest=True):
    '''
    Convert (l, b) to pixel index.
    '''

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

class BayestarQuery(DustMap):
    def __init__(self,
                 map_fname=os.path.join(data_dir(), 'bayestar', 'bayestar.h5'),
                 max_samples=None):
        f = h5py.File(map_fname, 'r')

        # Load pixel information
        self._pixel_info = f['/pixel_info'][:]
        self._DM_bin_edges = f['/pixel_info'].attrs['DM_bin_edges']
        self._n_distances = len(self._DM_bin_edges)

        # Load reddening, GR diagnostic
        if max_samples == None:
            self._samples = f['/samples'][:]
        else:
            self._samples = f['/samples'][:,:max_samples,:]
        self._n_samples = self._samples.shape[1]
        self._best_fit = f['/best_fit'][:]
        self._GR = f['/GRDiagnostic'][:]

        f.close()

        # Remove NaNs from reliable distance estimates
        # for k in ['DM_reliable_min', 'DM_reliable_max']:
        #     idx = ~np.isfinite(self._pixel_info[k])
        #     self._pixel_info[k][idx] = -999.

        # Get healpix indices at each nside level
        sort_idx = np.argsort(self._pixel_info, order=['nside', 'healpix_index'])

        self._nside_levels = np.unique(self._pixel_info['nside'])
        self._hp_idx_sorted = []
        self._data_idx = []

        start_idx = 0

        for nside in self._nside_levels:
            end_idx = np.searchsorted(self._pixel_info['nside'], nside,
                                      side='right', sorter=sort_idx)

            idx = sort_idx[start_idx:end_idx]

            self._hp_idx_sorted.append(self._pixel_info['healpix_index'][idx])
            self._data_idx.append(idx)

            start_idx = end_idx

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

    @ensure_flat_galactic
    def query(self, coords, mode='random_sample'):
        valid_modes = ['random_sample', 'samples', 'median', 'mean']
        if mode not in valid_modes:
            raise ValueError(
                '"{}" is not a valid `mode`. Valid modes are:\n'
                '  {}'.format(mode, valid_modes)
            )

        # gal = coords.transform_to('galactic')
        gal = coords
        l = gal.l.deg
        b = gal.b.deg

        # Determine if distance has been requested
        has_dist = hasattr(gal.distance, 'kpc')
        d = gal.distance.kpc if has_dist else None

        # Ensure that l and b are arrays
        # is_array = hasattr(gal.l.deg, '__len__')

        # if not is_array:
        #     l = np.array([l])
        #     b = np.array([b])
        #     d = np.array([d])

        # Extract the correct angular pixel
        pix_idx = self._find_data_idx(l, b)

        # Extract
        if mode == 'random_sample':
            samp_idx = np.random.randint(0, self._n_samples, pix_idx.size)
        else:
            samp_idx = slice(None)

        samples = self._samples[pix_idx, samp_idx]
        samples[pix_idx == -1] = np.nan

        # Extract the correct distance bin (possibly using linear interpolation)
        if has_dist:
            dm = 5. * (np.log10(d) + 2.)
            bin_idx_ceil = np.searchsorted(self._DM_bin_edges, dm)

            ret = np.zeros(samples.shape[:-1], dtype='f4')

            # d < d(nearest distance slice)
            idx_near = (bin_idx_ceil == 0)
            if np.any(idx_near):
                a = 10.**(0.2 * (dm[idx_near] - self._DM_bin_edges[0]))
                if len(samples.shape) == 2:
                    ret[idx_near] = a * samples[idx_near, 0]
                elif len(samples.shape) == 3:
                    ret[idx_near] = a[:,None] * samples[idx_near, :, 0]
                else:
                    raise ValueError('samples.shape = {:d}'.format(samples.shape))

            # d > d(farthest distance slice)
            idx_far = (bin_idx_ceil == self._n_distances)
            if np.any(idx_far):
                # ret[idx_far] = samples[idx_far, samp_idx, -1]
                if len(samples.shape) == 2:
                    ret[idx_far] = samples[idx_far, -1]
                elif len(samples.shape) == 3:
                    ret[idx_far] = samples[idx_far, :, -1]
                else:
                    raise ValueError('samples.shape = {:d}'.format(samples.shape))

            # d(nearest distance slice) < d < d(farthest distance slice)
            idx_btw = ~idx_near & ~idx_far
            if np.any(idx_btw):
                DM_ceil = self._DM_bin_edges[bin_idx_ceil[idx_btw]]
                DM_floor = self._DM_bin_edges[bin_idx_ceil[idx_btw]-1]
                a = (DM_ceil - dm[idx_btw]) / (DM_ceil - DM_floor)
                if len(samples.shape) == 2:
                    ret[idx_btw] = (
                        (1.-a) * samples[idx_btw, bin_idx_ceil[idx_btw]]
                        +    a * samples[idx_btw, bin_idx_ceil[idx_btw]-1]
                    )
                elif len(samples.shape) == 3:
                    ret[idx_btw] = (
                        (1.-a[:,None]) * samples[idx_btw, :, bin_idx_ceil[idx_btw]]
                        +    a[:,None] * samples[idx_btw, :, bin_idx_ceil[idx_btw]-1]
                    )
                else:
                    raise ValueError('samples.shape = {:d}'.format(samples.shape))
        else:
            ret = samples

        # Reduce the samples in the requested manner
        if mode == 'median':
            ret = np.median(ret, axis=1)
        elif mode == 'mean':
            ret = np.mean(ret, axis=1)

        # Transform back to scalar response if user supplied scalar coordinates
        # if not is_array:
        #     return ret[0]

        return ret


def fetch():
    """
    Download the Bayestar dust map of Green, Schlafly, Finkbeiner et al. (2015).
    """
    doi = '10.7910/DVN/40C44C'
    requirements = {'contentType': 'application/x-hdf'}
    local_fname = os.path.join(std_paths.data_dir(), 'bayestar', 'bayestar.h5')
    fetch_utils.dataverse_download_doi(
        doi,
        local_fname,
        file_requirements=requirements)
