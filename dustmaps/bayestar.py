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
import astropy.units as units
import h5py
import healpy as hp

from .std_paths import *
from .map_base import DustMap, ensure_flat_galactic
from . import fetch_utils


def lb2pix(nside, l, b, nest=True):
    """
    Converts Galactic (l, b) to HEALPix pixel index.

    Args:
        nside (int): The HEALPix `nside` parameter.
        l (float, or array of floats): Galactic longitude, in degrees.
        b (float, or array of floats): Galactic latitude, in degrees.
        nest (Optional[bool]): If `True` (the default), nested pixel ordering
            will be used. If `False`, ring ordering will be used.

    Returns:
        The HEALPix pixel index or indices. Has the same shape as the input `l`
        and `b`.
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


class BayestarQuery(DustMap):
    """
    Queries the Bayestar 3D dust maps, including Green, Schlafly & Finkbeiner
    (2015). The maps cover the Pan-STARRS 1 footprint, Dec > -30 deg, amounting
    to three-quarters of the sky.
    """

    def __init__(self, map_fname=None, max_samples=None):
        """
        Args:
            map_fname (Optional[str]): Filename of the Bayestar map. Defaults to
                `None`, meaning that the default location is used.
            max_samples (Optional[int]): Maximum number of samples of the map to
                load. Use a lower number in order to decrease memory usage.
                Defaults to `None`, meaning that all samples will be loaded.
        """

        if map_fname is None:
            map_fname = os.path.join(data_dir(), 'bayestar', 'bayestar.h5')

        with h5py.File(map_fname, 'r') as f:
            # Load pixel information
            self._pixel_info = f['/pixel_info'][:]
            self._DM_bin_edges = f['/pixel_info'].attrs['DM_bin_edges']
            self._n_distances = len(self._DM_bin_edges)
            self._n_pix = self._pixel_info.size

            # Load reddening, GR diagnostic
            if max_samples == None:
                self._samples = f['/samples'][:]
            else:
                self._samples = f['/samples'][:,:max_samples,:]

            self._n_samples = self._samples.shape[1]
            self._best_fit = f['/best_fit'][:]
            self._GR = f['/GRDiagnostic'][:]

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
        """
        Returns E(B-V) at the requested coordinates. There are several different
        query modes, which handle the probabilistic nature of the map
        differently.

        Args:
            coords (``astropy.coordinates.SkyCoord``): The coordinates to query.
            mode (Optional[str]): Five different query modes are available:
                'random_sample', 'random_sample_per_pix' 'samples', 'median' and
                'mean'. The ``mode`` determines how the output will reflect the
                probabilistic nature of the Bayestar dust maps.

        Returns:
            Reddening at the specified coordinates, in mags of E(B-V). The
            shape of the output depends on the ``mode``, and on whether
            ``coords`` contains distances.

            If ``coords`` does not specify distance(s), then the shape of the
            output begins with `coords.shape`. If ``coords`` does specify
            distance(s), then the shape of the output begins with
            ``coords.shape + ([number of distance bins],)``.

            If ``mode`` is 'random_sample', then at each coordinate/distance, a
            random sample of reddening is given.

            If ``mode`` is 'random_sample_per_pix', then the sample chosen for
            each angular pixel of the map will be consistent. For example, if
            two query coordinates lie in the same map pixel, then the same
            random sample will be chosen from the map for both query
            coordinates.

            If ``mode`` is 'median', then at each coordinate/distance, the
            median reddening is returned.

            If ``mode`` is 'mean', then at each coordinate/distance, the mean
            reddening is returned.

            Finally, if ``mode`` is 'samples', then all at each
            coordinate/distance, all samples are returned.
        """

        # Check that the query mode is supported
        valid_modes = [
            'random_sample',
            'random_sample_per_pix',
            'samples',
            'median',
            'mean']

        if mode not in valid_modes:
            raise ValueError(
                '"{}" is not a valid `mode`. Valid modes are:\n'
                '  {}'.format(mode, valid_modes)
            )

        n_coords_ret = coords.shape[0]

        # Determine if distance has been requested
        has_dist = hasattr(coords.distance, 'kpc')
        d = coords.distance.kpc if has_dist else None

        # Extract the correct angular pixel(s)
        pix_idx = self._find_data_idx(coords.l.deg, coords.b.deg)
        in_bounds_idx = (pix_idx != -1)

        # Extract the correct samples
        if mode == 'random_sample':
            samp_idx = np.random.randint(0, self._n_samples, pix_idx.size)
            n_samp_ret = 1
        elif mode == 'random_sample_per_pix':
            samp_idx = np.random.randint(0, self._n_samples, self._n_pix)[pix_idx]
            n_samp_ret = 1
        else:
            samp_idx = slice(None)
            n_samp_ret = self._n_samples

        # samples = self._samples[pix_idx, samp_idx]
        # samples[pix_idx == -1] = np.nan

        # Extract the correct distance bin (possibly using linear interpolation)
        if has_dist:
            dm = 5. * (np.log10(d) + 2.)
            bin_idx_ceil = np.searchsorted(self._DM_bin_edges, dm)

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
                        * self._samples[pix_idx[idx_near], samp_idx, 0])
                else:
                    # print('idx_near: {} true'.format(np.sum(idx_near)))
                    # print('ret[idx_near].shape = {}'.format(ret[idx_near].shape))
                    # print('self._samples.shape = {}'.format(self._samples.shape))
                    # print('pix_idx[idx_near].shape = {}'.format(pix_idx[idx_near].shape))

                    ret[idx_near] = (
                        a * self._samples[pix_idx[idx_near], samp_idx[idx_near], 0])

            # d > d(farthest distance slice)
            idx_far = (bin_idx_ceil == self._n_distances) & in_bounds_idx
            if np.any(idx_far):
                # print('idx_far: {} true'.format(np.sum(idx_far)))
                # print('pix_idx[idx_far].shape = {}'.format(pix_idx[idx_far].shape))
                # print('ret[idx_far].shape = {}'.format(ret[idx_far].shape))
                # print('self._samples.shape = {}'.format(self._samples.shape))
                if isinstance(samp_idx, slice):
                    ret[idx_far] = self._samples[pix_idx[idx_far], samp_idx, -1]
                else:
                    ret[idx_far] = self._samples[pix_idx[idx_far], samp_idx[idx_far], -1]

            # d(nearest distance slice) < d < d(farthest distance slice)
            idx_btw = ~idx_near & ~idx_far & in_bounds_idx
            if np.any(idx_btw):
                DM_ceil = self._DM_bin_edges[bin_idx_ceil[idx_btw]]
                DM_floor = self._DM_bin_edges[bin_idx_ceil[idx_btw]-1]
                a = (DM_ceil - dm[idx_btw]) / (DM_ceil - DM_floor)
                if isinstance(samp_idx, slice):
                    ret[idx_btw] = (
                        (1.-a[:,None])
                        * self._samples[pix_idx[idx_btw], samp_idx, bin_idx_ceil[idx_btw]]
                        + a[:,None]
                        * self._samples[pix_idx[idx_btw], samp_idx, bin_idx_ceil[idx_btw]-1]
                    )
                else:
                    ret[idx_btw] = (
                        (1.-a) * self._samples[pix_idx[idx_btw], samp_idx[idx_btw], bin_idx_ceil[idx_btw]]
                        +    a * self._samples[pix_idx[idx_btw], samp_idx[idx_btw], bin_idx_ceil[idx_btw]-1]
                    )
        else:
            ret = self._samples[pix_idx, samp_idx, :]
            ret[~in_bounds_idx] = np.nan

        # Reduce the samples in the requested manner
        if mode == 'median':
            ret = np.median(ret, axis=1)
        elif mode == 'mean':
            ret = np.mean(ret, axis=1)
        elif mode == 'samples':
            # Swap sample and distance axes to be consistent with other 3D dust
            # maps. The output shape will be (pixel, distance, sample).
            if not has_dist:
                np.swapaxes(ret, 1, 2)

        return ret

    @property
    def distances(self):
        """
        Returns the distance bins that the map uses. The return type is
        ``astropy.units.Quantity``, which stores unit-full quantities.
        """
        d = 10.**(0.2*self._DM_bin_edges - 2.)
        return d * units.kpc


def fetch():
    """
    Downloads the Bayestar dust map of Green, Schlafly, Finkbeiner et al. (2015).
    """
    doi = '10.7910/DVN/40C44C'
    requirements = {'contentType': 'application/x-hdf'}
    local_fname = os.path.join(data_dir(), 'bayestar', 'bayestar.h5')
    fetch_utils.dataverse_download_doi(
        doi,
        local_fname,
        file_requirements=requirements)
