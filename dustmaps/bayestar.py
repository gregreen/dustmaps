#!/usr/bin/env python
#
# bayestar.py
# Reads the Bayestar dust reddening maps, described in
# Green, Schlafly, Finkbeiner et al. (2015, 2018).
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


class BayestarQuery(DustMap):
    """
    Queries the Bayestar 3D dust maps (Green, Schlafly, Finkbeiner et al. 2015,
    2018). The maps cover the Pan-STARRS 1 footprint (dec > -30 deg) amounting
    to three-quarters of the sky.
    """

    def __init__(self, map_fname=None, max_samples=None, version='bayestar2019'):
        """
        Args:
            map_fname (Optional[:obj:`str`]): Filename of the Bayestar map. Defaults to
                :obj:`None`, meaning that the default location is used.
            max_samples (Optional[:obj:`int`]): Maximum number of samples of the map to
                load. Use a lower number in order to decrease memory usage.
                Defaults to :obj:`None`, meaning that all samples will be loaded.
            version (Optional[:obj:`str`]): The map version to download. Valid versions
                are :obj:`'bayestar2019'` (Green, Schlafly, Finkbeiner et al. 2019),
                :obj:`'bayestar2017'` (Green, Schlafly, Finkbeiner et al. 2018)
                and :obj:`'bayestar2015'` (Green, Schlafly, Finkbeiner et al. 2015).
                Defaults to :obj:`'bayestar2015'`.
        """

        if map_fname is None:
            map_fname = os.path.join(data_dir(), 'bayestar', '{}.h5'.format(version))

        t_start = time()
        
        with h5py.File(map_fname, 'r') as f:
            # Load pixel information
            print('Loading pixel_info ...')
            self._pixel_info = f['/pixel_info'][:]
            self._DM_bin_edges = f['/pixel_info'].attrs['DM_bin_edges']
            self._n_distances = len(self._DM_bin_edges)
            self._n_pix = self._pixel_info.size
            
            t_pix_info = time()

            # Load reddening
            print('Loading samples ...')
            if max_samples == None:
                self._samples = f['/samples'][:]
            else:
                self._samples = f['/samples'][:,:max_samples,:]
            
            t_samples = time()

            self._n_samples = self._samples.shape[1]
            print('Loading best_fit ...')
            self._best_fit = f['/best_fit'][:]
            
            t_best = time()

        # Reshape best fit
        s = self._best_fit.shape
        self._best_fit.shape = (s[0], 1, s[1])  # (pixels, samples=1, distances)

        # Replace NaNs in reliable distance estimates with +-infinity
        print('Replacing NaNs in reliable distance estimates ...')
        for k,v in [('DM_reliable_min',np.inf), ('DM_reliable_max',-np.inf)]:
            idx = ~np.isfinite(self._pixel_info[k])
            self._pixel_info[k][idx] = v

        t_nan = time()
        
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
        print('   samples: {: >7.3f} s'.format(t_samples-t_pix_info))
        print('      best: {: >7.3f} s'.format(t_best-t_samples))
        print('       nan: {: >7.3f} s'.format(t_nan-t_best))
        print('      sort: {: >7.3f} s'.format(t_sort-t_nan))
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

    def _raise_on_mode(self, mode):
        """
        Checks that the provided query mode is one of the accepted values. If
        not, raises a :obj:`ValueError`.
        """
        valid_modes = [
            'random_sample',
            'random_sample_per_pix',
            'samples',
            'median',
            'mean',
            'best',
            'percentile']

        if mode not in valid_modes:
            raise ValueError(
                '"{}" is not a valid `mode`. Valid modes are:\n'
                '  {}'.format(mode, valid_modes)
            )

    def _interpret_percentile(self, mode, pct):
        if mode == 'percentile':
            if pct is None:
                raise ValueError(
                    '"percentile" mode requires an additional keyword '
                    'argument: "pct"')
            if (type(pct) in (list,tuple)) or isinstance(pct, np.ndarray):
                try:
                    pct = np.array(pct, dtype='f8')
                except ValueError as err:
                    raise ValueError(
                        'Invalid "pct" specification. Must be number or '
                        'list/array of numbers.')
                if np.any((pct < 0) | (pct > 100)):
                    raise ValueError('"pct" must be between 0 and 100.')
                scalar_pct = False
            else:
                try:
                    pct = float(pct)
                except ValueError as err:
                    raise ValueError(
                        'Invalid "pct" specification. Must be number or '
                        'list/array of numbers.')
                if (pct < 0) or (pct > 100):
                    raise ValueError('"pct" must be between 0 and 100.')
                scalar_pct = True

            return pct, scalar_pct
        else:
            return None, None

    def get_query_size(self, coords, mode='random_sample',
                       return_flags=False, pct=None):
        # Check that the query mode is supported
        self._raise_on_mode(mode)

        # Validate percentile specification
        pct, scalar_pct = self._interpret_percentile(mode, pct)

        n_coords = np.prod(coords.shape, dtype=int)

        if mode == 'samples':
            n_samples = self._n_samples
        elif mode == 'percentile':
            if scalar_pct:
                n_samples = 1
            else:
                n_samples = len(pct)
        else:
            n_samples = 1

        if hasattr(coords.distance, 'kpc'):
            n_dists = 1
        else:
            n_dists = self._n_distances

        return n_coords * n_samples * n_dists

    @ensure_flat_galactic
    def query(self, coords, mode='random_sample', return_flags=False, pct=None):
        """
        Returns reddening at the requested coordinates. There are several
        different query modes, which handle the probabilistic nature of the map
        differently.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.
            mode (Optional[:obj:`str`]): Seven different query modes are available:
                'random_sample', 'random_sample_per_pix' 'samples', 'median',
                'mean', 'best' and 'percentile'. The :obj:`mode` determines how the
                output will reflect the probabilistic nature of the Bayestar
                dust maps.
            return_flags (Optional[:obj:`bool`]): If :obj:`True`, then QA flags will be
                returned in a second numpy structured array. That is, the query
                will return :obj:`ret`, :obj:'flags`, where :obj:`ret` is the normal return
                value, containing reddening. Defaults to :obj:`False`.
            pct (Optional[:obj:`float` or list/array of :obj:`float`]): If the mode is
                :obj:`percentile`, then :obj:`pct` specifies which percentile(s) is
                (are) returned.

        Returns:
            Reddening at the specified coordinates, in magnitudes of reddening.

            The conversion to E(B-V) (or other reddening units) depends on
            whether :obj:`version='bayestar2019'` (the default), :obj:`'bayestar2017'`
            or :obj:`'bayestar2015'` was selected when the :obj:`BayestarQuery` object
            was created. To convert Bayestar2019 to Pan-STARRS 1 extinctions,
            multiply by the coefficients given in Table 1 of Green et al. (2019).
            For Bayestar2017, use the coefficients given in Table 1 of Green et al.
            (2018). Conversion to extinction in non-PS1 passbands depends on the
            choice of extinction law. To convert Bayestar2015 to extinction in
            various passbands, multiply by the coefficients in Table 6 of
            Schlafly & Finkbeiner (2011). See Green et al. (2015, 2018) for more
            detailed discussion of how to convert the Bayestar dust maps into
            reddenings or extinctions in different passbands.

            The shape of the output depends on the :obj:`mode`, and on whether
            :obj:`coords` contains distances.

            If :obj:`coords` does not specify distance(s), then the shape of the
            output begins with :obj:`coords.shape`. If :obj:`coords` does specify
            distance(s), then the shape of the output begins with
            :obj:`coords.shape + ([number of distance bins],)`.

            If :obj:`mode` is :obj:`'random_sample'`, then at each
            coordinate/distance, a random sample of reddening is given.

            If :obj:`mode` is :obj:`'random_sample_per_pix'`, then the sample chosen
            for each angular pixel of the map will be consistent. For example,
            if two query coordinates lie in the same map pixel, then the same
            random sample will be chosen from the map for both query
            coordinates.

            If :obj:`mode` is :obj:`'median'`, then at each coordinate/distance, the
            median reddening is returned.

            If :obj:`mode` is :obj:`'mean'`, then at each coordinate/distance, the
            mean reddening is returned.

            If :obj:`mode` is :obj:`'best'`, then at each coordinate/distance, the
            maximum posterior density reddening is returned (the "best fit").

            If :obj:`mode` is :obj:`'percentile'`, then an additional keyword
            argument, :obj:`pct`, must be specified. At each coordinate/distance,
            the requested percentiles (in :obj:`pct`) will be returned. If :obj:`pct`
            is a list/array, then the last axis of the output will correspond to
            different percentiles.

            Finally, if :obj:`mode` is :obj:`'samples'`, then at each
            coordinate/distance, all samples are returned. The last axis of the
            output will correspond to different samples.

            If :obj:`return_flags` is :obj:`True`, then in addition to reddening, a
            structured array containing QA flags will be returned. If the input
            coordinates include distances, the QA flags will be :obj:`"converged"`
            (whether or not the line-of-sight fit converged in a given pixel)
            and :obj:`"reliable_dist"` (whether or not the requested distance is
            within the range considered reliable, based on the inferred
            stellar distances). If the input coordinates do not include
            distances, then instead of :obj:`"reliable_dist"`, the flags will
            include :obj:`"min_reliable_distmod"` and :obj:`"max_reliable_distmod"`,
            the minimum and maximum reliable distance moduli in the given pixel.
        """

        # Check that the query mode is supported
        self._raise_on_mode(mode)

        # Validate percentile specification
        pct, scalar_pct = self._interpret_percentile(mode, pct)

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

        # Extract the correct samples
        if mode == 'random_sample':
            # A different sample in each queried coordinate
            samp_idx = np.random.randint(0, self._n_samples, pix_idx.size)
            n_samp_ret = 1
        elif mode == 'random_sample_per_pix':
            # Choose same sample in all coordinates that fall in same angular
            # HEALPix pixel
            samp_idx = np.random.randint(0, self._n_samples, self._n_pix)[pix_idx]
            n_samp_ret = 1
        elif mode == 'best':
            samp_idx = slice(None)
            n_samp_ret = 1
        else:
            # Return all samples in each queried coordinate
            samp_idx = slice(None)
            n_samp_ret = self._n_samples

        # t2 = time.time()

        if mode == 'best':
            val = self._best_fit
        else:
            val = self._samples

        # Create empty array to store flags
        if return_flags:
            if has_dist:
                # If distances are provided in query, return only covergence and
                # whether or not this distance is reliable
                dtype = [('converged', 'bool'),
                         ('reliable_dist', 'bool')]
                # shape = (n_coords_ret)
            else:
                # Return convergence and reliable distance ranges
                dtype = [('converged', 'bool'),
                         ('min_reliable_distmod', 'f4'),
                         ('max_reliable_distmod', 'f4')]
            flags = np.empty(n_coords_ret, dtype=dtype)
        # samples = self._samples[pix_idx, samp_idx]
        # samples[pix_idx == -1] = np.nan

        # t3 = time.time()

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

            # Flag: distance in reliable range?
            if return_flags:
                dm_min = self._pixel_info['DM_reliable_min'][pix_idx]
                dm_max = self._pixel_info['DM_reliable_max'][pix_idx]
                flags['reliable_dist'] = (
                    (dm >= dm_min) &
                    (dm <= dm_max) &
                    np.isfinite(dm_min) &
                    np.isfinite(dm_max))
                flags['reliable_dist'][~in_bounds_idx] = False
        else:   # No distances provided
            ret = val[pix_idx, samp_idx, :]   # Return all distances
            ret[~in_bounds_idx] = np.nan

            # Flag: reliable distance bounds
            if return_flags:
                dm_min = self._pixel_info['DM_reliable_min'][pix_idx]
                dm_max = self._pixel_info['DM_reliable_max'][pix_idx]

                flags['min_reliable_distmod'] = dm_min
                flags['max_reliable_distmod'] = dm_max
                flags['min_reliable_distmod'][~in_bounds_idx] = np.nan
                flags['max_reliable_distmod'][~in_bounds_idx] = np.nan

        # t4 = time.time()

        # Flag: convergence
        if return_flags:
            flags['converged'] = (
                self._pixel_info['converged'][pix_idx].astype(np.bool))
            flags['converged'][~in_bounds_idx] = False

        # t5 = time.time()

        # Reduce the samples in the requested manner
        if mode == 'median':
            ret = np.median(ret, axis=1)
        elif mode == 'mean':
            ret = np.mean(ret, axis=1)
        elif mode == 'percentile':
            ret = np.nanpercentile(ret, pct, axis=1)
            if not scalar_pct:
                # (percentile, pixel) -> (pixel, percentile)
                # (pctile, pixel, distance) -> (pixel, distance, pctile)
                ret = np.moveaxis(ret, 0, -1)
        elif mode == 'best':
            # Remove "samples" axis
            s = ret.shape
            ret.shape = s[:1] + s[2:]
        elif mode == 'samples':
            # Swap sample and distance axes to be consistent with other 3D dust
            # maps. The output shape will be (pixel, distance, sample).
            if not has_dist:
                np.swapaxes(ret, 1, 2)

        # t6 = time.time()
        #
        # print('')
        # print('time inside bayestar.query: {:.4f} s'.format(t6-t0))
        # print('{: >7.4f} s : {: >6.4f} s : _find_data_idx'.format(t1-t0, t1-t0))
        # print('{: >7.4f} s : {: >6.4f} s : sample slice spec'.format(t2-t0, t2-t1))
        # print('{: >7.4f} s : {: >6.4f} s : create empty return flag array'.format(t3-t0, t3-t2))
        # print('{: >7.4f} s : {: >6.4f} s : extract results'.format(t4-t0, t4-t3))
        # print('{: >7.4f} s : {: >6.4f} s : convergence flag'.format(t5-t0, t5-t4))
        # print('{: >7.4f} s : {: >6.4f} s : reduce'.format(t6-t0, t6-t5))
        # print('')

        if return_flags:
            return ret, flags

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


def fetch(version='bayestar2019'):
    """
    Downloads the specified version of the Bayestar dust map.

    Args:
        version (Optional[:obj:`str`]): The map version to download. Valid versions are
            :obj:`'bayestar2019'` (Green, Schlafly, Finkbeiner et al. 2019),
            :obj:`'bayestar2017'` (Green, Schlafly, Finkbeiner et al. 2018) and
            :obj:`'bayestar2015'` (Green, Schlafly, Finkbeiner et al. 2015). Defaults
            to :obj:`'bayestar2019'`.

    Raises:
        :obj:`ValueError`: The requested version of the map does not exist.

        :obj:`DownloadError`: Either no matching file was found under the given DOI, or
            the MD5 sum of the file was not as expected.

        :obj:`requests.exceptions.HTTPError`: The given DOI does not exist, or there
            was a problem connecting to the Dataverse.
    """

    doi = {
        'bayestar2015': '10.7910/DVN/40C44C',
        'bayestar2017': '10.7910/DVN/LCYHJG',
        'bayestar2019': '10.7910/DVN/2EJ9TX'
    }

    # Raise an error if the specified version of the map does not exist
    try:
        doi = doi[version]
    except KeyError as err:
        raise ValueError('Version "{}" does not exist. Valid versions are: {}'.format(
            version,
            ', '.join(['"{}"'.format(k) for k in doi.keys()])
        ))

    requirements = {
        'bayestar2015': {'contentType': 'application/x-hdf'},
        'bayestar2017': {'filename': 'bayestar2017.h5'},
        'bayestar2019': {'filename': 'bayestar2019.h5'}
    }[version]

    local_fname = os.path.join(data_dir(), 'bayestar', '{}.h5'.format(version))

    # Download the data
    fetch_utils.dataverse_download_doi(
        doi,
        local_fname,
        file_requirements=requirements)


class BayestarWebQuery(WebDustMap):
    """
    Remote query over the web for the Bayestar 3D dust maps (Green,
    Schlafly, Finkbeiner et al. 2015, 2018, 2019). The maps cover the
    Pan-STARRS 1 footprint (dec > -30 deg) amounting to three-quarters of
    the sky.

    This query object does not require a local version of the data, but rather
    an internet connection to contact the web API. The query functions have the
    same inputs and outputs as their counterparts in :obj:`BayestarQuery`.
    """

    def __init__(self, api_url=None, version='bayestar2019'):
        """
        Args:
            version (Optional[:obj:`str`]): The map version to download. Valid versions
                are :obj:`'bayestar2019'` (Green, Schlafly, Finkbeiner et al. 2019),
                :obj:`'bayestar2017'` (Green, Schlafly, Finkbeiner et al. 2018)
                and :obj:`'bayestar2015'` (Green, Schlafly, Finkbeiner et al. 2015).
                Defaults to :obj:`'bayestar2019'`.
        """
        super(BayestarWebQuery, self).__init__(
            api_url=api_url,
            map_name=version)
