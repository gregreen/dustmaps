#!/usr/bin/env python
#
# bayestar.py
# Reads the Bayestar dust reddening maps, described in
# Green, Schlafly, Finkbeiner et al. (2015, 2018).
#
# Copyright (C) 2016-2019  Gregory M. Green
#
# dustmaps is free software: you can redistribute it and/or modify
# it under the terms of either:
#
# - The GNU General Public License as published by the Free Software Foundation,
#   either version 2 of the License, or (at your option) any later version, or
# - The 2-Clause BSD License (also known as the Simplified BSD License).
#
# You should have received copies of the GNU General Public License
# and the BSD License along with this program.
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


# When a Bayestar map is downloaded, it is repacked into an HDF5 file in which
# each chunk holds a block of adjacent HEALPix pixels, together with all of
# their samples and distance bins. A query for a random coordinate then reads
# only the blocks containing the pixels it asks for, rather than the whole map.
#
# Larger blocks compress slightly better, but a queried pixel has to read and
# decompress the entire block that contains it, so the cost of a query grows in
# proportion to the block size. Sixteen pixels is small enough to keep that
# cost low, while still compressing within a few percent of much larger blocks.
CHUNK_PIXELS = 16

# Deflate compression level used when repacking. zlib decompresses at a speed
# that is nearly independent of the level, so a higher level gives both a
# smaller file and faster reads; it only costs the time to write the file.
COMPRESSION_OPTS = 7


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

    def __init__(self, map_fname=None, max_samples=None, version='bayestar2019',
                 memmap=True):
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
                Defaults to :obj:`'bayestar2019'`.
            memmap (Optional[:obj:`bool`]): If ``True`` (the default), the
                reddening samples are not read into memory, and each query
                reads from the file only the pixels that it needs. This keeps
                memory usage low, and makes queries for a modest number of
                coordinates fast. If ``False``, the whole map is read into
                memory, which is faster if many queries will be made, or if a
                distance slice covering most of the sky is wanted. Memory
                mapping leaves the file open, so the file cannot be modified
                or deleted while the map is being queried.
        """

        if map_fname is None:
            map_fname = os.path.join(data_dir(), 'bayestar', '{}.h5'.format(version))

        t_start = time()
        
        self._memmap = memmap
        f = h5py.File(map_fname, 'r')

        # Load pixel information
        print('Loading pixel_info ...')
        self._pixel_info = f['/pixel_info'][:]
        self._DM_bin_edges = f['/pixel_info'].attrs['DM_bin_edges']
        self._n_distances = len(self._DM_bin_edges)
        self._n_pix = self._pixel_info.size

        t_pix_info = time()

        # Number of samples that will be returned by a query. When memory
        # mapping, the samples are truncated when they are read from the file,
        # rather than here.
        n_samples = f['/samples'].shape[1]
        self._n_samples = (n_samples if max_samples is None
                           else min(max_samples, n_samples))

        if memmap:
            if not f.attrs.get('repacked', False):
                print('Warning: this map has not been repacked, so '
                      'memory-mapped queries will be slow. Run '
                      '`bayestar.fetch(version=\'{}\')` to repack it.'.format(
                          version))

            # Keep the file open, and read the pixels that each query asks for
            self._f = f
            self._samples = f['/samples']
            self._best_fit = f['/best_fit']
        else:
            # Load reddening
            print('Loading samples ...')
            self._samples = f['/samples'][:,:self._n_samples,:]

            print('Loading best_fit ...')
            self._best_fit = f['/best_fit'][:]

            f.close()

            # Reshape best fit
            s = self._best_fit.shape
            self._best_fit = np.reshape(
                self._best_fit, (s[0], 1, s[1]))  # (pixels, samples=1, distances)

        t_reddening = time()

        # Replace NaNs in reliable distance estimates with +-infinity
        print('Replacing NaNs in reliable distance estimates ...')
        for k,v in [('DM_reliable_min',np.inf), ('DM_reliable_max',-np.inf)]:
            idx = ~np.isfinite(self._pixel_info[k])
            self._pixel_info[k][idx] = v

        t_nan = time()
        
        # The pixels are stored in order of (nside, healpix_index). This is
        # what allows a pixel to be located by binary search, and a query to
        # read only the pixels that it asks for, so check that it holds rather
        # than sorting the pixels ourselves.
        print('Checking that pixel_info is sorted ...')
        nside = self._pixel_info['nside']
        if np.any(nside[1:] < nside[:-1]):
            raise ValueError(
                'Pixels are not sorted by nside. The HDF5 file may be corrupt, '
                'or may not be one of the Bayestar maps.')

        t_check = time()

        self._nside_levels = np.unique(nside)
        self._hp_idx_sorted = []
        self._level_start = []

        start_idx = 0

        print('Extracting hp_idx_sorted at each nside ...')
        for nside_level in self._nside_levels:
            print('  nside = {}'.format(nside_level))
            end_idx = np.searchsorted(nside, nside_level, side='right')

            hp_idx = self._pixel_info['healpix_index'][start_idx:end_idx]
            if np.any(hp_idx[1:] <= hp_idx[:-1]):
                raise ValueError(
                    'Pixels at nside = {} are not sorted by healpix_index. '
                    'The HDF5 file may be corrupt, or may not be one of the '
                    'Bayestar maps.'.format(nside_level))

            self._hp_idx_sorted.append(hp_idx)
            self._level_start.append(start_idx)

            start_idx = end_idx

        t_finish = time()

        print('t = {:.3f} s'.format(t_finish - t_start))
        print('    pix_info: {: >7.3f} s'.format(t_pix_info-t_start))
        print('   reddening: {: >7.3f} s'.format(t_reddening-t_pix_info))
        print('         nan: {: >7.3f} s'.format(t_nan-t_reddening))
        print('       check: {: >7.3f} s'.format(t_check-t_nan))
        print('         idx: {: >7.3f} s'.format(t_finish-t_check))

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
                pix_idx[match_idx] = self._level_start[k] + idx

        return pix_idx

    def _gather_rows(self, dset, pix_idx, in_bounds_idx, is_best_fit=False):
        """
        Reads the rows of a memory-mapped dataset that a query is about to use
        into a small array, and returns that array along with the row within it
        of each queried coordinate.

        Args:
            dset (:obj:`h5py.Dataset`): The dataset to read from.
            pix_idx (:obj:`np.ndarray`): Row of the dataset for each queried
                coordinate. Coordinates outside the map have a row of ``-1``.
            in_bounds_idx (:obj:`np.ndarray`): Boolean array that is ``True``
                for coordinates that fall inside the map.
            is_best_fit (Optional[:obj:`bool`]): If ``True``, the dataset holds
                only the best-fit reddening, and is reshaped to add a (single)
                sample axis. Defaults to ``False``.

        Returns:
            The reddening at the requested pixels, and the row of each queried
            coordinate within it.
        """
        # h5py requires the rows to be in increasing order. This way, each
        # chunk that is needed is also read exactly once.
        rows = np.unique(pix_idx[in_bounds_idx])

        if rows.size == 0:
            # No coordinate is inside the map. A dummy row keeps the indexing
            # below valid, and all of its values are masked out by the caller.
            rows = np.zeros(1, dtype='i8')

        if is_best_fit:
            data = dset[rows]
            data = np.reshape(data, (data.shape[0], 1, data.shape[1]))
        else:
            data = dset[rows, :self._n_samples, :]

        # Coordinates outside the map (which have a row of -1) are pointed at
        # row 0, since they are masked out by the caller
        sel = np.searchsorted(rows, np.maximum(pix_idx, 0))

        return data, sel

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
            is_best_fit = True
        else:
            val = self._samples
            is_best_fit = False

        # When memory mapping, only the pixels that are asked for are read from
        # the file, into a small array
        if self._memmap:
            val, sel = self._gather_rows(val, pix_idx, in_bounds_idx,
                                         is_best_fit=is_best_fit)
        else:
            sel = pix_idx

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
                        * val[sel[idx_near], samp_idx, 0])
                else:
                    # print('idx_near: {} true'.format(np.sum(idx_near)))
                    # print('ret[idx_near].shape = {}'.format(ret[idx_near].shape))
                    # print('val.shape = {}'.format(val.shape))
                    # print('pix_idx[idx_near].shape = {}'.format(pix_idx[idx_near].shape))

                    ret[idx_near] = (
                        a * val[sel[idx_near], samp_idx[idx_near], 0])

            # d > d(farthest distance slice)
            idx_far = (bin_idx_ceil == self._n_distances) & in_bounds_idx
            if np.any(idx_far):
                # print('idx_far: {} true'.format(np.sum(idx_far)))
                # print('pix_idx[idx_far].shape = {}'.format(pix_idx[idx_far].shape))
                # print('ret[idx_far].shape = {}'.format(ret[idx_far].shape))
                # print('val.shape = {}'.format(val.shape))
                if isinstance(samp_idx, slice):
                    ret[idx_far] = val[sel[idx_far], samp_idx, -1]
                else:
                    ret[idx_far] = val[sel[idx_far], samp_idx[idx_far], -1]

            # d(nearest distance slice) < d < d(farthest distance slice)
            idx_btw = ~idx_near & ~idx_far & in_bounds_idx
            if np.any(idx_btw):
                DM_ceil = self._DM_bin_edges[bin_idx_ceil[idx_btw]]
                DM_floor = self._DM_bin_edges[bin_idx_ceil[idx_btw]-1]
                a = (DM_ceil - dm[idx_btw]) / (DM_ceil - DM_floor)
                if isinstance(samp_idx, slice):
                    ret[idx_btw] = (
                        (1.-a[:,None])
                        * val[sel[idx_btw], samp_idx, bin_idx_ceil[idx_btw]]
                        + a[:,None]
                        * val[sel[idx_btw], samp_idx, bin_idx_ceil[idx_btw]-1]
                    )
                else:
                    ret[idx_btw] = (
                        (1.-a) * val[sel[idx_btw], samp_idx[idx_btw], bin_idx_ceil[idx_btw]]
                        +    a * val[sel[idx_btw], samp_idx[idx_btw], bin_idx_ceil[idx_btw]-1]
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
            ret = val[sel, samp_idx, :]   # Return all distances
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
                self._pixel_info['converged'][pix_idx].astype(bool))
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
            ret = np.reshape(ret, s[:1] + s[2:])
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


def h5_repack(h5_in, h5_out, chunk_pixels=CHUNK_PIXELS,
              compression_opts=COMPRESSION_OPTS):
    """
    Converts an HDF5 file containing a Bayestar map into the layout that
    :obj:`BayestarQuery` expects, in which each chunk of the file holds a block
    of adjacent pixels, with all of their samples and distance bins.

    Every dataset of the input file is kept, along with the attributes that
    describe the map, and the output file is marked with a ``repacked``
    attribute recording the number of pixels per chunk.

    The pixels in the original file are already ordered by
    ``(nside, healpix_index)``, and are written out in the same order.

    Args:
        h5_in (:obj:`str`): Filename of the original HDF5 file.
        h5_out (:obj:`str`): Filename to write the repacked file to.
        chunk_pixels (Optional[:obj:`int`]): Number of pixels stored in each
            chunk. Defaults to :obj:`CHUNK_PIXELS`.
        compression_opts (Optional[:obj:`int`]): Deflate compression level,
            between 1 and 9. Defaults to :obj:`COMPRESSION_OPTS`.
    """

    # Number of pixels read and written at a time. This only affects how much
    # memory repacking uses.
    block_pixels = 100000

    with h5py.File(h5_in, 'r') as f_in, h5py.File(h5_out, 'w') as f_out:
        # Keep the attributes that describe the map, and record that this file
        # has been repacked, so that it can be identified without inspecting
        # the layout of the chunks
        for key in f_in.attrs:
            f_out.attrs[key] = f_in.attrs[key]
        f_out.attrs['repacked'] = True
        f_out.attrs['chunk_pixels'] = chunk_pixels

        # The pixel information is small, and is copied as it is
        f_in.copy('pixel_info', f_out)

        # GRDiagnostic is not used by `dustmaps`, but is kept so that the
        # repacked file is equivalent to the one that was downloaded. Not every
        # version of the map has it.
        for name in ('samples', 'best_fit', 'GRDiagnostic'):
            if name not in f_in:
                continue

            dset_in = f_in[name]
            n_pix = dset_in.shape[0]

            dset_out = f_out.create_dataset(
                name,
                shape=dset_in.shape,
                dtype=dset_in.dtype,
                chunks=(chunk_pixels,) + dset_in.shape[1:],
                compression='gzip',
                compression_opts=compression_opts)

            for start in range(0, n_pix, block_pixels):
                stop = min(start + block_pixels, n_pix)
                dset_out[start:stop] = dset_in[start:stop]


def h5_is_repacked(h5_fname):
    """
    Returns ``True`` if the given HDF5 file holds a Bayestar map that has been
    repacked by :obj:`h5_repack`, in which each chunk contains a block of
    adjacent pixels.

    Repacked files are marked with a ``repacked`` attribute, rather than being
    recognized by the layout of their chunks. That way, changing
    :obj:`CHUNK_PIXELS` does not make every existing file look as though it
    needs to be repacked again.

    Args:
        h5_fname (:obj:`str`): Filename of the HDF5 file.
    """
    try:
        with h5py.File(h5_fname, 'r') as f:
            repacked = f.attrs.get('repacked', False)
            chunks = f['samples'].chunks
    except (IOError, KeyError):
        return False

    return bool(repacked) and chunks is not None


# Expected size (in Bytes) and datasets of a repacked map file, for each version
# of the map. The size is a rough guide (see `fetch_utils.h5_file_exists`); the
# datasets are the check that matters.
REPACKED_SIZES = {
    'bayestar2015': 4613454007,
    'bayestar2017': 5079647547,
    'bayestar2019': 665933674
}

REPACKED_DSETS = {
    'bayestar2015': {
        'samples': (2437292, 20, 31),
        'best_fit': (2437292, 31),
        'GRDiagnostic': (2437292, 31),
        'pixel_info': (2437292,)
    },
    'bayestar2017': {
        'samples': (3420905, 18, 31),
        'best_fit': (3420905, 31),
        'GRDiagnostic': (3420905, 31),
        'pixel_info': (3420905,)
    },
    'bayestar2019': {
        'samples': (4214070, 5, 120),
        'best_fit': (4214070, 120),
        'pixel_info': (4214070,)
    }
}


def fetch(version='bayestar2019', clobber=False):
    """
    Downloads the specified version of the Bayestar dust map.

    The downloaded file is repacked into the layout that :obj:`BayestarQuery`
    expects, which allows a query for random coordinates to read only the small
    part of the file that it needs. The original file is then deleted, since
    the repacked file is equivalent to it.

    Args:
        version (Optional[:obj:`str`]): The map version to download. Valid versions are
            :obj:`'bayestar2019'` (Green, Schlafly, Finkbeiner et al. 2019),
            :obj:`'bayestar2017'` (Green, Schlafly, Finkbeiner et al. 2018) and
            :obj:`'bayestar2015'` (Green, Schlafly, Finkbeiner et al. 2015). Defaults
            to :obj:`'bayestar2019'`.
        clobber (Optional[:obj:`bool`]): If ``True``, any existing file will be
            overwritten, even if it appears to match. If ``False`` (the default),
            :obj:`fetch()` will attempt to determine if the dataset already exists.
            This determination is not 100% robust against data corruption.

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

    # Expected size and datasets of the repacked file, for the checks below
    h5_size = REPACKED_SIZES[version]
    h5_dsets = REPACKED_DSETS[version]

    map_fname = os.path.join(data_dir(), 'bayestar', '{}.h5'.format(version))

    # Check if a repacked file already exists. The size check alone is not
    # enough, because the file as published is similar in size to the repacked
    # file, so the layout has to be checked as well.
    if not clobber:
        if (fetch_utils.h5_file_exists(map_fname, h5_size, dsets=h5_dsets)
                and h5_is_repacked(map_fname)):
            print('File appears to exist already. Call `fetch(clobber=True)` '
                  'to force overwriting of existing file.')
            return

    # The original, as-published file is the input to the repacking. If such a
    # file is already on disk -- for example, one downloaded by an older
    # version of `dustmaps` -- then repack it, rather than download it again.
    # The file is checked over before it is reused, so that a truncated
    # download is replaced, rather than repacked.
    if (os.path.isfile(map_fname)
            and not h5_is_repacked(map_fname)
            and fetch_utils.h5_file_exists(map_fname, dsets=h5_dsets)):
        orig_fname = map_fname
    else:
        orig_fname = map_fname + '.orig'

        fetch_utils.dataverse_download_doi(
            doi,
            orig_fname,
            file_requirements=requirements)

    # Convert to the layout that BayestarQuery expects. The repacked file is
    # written to a temporary name, so that the original is left in place if
    # anything goes wrong.
    print('Repacking file...')
    h5_repack(orig_fname, map_fname + '.repacked')
    os.replace(map_fname + '.repacked', map_fname)

    # Cleanup
    if orig_fname != map_fname:
        print('Removing original file...')
        os.remove(orig_fname)


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
