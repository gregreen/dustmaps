#!/usr/bin/env python
#
# decaps.py
# Reads the DECaPS dust reddening maps, described in
# Zucker, Saydjari, & Speagle et al. 2025.
#
# Copyright (C) 2025  Catherine Zucker, Andrew Saydjari, and Gregory M. Green
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
from .map_base import DustMap, ensure_flat_galactic
from . import fetch_utils
import warnings

from time import time
from tqdm import tqdm


# The layout that `DECaPSQuery` expects: each chunk of the file holds a block of
# adjacent pixels, with all of their samples and distance bins. The published
# file instead stores a whole slab of 200842 pixels for one (sample, distance)
# pair at a time -- a layout that suits reading the map from end to end, which is
# how the original code loaded it, but that makes a query for a random coordinate
# read about 240 MB, and a 256-coordinate query read the whole 62 GB file.
#
# Sixteen pixels per chunk was measured to be the best trade-off. It also shrinks
# the file by 46%: the published 33.84 GB becomes 18.37 GB, because the extinction
# along a line of sight is smoother than it is from pixel to pixel, so encoding 16
# pixels x 5 samples x 120 distances together compresses better than encoding runs
# of single values does. The gain is almost all in the samples (0.491 -> 0.212 kB
# per pixel for `samples`, 0.145 -> 0.116 for `mean`). Sixty-four pixels per chunk
# save only a little more, while making queries slower, so the published
# compression settings were kept (deflate, with byte shuffling), except that the
# level is raised from 3 to 7, to match the rest of `dustmaps`.
#
# Repacking the published file on a laptop took 45 minutes (deflate is slow), and
# `scripts/repack_decaps.py` does the same job for the authors of the map.
CHUNK_PIXELS = 16
COMPRESSION_OPTS = 7
SHUFFLE = True


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


class DECaPSQuery(DustMap):
    """
    Queries the DECaPS 3D dust maps (Zucker, Saydjari, & Speagle et al. 2025).
    By default, the file is memory-mapped (i.e., `memmap=True`), which makes
    loading nearly instantaneous and keeps RAM usage low, but makes individual
    queries slower. Passing `memmap=False` loads the entire file into RAM: this
    comes with a large startup cost, but makes downstream queries significantly
    faster. The map covers the southern Galactic plane (239 < l < 6, |b| < 10),
    amounting to 6% of the sky. When combined with the BayestarQuery, this
    DECaPS query enables reddening estimates over the entire Galactic plane
    |b| < 10.
    """

    def __init__(self, map_fname=None, max_samples=None, mean_only=False,
                 memmap=True):
        """
        Args:
            map_fname (Optional[str]): Filename of the DECaPS map. Defaults to None, 
                meaning that the default location is used.
            max_samples (Optional[:obj:`int`]): Maximum number of samples of the map to
                load. Use a lower number in order to decrease memory usage.
                Defaults to :obj:`None`, meaning that all samples will be loaded.
            mean_only (Optional[bool]): If True, only the mean map file is available 
                to query (no samples). For users who did not download the larger 
                mean_and_samples file, this mode is required. However, you will not be 
                able to query in any other mode ('random_sample', 'random_sample_per_pix', 
                'samples', 'median', or 'percentile'). Defaults to False.
            memmap (Optional[bool]): If True, memory-map the file rather than
                loading it into RAM. Memory-mapping makes loading very fast and
                keeps memory usage low, which is ideal if you are querying a
                modest number of random coordinates. If you intend to perform
                large queries (e.g., generating a reddening map over a large
                region), then pass `memmap=False` to load the map into RAM,
                which is much faster per query, but incurs a large startup cost
                in both time and RAM. Defaults to True.
        """

        self._mean_only = mean_only
        self._memmap = memmap

        if map_fname is None:
            if self._mean_only:
                map_fname = os.path.join(data_dir(), 'decaps', 'decaps_mean.h5')
                if not os.path.isfile(map_fname):
                    map_fname = os.path.join(data_dir(), 'decaps', 'decaps_mean_and_samples.h5')
            else:
                map_fname = os.path.join(data_dir(), 'decaps', 'decaps_mean_and_samples.h5')
                if not os.path.isfile(map_fname):
                    raise ValueError(
                        "The file containing both the mean and samples was not found at the default location on disk. "
                        "Please confirm you have downloaded 'decaps_mean_and_samples.h5'."
                    )

        if self._memmap:
            print("Memory-mapping the DECaPS map. Loading is nearly instantaneous, "
                  "and RAM usage is low, but queries will read from disk. If you "
                  "are going to perform many large queries, then consider using "
                  "`memmap=False` instead, which loads the map into RAM.")

        f = h5py.File(map_fname, 'r')

        self._DM_bin_edges = f['/pixel_info'].attrs['DM_bin_edges']
        self._n_distances = len(self._DM_bin_edges)
        self._n_pix = f['/pixel_info/healpix_index'].size
        self._chunk_size = 100000

        # The pixel info is small, so it is always loaded into RAM, even when
        # memory-mapping the maps themselves. The pixel index is converted from
        # the uint64 it is stored as to int64, because comparing it against the
        # int64 that healpy returns otherwise promotes the whole array to
        # float64 on every query, which costs about 50 ms on a map this size.
        self._nside = f['/pixel_info'].attrs['nside']
        self._hp_idx_sorted = f['/pixel_info/healpix_index'][:].astype('i8')
        self._pixel_info = {name: ds[()] for name, ds in f['pixel_info'].items()}

        # Queries rely on the pixels being sorted by HEALPix index.
        if np.any(self._hp_idx_sorted[1:] <= self._hp_idx_sorted[:-1]):
            raise ValueError(
                'Pixels are not sorted by healpix_index. The DECaPS map must be '
                'sorted in order for it to be queried.')

        # The mean map is expected to carry an empty samples axis, as the maps
        # published on the Dataverse do
        if f['/mean'].ndim != 3 or f['/mean'].shape[1] != 1:
            raise ValueError(
                'Expected the mean map to have shape '
                '(n_pixels, 1, n_distances), but it has shape {}.'.format(
                    f['/mean'].shape))

        if self._memmap:
            if not h5_is_repacked(map_fname):
                print('Warning: this map has not been repacked, so '
                      'memory-mapped queries will read far more of the file '
                      'than they need to. Run `decaps.fetch()` to repack it.')

            # Keep the file open, and read the map from disk as it is queried.
            self._f = f
            self._mean = f['/mean']
            if not self._mean_only:
                self._samples = f['/samples']
                self._n_samples = self._samples.shape[1]
                if max_samples is not None:
                    self._n_samples = min(max_samples, self._n_samples)
        else:
            # Report the memory that the arrays themselves need, which is not the
            # same as the size of the file: this map has 51.4 million pixels, so
            # a single sample of each of them is already 12.3 GB, and all five
            # are 62 GB, however small the file has been compressed to.
            itemsize = np.dtype(f['/mean'].dtype).itemsize
            n_bytes = self._n_pix * self._n_distances * itemsize
            if not self._mean_only:
                n_samples = (f['/samples'].shape[1] if max_samples is None
                             else min(max_samples, f['/samples'].shape[1]))
                n_bytes += self._n_pix * n_samples * self._n_distances * itemsize

            print('You are about to read {:.1f} GB into RAM, which is a lot for a '
                  'map that can be memory-mapped instead (leave out '
                  '`memmap=False`).'.format(n_bytes / 1e9))

            print("Allocating memory for mean map...")
            mean_dataset = f['/mean']
            self._mean = np.empty_like(mean_dataset)

            print("Loading mean map in chunks...")
            for i in tqdm(range(0, mean_dataset.shape[0], self._chunk_size), desc="Mean Map Loading"):
                self._mean[i:i + self._chunk_size] = mean_dataset[i:i + self._chunk_size]

            if not self._mean_only:
                samples_dataset = f['/samples']
                print("Allocating memory for samples...")
                
                # If max_samples is provided, slice the samples to load only the desired number of samples
                if max_samples is not None:
                    self._samples = np.empty((samples_dataset.shape[0], max_samples, samples_dataset.shape[2]))
                    print(f"Loading only {max_samples} samples...")
                else:
                    self._samples = np.empty_like(samples_dataset)

                print("Loading samples in chunks...")
                for i in tqdm(range(0, samples_dataset.shape[0], self._chunk_size), desc="Samples Loading"):
                    if max_samples is not None:
                        self._samples[i:i + self._chunk_size] = samples_dataset[i:i + self._chunk_size, :max_samples, :]
                    else:
                        self._samples[i:i + self._chunk_size] = samples_dataset[i:i + self._chunk_size]

                self._n_samples = self._samples.shape[1]  # The number of samples loaded

            print("Data loading complete!")

            f.close()

    def _gather_rows(self, dset, pix_idx, in_bounds_idx, is_mean=False):
        """
        Gathers the rows of the map that are needed to answer a query. Only rows
        that are actually referenced are read from disk, and each row is read
        exactly once, no matter how many coordinates reference it.

        Args:
            dset: The dataset (or array) to read from.
            pix_idx (:obj:`ndarray` of `int`): The row of the map that each
                coordinate falls in (:obj:`-1` for coordinates that fall outside
                the footprint of the map).
            in_bounds_idx (:obj:`ndarray` of `bool`): Which of the coordinates
                given in :obj:`pix_idx` fall within the footprint of the map.
            is_mean (Optional[:obj:`bool`]): Whether the requested dataset is the
                mean map, which lacks the samples axis. Defaults to False.

        Returns:
            A tuple `(data, sel)`, where `data` contains the requested rows of
            the map, and `sel` gives, for each coordinate, the row of `data` that
            it falls in.
        """

        rows = np.unique(pix_idx[in_bounds_idx])

        # h5py will not accept an empty list of rows, so use a dummy row. The
        # results for these rows are discarded.
        if rows.size == 0:
            rows = np.zeros(1, dtype='i8')

        if is_mean:
            # The mean map is stored with an empty samples axis, i.e. the shape
            # is (n_pixels, 1, n_distances), so that it can be handed to the
            # rest of the query in the same way as the samples themselves.
            data = dset[rows]
        else:
            data = dset[rows, :self._n_samples, :]

        # For coordinates outside the map, `pix_idx` is -1; these are clamped to
        # zero here, and their results are discarded downstream.
        sel = np.searchsorted(rows, np.maximum(pix_idx, 0))

        return data, sel


    def _find_data_idx(self, l, b):
    
        pix_idx = np.empty(l.shape, dtype='i8')
        pix_idx[:] = -1

        # Search at each nside
        ipix = lb2pix(self._nside, l, b, nest=True)

        # Find the insertion points of the query pixels in the large, ordered pixel list
        idx = np.searchsorted(self._hp_idx_sorted, ipix, side='left')

        # Determine which insertion points are beyond the edge of the pixel list
        in_bounds = (idx < self._hp_idx_sorted.size)

        # Determine which query pixels are correctly placed
        idx[~in_bounds] = -1
        match_idx = (self._hp_idx_sorted[idx] == ipix)
        match_idx[~in_bounds] = False
        idx = idx[match_idx]

        if np.any(match_idx):
            # The pixels are sorted, and the map is read in that order, so the
            # row of the map that a pixel sits in is the pixel's own position in
            # this list
            pix_idx[match_idx] = idx

        return pix_idx
        
    def _raise_on_mode(self, mode):
    	"""
    	Checks that the provided query mode is one of the accepted values. If
    	not, raises a :obj:`ValueError`.
    	"""
    
    	if self._mean_only:
        	valid_modes = ['mean']
    	else:      
        	valid_modes = [
            	'random_sample',
            	'random_sample_per_pix',
            	'samples',
            	'median',
            	'mean',
            	'percentile'
        	]
    
    	if mode not in valid_modes:
        	raise ValueError(
        	    '"{}" is not a valid `mode`. Valid modes are:\n'
            	'  {}'.format(mode, valid_modes))
        
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
    def query(self, coords, mode='mean', return_flags=False, pct=None):
        """
        Returns reddening at the requested coordinates. There are several
        different query modes, which handle the probabilistic nature of the map
        differently.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.
            mode (Optional[:obj:`str`]): Six different query modes are available:
                'random_sample', 'random_sample_per_pix' 'samples', 'median',
                'mean', and 'percentile'. The :obj:`mode` determines how the
                output will reflect the probabilistic nature of the DECaPS
                dust maps.
            return_flags (Optional[:obj:`bool`]): If :obj:`True`, then QA flags will be
                returned in a second numpy structured array. That is, the query
                will return :obj:`ret`, :obj:`flags`, where :obj:`ret` is the normal return
                value, containing reddening. Defaults to :obj:`False`.
            pct (Optional[:obj:`float` or list/array of :obj:`float`]): If the mode is
                :obj:`percentile`, then :obj:`pct` specifies which percentile(s) is
                (are) returned.

        Returns:
            Reddening at the specified coordinates, in units of magnitudes of E(B-V). 
            Note that this reddening is different from the Bayestar19 query, which 
            returns reddening in an arbitrary unit and must be converted to E(B-V)

            The shape of the output depends on the :obj:`mode`, and on whether
            :obj:`coords` contains distances.

            If :obj:`coords` does not specify distance(s), then the shape of the
            output begins with :obj:`coords.shape`. If :obj:`coords` does specify
            distance(s), then the shape of the output begins with
            :obj:`coords.shape + ([number of distance bins],)`.
            
            If :obj:`mode` is :obj:`'mean'`, then at each coordinate/distance, the
            mean reddening is returned.  

            If :obj:`mode` is :obj:`'random_sample'`, then at each
            coordinate/distance, a random sample of reddening is given.

            If :obj:`mode` is :obj:`'random_sample_per_pix'`, then the sample chosen
            for each angular pixel of the map will be consistent. For example,
            if two query coordinates lie in the same map pixel, then the same
            random sample will be chosen from the map for both query
            coordinates.

            If :obj:`mode` is :obj:`'median'`, then at each coordinate/distance, the
            median reddening is returned.

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
            (whether or not the line-of-sight fit converged in a given pixel),
            :obj:`"infilled"`(whether or not the pixel was infilled), 
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
        pix_idx = self._find_data_idx(coords.l.deg, coords.b.deg)
        in_bounds_idx = (pix_idx != -1)

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
        elif mode == 'mean':
            samp_idx = slice(None)
            n_samp_ret = 1
        else:
            # Return all samples in each queried coordinate
            samp_idx = slice(None)
            n_samp_ret = self._n_samples

        # Read the rows of the map that the query touches. This reads each row
        # from disk at most once, and only reads rows that are actually needed.
        if mode == 'mean':
            val, sel = self._gather_rows(
                self._mean, pix_idx, in_bounds_idx, is_mean=True)
        else:
            val, sel = self._gather_rows(
                self._samples, pix_idx, in_bounds_idx, is_mean=False)

        # Create empty array to store flags
        if return_flags:
            if has_dist:
                # If distances are provided in query, return only covergence and
                # whether or not this distance is reliable
                dtype = [('converged', 'bool'),
                         ('infilled', 'bool'),
                         ('reliable_dist', 'bool')]
                # shape = (n_coords_ret)
            else:
                # Return convergence and reliable distance ranges
                dtype = [('converged', 'bool'),
                		 ('infilled', 'bool'),
                         ('min_reliable_distmod', 'f4'),
                         ('max_reliable_distmod', 'f4')]
            flags = np.empty(n_coords_ret, dtype=dtype)


        # Extract the correct distance bin (possibly using linear interpolation)
        if has_dist: # Distance has been provided
            # Determine ceiling bin index for each coordinate
            dm = 5. * (np.log10(d) + 2.)
            bin_idx_ceil = np.searchsorted(self._DM_bin_edges, dm)

            # Create NaN-filled return arrays
            if isinstance(samp_idx, slice):
                ret = np.full((n_coords_ret, n_samp_ret), np.nan, dtype='f2')
            else:
                ret = np.full((n_coords_ret,), np.nan, dtype='f2')

            # d < d(nearest distance slice)
            idx_near = (bin_idx_ceil == 0) & in_bounds_idx
            if np.any(idx_near):
                a = 10.**(0.2 * (dm[idx_near] - self._DM_bin_edges[0]))
                if isinstance(samp_idx, slice):
                    ret[idx_near] = (
                        a[:,None]
                        * val[sel[idx_near], samp_idx, 0])
                else:

                    ret[idx_near] = (
                        a * val[sel[idx_near], samp_idx[idx_near], 0])

            # d > d(farthest distance slice)
            idx_far = (bin_idx_ceil == self._n_distances) & in_bounds_idx
            if np.any(idx_far):

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

        # Flag: convergence
        if return_flags:
        
            flags['converged'] = (
                self._pixel_info['converged'][pix_idx].astype(bool))
            flags['converged'][~in_bounds_idx] = False
            
            flags['infilled'] = (
                self._pixel_info['infilled'][pix_idx].astype(bool))
            flags['infilled'][~in_bounds_idx] = False

        # Reduce the samples in the requested manner
        if mode == 'median':
            ret = np.median(ret, axis=1)
        elif mode == 'mean':
            # Remove "samples" axis
            s = ret.shape
            ret = np.reshape(ret, s[:1] + s[2:])
        elif mode == 'percentile':
            ret = np.nanpercentile(ret, pct, axis=1)
            if not scalar_pct:
                # (percentile, pixel) -> (pixel, percentile)
                # (pctile, pixel, distance) -> (pixel, distance, pctile)
                ret = np.moveaxis(ret, 0, -1)
        elif mode == 'samples':
            # Swap sample and distance axes to be consistent with other 3D dust
            # maps. The output shape will be (pixel, distance, sample).
            if not has_dist:
                np.swapaxes(ret, 1, 2)
        	
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

class DECaPSQueryLite(DustMap):
    """
    Queries the DECaPS 3D dust maps (Zucker, Saydjari, & Speagle et al. 2025) 
    **with memory mapping**, so the entire file does NOT need to be loaded into RAM.   
    This is the query you should use for smaller queries, since overhead is smaller. 
    The map covers the southern Galactic plane (239 < l < 6, |b| < 10), amounting
    to 6% of the sky. When combined with the BayestarQuery, this DECaPS query enables
    reddening estimates over the entire Galactic plane |b| < 10.
    """

    def __init__(self, map_fname=None, mean_only=False, contiguous=False):
        """
        Args:
            map_fname (Optional[str]): Filename of the DECaPS map. Defaults to None, 
                meaning that the default location is used.
            mean_only (Optional[bool]): If True, only the mean map file is available 
                to query (no samples). For users who did not download the larger 
                mean_and_samples file, this mode is required. However, you will not be 
                able to query in any other mode ('random_sample', 'random_sample_per_pix', 
                'samples', 'median', or 'percentile'). Defaults to False.
            contiguous (Optional[bool]): If you are querying a dense grid of points in a localized
            	area of sky, set contiguous=True to enable more efficient read into memory.
            	 Defaults to False.
        """
        

        self._mean_only = mean_only
        self._contiguous = contiguous

        if self._mean_only:
            if map_fname is None:
                map_fname = os.path.join(data_dir(), 'decaps', 'decaps_mean.h5')
            if not os.path.isfile(map_fname):
                map_fname = os.path.join(data_dir(), 'decaps', 'decaps_mean_and_samples.h5')
        else:
            if map_fname is None:
                map_fname = os.path.join(data_dir(), 'decaps', 'decaps_mean_and_samples.h5')
                if not os.path.isfile(map_fname):
                    raise ValueError(
                        "The file containing both the mean and samples was not found at the default location on disk. "
                        "Please confirm you have downloaded 'decaps_mean_and_samples.h5'."
                    )

        self._fhandle = h5py.File(map_fname, 'r')
        print("Loading meta pixel info...")
        self._DM_bin_edges = self._fhandle['pixel_info'].attrs['DM_bin_edges']
        self._n_distances = len(self._DM_bin_edges)
        self._n_pix = self._fhandle['pixel_info/healpix_index'].size
        self._n_samples = self._fhandle['pixel_info'].attrs['n_samples']

        # TODO update file structure
        self._nside = self._fhandle['pixel_info'].attrs['nside']
        self._hp_idx_sorted = self._fhandle['pixel_info/healpix_index'][:]
        self._data_idx = np.arange(len(self._hp_idx_sorted))
        print("Meta pixel info loaded!")

    def _find_data_idx(self, l, b):
        pix_idx = np.empty(l.shape, dtype='i8')
        pix_idx[:] = -1

        # Search at each nside
        ipix = lb2pix(self._nside, l, b, nest=True)

        # Find the insertion points of the query pixels in the large, ordered pixel list
        idx = np.searchsorted(self._hp_idx_sorted, ipix, side='left')

        # Determine which insertion points are beyond the edge of the pixel list
        in_bounds = (idx < self._hp_idx_sorted.size)

        # Determine which query pixels are correctly placed
        idx[~in_bounds] = -1
        match_idx = (self._hp_idx_sorted[idx] == ipix)
        match_idx[~in_bounds] = False
        idx = idx[match_idx]

        if np.any(match_idx):
            pix_idx[match_idx] = self._data_idx[idx]

        return pix_idx

    def _raise_on_mode(self, mode):
    
    	"""
    	Checks that the provided query mode is one of the accepted values. If
    	not, raises a :obj:`ValueError`.
    	"""
    
    	if self._mean_only:
        	valid_modes = ['mean']
    	else:      
        	valid_modes = [
            	'random_sample',
            	'random_sample_per_pix',
            	'samples',
            	'median',
            	'mean',
            	'percentile'
        	]
    
    	if mode not in valid_modes:
        	raise ValueError(
        	    '"{}" is not a valid `mode`. Valid modes are:\n'
            	'  {}'.format(mode, valid_modes))
        
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
    def query(self, coords, mode='mean', return_flags=False, pct=None):
        """
        Returns reddening at the requested coordinates. There are several
        different query modes, which handle the probabilistic nature of the map
        differently.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.
            mode (Optional[:obj:`str`]): Six different query modes are available:
                'random_sample', 'random_sample_per_pix' 'samples', 'median',
                'mean', and 'percentile'. The :obj:`mode` determines how the
                output will reflect the probabilistic nature of the DECaPS
                dust maps.
            return_flags (Optional[:obj:`bool`]): If :obj:`True`, then QA flags will be
                returned in a second numpy structured array. That is, the query
                will return :obj:`ret`, :obj:`flags`, where :obj:`ret` is the normal return
                value, containing reddening. Defaults to :obj:`False`.
            pct (Optional[:obj:`float` or list/array of :obj:`float`]): If the mode is
                :obj:`percentile`, then :obj:`pct` specifies which percentile(s) is
                (are) returned.

        Returns:
            Reddening at the specified coordinates, in units of magnitudes of E(B-V). 
            Note that this reddening is different from the Bayestar19 query, which 
            returns reddening in an arbitrary unit and must be converted to E(B-V)

            The shape of the output depends on the :obj:`mode`, and on whether
            :obj:`coords` contains distances.

            If :obj:`coords` does not specify distance(s), then the shape of the
            output begins with :obj:`coords.shape`. If :obj:`coords` does specify
            distance(s), then the shape of the output begins with
            :obj:`coords.shape + ([number of distance bins],)`.
            
            If :obj:`mode` is :obj:`'mean'`, then at each coordinate/distance, the
            mean reddening is returned.  

            If :obj:`mode` is :obj:`'random_sample'`, then at each
            coordinate/distance, a random sample of reddening is given.

            If :obj:`mode` is :obj:`'random_sample_per_pix'`, then the sample chosen
            for each angular pixel of the map will be consistent. For example,
            if two query coordinates lie in the same map pixel, then the same
            random sample will be chosen from the map for both query
            coordinates.

            If :obj:`mode` is :obj:`'median'`, then at each coordinate/distance, the
            median reddening is returned.

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
            coordinates include distances, the QA flags will be :obj:`"infilled"`
            (whether or not the pixel was infilled), :obj:`"converged"`
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
        pix_idx = self._find_data_idx(coords.l.deg, coords.b.deg)
        in_bounds_idx = (pix_idx != -1)

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
        elif mode == 'mean':
            samp_idx = slice(None)
            n_samp_ret = 1
        else:
            # Return all samples in each queried coordinate
            samp_idx = slice(None)
            n_samp_ret = self._n_samples

        if mode == 'mean':
            modekey = "mean"
        else:
            modekey = "samples"

        if self._contiguous:
            minIndx = np.min(pix_idx)
            maxIndx = np.max(pix_idx)
            
            if modekey == "mean":
                self._datacache = self._fhandle[modekey][minIndx:maxIndx+1,:]
            else:
                self._datacache = self._fhandle[modekey][minIndx:maxIndx+1, :, :]
            def data_handle(pix_idx, samp_slice, dist_slice):
                # Handle both single indices and array indexing
                if isinstance(pix_idx, (list, np.ndarray)):
                    adjusted_idx = pix_idx - minIndx
                    return self._datacache[adjusted_idx, samp_slice, dist_slice]
                else:
                    return self._datacache[pix_idx - minIndx, samp_slice, dist_slice]
            self._DM_reliable_min_cache = self._fhandle["pixel_info/DM_reliable_min"][minIndx:maxIndx+1]
            def DM_reliable_min_handle(pix_idx):
                return self._DM_reliable_min_cache[pix_idx-minIndx]
            self._DM_reliable_max_cache = self._fhandle["pixel_info/DM_reliable_max"][minIndx:maxIndx+1]
            def DM_reliable_max_handle(pix_idx):
                return self._DM_reliable_max_cache[pix_idx-minIndx]
            self._converged_cache = self._fhandle["pixel_info/converged"][minIndx:maxIndx+1]
            def converged_handle(pix_idx):
                return self._converged_cache[pix_idx-minIndx]
            self._infilled_cache = self._fhandle["pixel_info/infilled"][minIndx:maxIndx+1]
            def infilled_handle(pix_idx):
                return self._infilled_cache[pix_idx-minIndx]
        else:
            data_handle = self._fhandle[modekey]
            DM_reliable_min_handle = self._fhandle["pixel_info/DM_reliable_min"]
            DM_reliable_max_handle = self._fhandle["pixel_info/DM_reliable_max"]
            converged_handle = self._fhandle["pixel_info/converged"]
            infilled_handle = self._fhandle["pixel_info/infilled"]

        # Create empty array to store flags
        if return_flags:
            if has_dist:
                # If distances are provided in query, return only covergence and
                # whether or not this distance is reliable
                dtype = [('converged', 'bool'),
                         ('infilled', 'bool'),
                         ('reliable_dist', 'bool')]
                # shape = (n_coords_ret)
            else:
                # Return convergence and reliable distance ranges
                dtype = [('converged', 'bool'),
                		 ('infilled', 'bool'),
                         ('min_reliable_distmod', 'f4'),
                         ('max_reliable_distmod', 'f4')]
            flags = np.empty(n_coords_ret, dtype=dtype)


        # Extract the correct distance bin (possibly using linear interpolation)
        if has_dist: # Distance has been provided
            # Determine ceiling bin index for each coordinate
            dm = 5. * (np.log10(d) + 2.)
            bin_idx_ceil = np.searchsorted(self._DM_bin_edges, dm)

            # Create NaN-filled return arrays
            if isinstance(samp_idx, slice):
                ret = np.full((n_coords_ret, n_samp_ret), np.nan, dtype='f2')
            else:
                ret = np.full((n_coords_ret,), np.nan, dtype='f2')

            # d < d(nearest distance slice)
            idx_near = (bin_idx_ceil == 0) & in_bounds_idx
            if np.any(idx_near):
                dataindx2use = np.where(idx_near)[0]
                a = 10.**(0.2 * (dm[idx_near] - self._DM_bin_edges[0]))
                if isinstance(samp_idx, slice):
                    if self._contiguous:
                        ret[idx_near] = (
                            a[:,None]
                            * data_handle(pix_idx[idx_near], samp_idx, 0))
                    else:
                        for idx, pix in enumerate(pix_idx[idx_near]):
                            loc_indx = dataindx2use[idx]
                            ret[loc_indx] = (
                                a[idx]
                                * data_handle[pix, samp_idx, 0])
                else:
                    if self._contiguous:
                        ret[idx_near] = (
                            a * data_handle(pix_idx[idx_near], samp_idx, 0))
                    else:
                        for idx, pix in enumerate(pix_idx[idx_near]):
                            loc_indx = dataindx2use[idx]
                            ret[loc_indx] = (
                                a[idx]
                                * data_handle[pix, samp_idx[loc_indx], 0])

            # d > d(farthest distance slice)
            idx_far = (bin_idx_ceil == self._n_distances) & in_bounds_idx
            if np.any(idx_far):
                dataindx2use = np.where(idx_far)[0]
                if isinstance(samp_idx, slice):
                    if self._contiguous:
                        ret[idx_far] = (data_handle(pix_idx[idx_far], samp_idx, -1))
                    else:
                        for idx, pix in enumerate(pix_idx[idx_far]):
                            loc_indx = dataindx2use[idx]
                            ret[loc_indx] = (data_handle[pix, samp_idx, -1])
                else:
                    if self._contiguous:
                        ret[idx_far] = (data_handle(pix_idx[idx_far], samp_idx, -1))
                    else:
                        for idx, pix in enumerate(pix_idx[idx_far]):
                            loc_indx = dataindx2use[idx]
                            ret[loc_indx] = (data_handle[pix, samp_idx[loc_indx], -1])
            # d(nearest distance slice) < d < d(farthest distance slice)
            idx_btw = ~idx_near & ~idx_far & in_bounds_idx
            if np.any(idx_btw):
                dataindx2use = np.where(idx_btw)[0]
                DM_ceil = self._DM_bin_edges[bin_idx_ceil[idx_btw]]
                DM_floor = self._DM_bin_edges[bin_idx_ceil[idx_btw]-1]
                a = (DM_ceil - dm[idx_btw]) / (DM_ceil - DM_floor)
                if isinstance(samp_idx, slice):
                    if self._contiguous:
                        ret[idx_btw] = (
                            (1.-a[:,None])
                            * data_handle(pix_idx[idx_btw], samp_idx, bin_idx_ceil[idx_btw])
                            + a[:,None]
                            * data_handle(pix_idx[idx_btw], samp_idx, bin_idx_ceil[idx_btw]-1)
                        )
                    else:
                        for idx, pix in enumerate(pix_idx[idx_btw]):
                            loc_indx = dataindx2use[idx]
                            ret[loc_indx] = (
                                (1.-a[idx])
                                * data_handle[pix, samp_idx, bin_idx_ceil[loc_indx]]
                                + a[idx]
                                * data_handle[pix, samp_idx, bin_idx_ceil[loc_indx]-1]
                            )
                else:
                    if self._contiguous:
                        ret[idx_btw] = (
                            (1.-a[:,None])
                            * data_handle(pix_idx[idx_btw], samp_idx, bin_idx_ceil[idx_btw])
                            + a[:,None]
                            * data_handle(pix_idx[idx_btw], samp_idx, bin_idx_ceil[idx_btw]-1)
                        )
                    else:
                        for idx, pix in enumerate(pix_idx[idx_btw]):
                            loc_indx = dataindx2use[idx]
                            ret[loc_indx] = (
                                (1.-a[idx]) * data_handle[pix, samp_idx[loc_indx], bin_idx_ceil[loc_indx]]
                                + a[idx] * data_handle[pix, samp_idx[loc_indx], bin_idx_ceil[loc_indx]-1]
                            )
            # Flag: distance in reliable range?
            if return_flags:
                if self._contiguous:
                    dm_min = DM_reliable_min_handle(pix_idx)
                    dm_max = DM_reliable_max_handle(pix_idx)
                else:
                    dm_min = np.empty(n_coords_ret)
                    dm_max = np.empty(n_coords_ret)
                    for idx, pix in enumerate(pix_idx):
                        dm_min[idx] = DM_reliable_min_handle[pix]
                        dm_max[idx] = DM_reliable_max_handle[pix]
                flags['reliable_dist'] = (
                    (dm >= dm_min) &
                    (dm <= dm_max) &
                    np.isfinite(dm_min) &
                    np.isfinite(dm_max))
                flags['reliable_dist'][~in_bounds_idx] = False
        else:   # No distances provided
            if isinstance(samp_idx, slice):
                ret = np.empty((n_coords_ret, n_samp_ret, self._n_distances), dtype='f2')
            else:
                ret = np.empty((n_coords_ret, self._n_distances), dtype='f2')
            if self._contiguous:
                ret = data_handle(pix_idx, samp_idx, slice(None))
            else:
                if isinstance(samp_idx, slice):
                    for idx, pix in enumerate(pix_idx):
                        ret[idx] = data_handle[pix, samp_idx, slice(None)]
                else:
                    for idx, pix in enumerate(pix_idx):
                        ret[idx] = data_handle[pix, samp_idx[idx], slice(None)]
            ret[~in_bounds_idx, ...] = np.nan

            # Flag: reliable distance bounds
            if return_flags:
                if self._contiguous:
                    flags['min_reliable_distmod'] = DM_reliable_min_handle(pix_idx)
                    flags['max_reliable_distmod'] = DM_reliable_max_handle(pix_idx)
                else:
                    for idx, pix in enumerate(pix_idx):
                        flags['min_reliable_distmod'][idx] = DM_reliable_min_handle[pix]
                        flags['max_reliable_distmod'][idx] = DM_reliable_max_handle[pix]
                flags['min_reliable_distmod'][~in_bounds_idx] = np.nan
                flags['max_reliable_distmod'][~in_bounds_idx] = np.nan

        # Flag: convergence
        if return_flags:
            if self._contiguous:
                flags['converged'] = converged_handle( pix_idx).astype(bool)
                flags['infilled'] = infilled_handle(pix_idx).astype(bool)
            else:
                for idx, pix in enumerate(pix_idx):
                    flags['converged'][idx] = converged_handle[pix].astype(bool)
                    flags['infilled'][idx] = infilled_handle[pix].astype(bool)
            flags['converged'][~in_bounds_idx] = False
            flags['infilled'][~in_bounds_idx] = False

        # Reduce the samples in the requested manner
        if mode == 'median':
            ret = np.median(ret, axis=1)
        elif mode == 'mean':
            # Remove "samples" axis
            s = ret.shape
            ret = np.reshape(ret, s[:1] + s[2:])
        elif mode == 'percentile':
            ret = np.nanpercentile(ret, pct, axis=1)
            if not scalar_pct:
                # (percentile, pixel) -> (pixel, percentile)
                # (pctile, pixel, distance) -> (pixel, distance, pctile)
                ret = np.moveaxis(ret, 0, -1)
        elif mode == 'samples':
            # Swap sample and distance axes to be consistent with other 3D dust
            # maps. The output shape will be (pixel, distance, sample).
            if not has_dist:
                np.swapaxes(ret, 1, 2)
        	
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


def _prepare_download(local_fname, file_requirement, h5_size, h5_dsets,
                      clobber, silence_warnings, n_gb, hint=None):
    """
    Makes sure that there is something on disk for :obj:`h5_repack` to read,
    downloading it if necessary.

    Returns the name of the file to repack, or :obj:`None` if the repacked file
    is already in place and :obj:`clobber` is :obj:`False`.
    """
    # Check if a repacked file already exists. The size check alone is not
    # enough, because the file as published is similar in size to the repacked
    # file, so the layout has to be checked as well.
    if not clobber:
        if (h5_is_repacked(local_fname)
                and fetch_utils.h5_file_exists(local_fname, h5_size,
                                               dsets=h5_dsets)):
            print('File appears to exist already. Call `fetch(clobber=True)` '
                  'to force overwriting of existing file.')
            return None

    # The original, as-published file is the input to the repacking. If such a
    # file is already on disk -- for example, one downloaded by an older
    # version of `dustmaps` -- then repack it, rather than download it again.
    # The file is checked over before it is reused, so that a truncated
    # download is replaced, rather than repacked.
    if (os.path.isfile(local_fname)
            and not h5_is_repacked(local_fname)
            and fetch_utils.h5_file_exists(local_fname, dsets=h5_dsets)):
        print('Found an existing file to repack: {}'.format(local_fname))
        return local_fname

    if not silence_warnings:
        print('Warning: You are about to download a large file ({} GB).'
              .format(n_gb))
        if hint is not None:
            print(hint)
        print('Tip: To suppress this warning and skip confirmation in future '
              'runs, use silence_warnings=True.')
        response = input("Do you want to proceed? (Yes/No): ").strip().lower()
        if response != "yes":
            print("Download aborted.")
            return None

    orig_fname = local_fname + '.orig'

    print("Proceeding with the download...")

    # Download the data
    fetch_utils.dataverse_download_doi(
        '10.7910/DVN/J9JCKO',
        orig_fname,
        file_requirements={'filename': file_requirement}
    )

    return orig_fname


def h5_repack(h5_in, h5_out, chunk_pixels=CHUNK_PIXELS,
              compression_opts=COMPRESSION_OPTS, shuffle=SHUFFLE):
    """
    Converts an HDF5 file containing a DECaPS map into the layout that
    :obj:`DECaPSQuery` expects, in which each chunk of the file holds a block of
    adjacent pixels, with all of their samples and distance bins.

    Every dataset of the input file is kept, along with the attributes that
    describe the map, and the output file is marked with a ``repacked``
    attribute recording the number of pixels per chunk.

    The pixels in the original file are already ordered by ``healpix_index``,
    and are written out in the same order.

    Args:
        h5_in (:obj:`str`): Filename of the original HDF5 file.
        h5_out (:obj:`str`): Filename to write the repacked file to.
        chunk_pixels (Optional[:obj:`int`]): Number of pixels stored in each
            chunk. Defaults to :obj:`CHUNK_PIXELS`.
        compression_opts (Optional[:obj:`int`]): Deflate compression level,
            between 1 and 9. Defaults to :obj:`COMPRESSION_OPTS`.
        shuffle (Optional[:obj:`bool`]): Whether to shuffle the bytes of each
            value before compressing it. Defaults to :obj:`SHUFFLE`.
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

        # `samples` holds the reddening samples, and `mean` the mean reddening.
        # A file that was downloaded with `mean_only=True` has only the latter.
        for name in ('mean', 'samples'):
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
                compression_opts=compression_opts,
                shuffle=shuffle)

            # Keep the attributes that describe the dataset, such as its units
            for key in dset_in.attrs:
                dset_out.attrs[key] = dset_in.attrs[key]

            for start in range(0, n_pix, block_pixels):
                stop = min(start + block_pixels, n_pix)
                dset_out[start:stop] = dset_in[start:stop]


def h5_is_repacked(h5_fname):
    """
    Returns ``True`` if the given HDF5 file holds a DECaPS map that has been
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
            name = 'samples' if 'samples' in f else 'mean'
            chunks = f[name].chunks
    except (IOError, KeyError):
        return False

    return bool(repacked) and chunks is not None


# Expected size (in Bytes) and datasets of a repacked file, for each of the two
# files that `fetch()` can download. The size is a rough guide (see
# `fetch_utils.h5_file_exists`, which allows a relative tolerance of 30%); the
# datasets are the check that matters. No size is given for `mean_only`,
# because that file has not been inspected here, and so the shape of its `mean`
# dataset is not checked either.
REPACKED_SIZES = {
    # Measured for the file published on the Dataverse on 2026-10-09, before the
    # authors replaced it with a repacked one
    'mean_and_samples': 18370195128,
    'mean_only': None
}

REPACKED_DSETS = {
    'mean_and_samples': {
        'mean': (51415400, 1, 120),
        'samples': (51415400, 5, 120),
        'pixel_info': None
    },
    'mean_only': {
        'mean': None,
        'pixel_info': None
    }
}


def fetch(mean_only=False, silence_warnings=False, clobber=False):
    """
    Downloads the specified version of the DECaPS dust map.

    The downloaded file is repacked into the layout that :obj:`DECaPSQuery`
    expects, in which each chunk holds a block of adjacent pixels, with all of
    their samples and distance bins. The original file is then deleted, since
    the repacked file is equivalent to it, and is nearly a third smaller.
    
    Args:
        mean_only (Optional[bool]): If True, only the mean map (8 GB) will be downloaded 
            and available to query. If False (the default), both the mean and samples
            will be downloaded (33 GB) and available to query.
        silence_warnings (Optional[bool]): If True, suppresses all warnings and proceeds 
            without requiring user confirmation. Defaults to False.
        clobber (Optional[bool]): If True, overwrites any existing files. Defaults to False.

    Raises:
        DownloadError: Either no matching file was found under the given DOI, or
            the MD5 sum of the file was not as expected.
        requests.exceptions.HTTPError: The given DOI does not exist, or there
            was a problem connecting to the Dataverse.
    """

    if not mean_only:
        key = 'mean_and_samples'
        n_gb = 33
        hint = ('If you only want the mean map file (8 GB), use '
                'mean_only=True.')
        local_fname = os.path.join(data_dir(), 'decaps',
                                   'decaps_mean_and_samples.h5')
        file_requirement = 'decaps_mean_and_samples.h5'
    else:
        key = 'mean_only'
        n_gb = 8
        hint = None
        local_fname = os.path.join(data_dir(), 'decaps', 'decaps_mean.h5')
        file_requirement = 'decaps_mean.h5'

    orig_fname = _prepare_download(
        local_fname, file_requirement, REPACKED_SIZES[key],
        REPACKED_DSETS[key], clobber, silence_warnings, n_gb, hint)

    if orig_fname is None:
        return

    # Convert to the layout that DECaPSQuery expects. The repacked file is
    # written to a temporary name, so that the original is left in place if
    # anything goes wrong.
    print('Repacking file...')
    h5_repack(orig_fname, local_fname + '.repacked')
    os.replace(local_fname + '.repacked', local_fname)

    # Cleanup
    if orig_fname != local_fname:
        print('Removing original file...')
        os.remove(orig_fname)


def example_plot(fname, max_samples=None, memmap=True):
    """
    Example plot of the DECaPS dust map.

    Args:
        fname (:obj:`str`): Filename to write the plot to.
        max_samples (Optional[:obj:`int`]): Maximum number of samples to use.
            Defaults to :obj:`None`, meaning that all samples are used. Note
            that the plot queries over half a million coordinates, and that
            loading all of the samples requires about 62 GB of RAM (this map
            has 51.4 million pixels), so on a smaller machine a small number of
            samples (e.g., 1) should be passed instead.
        memmap (Optional[:obj:`bool`]): If :obj:`True` (the default), the
            reddening samples are memory-mapped. That is the default here, and
            not for the other maps, because loading them into memory needs
            about 62 GB: on a machine that cannot do that, reading them into
            memory is not merely slower, it is impossible. With a repacked file
            the query above reads only a few hundred MB, so memory mapping is
            fast as well as frugal.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import PowerNorm
    from matplotlib.ticker import FuncFormatter
    from astropy.coordinates import SkyCoord
    from astropy import units

    plt.rcParams.update({
        'font.size': 6,
        'axes.titlesize': 6,
        'axes.labelsize': 6,
        'xtick.labelsize': 5,
        'ytick.labelsize': 5,
        'legend.fontsize': 5,
    })

    q = DECaPSQuery(max_samples=max_samples, memmap=memmap)

    # Query a grid of coordinates covering the footprint of the map, which
    # spans 239 < l < 6 deg and |b| < 10 deg. Longitude is unwrapped here (so
    # that it runs from 239 to 366 deg), in order to avoid the discontinuity at
    # l = 0, and the x axis is reversed below, so that longitude increases to
    # the left, as is conventional.
    #
    # The grid is much coarser than the map's pixels (which are 0.007 deg
    # across), because each coordinate that falls in its own chunk costs one
    # read: at 1024 x 160 the whole plot takes a few minutes, while a grid fine
    # enough to resolve the pixels would take hours.
    l = np.linspace(239., 366., 1024)
    b = np.linspace(-10., 10., 160)

    l, b = np.meshgrid(l, b, indexing='ij')

    coords = SkyCoord(
        np.mod(l, 360.)*units.deg, b*units.deg,
        distance=1.0*units.kpc,
        frame='galactic'
    )

    modes = ['mean', 'median', 'random_sample_per_pix']

    # Set up the figure. It is shorter than the Bayestar example plot, because
    # the footprint of this map is a narrow strip of sky.
    fig,axes = plt.subplots(2,2, figsize=(6,3.2), constrained_layout=True)

    # Query the dust map for each mode. These three panels should look almost
    # identical; if they do not, something is wrong.
    for ax,m in zip(axes.flat, modes):
        E = q(coords, mode=m)
        im = ax.imshow(
            E.T,
            origin='lower',
            extent=[239., 366., -10., 10.],
            norm=PowerNorm(0.5, vmin=0, vmax=2.0)
        )
        ax.invert_xaxis()
        # Longitude is unwrapped above, so label the ticks in the range
        # 0 to 360 deg.
        ax.xaxis.set_major_formatter(
            FuncFormatter(lambda x, pos: '{:.0f}'.format(np.mod(x, 360.))))
        ax.set_xlabel(r'$\ell$ (deg)')
        ax.set_ylabel(r'$b$ (deg)')
        ax.set_title('mode = {}'.format(m))
        ax.set_aspect('equal')

    # One colorbar, shared by the three panels of the sky
    fig.colorbar(im, ax=axes.flat[:3], label=r'$E(B-V)$ (mag)')

    # Query selected lines of sight
    sightline_names = ['Vela', 'Carina', 'Sagittarius B2']
    coords = SkyCoord(
        [266.0, 287.6, 0.7]*units.deg,
        [-1.0, -0.6, -0.05]*units.deg,
        frame='galactic'
    )
    E = q(coords, mode='mean')
    ax = axes[1,1]
    for k, (name, e) in enumerate(zip(sightline_names, E)):
        ax.plot(q.distances, e, ls=['-', '--', ':'][k % 3], lw=1, label=name)
    ax.legend()
    ax.set_title('Selected sightlines')
    ax.set_xlabel(r'$r$ (kpc)')
    ax.set_ylabel(r'$E$ (mag)')
    ax.set_xlim(0, 10.0)
    ax.grid(True, alpha=0.1)

    # Save and close figure
    fig.savefig(fname, dpi=300)
    plt.close(fig)


def diagnostic(fname=None, n_coords=(16, 256, 4096), max_samples=1, seed=0,
               in_memory=False):
    """
    Prints the time taken to load and to query the DECaPS dust map, along with
    the memory that it used, for checking by eye before a release that nothing
    has become dramatically slower or more memory-hungry.

    ``memmap=True`` is measured by default, because this map is far too large to
    read into RAM: it has 51.4 million pixels, so a single sample of each of them
    already needs 12 GB for the samples plus 12 GB for the mean map, and all five
    samples need 62 GB. Pass ``in_memory=True`` to measure ``memmap=False`` as
    well, if the machine has the memory for it.

    Args:
        fname (Optional[:obj:`str`]): Filename of the DECaPS map. Defaults to
            :obj:`None`, meaning that the default location is used.
        n_coords (Optional[:obj:`list` or :obj:`tuple`]): Numbers of
            coordinates to query. Defaults to ``(16, 256, 4096)``. The
            coordinates are all drawn at once, so that the shorter queries use
            the first few of the longer ones.
        max_samples (Optional[:obj:`int`]): Maximum number of samples to load.
            Defaults to 1, because this map has 51.4 million pixels: loading a
            single sample of each of them needs about 12 GB of RAM, and all
            five need about 62 GB. Pass :obj:`None` to use all of the samples,
            if the machine can take it. The memory-mapped numbers below are the
            ones that matter for this map.
        seed (Optional[:obj:`int`]): Seed for the random coordinates, so that
            runs can be compared. This also seeds numpy's global random state,
            so that the ``random_sample`` mode gives the same answer on each
            run.
        in_memory (Optional[:obj:`bool`]): If :obj:`True`, the in-RAM case
            (``memmap=False``) is measured as well. Defaults to :obj:`False`,
            because loading even one sample of this map takes about 25 GB.
    """
    import resource
    import sys

    def peak_memory():
        """Peak memory used by this process so far, in MB."""
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # ru_maxrss is in bytes on macOS, and in kB everywhere else
        return peak / 1024.**2 if sys.platform == 'darwin' else peak / 1024.

    rng = np.random.RandomState(seed)
    np.random.seed(seed)

    # Coordinates drawn from within the footprint of the map (239 < l < 6,
    # |b| < 10 deg), with distances drawn uniformly in log between 0.1 and
    # 10 kpc
    n_max = max(n_coords)
    l = np.mod(rng.uniform(239., 366., n_max), 360.)
    b = rng.uniform(-10., 10., n_max)
    d = 10.**rng.uniform(-1., 1., n_max)
    coords = coordinates.SkyCoord(
        l*units.deg, b*units.deg, distance=d*units.kpc, frame='galactic')

    if fname is None:
        fname = os.path.join(data_dir(), 'decaps',
                             'decaps_mean_and_samples.h5')

    info = {}
    modes = (True, False) if in_memory else (True,)

    for memmap in modes:
        peak_before = peak_memory()

        t0 = time()
        q = DECaPSQuery(fname, max_samples=max_samples, memmap=memmap)
        t_load = time() - t0

        stats = dict((mode, []) for mode in ('mean', 'random_sample'))
        mean = None
        for mode in ('mean', 'random_sample'):
            for n in n_coords:
                coords_n = coords[:n]

                t0 = time()
                val = q.query(coords_n, mode=mode)
                t_cold = time() - t0

                # The shortest of several repeats, so that interference from
                # other work on the machine does not look like a slow-down
                t_warm = []
                for _ in range(3):
                    t0 = time()
                    q.query(coords_n, mode=mode)
                    t_warm.append(time() - t0)

                stats[mode].append((n, t_cold, min(t_warm)))

                if mode == 'mean' and n == n_max:
                    mean = val

        in_map = np.isfinite(mean)

        info[memmap] = dict(
            load=t_load,
            stats=stats,
            peak=peak_memory() - peak_before,
            fraction=np.mean(in_map),
            mean_val=np.mean(mean[in_map]) if np.any(in_map) else np.nan,
            mean=mean)

        del q

    print('')
    print('decaps {}'.format(os.path.basename(fname)))
    print('  file       {:>8.1f} MB  {}'.format(
        os.path.getsize(fname) / 1e6, fname))

    with h5py.File(fname, 'r') as f:
        for name in ('mean', 'samples'):
            if name in f:
                dset = f[name]
                print('  {:<10} chunks {}, {}'.format(
                    name, dset.chunks,
                    dset.compression if dset.compression is not None
                    else 'uncompressed'))

    for memmap in modes:
        i = info[memmap]
        print('')
        print('  memmap={}'.format(memmap))
        print('    {:<14} {:>9.2f} s'.format('load', i['load']))
        for mode in ('mean', 'random_sample'):
            print('    {}'.format(mode))
            for n, t_cold, t_warm in i['stats'][mode]:
                print('      {:>6} coords {:>9.4f} s cold,'
                      ' {:>9.4f} s warm'.format(n, t_cold, t_warm))
        print('    {:<14} {:>9.1f} MB'.format('peak increase', i['peak']))

    print('')
    print('  {:.1f}% of coordinates are in the map; mean {:.3f} mag'.format(
        100.*info[True]['fraction'], info[True]['mean_val']))
    if in_memory:
        print('  memmap=True and memmap=False agree: {}'.format(
            np.array_equal(info[True]['mean'], info[False]['mean'],
                           equal_nan=True)))
    else:
        with h5py.File(fname, 'r') as f:
            itemsize = np.dtype(f['mean'].dtype).itemsize
            n_pix, n_dist = f['mean'].shape[0], f['mean'].shape[2]
            n_samp = 1 if max_samples is None else min(
                max_samples, f['samples'].shape[1]) if 'samples' in f else 0
        needed = (n_pix * n_dist + n_pix * n_samp * n_dist) * itemsize
        print('  memmap=False was not measured: it needs about {:.0f} GB of '
              'RAM (pass `in_memory=True` to measure it anyway).'.format(
                  needed / 1e9))
