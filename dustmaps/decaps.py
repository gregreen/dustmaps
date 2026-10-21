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
from .rng_utils import get_rng
from . import fetch_utils
from . import h5_utils
import warnings

from time import time
from tqdm import tqdm


# The layout that `DECaPSQuery` expects: each chunk of the file holds a block of
# adjacent pixels, their samples, and a slice of the distance axis. The published
# file instead stores a whole slab of 200842 pixels for one (sample, distance)
# pair at a time -- a layout that suits reading the map from end to end, which is
# how the original code loaded it, but that makes a query for a random coordinate
# read about 240 MB, and a 256-coordinate query read the whole 62 GB file.
#
# Repacking also shrinks the file by 46%: the published 33.84 GB becomes 18.37 GB,
# because the extinction along a line of sight is smoother than it is from pixel
# to pixel, so encoding a block of pixels together compresses better than encoding
# runs of single values does. The gain is almost all in the samples (0.491 ->
# 0.212 kB per pixel for `samples`, 0.145 -> 0.116 for `mean`). The published
# compression settings were kept (deflate, with byte shuffling), except that the
# level is raised from 3 to 7, to match the rest of `dustmaps`.
#
# The shape of a chunk is a trade-off between the two ways of querying the map.
# A query for scattered sightlines wants chunks with few pixels, so that each one
# returns a small block, while a query over a contiguous area at one distance
# wants chunks that hold many pixels and few distance bins, so that the bins it
# does not need are not read along with it. Sixteen pixels with all 120 bins was
# the best point of that trade-off for a query that read every bin and picked two
# out of the array, which is what the first version of this reader did.
#
# A query with a distance now asks for the two bins it interpolates between (see
# `DECaPSQuery.query`), so the distance axis is worth chunking. Measured on a
# wedge of the map (l = 350-360 deg, |b| < 5, 1.95 million pixels), with the same
# bytes in a chunk either way:
#
#     pattern                            (16, 120)  (64, 30)  (128, 15)
#     contiguous area, one distance         3.23 s    1.27 s    1.66 s
#     scattered sightlines, one distance    1.35 s    1.12 s    1.04 s
#     scattered sightlines, all distances   1.10 s    2.26 s    3.03 s
#
# Four distance chunks is the best of the three on the first pattern, which is the
# one an example plot of the map exercises, and costs the least on the last. Eight
# is worse than four even on the pattern it was meant to help: a sweep at one
# distance then touches one chunk in eight, so its reads stride over pages of the
# file instead of staying inside them.
#
# Repacking the published file on a laptop takes about 45 minutes (deflate is
# slow), and `scripts/repack_decaps.py` does the same job for the authors of the
# map.
CHUNK_PIXELS = 64
CHUNK_DISTANCES = 30
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

        # The pixel index has just been loaded, so leave it out of here: a
        # second copy of it would be another 411 MB of RAM on this map.
        self._pixel_info = {name: ds[()]
                            for name, ds in f['pixel_info'].items()
                            if name != 'healpix_index'}

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
    def query(self, coords, mode='mean', return_flags=False, pct=None,
              seed_or_rng=None):
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
            seed_or_rng (Optional[:obj:`int` or :obj:`numpy.random.Generator`]): A
                seed for, or an existing generator of, the random numbers drawn
                by the 'random_sample' and 'random_sample_per_pix' modes. An
                integer makes each call draw the same samples, so that a result
                can be repeated; a generator is drawn from as it stands, so that
                successive calls give different samples. Defaults to
                :obj:`None`, meaning a fresh, unpredictable generator.

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

        # Extract the correct samples. The two random modes draw from a
        # generator; the modes that do not draw anything ignore it.
        if mode in ('random_sample', 'random_sample_per_pix'):
            rng = get_rng(seed_or_rng)

        if mode == 'random_sample':
            # A different sample in each queried coordinate
            samp_idx = rng.integers(0, self._n_samples, pix_idx.size)
            n_samp_ret = 1
        elif mode == 'random_sample_per_pix':
            # Choose same sample in all coordinates that fall in same angular
            # HEALPix pixel
            samp_idx = rng.integers(0, self._n_samples,
                                    self._n_pix)[pix_idx]
            n_samp_ret = 1
        elif mode == 'mean':
            samp_idx = slice(None)
            n_samp_ret = 1
        else:
            # Return all samples in each queried coordinate
            samp_idx = slice(None)
            n_samp_ret = self._n_samples

        # The rows of the map that the query touches are read in the loop below,
        # one group at a time, and reduced into `ret` as they arrive. The mean
        # map carries an empty samples axis; the samples map is read only for
        # the samples that the query mode asked for.
        if mode == 'mean':
            dset = self._mean
            columns = (slice(0, 1),)
        else:
            dset = self._samples
            columns = (slice(0, self._n_samples),)

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
            # Return all distances, for the samples that were requested. As
            # before, the results take their precision from the map itself.
            if isinstance(samp_idx, slice):
                ret = np.full((n_coords_ret, n_samp_ret, self._n_distances),
                              np.nan, dtype=dset.dtype)
            else:
                ret = np.full((n_coords_ret, self._n_distances), np.nan,
                              dtype=dset.dtype)

            # Flag: reliable distance bounds
            if return_flags:
                dm_min = self._pixel_info['DM_reliable_min'][pix_idx]
                dm_max = self._pixel_info['DM_reliable_max'][pix_idx]

                flags['min_reliable_distmod'] = dm_min
                flags['max_reliable_distmod'] = dm_max
                flags['min_reliable_distmod'][~in_bounds_idx] = np.nan
                flags['max_reliable_distmod'][~in_bounds_idx] = np.nan

        # With a distance given, only the two bins either side of it are worth
        # reading, and asking for them coordinate by coordinate keeps the read to
        # the chunks that hold them. The bins are given in order, so the lower
        # one comes first, and a coordinate beyond either end of the axis asks
        # for its own bin twice.
        if has_dist:
            bins = np.empty((n_coords_ret, 2), dtype='i4')
            bins[:, 0] = np.clip(bin_idx_ceil - 1, 0, self._n_distances - 1)
            bins[:, 1] = np.clip(bin_idx_ceil, 0, self._n_distances - 1)
            distance_selection = (bins,)
        else:
            distance_selection = (slice(None),)

        # Read the rows of the map that the query touches, one group at a time,
        # and reduce them into the results as they arrive. Because the data is
        # never accumulated, a query for the whole map needs no more memory than
        # a query for a single coordinate. Coordinates outside the map keep the
        # NaN values that the results were filled with.
        for coords_here, val in h5_utils.iter_nd_points(
                dset, (pix_idx,) + columns + distance_selection, in_bounds_idx):
            if not has_dist:
                if isinstance(samp_idx, slice):
                    ret[coords_here] = val
                else:
                    ret[coords_here] = val[np.arange(coords_here.size),
                                           samp_idx[coords_here], :]
                continue

            dm_here = dm[coords_here]
            bin_here = bin_idx_ceil[coords_here]

            if isinstance(samp_idx, slice):
                samp_here = samp_idx
            else:
                samp_here = samp_idx[coords_here]

            # The two bins that each coordinate asked for
            low = val[..., 0]
            high = val[..., 1]

            # d < d(nearest distance slice)
            idx_near = (bin_here == 0)
            if np.any(idx_near):
                a = 10.**(0.2 * (dm_here[idx_near] - self._DM_bin_edges[0]))
                if isinstance(samp_idx, slice):
                    ret[coords_here[idx_near]] = (
                        a[:,None] * low[idx_near])
                else:
                    ret[coords_here[idx_near]] = (
                        a * low[idx_near, samp_here[idx_near]])

            # d > d(farthest distance slice)
            idx_far = (bin_here == self._n_distances)
            if np.any(idx_far):
                if isinstance(samp_idx, slice):
                    ret[coords_here[idx_far]] = high[idx_far]
                else:
                    ret[coords_here[idx_far]] = high[idx_far,
                                                     samp_here[idx_far]]

            # d(nearest distance slice) < d < d(farthest distance slice)
            idx_btw = ~idx_near & ~idx_far
            if np.any(idx_btw):
                DM_ceil = self._DM_bin_edges[bin_here[idx_btw]]
                DM_floor = self._DM_bin_edges[bin_here[idx_btw]-1]
                a = (DM_ceil - dm_here[idx_btw]) / (DM_ceil - DM_floor)
                if isinstance(samp_idx, slice):
                    ret[coords_here[idx_btw]] = (
                        (1.-a[:,None]) * high[idx_btw]
                        + a[:,None] * low[idx_btw]
                    )
                else:
                    ret[coords_here[idx_btw]] = (
                        (1.-a) * high[idx_btw, samp_here[idx_btw]]
                        +    a * low[idx_btw, samp_here[idx_btw]]
                    )

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
    def query(self, coords, mode='mean', return_flags=False, pct=None,
              seed_or_rng=None):
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
            seed_or_rng (Optional[:obj:`int` or :obj:`numpy.random.Generator`]): A
                seed for, or an existing generator of, the random numbers drawn
                by the 'random_sample' and 'random_sample_per_pix' modes. An
                integer makes each call draw the same samples, so that a result
                can be repeated; a generator is drawn from as it stands, so that
                successive calls give different samples. Defaults to
                :obj:`None`, meaning a fresh, unpredictable generator.

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

        # Extract the correct samples. The two random modes draw from a
        # generator; the modes that do not draw anything ignore it.
        if mode in ('random_sample', 'random_sample_per_pix'):
            rng = get_rng(seed_or_rng)

        if mode == 'random_sample':
            # A different sample in each queried coordinate
            samp_idx = rng.integers(0, self._n_samples, pix_idx.size)
            n_samp_ret = 1
        elif mode == 'random_sample_per_pix':
            # Choose same sample in all coordinates that fall in same angular
            # HEALPix pixel
            samp_idx = rng.integers(0, self._n_samples,
                                    self._n_pix)[pix_idx]
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


def _md5_matches(fname, md5sum):
    """
    Returns ``True`` if a file matches the MD5 sum that is expected for it.

    A sum of :obj:`None` means that the file has not been inspected here, so
    there is nothing to check it against, and it counts as a match.
    """
    if md5sum is None:
        return True

    return fetch_utils.check_md5sum(fname, md5sum)


def _prepare_download(local_fname, file_requirement, h5_size, h5_dsets,
                      published_md5, clobber, silence_warnings, n_gb, hint=None):
    """
    Makes sure that there is something on disk for :obj:`h5_repack` to read,
    downloading it if necessary.

    A repacked file that is already in place is accepted if it carries the
    ``repacked`` attribute, holds the datasets that are expected with the right
    shapes, and is of about the right size. A file in the published layout is
    checked the same way, and also against the MD5 sum that the Dataverse gives
    for it, before it is handed over to the repacker.

    Returns the name of the file to repack, or :obj:`None` if the repacked file
    is already in place and :obj:`clobber` is :obj:`False`.
    """
    # Check if a repacked file already exists. What marks it out is the
    # `repacked` attribute, the datasets inside it and their shapes, and its
    # size, which is checked to within 30%. No checksum is imposed on it: a file
    # that a user repacked, or that the authors of a map uploaded, need not be
    # byte for byte the file that `dustmaps` would write, since a different
    # version of HDF5 can lay the same data out in different bytes.
    if not clobber:
        if (h5_is_repacked(local_fname)
                and fetch_utils.h5_file_exists(local_fname, h5_size,
                                               dsets=h5_dsets)):
            print('File appears to exist already. Call `fetch(clobber=True)` '
                  'to force overwriting of existing file.')
            return None

    # The original, as-published file is the input to the repacking. If such a
    # file is already on disk -- for example, one downloaded by an older
    # version of `dustmaps` -- then repack it, rather than download it again,
    # but only after checking it against the MD5 sum that the Dataverse gives
    # for it. That sum is of the bytes the Dataverse serves, so it can be relied
    # on. Checking it reads the whole file, which takes about 75 s for 34 GB,
    # and it is worth every one of them: repacking takes 45 minutes, and a
    # download that was interrupted halfway through would waste all of them, and
    # leave behind an 18 GB file built from a corrupt one, which nothing
    # downstream would notice.
    if (os.path.isfile(local_fname)
            and not h5_is_repacked(local_fname)
            and fetch_utils.h5_file_exists(local_fname, dsets=h5_dsets)):
        print('Checking the checksum of the file to repack: {}'
              .format(local_fname))
        if _md5_matches(local_fname, published_md5):
            print('Found an existing file to repack: {}'.format(local_fname))
            return local_fname

        print('Warning: {} does not have the MD5 sum that the Dataverse gives '
              'for a published file. It may be incomplete, or corrupted, or of '
              'a different version, so it will be downloaded again.'.format(
                  local_fname))

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
              compression_opts=COMPRESSION_OPTS, shuffle=SHUFFLE,
              max_memory=256 << 20, chunk_distances=CHUNK_DISTANCES):
    """
    Converts an HDF5 file containing a DECaPS map into the layout that
    :obj:`DECaPSQuery` expects, in which each chunk of the file holds a block of
    adjacent pixels, with their samples and a slice of the distance axis.

    Every dataset of the input file is kept, along with the attributes that
    describe the map, and the output file is marked with a ``repacked``
    attribute recording the number of pixels per chunk. A dataset is written a
    block of pixels at a time, so a map larger than memory can be repacked, and
    each reports how far along it is with a progress bar. The number of pixels
    in a block follows from `max_memory`.

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
        max_memory (Optional[:obj:`int`]): Bytes of a map to hold in memory at
            a time, which sets how many pixels are written in one block. The
            larger the block, the fewer, larger writes a dataset is made of.
            Defaults to 256 MiB.
        chunk_distances (Optional[:obj:`int`]): Number of distance bins stored
            in each chunk, or :obj:`None` to keep the whole distance axis in
            one chunk. Defaults to :obj:`CHUNK_DISTANCES`.
    """

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

            # The last axis is the distance axis, which a query with a distance
            # reads a slice of. The axes between it and the pixels -- the
            # samples, and nothing else in practice -- are kept whole.
            if chunk_distances is None:
                chunks = (chunk_pixels,) + dset_in.shape[1:]
            else:
                chunks = ((chunk_pixels,) + dset_in.shape[1:-1]
                          + (min(chunk_distances, dset_in.shape[-1]),))

            dset_out = f_out.create_dataset(
                name,
                shape=dset_in.shape,
                dtype=dset_in.dtype,
                chunks=chunks,
                compression='gzip',
                compression_opts=compression_opts,
                shuffle=shuffle)

            # Keep the attributes that describe the dataset, such as its units
            for key in dset_in.attrs:
                dset_out.attrs[key] = dset_in.attrs[key]

            # The pixels of a dataset, as many as the memory allows for one
            # block of them: a whole pixel is all of its samples and distance
            # bins, or just its distance bins for the mean map
            row_bytes = int(np.prod(dset_in.shape[1:])) * dset_in.dtype.itemsize
            block_pixels = max(1, int(max_memory) // max(1, row_bytes))

            for start in tqdm(range(0, n_pix, block_pixels),
                              desc="Repacking {}".format(name)):
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

# The MD5 sum that the Dataverse gives for the file as published, used to check
# a file that is already on disk before it is repacked -- here the sum of the
# file published on 2026-10-09, which is also the sum that
# `dataverse_download_doi` checks a download against. It is the sum of the whole
# file, because the check has to come before the repacking, which takes 45
# minutes: a download that was interrupted halfway through would waste all of
# them, and leave behind an 18 GB file built from a corrupt one. Reading the
# whole file to check it takes about 75 s, which is nothing next to that. The
# sum can be relied on, because it is of bytes that the Dataverse serves as they
# are. ``None`` means that the file has not been inspected here, so there is
# nothing to check it against.
#
# No checksum is imposed on a *repacked* file, whether it was downloaded from
# the authors or written here: two HDF5 files can hold the same data and differ
# in their bytes, so a file that a user repacked would be rejected for no good
# reason. Size and layout are what is checked for one.
# `fetch_utils.get_sparse_md5sum` is the tool for checking a huge file cheaply,
# if a case turns up where its bytes can be known in advance.
PUBLISHED_MD5 = {
    'mean_and_samples': '7ccf070772769c8225bc1568b32de76c',
    'mean_only': None
}


def fetch(mean_only=False, silence_warnings=False, clobber=False,
          max_memory=256 << 20):
    """
    Downloads the specified version of the DECaPS dust map.

    The downloaded file is repacked into the layout that :obj:`DECaPSQuery`
    expects, in which each chunk holds a block of adjacent pixels, with all of
    their samples and distance bins. The original file is then deleted, since
    the repacked file is equivalent to it, and is nearly a third smaller.
    A file that is already on disk is checked over before it is reused, so
    that a download that was interrupted, or a file of the wrong version, is
    fetched again rather than used. The file as published, before it is
    repacked, has to match the MD5 sum that the Dataverse gives for it. A
    repacked file is accepted on its ``repacked`` attribute, its datasets and
    their shapes, and its size: no checksum is imposed on it, since a file
    repacked elsewhere need not be byte for byte the one that `dustmaps` would
    write.
    
    Args:
        mean_only (Optional[bool]): If True, only the mean map (8 GB) will be downloaded 
            and available to query. If False (the default), both the mean and samples
            will be downloaded (33 GB) and available to query.
        silence_warnings (Optional[bool]): If True, suppresses all warnings and proceeds 
            without requiring user confirmation. Defaults to False.
        clobber (Optional[bool]): If True, overwrites any existing files. Defaults to False.
        max_memory (Optional[:obj:`int`]): Bytes of the map to hold in memory
            at a time while repacking, which sets how many pixels are written
            in one block. The larger the block, the fewer, larger writes the
            repacked file is made of. Defaults to 256 MiB.

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
        local_fname, file_requirement,
        h5_size=REPACKED_SIZES[key],
        h5_dsets=REPACKED_DSETS[key],
        published_md5=PUBLISHED_MD5[key],
        clobber=clobber, silence_warnings=silence_warnings, n_gb=n_gb,
        hint=hint)

    if orig_fname is None:
        return

    # Convert to the layout that DECaPSQuery expects. The repacked file is
    # written to a temporary name, so that the original is left in place if
    # anything goes wrong.
    print('Repacking file...')
    h5_repack(orig_fname, local_fname + '.repacked', max_memory=max_memory)
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
            that the plot queries over 160000 coordinates, and that loading all
            of the samples requires about 62 GB of RAM (this map has 51.4
            million pixels), so on a smaller machine a small number of samples
            (e.g., 1) should be passed instead.
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
        seed (Optional[:obj:`int`]): Seed for the random coordinates, and for
            the map queries, so that runs can be compared. The
            ``random_sample`` mode gives the same answer on each run. It has to
            be an integer, rather than a generator: a generator would be drawn
            from by each timed query in turn, so the ``memmap=True`` and
            ``memmap=False`` passes would no longer draw the same samples.
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

    rng = np.random.default_rng(seed)

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
                val = q.query(coords_n, mode=mode, seed_or_rng=seed)
                t_cold = time() - t0

                # The shortest of several repeats, so that interference from
                # other work on the machine does not look like a slow-down
                t_warm = []
                for _ in range(3):
                    t0 = time()
                    q.query(coords_n, mode=mode, seed_or_rng=seed)
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
