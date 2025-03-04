#!/usr/bin/env python
#
# bayestar.py
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
    Queries the DECaPS 3D dust maps (Zucker, Saydjari, & Speagle et al. 2025)
    **without memory mapping**, requiring the entire file be loaded into RAM. This will 
    be significantly faster downstream if you need to execute many large queries, but you will
    pay a startup cost up front loading many GBs into RAM. The map covers the southern 
    Galactic plane (239 < l < 6, |b| < 10), amounting to 6% of the sky. When combined 
    with the BayestarQuery, this DECaPS query enables reddening estimates over the 
    entire Galactic plane |b| < 10.
    """

    def __init__(self, map_fname=None, max_samples=None, mean_only=False):
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
        """

        self._mean_only = mean_only

        if self._mean_only:
            print("You are about to load 7 GB into RAM. If you are not performing intensive queries, consider using the DECaPSQueryLite class.")

            if map_fname is None:
                map_fname = os.path.join(data_dir(), 'decaps', 'decaps_mean.h5')
            if not os.path.isfile(map_fname):
                map_fname = os.path.join(data_dir(), 'decaps', 'decaps_mean_and_samples.h5')
        else:
            print("You are about to load 30 GB into RAM. If you are not performing intensive queries, consider using the DECaPSQueryLite class.")

            if map_fname is None:
                map_fname = os.path.join(data_dir(), 'decaps', 'decaps_mean_and_samples.h5')
                if not os.path.isfile(map_fname):
                    raise ValueError(
                        "The file containing both the mean and samples was not found at the default location on disk. "
                        "Please confirm you have downloaded 'decaps_mean_and_samples.h5'."
                    )

        with h5py.File(map_fname, 'r') as f:
            print("Loading pixel info...")
            self._DM_bin_edges = f['/pixel_info'].attrs['DM_bin_edges']
            self._n_distances = len(self._DM_bin_edges)
            self._n_pix = f['/pixel_info/healpix_index'].size
            self._chunk_size = 100000

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


            self._nside = f['/pixel_info'].attrs['nside']
            self._hp_idx_sorted = f['/pixel_info/healpix_index'][:]
            self._data_idx = np.arange(len(self._hp_idx_sorted))
            self._pixel_info = {name: ds[()] for name, ds in f['pixel_info'].items()}


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
            mode (Optional[:obj:`str`]): Seven different query modes are available:
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

        if mode == 'mean':
        	val = self._mean
        else:
            val = self._samples

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
                         ('min_reliable_distmod', 'f2'),
                         ('max_reliable_distmod', 'f2')]
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
                        * val[pix_idx[idx_near], samp_idx, 0])
                else:

                    ret[idx_near] = (
                        a * val[pix_idx[idx_near], samp_idx[idx_near], 0])

            # d > d(farthest distance slice)
            idx_far = (bin_idx_ceil == self._n_distances) & in_bounds_idx
            if np.any(idx_far):

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
            ret.shape = s[:1] + s[2:]
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
            mode (Optional[:obj:`str`]): Seven different query modes are available:
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
                         ('min_reliable_distmod', 'f2'),
                         ('max_reliable_distmod', 'f2')]
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
            ret.shape = s[:1] + s[2:]
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


def fetch(mean_only=False, silence_warnings=False, clobber=False):
    """
    Downloads the specified version of the DECaPS dust map.
    
    Args:
        mean_only (Optional[bool]): If True, only the mean map (7 GB) will be downloaded 
            and available to query. If False (the default), both the mean and samples
            will be downloaded (30 GB) and available to query.
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
        local_fname = os.path.join(data_dir(), 'decaps', 'decaps_mean_and_samples.h5')
        
        if os.path.isfile(local_fname) and not clobber:
            print(f"File '{local_fname}' already exists. Skipping download. To overwrite, use clobber=True.")
            return

        if not silence_warnings:
            print(
                "Warning: You are about to download a large file (30 GB), containing the mean map and samples.\n"
                "If only want to download the mean map file (7 GB), use mean_only=True."
            )
            print("Tip: To suppress this warning and skip confirmation in future runs, use silence_warnings=True.")
            response = input("Do you want to proceed? (Yes/No): ").strip().lower()
            if response != "yes":
                print("Download aborted.")
                return
        
        print("Proceeding with the download...")
        doi = '10.7910/DVN/J9JCKO'

        # Download the data
        fetch_utils.dataverse_download_doi(
            doi,
            local_fname,
            file_requirements={'filename': 'decaps_mean_and_samples.h5'}
        )
    
    else:
        local_fname = os.path.join(data_dir(), 'decaps', 'decaps_mean.h5')

        if os.path.isfile(local_fname) and not clobber:
            print(f"File '{local_fname}' already exists. Skipping download. To overwrite, use clobber=True.")
            return

        if not silence_warnings:
            print("Warning: You are about to download a large file (7 GB) containing the mean map.")
            print("Tip: To suppress this warning and skip confirmation in future runs, use silence_warnings=True.")
            response = input("Do you want to proceed? (Yes/No): ").strip().lower()
            if response != "yes":
                print("Download aborted.")
                return
        
        print("Proceeding with the download...")
        doi = '10.7910/DVN/J9JCKO'

        # Download the data
        fetch_utils.dataverse_download_doi(
            doi,
            local_fname,
            file_requirements={'filename': 'decaps_mean.h5'}
        )
