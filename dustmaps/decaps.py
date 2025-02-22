#!/usr/bin/env python
#
# bayestar.py
# Reads the DECaPS dust reddening maps, described in
# Zucker, Saydjari, & Speagle et al. 2015.
#
# Copyright (C) 2025  Catherine Zucker and Gregory M. Green
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
    Queries the DECaPS 3D dust maps (Zucker, Saydjari, & Speagle et al. 2025).
    The map covers the southern Galactic plane (238 < l < 6, |b| < 10), amounting
    to 6% of the sky. When combined with the BayestarQuery, this DECaPS query enables
    reddening estimates over the entire Galactic plane |b| < 10.
    """

    def __init__(self, map_fname=None, mmap=True, mean_only=False):
        """
        Args:
            map_fname (Optional[str]): Filename of the DECaPS map. Defaults to None, 
                meaning that the default location is used.

            mmap (Optional[bool]): Whether to implement memory mapping and load only 
                the pixels being accessed for a given query. This saves substantial 
                read-in time and RAM upfront but may be slower for larger or complex 
                queries downstream. Defaults to True.

            mean_only (Optional[bool]): If True, only the mean map file is available 
                to query (no samples). For users who did not download the larger 
                mean_and_samples file, this mode is required. However, you will not be 
                able to query in any other mode ('random_sample', 'random_sample_per_pix', 
                'samples', 'median', or 'percentile'). Defaults to False.
        """

        self._mean_only = mean_only
        self._mmap = mmap

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

        if self._mmap:
            f = h5py.File(map_fname, 'r')
            self._pixel_info = f['/pixel_info'][:]
            self._DM_bin_edges = f['/pixel_info'].attrs['DM_bin_edges']
            self._n_distances = len(self._DM_bin_edges)
            self._n_pix = self._pixel_info.size
            
            self._mean = f['/mean']
            
            if not self._mean_only: 
                self._samples = f['/samples']
                self._n_samples = f['/samples'].shape[1]
            
        else:
            with h5py.File(map_fname, 'r') as f:
                print("Loading pixel info...")
                self._pixel_info = f['/pixel_info'][:]
                self._DM_bin_edges = f['/pixel_info'].attrs['DM_bin_edges']
                self._n_distances = len(self._DM_bin_edges)
                self._n_pix = self._pixel_info.size
                self._chunk_size = 100000

                print("Pre-allocating memory for mean map...")
                mean_dataset = f['/mean']
                self._mean = np.empty_like(mean_dataset)

                print("Loading mean map in chunks...")
                for i in tqdm(range(0, mean_dataset.shape[0], self._chunk_size), desc="Mean Map Loading"):
                    self._mean[i:i + self._chunk_size] = mean_dataset[i:i + self._chunk_size]

                if not self._mean_only:
                    samples_dataset = f['/samples']
                    print("Pre-allocating memory for samples...")
                    self._samples = np.empty_like(samples_dataset)

                    print("Loading samples in chunks...")
                    for i in tqdm(range(0, samples_dataset.shape[0], self._chunk_size), desc="Samples Loading"):
                        self._samples[i:i + self._chunk_size] = samples_dataset[i:i + self._chunk_size]

                    self._n_samples = f['/samples'].shape[1]  # Fix indentation here

                print("Data loading complete!")


        self._nside_levels = np.unique(self._pixel_info['nside'])
        self._hp_idx_sorted = [self._pixel_info['healpix_index']]
        self._data_idx = [np.arange(len(self._pixel_info['healpix_index']))]


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
                will return :obj:`ret`, :obj:'flags`, where :obj:`ret` is the normal return
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
                        * val[pix_idx[idx_near]][:, samp_idx, 0])
                else:
                    ret[idx_near] = (
                        a * val[pix_idx[idx_near]][:,samp_idx[idx_near], 0])

            # d > d(farthest distance slice)
            idx_far = (bin_idx_ceil == self._n_distances) & in_bounds_idx
            if np.any(idx_far):

                if isinstance(samp_idx, slice):
                    ret[idx_far] = val[pix_idx[idx_far]][:,samp_idx, -1]
                else:
                    ret[idx_far] = val[pix_idx[idx_far]][:,samp_idx[idx_far], -1]

            # d(nearest distance slice) < d < d(farthest distance slice)
            idx_btw = ~idx_near & ~idx_far & in_bounds_idx
            if np.any(idx_btw):
                DM_ceil = self._DM_bin_edges[bin_idx_ceil[idx_btw]]
                DM_floor = self._DM_bin_edges[bin_idx_ceil[idx_btw]-1]
                a = (DM_ceil - dm[idx_btw]) / (DM_ceil - DM_floor)
                if isinstance(samp_idx, slice):
                    ret[idx_btw] = (
                        (1.-a[:,None])
                        * val[pix_idx[idx_btw]][:,samp_idx, bin_idx_ceil[idx_btw]]
                        + a[:,None]
                        * val[pix_idx[idx_btw]][:,samp_idx, bin_idx_ceil[idx_btw]-1]
                    )
                else:
                    ret[idx_btw] = (
                        (1.-a) * val[pix_idx[idx_btw]][:,samp_idx[idx_btw], bin_idx_ceil[idx_btw]]
                        +    a * val[pix_idx[idx_btw]][:,samp_idx[idx_btw], bin_idx_ceil[idx_btw]-1]
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
            ret = val[pix_idx][:,samp_idx, :]   # Return all distances
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
