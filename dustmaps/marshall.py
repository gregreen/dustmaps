#!/usr/bin/env python
#
# marshall.py
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
            map_fname (Optional[:obj:`str`]): Filename at which the map is stored.
                Defaults to ``None``, meaning that the default filename is used.
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
        self._n_dists = np.sum(np.isfinite(self._dist), axis=2)

        # idx = ~np.isfinite(self._dist)
        # if np.any(idx):
        #     self._dist[idx] = np.inf

        self._l_bounds = (-100., 100.) # min,max Galactic longitude, in deg
        self._b_bounds = (-10., 10.)    # min,max Galactic latitude, in deg
        self._inv_pix_scale = 4.        # 1 / (pixel scale, in deg)

    def _gal2idx(self, gal):
        """
        Converts from Galactic coordinates to pixel indices.

        Args:
            gal (:obj:`astropy.coordinates.SkyCoord`): Galactic coordinates. Must
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
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.
                Must contain distances.
            return_sigma (Optional[:obj:`bool`]): If ``True``, returns the uncertainty in
                extinction as well. Defaults to ``False``.

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
        dist_idx_ceil = np.sum(d[:,None] > self._dist[j, k, :], axis=1)

        # Initialize return arrays
        A_ret = np.full(coords.shape, np.nan, dtype='f4')
        if return_sigma:
            sigma_ret = np.full(coords.shape, np.nan, dtype='f4')

        # d < d(nearest distance slice)
        idx_near = (dist_idx_ceil == 0) & mask_idx

        if np.any(idx_near):
            a = d[idx_near] / self._dist[j[idx_near], k[idx_near], 0]
            A_ret[idx_near] = a * self._A[j[idx_near], k[idx_near], 0]

            if return_sigma:
                sigma_ret[idx_near] = (
                    self._sigma_A[
                        j[idx_near],
                        k[idx_near],
                        0
                    ])

        # d > d(farthest distance slice)
        idx_far = (dist_idx_ceil == self._n_dists[j,k]) & mask_idx

        if np.any(idx_far):
            A_ret[idx_far] = (
                self._A[
                    j[idx_far],
                    k[idx_far],
                    self._n_dists[j[idx_far],k[idx_far]]-1
                ])

            if return_sigma:
                sigma_ret[idx_far] = (
                    self._sigma_A[
                        j[idx_far],
                        k[idx_far],
                        self._n_dists[j[idx_far],k[idx_far]]-1
                    ])

        # d(nearest distance slice) < d < d(farthest distance slice)
        idx_btw = (~idx_near & ~idx_far) & mask_idx

        if np.any(idx_btw):
            d_ceil = (
                self._dist[
                    j[idx_btw],
                    k[idx_btw],
                    dist_idx_ceil[idx_btw]
                ])
            d_floor = (
                self._dist[
                    j[idx_btw],
                    k[idx_btw],
                    dist_idx_ceil[idx_btw]-1
                ])
            a = (d_ceil - d[idx_btw]) / (d_ceil - d_floor)

            A_ret[idx_btw] = (
                (1.-a) * self._A[j[idx_btw], k[idx_btw], dist_idx_ceil[idx_btw]]
                + a * self._A[j[idx_btw], k[idx_btw], dist_idx_ceil[idx_btw]-1]
            )

            if return_sigma:
                w0 = (1.-a)**2
                w1 = a**2
                norm = 1. / (w0 + w1)
                w0 *= norm
                w1 *= norm
                sigma_ret[idx_btw] = np.sqrt(
                    w0 * self._sigma_A[j[idx_btw], k[idx_btw], dist_idx_ceil[idx_btw]]**2
                    + w1 * self._sigma_A[j[idx_btw], k[idx_btw], dist_idx_ceil[idx_btw]-1]**2
                )

        if return_sigma:
            return A_ret, sigma_ret

        return A_ret


def dat2hdf5(table_dir):
    """
    Convert the Marshall et al. (2006) map from \*.dat.gz to \*.hdf5.
    """

    import astropy.io.ascii as ascii
    import gzip
    from contextlib import closing

    readme_fname = os.path.join(table_dir, 'ReadMe')
    table_fname = os.path.join(table_dir, 'table1.dat.gz')
    h5_fname = os.path.join(table_dir, 'marshall.h5')

    # Extract the gzipped table
    with gzip.open(table_fname, 'rb') as f:
        # Read in the table using astropy's CDS table reader
        r = ascii.get_reader(ascii.Cds, readme=readme_fname)
        r.data.table_name = 'table1.dat' # Hack to deal with bug in CDS reader.
        table = r.read(f)
        print(table)

    # Reorder table entries according to Galactic (l, b)
    l = coordinates.Longitude(
        table['GLON'][:],
        wrap_angle=180.*units.deg)
    b = table['GLAT'][:]

    sort_idx = np.lexsort((b, l))

    l = l[sort_idx].astype('f4')
    b = b[sort_idx].astype('f4')
    l.shape = (801, 81)
    b.shape = (801, 81)

    # Extract arrays from the table
    chi2_all = np.reshape((table['x2all'][sort_idx]).astype('f4'), (801,81))
    chi2_giants = np.reshape((table['x2gts'][sort_idx]).astype('f4'), (801,81))

    A = np.empty((801*81,33), dtype='f4')
    sigma_A = np.empty((801*81,33), dtype='f4')
    dist = np.empty((801*81,33), dtype='f4')
    sigma_dist = np.empty((801*81,33), dtype='f4')

    for k in range(33):
        A[:,k] = table['ext{:d}'.format(k+1)][sort_idx]
        sigma_A[:,k] = table['e_ext{:d}'.format(k+1)][sort_idx]
        dist[:,k] = table['r{:d}'.format(k+1)][sort_idx]
        sigma_dist[:,k] = table['e_r{:d}'.format(k+1)][sort_idx]

    A.shape = (801,81,33)
    sigma_A.shape = (801,81,33)
    dist.shape = (801,81,33)
    sigma_dist.shape = (801,81,33)

    # Construct the HDF5 file
    h5_fname = os.path.join(table_dir, 'marshall.h5')
    filter_kwargs = dict(
        chunks=True,
        compression='gzip',
        compression_opts=3,
        # scaleoffset=4
    )

    with h5py.File(h5_fname, 'w') as f:
        dset = f.create_dataset('A', data=A, **filter_kwargs)
        dset.attrs['description'] = 'Extinction of each bin'
        dset.attrs['band'] = 'Ks (2MASS)'
        dset.attrs['units'] = 'mag'

        dset = f.create_dataset('sigma_A', data=sigma_A, **filter_kwargs)
        dset.attrs['description'] = 'Extinction uncertainty of each bin'
        dset.attrs['band'] = 'Ks (2MASS)'
        dset.attrs['units'] = 'mag'

        dset = f.create_dataset('dist', data=dist, **filter_kwargs)
        dset.attrs['description'] = 'Distance of each bin'
        dset.attrs['units'] = 'kpc'

        dset = f.create_dataset('sigma_dist', data=sigma_dist, **filter_kwargs)
        dset.attrs['description'] = 'Distance uncertainty of each bin'
        dset.attrs['units'] = 'kpc'

        dset = f.create_dataset('chi2_all', data=chi2_all, **filter_kwargs)
        dset.attrs['description'] = 'Chi^2, based on all the stars'
        dset.attrs['units'] = 'unitless'

        dset = f.create_dataset('chi2_giants', data=chi2_giants, **filter_kwargs)
        dset.attrs['description'] = 'Chi^2, based on giants only'
        dset.attrs['units'] = 'unitless'

        # filter_kwargs.pop('scaleoffset')

        dset = f.create_dataset('l', data=l, **filter_kwargs)
        dset.attrs['description'] = 'Galactic longitude'
        dset.attrs['units'] = 'deg'

        dset = f.create_dataset('b', data=b, **filter_kwargs)
        dset.attrs['description'] = 'Galactic latitude'
        dset.attrs['units'] = 'deg'


def fetch(clobber=False):
    """
    Downloads the Marshall et al. (2006) dust map, which is based on 2MASS
    stellar photometry.

    Args:
        clobber (Optional[:obj:`bool`]): If ``True``, any existing file will be
            overwritten, even if it appears to match. If ``False`` (the
            default), :obj:`fetch()` will attempt to determine if the dataset
            already exists. This determination is not 100\% robust against data
            corruption.
    """

    table_dir = os.path.join(data_dir(), 'marshall')

    # Check if file already exists
    if not clobber:
        h5_fname = os.path.join(table_dir, 'marshall.h5')
        h5_size = 5033290 # Guess, in Bytes
        h5_dsets = {
            'l': (801, 81),
            'b': (801, 81),
            'chi2_all': (801, 81),
            'chi2_giants': (801, 81),
            'A': (801, 81, 33),
            'sigma_A': (801, 81, 33),
            'dist': (801, 81, 33),
            'sigma_dist': (801, 81, 33)
        }
        if fetch_utils.h5_file_exists(h5_fname, h5_size, dsets=h5_dsets):
            print('File appears to exist already. Call ``fetch(clobber=True)`` '
                  'to force overwriting of existing file.')
            return

    # Download the ASCII table
    url = 'ftp://cdsarc.u-strasbg.fr/pub/cats/J/A%2BA/453/635/table1.dat.gz'
    md5 = '637b95b025517a8b9757b6465b632285'
    table_fname = os.path.join(table_dir, 'table1.dat.gz')
    fetch_utils.download_and_verify(url, md5, fname=table_fname)

    # Download the README
    url = 'ftp://cdsarc.u-strasbg.fr/pub/cats/J/A%2BA/453/635/ReadMe'
    md5 = '3b7c1296b181b3d77106ab50193dc7ee'
    readme_fname = os.path.join(table_dir, 'ReadMe')
    fetch_utils.download_and_verify(url, md5, fname=readme_fname)

    # Convert from ASCII table to HDF5
    dat2hdf5(table_dir)

    # Cleanup
    print('Cleaning up ...')
    os.remove(table_fname)
    os.remove(readme_fname)
