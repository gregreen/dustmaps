#!/usr/bin/env python
#
# chen2014.py
# Reads the Chen et al. (2014) dust map, based on stellar photometry from the
# Xuyi Schmidt Telescope Photometric Survey of the Galactic Anticentre.
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
import h5py
import os

import astropy.coordinates as coordinates
import astropy.units as units

from .std_paths import *
from .map_base import DustMap, ensure_flat_galactic
from .unstructured_map import UnstructuredDustMap
from . import fetch_utils


class Chen2014Query(UnstructuredDustMap):
    """
    The 3D dust map of Chen et al. (2014), based on stellar photometry from the
    Xuyi Schmidt Telescope Photometric Survey of the Galactic Anticentre. The
    map covers 140 deg < l < 240 deg, -60 deg < b < 40 deg.
    """

    def __init__(self, map_fname=None):
        """
        Args:
            map_fname (Optional[:obj:`str`]): Filename at which the map is stored.
                Defaults to ``None``, meaning that the default filename is used.
        """
        if map_fname is None:
            map_fname = os.path.join(data_dir(), 'chen2014', 'chen2014.h5')

        with h5py.File(map_fname, 'r') as f:
            self._dists = f['dists'][:]
            self._lb = f['pix_lb'][:]
            self._A = f['A_r'][:]
            self._sigma_A = f['A_r_err'][:]

        # Have to filter out zero pixels
        # idx = ~np.all(self._A < 1.e-5, axis=1)
        # self._lb = self._lb[idx]
        # self._A = self._A[idx]
        # self._sigma_A = self._sigma_A[idx]

        self._n_dists = self._dists.size

        # Don't query more than this angular distance from any point
        max_pix_scale = 0.5 * units.deg

        # Tesselate the sphere
        coords = coordinates.SkyCoord(
            self._lb[:,0],
            self._lb[:,1],
            unit='deg',
            frame='galactic')

        super(Chen2014Query, self).__init__(coords, max_pix_scale, metric_p=2)

    @ensure_flat_galactic
    def query(self, coords, return_sigma=False):
        """
        Returns r-band extinction, A_r, at the given coordinates. Can also
        return uncertainties.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.
            return_sigma (Optional[:obj:`bool`]): If ``True``, returns the uncertainty in
                extinction as well. Defaults to ``False``.

        Returns:
            Extinction in the r-band at the specified coordinates, in mags.
            The shape of the output depends on whether :obj:`coords` contains
            distances.

            If :obj:`coords` does not specify distance(s), then the shape of the
            output begins with :obj:`coords.shape`. If :obj:`coords` does specify
            distance(s), then the shape of the output begins with
            ``coords.shape + ([number of distance bins],)``.
        """
        n_coords_ret = coords.shape[0]

        # Determine if distance has been requested
        has_dist = hasattr(coords.distance, 'kpc')
        d = coords.distance.kpc if has_dist else None

        # Convert coordinates to pixel indices
        pix_idx = self._coords2idx(coords)

        # Determine which coordinates are out of bounds
        mask_idx = (pix_idx == self._n_pix)
        if np.any(mask_idx):
            pix_idx[mask_idx] = 0

        # Which distances to extract
        if has_dist:
            d = coords.distance.kpc
            dist_idx_ceil = np.searchsorted(self._dists, d)

            ret = np.empty((n_coords_ret,), dtype='f8')
            if return_sigma:
                sigma_ret = np.empty((n_coords_ret,), dtype='f8')

            # d < d(nearest distance slice)
            idx_near = (dist_idx_ceil == 0) & ~mask_idx
            print('d < d(nearest): {:d}'.format(np.sum(idx_near)))
            if np.any(idx_near):
                a = d[idx_near] / self._dists[0]
                ret[idx_near] = a[:] * self._A[pix_idx[idx_near], 0]
                if return_sigma:
                    sigma_ret[idx_near] = a[:] * self._sigma_A[pix_idx[idx_near], 0]

            # d > d(farthest distance slice)
            idx_far = (dist_idx_ceil == self._n_dists) & ~mask_idx
            print('d > d(farthest): {:d}'.format(np.sum(idx_far)))
            if np.any(idx_far):
                ret[idx_far] = self._A[pix_idx[idx_far], -1]
                if return_sigma:
                    sigma_ret[idx_far] = self._sigma_A[pix_idx[idx_far], -1]

            # d(nearest distance slice) < d < d(farthest distance slice)
            idx_btw = ~idx_near & ~idx_far & ~mask_idx
            print('d(nearest) < d < d(farthest): {:d}'.format(np.sum(idx_btw)))
            if np.any(idx_btw):
                d_ceil = self._dists[dist_idx_ceil[idx_btw]]
                d_floor = self._dists[dist_idx_ceil[idx_btw]-1]
                a = (d_ceil - d[idx_btw]) / (d_ceil - d_floor)
                ret[idx_btw] = (
                    (1.-a[:]) * self._A[pix_idx[idx_btw], dist_idx_ceil[idx_btw]]
                    +    a[:] * self._A[pix_idx[idx_btw], dist_idx_ceil[idx_btw]-1])
                if return_sigma:
                    w0 = (1.-a)**2
                    w1 = a**2
                    norm = 1. / (w0 + w1)
                    w0 *= norm
                    w1 *= norm
                    sigma_ret[idx_btw] = np.sqrt(
                        w0 * self._sigma_A[pix_idx[idx_btw], dist_idx_ceil[idx_btw]]**2
                        + w1 * self._sigma_A[pix_idx[idx_btw], dist_idx_ceil[idx_btw]-1]**2
                    )
        else:
            # TODO: Harmonize order of distances & samples with Bayestar.
            ret = self._A[pix_idx, :]
            if return_sigma:
                sigma_ret = self._sigma_A[pix_idx, :]

        if np.any(mask_idx):
            ret[mask_idx] = np.nan
            if return_sigma:
                sigma_ret[mask_idx] = np.nan

        if return_sigma:
            return ret, sigma_ret

        return ret

    @property
    def distances(self):
        """
        Returns the distance bins that the map uses. The return type is
        :obj:`astropy.units.Quantity`, which stores unit-full quantities.
        """
        return self._dists * units.kpc


def ascii2h5(dat_fname, h5_fname):
    """
    Converts from the original ASCII format of the Chen+ (2014) 3D dust map to
    the HDF5 format.

    Args:
        dat_fname (:obj:`str`): Filename of the original ASCII .dat file.
        h5_fname (:obj:`str`): Output filename to write the resulting HDF5 file to.
    """
    table = np.loadtxt(dat_fname, skiprows=1, dtype='f4')

    filter_kwargs = dict(
        chunks=True,
        compression='gzip',
        compression_opts=3)

    # Filter out pixels with all zeros
    idx = ~np.all(table[:,2:32] < 1.e-5, axis=1)

    with h5py.File(h5_fname, 'w') as f:
        d = np.arange(0., 4.351, 0.15).astype('f4')

        dset = f.create_dataset('dists', data=d, **filter_kwargs)
        dset.attrs['description'] = 'Distances at which extinction is measured'
        dset.attrs['units'] = 'kpc'

        dset = f.create_dataset('pix_lb', data=table[idx,0:2], **filter_kwargs)
        dset.attrs['description'] = 'Galactic (l, b) of each pixel'
        dset.attrs['units'] = 'deg'

        dset = f.create_dataset('A_r', data=table[idx,2:32], **filter_kwargs)
        dset.attrs['description'] = 'Extinction'
        dset.attrs['shape'] = '(pixel, distance)'
        dset.attrs['band'] = 'r'
        dset.attrs['units'] = 'mag'

        dset = f.create_dataset('A_r_err', data=table[idx,32:], **filter_kwargs)
        dset.attrs['description'] = 'Gaussian uncertainty in extinction'
        dset.attrs['shape'] = '(pixel, distance)'
        dset.attrs['band'] = 'r'
        dset.attrs['units'] = 'mag'


def fetch(clobber=False):
    """
    Downloads the Chen et al. (2014) dust map.

    Args:
        clobber (Optional[:obj:`bool`]): If ``True``, any existing file will be
            overwritten, even if it appears to match. If ``False`` (the
            default), :obj:`fetch()` will attempt to determine if the dataset
            already exists. This determination is not 100\% robust against data
            corruption.
    """

    dest_dir = fname_pattern = os.path.join(data_dir(), 'chen2014')
    url = 'http://lamost973.pku.edu.cn/site/Photometric-Extinctions-and-Distances/table2.dat'
    dat_fname = os.path.join(dest_dir, 'chen2014.dat')
    h5_fname = os.path.join(dest_dir, 'chen2014.h5')
    md5 = 'f8a2bc46d411c57ca4c76dc344e291f1'

    # Check if file already exists
    if not clobber:
        h5_size = 52768768 # Guess, in Bytes
        h5_dsets = {
            'dists': (30,),
            'pix_lb': (557398, 2),
            'A_r': (557398, 30),
            'A_r_err': (557398, 30)
        }
        if fetch_utils.h5_file_exists(h5_fname, h5_size, dsets=h5_dsets):
            print('File appears to exist already. Call `fetch(clobber=True)` '
                  'to force overwriting of existing file.')
            return

    # Download the table
    print('Downloading {}'.format(url))
    fetch_utils.download_and_verify(url, md5, fname=dat_fname)

    # Convert from ASCII to HDF5 format
    print('Repacking files...')
    ascii2h5(dat_fname, h5_fname)

    # Cleanup
    print('Removing original file...')
    os.remove(dat_fname)
