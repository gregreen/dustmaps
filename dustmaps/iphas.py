#!/usr/bin/env python
#
# iphas.py
# Reads the IPHAS 3D dust map of Sale et al. (2014).
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


class IPHASQuery(UnstructuredDustMap):
    """
    The 3D dust map of Sale et al. (2014), based on IPHAS imaging in the
    Galactic plane. The map covers 30 deg < l < 115 deg, -5 deg < b < 5 deg.
    """

    def __init__(self, map_fname=None):
        """
        Args:
            map_fname (Optional[:obj:`str`]): Filename at which the map is stored.
                Defaults to ``None``, meaning that the default filename is used.
        """
        if map_fname is None:
            map_fname = os.path.join(data_dir(), 'iphas', 'iphas.h5')

        with h5py.File(map_fname, 'r') as f:
            self._data = f['samples'][:]

        self._n_pix = self._data.size
        self._n_dists = self._data['A0'].shape[1]
        self._n_samples = self._data['A0'].shape[2]

        # All the distance bins are the same
        self._dists = self._data['dist'][0]

        # Don't query more than this angular distance from any point
        max_pix_scale = 0.5 * units.deg

        # Tesselate the sphere
        coords = coordinates.SkyCoord(
            self._data['l'],
            self._data['b'],
            unit='deg',
            frame='galactic')

        super(IPHASQuery, self).__init__(coords, max_pix_scale, metric_p=2)

    @ensure_flat_galactic
    def query(self, coords, mode='random_sample'):
        """
        Returns A0 at the given coordinates. There are several different query
        modes, which handle the probabilistic nature of the map differently.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.
            mode (Optional[:obj:`str`]): Five different query modes are available:
                ``'random_sample'``, ``'random_sample_per_pix'``, ``'samples'``,
                ``'median'`` and ``'mean'``. The ``mode`` determines how the output
                will reflect the probabilistic nature of the IPHAS dust map.

        Returns:
            Monochromatic extinction, A0, at the specified coordinates, in mags.
            The shape of the output depends on the ``mode``, and on whether
            ``coords`` contains distances.

            If ``coords`` does not specify distance(s), then the shape of the
            output begins with `coords.shape`. If `coords` does specify
            distance(s), then the shape of the output begins with
            ``coords.shape + ([number of distance bins],)``.

            If ``mode`` is ``'random_sample'``, then at each coordinate/distance, a
            random sample of reddening is given.

            If ``mode`` is ``'random_sample_per_pix'``, then the sample chosen for
            each angular pixel of the map will be consistent. For example, if
            two query coordinates lie in the same map pixel, then the same
            random sample will be chosen from the map for both query
            coordinates.

            If ``mode`` is ``'median'``, then at each coordinate/distance, the
            median reddening is returned.

            If ``mode`` is ``'mean'``, then at each coordinate/distance, the mean
            reddening is returned.

            Finally, if ``mode`` is ``'samples'``, then all at each
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
                '  {}'.format(mode, valid_modes))

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

        # Which samples to extract
        if mode == 'random_sample':
            samp_idx = np.random.randint(0, self._n_samples, pix_idx.size)
            n_samp_ret = 1
        elif mode == 'random_sample_per_pix':
            samp_idx = np.random.randint(0, self._n_samples, self._n_pix)[pix_idx]
            n_sample_ret = 1
        else:
            samp_idx = slice(None)
            n_samp_ret = self._n_samples

        # Which distances to extract
        if has_dist:
            d = coords.distance.pc
            dist_idx_ceil = np.searchsorted(self._dists, d)

            if isinstance(samp_idx, slice):
                ret = np.empty((n_coords_ret, n_samp_ret), dtype='f4')
            else:
                ret = np.empty((n_coords_ret,), dtype='f4')

            # d < d(nearest distance slice)
            idx_near = (dist_idx_ceil == 0)
            if np.any(idx_near):
                a = d[idx_near] / self._dists[0]
                if isinstance(samp_idx, slice):
                    ret[idx_near] = a[:,None] * self._data['A0'][pix_idx[idx_near], 0, samp_idx]
                else:
                    ret[idx_near] = a[:] * self._data['A0'][pix_idx[idx_near], 0, samp_idx[idx_near]]

            # d > d(farthest distance slice)
            idx_far = (dist_idx_ceil == self._n_dists)
            if np.any(idx_far):
                if isinstance(samp_idx, slice):
                    ret[idx_far] = self._data['A0'][pix_idx[idx_far], -1, samp_idx]
                else:
                    ret[idx_far] = self._data['A0'][pix_idx[idx_far], -1, samp_idx[idx_far]]

            # d(nearest distance slice) < d < d(farthest distance slice)
            idx_btw = ~idx_near & ~idx_far

            if np.any(idx_btw):
                d_ceil = self._dists[dist_idx_ceil[idx_btw]]
                d_floor = self._dists[dist_idx_ceil[idx_btw]-1]
                a = (d_ceil - d[idx_btw]) / (d_ceil - d_floor)
                if isinstance(samp_idx, slice):
                    ret[idx_btw] = (
                        (1.-a[:,None]) * self._data['A0'][pix_idx[idx_btw], dist_idx_ceil[idx_btw], samp_idx]
                        +    a[:,None] * self._data['A0'][pix_idx[idx_btw], dist_idx_ceil[idx_btw]-1, samp_idx])
                else:
                    ret[idx_btw] = (
                        (1.-a[:]) * self._data['A0'][pix_idx[idx_btw], dist_idx_ceil[idx_btw], samp_idx[idx_btw]]
                        +    a[:] * self._data['A0'][pix_idx[idx_btw], dist_idx_ceil[idx_btw]-1, samp_idx[idx_btw]])
        else:
            # TODO: Harmonize order of distances & samples with Bayestar.
            ret = self._data['A0'][pix_idx, :, samp_idx]

        # Reduce the samples in the requested manner
        samp_axis = 1 if has_dist else 2

        if mode == 'median':
            ret = np.median(ret, axis=samp_axis)
        elif mode == 'mean':
            ret = np.mean(ret, axis=samp_axis)

        if np.any(mask_idx):
            ret[mask_idx] = np.nan

        return ret

    @property
    def distances(self):
        """
        Returns the distance bins that the map uses. The return type is
        :obj:`astropy.units.Quantity`, which stores unit-full quantities.
        """
        return self._dists * units.kpc






def ascii2h5(dirname, output_fname):
    """
    Converts from a directory of tarballed ASCII ".samp" files to a single
    HDF5 file. Essentially, converts from the original release format to a
    single HDF5 file.
    """

    import tarfile
    import sys
    from glob import glob
    from contextlib import closing

    # The datatype that will be used to store extinction, A0
    A0_dtype = 'float16'

    def load_samp_file(f, fname):
        # Parse filename
        fname_chunks = os.path.split(fname)[1].split('_')
        l = float(fname_chunks[0])
        b = float(fname_chunks[1])

        # Load ASCII data
        data_raw = np.loadtxt(f, dtype='float64')

        n_samples = data_raw.shape[1] - 1
        n_dists = data_raw.shape[0]

        # Construct output
        dtype = [
            ('dist', 'int32'),
            ('A0', A0_dtype, (n_samples,))]

        data = np.empty(n_dists, dtype=dtype)
        data['dist'][:] = data_raw[:,0]
        data['A0'][:,:] = data_raw[:,1:]

        return (l,b), data

    def process_tarball(tarball_fname):
        # Write to the progress bar
        print('.', end='')
        sys.stdout.flush()

        with closing(tarfile.open(tarball_fname, mode='r:gz')) as f_tar:
            fnames = f_tar.getnames()

            f = f_tar.extractfile(fnames[0])
            (l,b), data = load_samp_file(f, fnames[0])
            n_dists, n_samples = data['A0'].shape
            n_coords = len(fnames)

            dtype = [
                ('l', 'float32'),
                ('b', 'float32'),
                ('dist', 'int32', (n_dists,)),
                ('A0', A0_dtype, (n_dists, n_samples))]
            data_combined = np.empty(n_coords, dtype=dtype)

            for k,fn in enumerate(fnames):
                # print('File {: >4d} of {:d}'.format(k+1, n_coords))
                f = f_tar.extractfile(fn)
                (l,b), data = load_samp_file(f, fn)
                data_combined['l'][k] = l
                data_combined['b'][k] = b
                data_combined['dist'][k] = data['dist']
                data_combined['A0'][k] = data['A0']

        return data_combined

    def save_data(data, fname):
        with closing(h5py.File(fname, 'w')) as f:
            f.create_dataset(
                'samples',
                data=data,
                chunks=True,
                compression='gzip',
                compression_opts=3)

    print('Progress: ', end='')
    sys.stdout.flush()

    tar_fname_list = glob(os.path.join(dirname, 'A_samp_*.tar.gz'))
    d = np.hstack([process_tarball(fn) for fn in tar_fname_list])

    print('+', end='')
    sys.stdout.flush()

    save_data(d, output_fname)

    print('')


def fetch(clobber=False):
    """
    Downloads the IPHAS 3D dust map of Sale et al. (2014).

    Args:
        clobber (Optional[bool]): If ``True``, any existing file will be
            overwritten, even if it appears to match. If ``False`` (the
            default), ``fetch()`` will attempt to determine if the dataset
            already exists. This determination is not 100\% robust against data
            corruption.
    """

    dest_dir = fname_pattern = os.path.join(data_dir(), 'iphas')
    url_pattern = 'http://www.iphas.org/data/extinction/A_samp_{:03d}.tar.gz'
    fname_pattern = os.path.join(dest_dir, 'A_samp_') + '{:03d}.tar.gz'

    # Check if file already exists
    if not clobber:
        h5_fname = os.path.join(dest_dir, 'iphas.h5')
        h5_size = 227817543 # Guess, in Bytes
        h5_dsets = {
            'samples': (61130,)
        }
        if fetch_utils.h5_file_exists(h5_fname, h5_size, dsets=h5_dsets):
            print('File appears to exist already. Call `fetch(clobber=True)` '
                  'to force overwriting of existing file.')
            return

    # Expected MD5 sums of .samp files
    file_md5sum = {
        30:  'dd531e397622bc97d4ff92b6c7863ade',
        40:  'b0f925eb3e46b77876e4054a26ad5b52',
        50:  'ea3b9500f0419d66dd92d9f9c127c2b5',
        60:  'cccf136f4e2306a6038e8093499216fd',
        70:  'a05fe2f815086686056c18087cc5410b',
        80:  '799bf618c8827b3d7250c884ec66ec49',
        90:  'd2a302d917da768bacf6ea74cb9dcfad',
        100: '2c75e31ad9320818556c4c9964b6af65',
        110: '742ea8de6f5f8a7e549f6c56b0088789',
        120: '9beabfa2c9634f953adadb5016eab072',
        130: '7cd7313f466eb60e8318d0f1bd32e035',
        140: 'fb6d09e4d939081b891e245c30b791f1',
        150: '8e9b6dc1561183aeadc64f41c85a64a8',
        160: '8a35828457b7b1d53d06998114553674',
        170: '7ffb29ec23e2f625dcfaaa84c293821d',
        180: 'c737da479d132b88483d6ddab5b25fc8',
        190: '9bc5fc7f7ba55f36a167473bb3679601',
        200: '7d8ffc4aa2f7c7026d8aa3ffb670d48e',
        210: 'e31b04964b7970b81fc90c120b4ebc24'
    }

    # Download the .samp files
    for key in file_md5sum:
        url = url_pattern.format(key)
        print('Downloading {}'.format(url))

        fetch_utils.download_and_verify(
            url,
            file_md5sum[key],
            fname_pattern.format(key))

    # Convert from ASCII to HDF5 format
    print('Repacking files...')
    ascii2h5(dest_dir, os.path.join(dest_dir, 'iphas.h5'))

    # Cleanup
    print('Removing original files...')
    for key in file_md5sum:
        os.remove(fname_pattern.format(key))
