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
from scipy.spatial import cKDTree as KDTree
import h5py
import os.path

import astropy.coordinates as coordinates
import astropy.units as units

from contextlib import closing

from .std_paths import *
from .map_base import DustMap, ensure_flat_galactic


class IPHASQuery(DustMap):
    """
    The 3D dust maps of Sale et al. (2014), based on IPHAS imaging in the
    Galactic plane. The map covers 30 deg < l < 115 deg, -5 deg < b < 5 deg.
    """

    def __init__(self, map_fname=None):
        """
        Args:
            map_fname (Optional[str]): Filename at which the map is stored.
                Defaults to `None`, meaning that the default filename is used.
        """
        if map_fname is None:
            map_fname = os.path.join(data_dir(), 'iphas', 'iphas.h5')

        with closing(h5py.File(map_fname, 'r')) as f:
            self._data = f['samples'][:]

        self._n_pix = self._data.size
        self._n_dists = self._data['A0'].shape[1]
        self._n_samples = self._data['A0'].shape[2]

        # All the distance bins are the same
        self._dists = self._data['dist'][0]

        # Tesselate (l,b)-space
        self._lb = np.empty((self._data.size, 2), dtype='f8')
        self._lb[:,0] = self._data['l']
        self._lb[:,1] = self._data['b']
        self._kd = KDTree(self._lb)

        # Don't query more than this many degrees from any point
        self._border = 0.5

    def _gal2idx(self, gal):
        """
        Converts from Galactic coordinates to pixel indices.

        Args:
            gal (``astropy.coordinates.SkyCoord``): Galactic coordinates.

        Returns:
            Pixel indices of the coordinates, with the same shape as the input
            coordinates.
        """
        x = np.empty((gal.shape[0], 2), dtype='f8')
        x[:,0] = gal.l.deg
        x[:,1] = gal.b.deg
        idx = self._kd.query(x, p=1, distance_upper_bound=self._border)
        return idx[1]

    @ensure_flat_galactic
    def query(self, coords, mode='random_sample'):
        """
        Returns A0 at the given coordinates. There are several different query
        modes, which handle the probabilistic nature of the map differently.

        Args:
            coords (`astropy.coordinates.SkyCoord`): The coordinates to query.
            mode (Optional[str]): Four different query modes are available:
                'random_sample', 'samples', 'median' and 'mean'. The `mode`
                determines how the output will reflect the probabilistic nature
                of the Bayestar dust maps.

        Returns:
            Reddening at the specified coordinates, in mags of extinction, A0.
            The shape of the output depends on the `mode`, and on whether
            `coords` contains distances.

            If `coords` does not specify distance(s), then the shape of the
            output begins with `coords.shape`. If `coords` does specify
            distance(s), then the shape of the output begins with
            `coords.shape + ([number of distance bins],)`.

            If `mode` is 'random_sample', then at each coordinate/distance, a
            random sample of reddening is given.

            If `mode` is 'median', then at each coordinate/distance, the median
            reddening is returned.

            If `mode` is 'mean', then at each coordinate/distance, the mean
            reddening is returned.

            Finally, if `mode` is 'samples', then all at each
            coordinate/distance, all samples are returned.
        """

        # Check that the query mode is supported
        valid_modes = ['random_sample', 'samples', 'median', 'mean']
        if mode not in valid_modes:
            raise ValueError(
                '"{}" is not a valid `mode`. Valid modes are:\n'
                '  {}'.format(mode, valid_modes))

        n_coords_ret = coords.shape[0]

        # Determine if distance has been requested
        has_dist = hasattr(coords.distance, 'kpc')
        d = coords.distance.kpc if has_dist else None

        # Convert coordinates to pixel indices
        pix_idx = self._gal2idx(coords)

        # Determine which coordinates are out of bounds
        mask_idx = (pix_idx == self._n_pix)
        if np.any(mask_idx):
            pix_idx[mask_idx] = 0

        # Which samples to extract
        if mode == 'random_sample':
            samp_idx = np.random.randint(0, self._n_samples, pix_idx.size)
            n_samp_ret = 1
        else:
            samp_idx = slice(None)
            n_samp_ret = self._n_samples

        # Which distances to extract
        if has_dist:
            d = coords.distance.pc
            dist_idx_ceil = np.searchsorted(self._dists, d)

            if isinstance(samp_idx, slice):
                ret = np.empty((n_coords_ret, n_samp_ret), dtype='f8')
            else:
                ret = np.empty((n_coords_ret,), dtype='f8')

            # d < d(nearest distance slice)
            idx_near = (dist_idx_ceil == 0)
            if np.any(idx_near):
                a = d[idx_near] / self._dists[0]
                if isinstance(samp_idx, slice):
                    ret[idx_near] = a[:,None] * self._data['A0'][pix_idx[idx_near], 0, samp_idx]
                else:
                    ret[idx_near] = a[:] * self._data['A0'][pix_idx[idx_near], 0, samp_idx]

            # d > d(farthest distance slice)
            idx_far = (dist_idx_ceil == self._n_dists)
            if np.any(idx_far):
                ret[idx_far] = self._data['A0'][pix_idx[idx_far], -1, samp_idx]

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
                        (1.-a[:]) * self._data['A0'][pix_idx[idx_btw], dist_idx_ceil[idx_btw], samp_idx]
                        +    a[:] * self._data['A0'][pix_idx[idx_btw], dist_idx_ceil[idx_btw]-1, samp_idx])
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
        ``astropy.units.Quantity``, which stores unit-full quantities.
        """
        return self._dists * units.kpc






def ascii2h5(dirname, output_fname):
    """
    Converts from a directory of tarballed ASCII ".samp" files to a single
    HDF5 file. Essentially, converts from the original release format to a
    single HDF5 file.
    """

    import tarfile
    from glob import glob

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
                print('File {: >4d} of {:d}'.format(k+1, n_coords))
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

    tar_fname_list = glob(os.path.join(dirname, 'A_samp_*.tar.gz'))
    d = np.hstack([process_tarball(fn) for fn in tar_fname_list])
    save_data(d, output_fname)


def fetch():
    """
    Downloads the IPHAS 3D dust map of Sale et al. (2014).
    """
    raise NotImplementedError(
        'The automatic downloading of this map has not yet been implemented.')



def main():
    dirname = os.path.expanduser('~/Downloads/sale-iphas')
    output_fname = os.path.join(dirname, 'samp.h5')
    # ascii2h5(dirname, output_fname)

    iphas = IPHASQuery(map_fname=output_fname)

    from astropy.coordinates import SkyCoord
    c = SkyCoord(120.*units.deg, 30.*units.deg, distance=100.*units.pc, frame='galactic')

    A0 = iphas(c)
    print(A0)

    return 0


if __name__ == '__main__':
    main()
