#!/usr/bin/env python
#
# leike2020.py
# 
# The Leike, Glatzle & Ensslin (2020) dust map.
#
# Copyright (C) 2020  Gregory M. Green
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
import astropy.coordinates as coordinates
import astropy.units as units
from astropy.coordinates import Longitude

from .map_base import DustMap, ensure_flat_coords
from .std_paths import *
from . import fetch_utils


class Leike2020Query(DustMap):
    """
    A class for querying the Leike, Glatzle & Ensslin (2020) dust map.

    For details on how to use this map, see the original paper:
    https://ui.adsabs.harvard.edu/abs/2020A%26A...639A.138L/abstract.

    The data is deposited at Zenodo: https://doi.org/10.5281/zenodo.3993082.
    """

    def __init__(self, map_fname=None):
        """
        Args:
            map_fname (Optional[str]): Filename of the map. Defaults
                to :obj:`None`, meaning that the default location
                is used.
        """

        if map_fname is None:
            map_fname = os.path.join(
                data_dir(),
                'leike_2020',
                'mean_std.h5'
            )

        self._data = {}

        with h5py.File(map_fname) as f:
            self._data['mean'] = f['mean'][:]
            self._data['std'] = f['std'][:]

        self._xyz0 = (-370., -370., -270.) # Lower edge of map, in pc
        self._shape = self._data['mean'].shape

    def _coords2idx(self, coords):
        c = coords.transform_to('galactic').represent_as('cartesian')
        
        idx = np.empty((3,) + c.shape, dtype='i4')
        mask = np.zeros(c.shape, dtype=np.bool)

        for i,x in enumerate((c.x, c.y, c.z)):
            idx[i,...] = np.floor(x.to('pc').value - self._xyz0[i])
            mask |= (idx[i] < 0) | (idx[i] >= self._shape[i])

        for i in range(3):
            idx[i, mask] = -1

        return idx, mask

    @ensure_flat_coords
    def query(self, coords, component='mean'):
        """
        Returns the extinction density (in e-foldings / kpc, in Gaia G-band)
        at the given coordinates.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): Coordinates at which
                to query the extinction. Must be 3D (i.e., include distance
                information).
            component (str): Which component to return. Allowable values are
                'mean' (for the mean extinction density) and 'std' (for the
                standard deviation of extinction density). Defaults to 'mean'.

        Returns:
            The extinction density, in units of e-foldings / pc, as either a
            numpy array or float, with the same shape as the input
            :obj:`coords`.
        """
        idx,mask = self._coords2idx(coords)

        v = self._data[component][idx[0], idx[1], idx[2]]
        
        if np.any(mask):
            # Set extinction to NaN for out-of-bounds (x, y, z)
            v[mask] = np.nan

        return v


def fetch(clobber=False, fetch_samples=False):
    """
    Downloads the 3D dust map of Leike & Ensslin (2020).

    Args:
        clobber (Optional[bool]): If ``True``, any existing file will be
            overwritten, even if it appears to match. If ``False`` (the
            default), ``fetch()`` will attempt to determine if the dataset
            already exists. This determination is not 100\% robust against data
            corruption.
        fetch_samples (Optional[bool]): If ``True``, the samples will also be
            downloaded. If ``False`` (the default), only the mean and standard
            deviation will be downloaded. The samples take up 14 GB, which is
            why the default is not to download them.
    """
    dest_dir = fname_pattern = os.path.join(data_dir(), 'leike_2020')

    file_spec = [
        ('mean_std.h5', '1ea998fdaef58f53da639356362223ba')
    ]
    if fetch_samples:
        file_spec += [
            ('samples.h5', '581f9ebc4775d37fd431fc6c0984dcf6')
        ]

    for fn,md5sum in file_spec:
        fname = os.path.join(dest_dir, fn)

        # Check if the file already exists
        if (not clobber) and fetch_utils.check_md5sum(fname, md5sum):
            print('File "{}" appears to exist already. Call '.format(fn)+
                  '`fetch(clobber=True)` to force overwriting of existing '+
                  'file.')
            continue
        
        # Download from the server
        url = 'https://zenodo.org/record/3993082/files/{}?download=1'.format(fn)
        fetch_utils.download_and_verify(url, md5sum, fname)

