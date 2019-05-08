#!/usr/bin/env python
#
# leike_ensslin_2019.py
# 
# The Leike & Ensslin (2019) dust map.
#
# Copyright (C) 2019  Gregory M. Green
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


class LeikeEnsslin2019Query(DustMap):
    """
    A class for querying the Leike & Ensslin (2019) dust map.
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
                'leike_ensslin_2019',
                'simple_cube.h5'
            )

        self._data = {}

        with h5py.File(map_fname) as f:
            self._data['mean'] = f['mean'][:]
            self._data['std'] = f['std'][:]

        self._shape = self._data['mean'].shape

    def _coords2idx(self, coords):
        c = coords.transform_to('galactic').represent_as('cartesian')
        
        idx = np.empty((3,) + c.shape, dtype='i4')
        mask = np.zeros(c.shape, dtype=np.bool)

        for i,x in enumerate((c.x, c.y, c.z)):
            idx[i,...] = np.floor(x.to('pc').value + 300) * 256/600.
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


def fetch(clobber=False):
    """
    Downloads the 3D dust map of Leike & Ensslin (2019).

    Args:
        clobber (Optional[bool]): If ``True``, any existing file will be
            overwritten, even if it appears to match. If ``False`` (the
            default), ``fetch()`` will attempt to determine if the dataset
            already exists. This determination is not 100\% robust against data
            corruption.
    """
    dest_dir = fname_pattern = os.path.join(data_dir(), 'leike_ensslin_2019')
    fname = os.path.join(dest_dir, 'simple_cube.h5')
    
    # Check if the FITS table already exists
    md5sum = 'f54e01c253453117e3770575bed35078'

    if (not clobber) and fetch_utils.check_md5sum(fname, md5sum):
        print('File appears to exist already. Call `fetch(clobber=True)` '
              'to force overwriting of existing file.')
        return

    # Download from the server
    url = 'https://zenodo.org/record/2577337/files/simple_cube.h5?download=1'
    fetch_utils.download_and_verify(url, md5sum, fname)
