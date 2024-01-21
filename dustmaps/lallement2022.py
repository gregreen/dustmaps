#!/usr/bin/env python
#
# leike2020.py
#
# The Lallement (2022) dust map.
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
from astropy.io import fits

from .map_base import DustMap, ensure_flat_coords
from .std_paths import *
from . import fetch_utils


class LallementQuery(DustMap):
    """
    A class for querying the Lallement (2022) dust map.

    For details on how to use this map, see the original paper:
    https://ui.adsabs.harvard.edu/abs/2022A%26A...661A.147L/abstract.

    The data is deposited at Vizier: https://cdsarc.cds.unistra.fr/ftp/J/A+A/661/A147/.
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
                'Lallement_2022',
                'cube_ext.fits.gz'
            )

        self._data = df = fits.getdata(map_fname).T
        header = fits.getheader(map_fname)
        sx, sy, sz = (header["SUN_POSX"], header["SUN_POSY"], header["SUN_POSZ"]) * units.pixel
        self.step = header["STEP"] * units.parsec / units.pixel
        # Lower edge of map, in pc
        xyz0 = (sx * -self.step, sy * -self.step, sz * -self.step)


        self._xyz0 = tuple((i.to(units.pc).value for i in xyz0)) # Lower edge of map, in pc
        self._shape = self._data.shape

    def _coords2idx(self, coords):
        c = coords.transform_to('galactic').represent_as('cartesian')

        idx = np.empty((3,) + c.shape, dtype='i4')
        mask = np.zeros(c.shape, dtype=np.bool)

        for i,x in enumerate((c.x, c.y, c.z)):
            idx[i,...] = np.floor((x.to('pc').value - self._xyz0[i]) / self.step.value)
            mask |= (idx[i] < 0) | (idx[i] >= self._shape[i])

        for i in range(3):
            idx[i, mask] = -1

        return idx, mask

    @ensure_flat_coords
    def query(self, coords):
        """
        Returns the extinction density (in e-foldings / kpc, in Gaia G-band)
        at the given coordinates.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): Coordinates at which
                to query the extinction. Must be 3D (i.e., include distance
                information).

        Returns:
            The extinction density, in units of e-foldings / pc, as either a
            numpy array or float, with the same shape as the input
            :obj:`coords`.
        """
        idx,mask = self._coords2idx(coords)

        v = self._data[idx[0], idx[1], idx[2]]

        if np.any(mask):
            # Set extinction to NaN for out-of-bounds (x, y, z)
            v[mask] = np.nan

        return v


def fetch(clobber=False):
    """
    Downloads the 3D dust map of Lallement (2022).

    Args:
        clobber (Optional[bool]): If ``True``, any existing file will be
            overwritten, even if it appears to match. If ``False`` (the
            default), ``fetch()`` will attempt to determine if the dataset
            already exists. This determination is not 100\% robust against data
            corruption.
    """
    dest_dir = fname_pattern = os.path.join(data_dir(), 'Lallement_2022')
    fname = os.path.join(dest_dir, 'cube_ext.fits.gz')

    # Check if the FITS table already exists
    md5sum = '5600c99aa869c6ad8971c5bb3b6aaea1'

    if (not clobber) and fetch_utils.check_md5sum(fname, md5sum):
        print('File appears to exist already. Call `fetch(clobber=True)` '
              'to force overwriting of existing file.')
        return

    # Download from the server
    url = 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/661/A147/cube_ext.fits.gz'
    fetch_utils.download_and_verify(url, md5sum, fname)


