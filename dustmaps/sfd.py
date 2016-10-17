#!/usr/bin/env python
#
# sfd.py
# Reads the Schlegel, Finkbeiner & Davis (1998; SFD) dust reddening map.
#
# Copyright (C) 2016  Gregory M. Green, Edward F. Schlafly
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
import numpy as np

import astropy.wcs as wcs
import astropy.io.fits as fits
from scipy.ndimage import map_coordinates

from .std_paths import *
from .map_base import DustMap, ensure_flat_galactic
from . import fetch_utils
from . import dustexceptions


class SFDQuery(DustMap):
    """
    Queries the Schlegel, Finbeiner & Davis (1998) dust reddening map.
    """

    def __init__(self, map_dir=None):
        """
        Args:
            map_dir (Optional[str]): The directory containing the SFD map.
                Defaults to `None`, which means that `dustmaps` will look in its
                default data directory.
        """

        if map_dir is None:
            map_dir = os.path.join(data_dir(), 'sfd')

        self._data = {}

        base_fname = os.path.join(map_dir, 'SFD_dust_4096')

        for pole in ['ngp', 'sgp']:
            fname = '{}_{}.fits'.format(base_fname, pole)
            try:
                with fits.open(fname) as hdulist:
                    self._data[pole] = [hdulist[0].data, wcs.WCS(hdulist[0].header)]
            except IOError as error:
                print(dustexceptions.data_missing_message('sfd',
                                                          "SFD'98"))
                raise error

    @ensure_flat_galactic
    def query(self, coords, order=1):
        """
        Returns E(B-V) at the specified location(s) on the sky. See Table 6 of
        Schlafly & Finkbeiner (2011) for instructions on how to convert this
        quantity to extinction in various passbands.

        Args:
            coords (`astropy.coordinates.SkyCoord`): The coordinates to query.
            order (Optional[int]): Interpolation order to use. Defaults to `1`,
                for linear interpolation.

        Returns:
            A float array containing the SFD E(B-V) at every input coordinate.
            The shape of the output will be the same as the shape of the
            coordinates stored by `coords`.
        """
        out = np.zeros(len(coords.l.deg), dtype='f4')

        for pole in ['ngp', 'sgp']:
            m = (coords.b.deg >= 0) if pole == 'ngp' else (coords.b.deg < 0)

            if np.any(m):
                data, w = self._data[pole]
                x, y = w.wcs_world2pix(coords.l.deg[m], coords.b.deg[m], 0)
                out[m] = map_coordinates(data, [y, x], order=order, mode='nearest')

        return out


def fetch():
    """
    Downloads the Schlegel, Finkbeiner & Davis (1998) dust map, placing it in
    the data directory for `dustmap`.
    """
    doi = '10.7910/DVN/EWCNL5'

    for pole in ['ngp', 'sgp']:
        requirements = {'filename': 'SFD_dust_4096_{}.fits'.format(pole)}
        local_fname = os.path.join(
            data_dir(),
            'sfd', 'SFD_dust_4096_{}.fits'.format(pole))
        print('Downloading SFD data file to {}'.format(local_fname))
        fetch_utils.dataverse_download_doi(
            doi,
            local_fname,
            file_requirements=requirements)
