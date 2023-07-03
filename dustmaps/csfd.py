#!/usr/bin/env python
#
# csfd.py
# Reads the "Corrected SFD" dust reddening map of Chiang (2023).
#
# Copyright (C) 2023  Gregory M. Green
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
import healpy as hp
import astropy.io.fits as fits
import astropy.units as units

from .std_paths import *
from .healpix_map import HEALPixQuery
from . import fetch_utils
from . import dustexceptions


class CSFDQuery(HEALPixQuery):
    """
    Queries the Corrected SFD dust map of Chiang (2023). This map is based
    on SFD, but contains a correction to remove contamination from
    large-scale structure (i.e., external galaxies).
    """

    def __init__(self, map_fname=None, mask_fname=None):
        """
        Args:
            map_fname (Optional[:obj:`str`]): Filename of the CSFD EBV map.
                Defaults to ```None``, meaning that the default location is
                used.
            mask_fname (Optional[:obj:`str`]): Filename of the CSFD mask map.
                Defaults to ```None``, meaning that the default location is
                used.
        """

        if map_fname is None:
            map_fname = os.path.join(data_dir(), 'csfd', 'csfd_ebv.fits')
        if mask_fname is None:
            mask_fname = os.path.join(data_dir(), 'csfd', 'mask.fits')

        try:
            with fits.open(map_fname) as hdulist:
                ebv_data = hdulist['xtension'].data[:]['T'].flatten()
            with fits.open(mask_fname) as hdulist:
                mask_data = hdulist['xtension'].data[:]['T'].flatten()
            mask_data = mask_data.astype('bool')
        except IOError as error:
            print(dustexceptions.data_missing_message('csfd',
                                                      'CSFD (Chiang 2023)'))
            raise error

        super(CSFDQuery, self).__init__(ebv_data, False, 'galactic',
                                        flags=mask_data)

    def query(self, coords, **kwargs):
        """
        Returns CSFD reddening on the same scale as SFD (similar to E(B-V)) at
        the specified location(s) on the sky.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to
                query.
            return_flags (Optional[:obj:`bool`]): If ``True``, then a
                boolean mask is returned as well, indicating where CSFD
                has been corrected for large-scale structure. Defaults to
                ``False``.

        Returns:
            A float array of the reddening, at the given coordinates. The
            shape of the output is the same as the shape of the input
            coordinate array, ``coords``. If ``return_flags`` is ``True``,
            a second array of the same shape (and boolean type) is returned,
            indicating whether the map has been corrected for large-scale
            structure at the given coordinates.
        """
        return super(CSFDQuery, self).query(coords, **kwargs)


def fetch():
    raise NotImplementedError()


