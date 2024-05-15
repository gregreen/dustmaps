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
        except IOError as error:
            print(dustexceptions.data_missing_message('csfd',
                                                      'CSFD (Chiang 2023)'))
            raise error

        super(CSFDQuery, self).__init__(ebv_data, False, 'galactic',
                                        flags=mask_data)

    def query(self, coords, **kwargs):
        """
        Returns CSFD reddening on the same scale as SFD (similar to E(B-V)) at
        the specified location(s) on the sky. Also optionally returns a
        bit mask, where the bits (ordered from least to most significant) have
        the following meanings::
        
            Bit 0: 'LSS_corr' - This bit is set in the footprint within which
                   the LSS is reconstructed, and CSFD = SFD - LSS (otherwise
                   CSFD = SFD).
            Bit 1: 'no_IRAS' - Set in the area with no IRAS data (DIRBE data
                   filled in SFD); LSS removal in CSFD is done using a 1 deg
                   smoothed LSS.
            Bit 2: 'cosmology' - Set in the area where both the LSS and CSFD
                   are most reliable for precision cosmology analyses.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to
                query.
            return_flags (Optional[:obj:`bool`]): If ``True``, then a
                bit mask is returned as well, indicating where CSFD
                has been corrected for large-scale structure, where IRAS data
                was used, and where the map is suitable for cosmology. See
                above description of bits. Defaults to ``False``.

        Returns:
            A float array of the reddening, at the given coordinates. The
            shape of the output is the same as the shape of the input
            coordinate array, ``coords``. If ``return_flags`` is ``True``,
            a second array (a bit mask) of the same shape is returned. See
            above description of the meaning of each bit.
        """
        return super(CSFDQuery, self).query(coords, **kwargs)


def fetch(clobber=False):
    """
    Downloads the Corrected SFD dust map of Chiang (2023).

    Args:
        clobber (Optional[bool]): If ``True``, any existing file will be
            overwritten, even if it appears to match. If ``False`` (the
            default), ``fetch()`` will attempt to determine if the dataset
            already exists. This determination is not 100\% robust against data
            corruption.
    """
    file_spec = [
        ('csfd_ebv.fits', '31cd2eec51bcb5f106af84a610ced53c'),
        ('mask.fits', '9142f5a5d184125836a68b6f48d1113f')
    ]
    for fn,md5sum in file_spec:
        fname = os.path.join(data_dir(), 'csfd', fn)
        # Download from Zenodo
        url = 'https://zenodo.org/record/8207175/files/{}'.format(fn)
        fetch_utils.download_and_verify(url, md5sum, fname, clobber=clobber)

