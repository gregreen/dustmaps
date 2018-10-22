#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# lenz2017.py
# Reads the Lenz, Hensley & Doré (2017) dust reddening map.
# http://arxiv.org/abs/1706.00011
#
# Copyright (C) 2018  Gregory M. Green
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
from .healpix_map import HEALPixFITSQuery
from . import fetch_utils
from . import dustexceptions


class Lenz2017Query(HEALPixFITSQuery):
    """
    Queries the Lenz, Hensley & Doré (2017) dust map:
    http://arxiv.org/abs/1706.00011
    """

    def __init__(self, map_fname=None):
        """
        Args:
            map_fname (Optional[:obj:`str`]): Filename for the Lenz map. Defaults to
                ``None``, meaning that the default location is used.
        """

        if map_fname is None:
            map_fname = os.path.join(
                data_dir(),
                'lenz2017',
                'ebv_lhd.hpx.fits')

        try:
            super(Lenz2017Query, self).__init__(
                map_fname, 'galactic',
                hdu=1,
                field='EBV')
        except IOError as error:
            print(dustexceptions.data_missing_message('lenz2017',
                                                      'Lenz et al. (2017)'))
            raise error

    def query(self, coords, **kwargs):
        """
        Returns E(B-V), in mags, at the specified location(s) on the sky.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.

        Returns:
            A float array of the reddening, in magnitudes of E(B-V), at the
            selected coordinates.
        """
        return super(Lenz2017Query, self).query(coords, **kwargs)


def fetch():
    """
    Downloads the Lenz, Hensley & Doré (2017) dust map, placing it in the
    default :obj:`dustmaps` data directory.
    """
    doi = '10.7910/DVN/AFJNWJ'
    fname = os.path.join(
        data_dir(),
        'lenz2017',
        'ebv_lhd.hpx.fits')
    fetch_utils.dataverse_download_doi(
        doi, fname,
        file_requirements={'filename': 'ebv_lhd.hpx.fits'})


def main():
    from astropy.coordinates import SkyCoord
    q = Lenz2017Query()
    c = SkyCoord([0., 180., 0.], [0., 0., 90.], frame='galactic', unit='deg')
    print(q(c))


if __name__ == '__main__':
    main()
