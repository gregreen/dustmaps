#!/usr/bin/env python
#
# chen2018.py
# Reads the 3D dust map of Chen et al. (2018).
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
import os

import astropy.coordinates as coordinates
import astropy.units as units
import astropy.io.fits as fits

from .equirectangular_map import EquirectangularDustMap
from .std_paths import *
from . import fetch_utils


class Chen2018Query(EquirectangularDustMap):
    """
    The 3D dust map of Chen et al. (2018), based on parallaxes from Gaia and
    photometry from Gaia and 2MASS imaging in the Galactic plane. The map
    covers |b| < 10 deg.
    """

    def __init__(self, map_fname=None, color='BR'):
        """
        Args:
            map_fname (Optional[:obj:`str`]): Filename at which the map is stored.
                Defaults to ``None``, meaning that the default filename is used.
        """
        if map_fname is None:
            map_fname = os.path.join(data_dir(), 'chen2018', 'chen2018.fits')

        with fits.open(map_fname) as f:
            d = f[1].data[:]
        
        lon0,lon1 = 0., 360.
        lat0,lat1 = -10., 10.
        dist0,dist1 = (0.2, 6.0) * units.kpc
        shape = (3600, 200, 30)

        pix_val = d['E'+color]
        pix_val.shape = shape

        super(Chen2018Query, self).__init__(
            pix_val,
            lon0, lon1,
            lat0, lat1,
            dist0=dist0, dist1=dist1,
            axis_order=('lon', 'lat', 'dist'),
            dist_interp='linear',
            frame='galactic'
        )


def fetch(clobber=False):
    """
    Downloads the 3D dust map of Chen et al. (2018).

    Args:
        clobber (Optional[bool]): If ``True``, any existing file will be
            overwritten, even if it appears to match. If ``False`` (the
            default), ``fetch()`` will attempt to determine if the dataset
            already exists. This determination is not 100\% robust against data
            corruption.
    """
    dest_dir = fname_pattern = os.path.join(data_dir(), 'chen2018')
    table_fname = os.path.join(dest_dir, 'chen2018.fits')
    
    # Check if the FITS table already exists
    table_md5sum = '043f2aa2064af607b56c040971f8d786'

    if (not clobber) and fetch_utils.check_md5sum(table_fname, table_md5sum):
        print('File appears to exist already. Call `fetch(clobber=True)` '
              'to force overwriting of existing file.')
        return

    # Download from the server
    url = 'http://paperdata.china-vo.org/diskec/cestar/table1.zip'
    archive_fname = os.path.join(dest_dir, 'table1.zip')
    archive_md5sum = '4acedf3f11ee8045e102a78b5b72036b'

    fetch_utils.download_and_verify(url, archive_md5sum, archive_fname)

    # Extract the FITS table
    print('Exracting FITS table from Zip archive ...')
    import zipfile

    readme_fname = os.path.join(dest_dir, 'readme.txt')

    with zipfile.ZipFile(archive_fname, 'r') as f:
        f.extract('table1.fits', path=dest_dir)
        f.extract('readme.txt', path=dest_dir)

    os.rename(os.path.join(dest_dir, 'table1.fits'), table_fname)

    # Delete the Zip archive
    print('Removing Zip archive ...')
    os.remove(archive_fname)
