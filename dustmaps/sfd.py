#!/usr/bin/env python
#
# sfd.py
# Reads the Schlegel, Finkbeiner & Davis (1998; SFD) dust reddening map.
#
# Copyright (C) 2016-2018  Gregory M. Green, Edward F. Schlafly
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
from .map_base import DustMap, WebDustMap, ensure_flat_galactic
from . import fetch_utils
from . import dustexceptions


class SFDBase(DustMap):
    """
    Queries maps stored in the same format as Schlegel, Finkbeiner & Davis (1998).
    """

    map_name = ''
    map_name_long = ''
    poles = ['ngp', 'sgp']

    def __init__(self, base_fname):
        """
        Args:
            base_fname (str): The map should be stored in two FITS files, named
                ``base_fname + '_' + X + '.fits'``, where ``X`` is ``'ngp'`` and
                ``'sgp'``.
        """
        self._data = {}

        for pole in self.poles:
            fname = '{}_{}.fits'.format(base_fname, pole)
            try:
                with fits.open(fname) as hdulist:
                    self._data[pole] = [hdulist[0].data, wcs.WCS(hdulist[0].header)]
            except IOError as error:
                print(dustexceptions.data_missing_message(self.map_name,
                                                          self.map_name_long))
                raise error

    @ensure_flat_galactic
    def query(self, coords, order=1):
        """
        Returns the map value at the specified location(s) on the sky.

        Args:
            coords (`astropy.coordinates.SkyCoord`): The coordinates to query.
            order (Optional[int]): Interpolation order to use. Defaults to `1`,
                for linear interpolation.

        Returns:
            A float array containing the map value at every input coordinate.
            The shape of the output will be the same as the shape of the
            coordinates stored by `coords`.
        """
        out = np.full(len(coords.l.deg), np.nan, dtype='f4')

        for pole in self.poles:
            m = (coords.b.deg >= 0) if pole == 'ngp' else (coords.b.deg < 0)

            if np.any(m):
                data, w = self._data[pole]
                x, y = w.wcs_world2pix(coords.l.deg[m], coords.b.deg[m], 0)
                out[m] = map_coordinates(data, [y, x], order=order, mode='nearest')

        return out


class SFDQuery(SFDBase):
    """
    Queries the Schlegel, Finkbeiner & Davis (1998) dust reddening map.
    """

    map_name = 'sfd'
    map_name_long = "SFD'98"

    def __init__(self, map_dir=None):
        """
        Args:
            map_dir (Optional[str]): The directory containing the SFD map.
                Defaults to `None`, which means that `dustmaps` will look in its
                default data directory.
        """
        
        if map_dir is None:
            map_dir = os.path.join(data_dir(), 'sfd')

        base_fname = os.path.join(map_dir, 'SFD_dust_4096')
        
        super(SFDQuery, self).__init__(base_fname)
    
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
        return super(SFDQuery, self).query(coords, order=order)
    

class SFDWebQuery(WebDustMap):
    """
    Remote query over the web for the Schlegel, Finkbeiner & Davis (1998) dust
    map.

    This query object does not require a local version of the data, but rather
    an internet connection to contact the web API. The query functions have the
    same inputs and outputs as their counterparts in ``SFDQuery``.
    """

    def __init__(self, api_url=None):
        super(SFDWebQuery, self).__init__(
            api_url=api_url,
            map_name='sfd')


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
