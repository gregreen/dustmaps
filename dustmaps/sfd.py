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

import os
import h5py
import numpy as np

import astropy.coordinates as coordinates
import astropy.wcs as pywcs
import astropy.io.fits as fits
from scipy.ndimage import map_coordinates

from std_paths import *

class SFDQuery(object):
    def __init__(self, map_dir=os.path.join(data_dir, 'sfd')):
        self._data = {}

        base_fname = os.path.join(map_dir, 'SFD_dust_4096')

        for pole in ['ngp', 'sgp']:
            fname = '{}_{}.fits'.format(base_fname, pole)
            with fits.open(fname) as hdulist:
                self._data[pole] = [hdulist[0].data, pywcs.WCS(hdulist[0].header)]

    def query(self, coords, order=1):
        gal = coords.transform_to('galactic')

        is_array = hasattr(gal.l.deg, '__len__')

        if is_array:
            out = np.zeros(len(gal.l.deg), dtype='f4')

        for pole in ['ngp', 'sgp']:
            m = (gal.b.deg >= 0) if pole == 'ngp' else (gal.b.deg < 0)

            if np.any(m):
                data, wcs = self._data[pole]

                if not is_array: # Support for 0-dimensional arrays (scalars). Otherwise it barfs on l[m], b[m]
                    x, y = wcs.wcs_world2pix(gal.l.deg, gal.b.deg, 0)
                    out = map_coordinates(data, [[y], [x]], order=order, mode='nearest')[0]
                    continue

                x, y = wcs.wcs_world2pix(gal.l.deg[m], gal.b.deg[m], 0)
                out[m] = map_coordinates(data, [y, x], order=order, mode='nearest')

        return out

    def __call__(self, *args, **kwargs):
        return self.query(*args, **kwargs)
