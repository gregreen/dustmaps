#!/usr/bin/env python
#
# sfd.py
# Reads the Planck Collaboration dust reddening map.
#
# Copyright (C) 2016  Gregory M. Green
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

from .std_paths import *
from .healpix_map import HEALPixFITSQuery
from . import fetch_utils
from . import dustexceptions


class PlanckQuery(HEALPixFITSQuery):
    def __init__(self,
                 map_fname=os.path.join(
                    data_dir(),
                    'planck',
                    'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'),
                 component='extragalactic'):
        if component.lower() in ('ebv', 'extragalactic'):
            field = 'EBV'
            self._scale = 1.
        elif component.lower() in ('tau', 'tau353', 'tau_353', 'optical depth'):
            field = 'TAU353'
            self._scale = 1.49e4
        elif component.lower() in ('radiance', 'r'):
            field = 'RADIANCE'
            self._scale = 5.4e5

        try:
            with fits.open(map_fname) as hdulist:
                super(PlanckQuery, self).__init__(
                    hdulist, 'galactic',
                    hdu='COMP-MAP',
                    field=field)
        except IOError as error:
            print(dustexceptions.data_missing_message('planck',
                                                      'Planck Collaboration'))
            raise error

    def query(self, *args, **kwargs):
        return self._scale * super(PlanckQuery, self).query(*args, **kwargs)


def fetch():
    """
    Download the Planck dust map.
    """
    url = 'http://pla.esac.esa.int/pla/aio/product-action?MAP.MAP_ID=HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
    md5 = '8d804f4e64e709f476a63f0dfed1fd11'
    fname = os.path.join(
        data_dir(),
        'planck',
        'HFI_CompMap_ThermalDustModel_2048_R1.20.fits')
    fetch_utils.download_and_verify(url, md5, fname=fname)


def main():
    from astropy.coordinates import SkyCoord
    q = PlanckQuery()
    c = SkyCoord([0., 180., 0.], [0., 0., 90.], frame='galactic', unit='deg')
    print(q(c))


if __name__ == '__main__':
    main()
