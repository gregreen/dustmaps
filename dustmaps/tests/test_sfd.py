#!/usr/bin/env python
#
# test_sfd.py
# Test query code for Schlegel, Finkbeiner & Davis (1998) dust reddening map.
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

import unittest

import numpy as np
import astropy.coordinates as coords
import ujson as json
import os
import time

from .. import sfd
from ..std_paths import *

class TestSFD(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        t0 = time.time()
        with open(os.path.join(test_dir, 'ned_output.json'), 'r') as f:
            self._test_data = json.load(f)
            self._sfd = sfd.SFDQuery()
        t1 = time.time()
        print 'Loaded test data in {:.5f} s.'.format(t1-t0)

    def _get_equ(self, d):
        return coords.SkyCoord(d['equ'][0], d['equ'][1], frame='icrs')

    def _get_gal(self, d):
        return coords.SkyCoord(
            d['gal'][0], d['gal'][1],
            frame='galactic', unit='deg'
        )

    def test_sfd_scalar(self):
        for d in self._test_data:
            c = self._get_equ(d)
            Av = 2.742 * self._sfd(c)
            # print (d['sf11_Av'] - Av) / (0.001 + 0.001 * d['sf11_Av'])
            np.testing.assert_allclose(d['sf11_Av'], Av, atol=0.001, rtol=0.001)

    def test_sfd_vector(self):
        ra = [d['equ'][0] for d in self._test_data]
        dec = [d['equ'][1] for d in self._test_data]
        sf11_Av = np.array([d['sf11_Av'] for d in self._test_data])
        c = coords.SkyCoord(ra, dec, frame='icrs')

        Av = 2.742 * self._sfd(c)
        
        np.testing.assert_allclose(sf11_Av, Av, atol=0.001, rtol=0.001)

if __name__ == '__main__':
    unittest.main()
