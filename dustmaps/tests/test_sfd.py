#!/usr/bin/env python
#
# test_sfd.py
# Test query code for Schlegel, Finkbeiner & Davis (1998) dust reddening map.
#
# Copyright (C) 2016  Gregory M. Green
#
# dustmaps is free software: you can redistribute it and/or modify
# it under the terms of either:
#
# - The GNU General Public License as published by the Free Software Foundation,
#   either version 2 of the License, or (at your option) any later version, or
# - The 2-Clause BSD License (also known as the Simplified BSD License).
#
# You should have received copies of the GNU General Public License
# and the BSD License along with this program.
#

from __future__ import print_function, division

import unittest

import numpy as np
import astropy.coordinates as coords

try:
    import ujson as json
except ImportError as error:
    import json

import os
import time

from .. import sfd
from ..std_paths import *

class TestSFD(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        t0 = time.time()

        # Test data comes from NED
        with open(os.path.join(test_dir, 'ned_output.json'), 'r') as f:
            self._test_data = json.load(f)

        # Set up SFD query object
        self._sfd = sfd.SFDQuery()

        t1 = time.time()
        print('Loaded SFD test data in {:.5f} s.'.format(t1-t0))

    def _get_equ(self, d):
        """
        Get Equatorial (ICRS) coordinates of test data point.
        """
        return coords.SkyCoord(d['equ'][0], d['equ'][1], frame='icrs')

    def _get_gal(self, d):
        """
        Get Galactic coordinates of test data point.
        """
        return coords.SkyCoord(
            d['gal'][0], d['gal'][1],
            frame='galactic', unit='deg'
        )

    def test_sfd_equ_scalar(self):
        """
        Test SFD query of individual ICRS coordinates.
        """
        # print 'Equatorial'
        # print '=========='

        for d in self._test_data:
            c = self._get_equ(d)
            Av = 2.742 * self._sfd(c)
            # c_gal = c.transform_to('galactic')
            # print '* (l, b) = ({: >16.8f} {: >16.8f})'.format(c_gal.l.deg, c_gal.b.deg)
            # print (d['sf11_Av'] - Av) / (0.001 + 0.001 * d['sf11_Av'])
            np.testing.assert_allclose(d['sf11_Av'], Av, atol=0.001, rtol=0.001)

    def test_sfd_gal_scalar(self):
        """
        Test SFD query of individual Galactic coordinates.
        """
        # print 'Galactic'
        # print '========'

        for d in self._test_data:
            c = self._get_gal(d)
            Av = 2.742 * self._sfd(c)
            # print '* (l, b) = ({: >16.8f} {: >16.8f})'.format(c.l.deg, c.b.deg)
            # print (d['sf11_Av'] - Av) / (0.001 + 0.001 * d['sf11_Av'])
            np.testing.assert_allclose(d['sf11_Av'], Av, atol=0.001, rtol=0.001)

    def test_sfd_equ_vector(self):
        """
        Test SFD query of multiple ICRS coordinates at once.
        """
        ra = [d['equ'][0] for d in self._test_data]
        dec = [d['equ'][1] for d in self._test_data]
        sf11_Av = np.array([d['sf11_Av'] for d in self._test_data])
        c = coords.SkyCoord(ra, dec, frame='icrs')

        Av = 2.742 * self._sfd(c)

        np.testing.assert_allclose(sf11_Av, Av, atol=0.001, rtol=0.001)

    def test_sfd_gal_vector(self):
        """
        Test SFD query of multiple Galactic coordinates at once.
        """
        l = [d['gal'][0] for d in self._test_data]
        b = [d['gal'][1] for d in self._test_data]
        sf11_Av = np.array([d['sf11_Av'] for d in self._test_data])
        c = coords.SkyCoord(l, b, frame='galactic', unit='deg')

        Av = 2.742 * self._sfd(c)

        np.testing.assert_allclose(sf11_Av, Av, atol=0.001, rtol=0.001)

    def test_shape(self):
        """
        Test that the output shapes are as expected with input coordinate arrays
        of different shapes.
        """

        for reps in range(10):
            # Draw random coordinates, with different shapes
            n_dim = np.random.randint(1,4)
            shape = np.random.randint(1,7, size=(n_dim,))

            ra = -180. + 360.*np.random.random(shape)
            dec = -90. + 180. * np.random.random(shape)
            c = coords.SkyCoord(ra, dec, frame='icrs', unit='deg')

            ebv_calc = self._sfd(c)

            np.testing.assert_equal(ebv_calc.shape, shape)

    def test_malformed_coords(self):
        """
        Test that SFD query errors with malformed input.
        """
        c = np.array([
            [d['equ'][0] for d in self._test_data],
            [d['equ'][1] for d in self._test_data]
        ])

        with self.assertRaises(TypeError):
            self._sfd(c)

if __name__ == '__main__':
    unittest.main()
