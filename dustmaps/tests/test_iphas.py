#!/usr/bin/env python
#
# test_iphas.py
# Test query code for the IPHAS dust extinction map of Sale et al. (2014).
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
import astropy.units as units
import os
import re
import time

from .. import iphas
from ..std_paths import *


class TestIPHAS(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        t0 = time.time()

        # Set up IPHAS query object
        self._iphas = iphas.IPHASQuery()

        t1 = time.time()
        print('Loaded IPHAS test data in {:.5f} s.'.format(t1-t0))

    def test_bounds(self):
        """
        Test that out-of-bounds coordinates return NaN reddening, and that
        in-bounds coordinates do not return NaN reddening.
        """

        for mode in (['random_sample', 'random_sample_per_pix',
                      'median', 'samples', 'mean']):
            # Draw random coordinates on the sphere
            n_pix = 10000
            u, v = np.random.random((2,n_pix))
            l = 360. * u
            b = 90. - np.degrees(np.arccos(2.*v - 1.))
            c = coords.SkyCoord(l, b, frame='galactic', unit='deg')

            A_calc = self._iphas(c, mode=mode)

            in_bounds = (l > 32.) & (l < 213.) & (b < 4.5) & (b > -4.5)
            out_of_bounds = (l < 28.) | (l > 217.) | (b > 7.) | (b < -7.)

            n_nan_in_bounds = np.sum(np.isnan(A_calc[in_bounds]))
            n_finite_out_of_bounds = np.sum(np.isfinite(A_calc[out_of_bounds]))

            self.assertTrue(n_nan_in_bounds == 0)
            self.assertTrue(n_finite_out_of_bounds == 0)

    def test_shape(self):
        """
        Test that the output shapes are as expected with input coordinate arrays
        of different shapes.
        """

        for mode in (['random_sample', 'random_sample_per_pix',
                      'median', 'mean', 'samples']):
            for include_dist in [False, True]:
                for reps in range(5):
                    # Draw random coordinates, with different shapes
                    n_dim = np.random.randint(1,4)
                    shape = np.random.randint(1,7, size=(n_dim,))

                    ra = (-180. + 360.*np.random.random(shape)) * units.deg
                    dec = (-90. + 180. * np.random.random(shape)) * units.deg
                    if include_dist:
                        dist = 5. * np.random.random(shape) * units.kpc
                    else:
                        dist = None
                    c = coords.SkyCoord(ra, dec, distance=dist, frame='icrs')

                    A_calc = self._iphas(c, mode=mode)

                    np.testing.assert_equal(A_calc.shape[:n_dim], shape)

                    extra_dims = 0
                    if mode == 'samples':
                        extra_dims += 1
                    if not include_dist:
                        extra_dims += 1

                    self.assertEqual(len(A_calc.shape), n_dim+extra_dims)



if __name__ == '__main__':
    unittest.main()
