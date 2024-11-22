#!/usr/bin/env python
#
# test_marshall.py
# Test query code for the dust extinction map of Marshall et al. (2006).
#
# Copyright (C) 2018  Gregory M. Green
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

from .. import marshall
from ..std_paths import *


class TestMarshall(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        t0 = time.time()

        # Set up IPHPAS query object
        self._marshall = marshall.MarshallQuery()

        t1 = time.time()
        print('Loaded Marshall+(2006) test data in {:.5f} s.'.format(t1-t0))

    def test_bounds(self):
        """
        Test that out-of-bounds coordinates return NaN reddening, and that
        in-bounds coordinates do not return NaN reddening.
        """

        for return_sigma in [False, True]:
            # Draw random coordinates on the sphere
            n_pix = 10000
            u, v = np.random.random((2,n_pix))
            l = 360. * u - 180.
            b = 90. - np.degrees(np.arccos(2.*v - 1.))
            d = 5. * np.random.random(l.shape)
            c = coords.SkyCoord(l*units.deg, b*units.deg,
                                distance=d*units.kpc, frame='galactic')

            res = self._marshall(c, return_sigma=return_sigma)

            if return_sigma:
                self.assertTrue(len(res) == 2)
                A, sigma = res
                np.testing.assert_equal(A.shape, sigma.shape)
            else:
                self.assertFalse(isinstance(res, tuple))
                A = res

            in_bounds = (l > -99.) & (l < 99.) & (b < 9.5) & (b > -9.5)
            out_of_bounds = (l < -101.) | (l > 101.) | (b > 10.5) | (b < -10.5)

            n_nan_in_bounds = np.sum(np.isnan(A[in_bounds]))
            n_finite_out_of_bounds = np.sum(np.isfinite(A[out_of_bounds]))

            self.assertTrue(n_nan_in_bounds == 0)
            self.assertTrue(n_finite_out_of_bounds == 0)

    def test_shape(self):
        """
        Test that the output shapes are as expected with input coordinate arrays
        of different shapes.
        """

        for return_sigma in [False, True]:
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

                    if not include_dist:
                        self.assertRaises(ValueError, self._marshall,
                                          c, return_sigma=return_sigma)
                        continue

                    res = self._marshall(c, return_sigma=return_sigma)

                    if return_sigma:
                        self.assertTrue(len(res) == 2)
                        A, sigma = res
                        np.testing.assert_equal(A.shape, sigma.shape)
                    else:
                        self.assertFalse(isinstance(res, tuple))
                        A = res

                    np.testing.assert_equal(A.shape, shape)



if __name__ == '__main__':
    unittest.main()
