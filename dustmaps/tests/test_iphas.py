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

    def _coords_in_map(self, n=500, seed=0):
        """Random coordinates inside the footprint of the map."""
        rng = np.random.default_rng(seed)
        l = rng.uniform(35., 110., n)
        b = rng.uniform(-4., 4., n)
        return coords.SkyCoord(l, b, unit='deg', frame='galactic')

    def test_the_same_seed_gives_the_same_random_samples(self):
        """
        The same seed draws the same random samples.
        """
        c = self._coords_in_map()
        for mode in ('random_sample', 'random_sample_per_pix'):
            first = self._iphas.query(c, mode=mode, seed_or_rng=2024)
            second = self._iphas.query(c, mode=mode, seed_or_rng=2024)
            np.testing.assert_equal(first, second)

    def test_different_seeds_give_different_random_samples(self):
        """
        Different seeds draw different samples. There are 500 coordinates, so
        agreement by chance is vanishingly unlikely.
        """
        c = self._coords_in_map()
        for mode in ('random_sample', 'random_sample_per_pix'):
            a = self._iphas.query(c, mode=mode, seed_or_rng=11)
            b = self._iphas.query(c, mode=mode, seed_or_rng=22)
            self.assertTrue(np.all(np.isfinite(a)))
            self.assertFalse(np.array_equal(a, b))

    def test_a_generator_is_drawn_from_as_it_stands(self):
        """
        A generator that is passed in is drawn from as it stands, so that
        successive calls give different samples, while two generators seeded
        alike agree with each other.
        """
        c = self._coords_in_map()

        rng = np.random.default_rng(7)
        first = self._iphas.query(c, mode='random_sample', seed_or_rng=rng)
        second = self._iphas.query(c, mode='random_sample', seed_or_rng=rng)
        self.assertFalse(np.array_equal(first, second))

        # A generator of the same seed starts from the same place
        again = self._iphas.query(c, mode='random_sample',
                                  seed_or_rng=np.random.default_rng(7))
        np.testing.assert_equal(first, again)

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
