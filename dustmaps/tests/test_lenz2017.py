#!/usr/bin/env python
#
# test_lenz2017.py
# Test query code for the Lenz, Hensley & Doré (2017) dust map.
#
# Copyright (C) 2026  Gregory M. Green
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
import time

import numpy as np
import astropy.coordinates as coords
import astropy.units as units

from .. import lenz2017
from . import is_mmap


class TestLenz2017(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        print('Loading Lenz2017 dust map ...')
        t0 = time.time()
        self._lenz = lenz2017.Lenz2017Query()
        t1 = time.time()
        print('Loaded Lenz2017 test data in {:.5f} s.'.format(t1-t0))

    def test_shape(self):
        """
        Test that the output shapes are as expected with input coordinate arrays
        of different shapes.
        """
        for reps in range(5):
            n_dim = np.random.randint(1, 4)
            shape = np.random.randint(1, 7, size=(n_dim,))
            ra = np.random.uniform(-180., 180., size=shape) * units.deg
            dec = np.random.uniform(-90., 90., size=shape) * units.deg
            c = coords.SkyCoord(ra, dec, frame='icrs')
            E = self._lenz(c)
            np.testing.assert_equal(E.shape, shape)

    def test_frame(self):
        """
        Test that the results are independent of the coordinate frame.
        """
        frames = ('icrs', 'galactic', 'fk5', 'fk4', 'barycentrictrueecliptic')
        shape = (100,)
        ra = np.random.uniform(-180., 180., size=shape) * units.deg
        dec = np.random.uniform(-90., 90., size=shape) * units.deg
        c = coords.SkyCoord(ra, dec, frame='icrs')
        E0 = self._lenz(c)
        for fr in frames:
            cc = c.transform_to(fr)
            E = self._lenz(cc)
            np.testing.assert_equal(E, E0)

    def test_mmap(self):
        """
        Test that the underlying pixel data is memory-mapped (not copied into
        RAM) by default, and read into memory when ``memmap=False``.
        """
        self.assertTrue(is_mmap(self._lenz._pix_val))
        self.assertFalse(is_mmap(
            lenz2017.Lenz2017Query(memmap=False)._pix_val))


    def test_memmap_values(self):
        """
        Test that ``memmap=True`` and ``memmap=False`` give identical values
        at a set of pseudorandom, fixed-seed coordinates.

        The Lenz et al. (2017) map only covers high Galactic latitudes, so the
        candidates are drawn at ``|b| > 60`` deg and filtered to those where the
        map is defined (finite); closer to the plane the map returns NaN.
        """
        rng = np.random.default_rng(60217)
        l = rng.uniform(0., 360., 512)
        b = rng.uniform(60., 90., 512) * rng.choice([-1., 1.], 512)
        c_all = coords.SkyCoord(l*units.deg, b*units.deg, frame='galactic')

        values_all = self._lenz(c_all)
        c = c_all[np.isfinite(values_all)][:128]
        self.assertGreaterEqual(len(c), 64)

        values_mmap = self._lenz(c)
        values_ram = lenz2017.Lenz2017Query(memmap=False)(c)
        self.assertTrue(np.all(np.isfinite(values_mmap)))
        np.testing.assert_equal(values_mmap, values_ram)


if __name__ == '__main__':
    unittest.main()
