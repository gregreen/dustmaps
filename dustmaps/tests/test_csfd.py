#!/usr/bin/env python
#
# test_csfd.py
# Test query code for the Corrected SFD (CSFD) dust map of Chiang (2023).
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

from .. import csfd
from . import is_mmap


class TestCSFD(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        print('Loading CSFD dust map ...')
        t0 = time.time()
        self._csfd = csfd.CSFDQuery()
        t1 = time.time()
        print('Loaded CSFD test data in {:.5f} s.'.format(t1-t0))

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
            E = self._csfd(c)
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
        E0 = self._csfd(c)
        for fr in frames:
            cc = c.transform_to(fr)
            E = self._csfd(cc)
            np.testing.assert_equal(E, E0)

    def test_mmap(self):
        """
        Test that the underlying pixel data is memory-mapped (not copied into
        RAM) by default, and read into memory when ``memmap=False``.
        """
        self.assertTrue(is_mmap(self._csfd._pix_val))
        self.assertTrue(is_mmap(self._csfd._flags))
        q = csfd.CSFDQuery(memmap=False)
        self.assertFalse(is_mmap(q._pix_val))
        self.assertFalse(is_mmap(q._flags))


    def test_memmap_values(self):
        """
        Test that ``memmap=True`` and ``memmap=False`` give identical values
        (both the reddening map and the mask/flags) at a set of pseudorandom,
        fixed-seed coordinates.
        """
        rng = np.random.default_rng(35)
        n = 128
        c = coords.SkyCoord(rng.uniform(-180., 180., n)*units.deg,
                            rng.uniform(-90., 90., n)*units.deg,
                            frame='icrs')
        q_ram = csfd.CSFDQuery(memmap=False)
        np.testing.assert_equal(self._csfd(c), q_ram(c))
        values_mmap, flags_mmap = self._csfd(c, return_flags=True)
        values_ram, flags_ram = q_ram(c, return_flags=True)
        np.testing.assert_equal(values_mmap, values_ram)
        np.testing.assert_equal(flags_mmap, flags_ram)


if __name__ == '__main__':
    unittest.main()
