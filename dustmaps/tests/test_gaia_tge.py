#!/usr/bin/env python
#
# test_gaia_tge.py
# Test query code for the Gaia Total Galactic Extinction (TGE) map.
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

import os
import shutil
import tempfile
import unittest

import numpy as np
import astropy.coordinates as coords
import astropy.units as units
import healpy as hp

from astropy.table import Table

from .. import gaia_tge
from ..std_paths import *


# A complete (if very coarse) HEALPix level: all 48 level-1 pixels, none of
# them flagged as optimum. This is needed to test loading a single, specified
# level, since that code path requires every pixel of the level to be present.
# (The uncertainty and min/max columns are just filler.)
_LEVEL1_ROWS = ''.join(
    '1,{},1,{:.2f},0.01,0.99,1.01,1,"False",0\n'.format(i, 1. + 0.01*i)
    for i in range(48))

# A miniature catalog, in the same format as the published CSV. On top of the
# complete level above, it has an optimum pixel at each of the levels 6-9, all
# covering the same place on the sky (NEST index 0, the pole), plus one pixel
# that is not optimum at any level. The A0 values increase with level, so that
# the resulting map reveals which level was used at each position.
SYNTHETIC_CSV = (
    '# %ECSV 1.0\n'
    '# ---\n'
    'solution_id,healpix_id,healpix_level,a0,a0_uncertainty,a0_min,a0_max,'
    'num_tracers_used,optimum_hpx_flag,status\n'
    '1,0,6,0.1,0.01,0.09,0.11,10,"True",0\n'
    '1,1,6,0.2,0.01,0.19,0.21,20,"True",0\n'
    '1,0,7,0.3,0.01,0.29,0.31,30,"True",0\n'
    '1,0,8,0.4,0.01,0.39,0.41,40,"True",0\n'
    '1,0,9,0.5,0.01,0.49,0.51,50,"True",0\n'
    '1,64,9,0.9,0.01,0.89,0.91,60,"False",0\n'
    + _LEVEL1_ROWS
)

COLUMNS = [
    'solution_id', 'healpix_id', 'healpix_level', 'a0', 'a0_uncertainty',
    'a0_min', 'a0_max', 'num_tracers_used', 'optimum_hpx_flag', 'status'
]


def nest_coords(nside, nest_idx):
    """
    An ICRS coordinate at the center of each of the given NEST-ordered HEALPix
    pixels. Note that the map is in NEST ordering, so this is the ordering that
    the query will see.
    """
    nest_idx = np.atleast_1d(nest_idx)
    theta, phi = hp.pix2ang(nside, nest_idx, nest=True)
    return coords.SkyCoord(phi, 0.5*np.pi - theta, frame='icrs', unit='rad')


class TestGaiaTGE(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self._tmp_dir = tempfile.mkdtemp()
        csv_fname = os.path.join(self._tmp_dir, 'catalog.csv')
        self._h5_fname = os.path.join(self._tmp_dir, 'catalog.h5')

        with open(csv_fname, 'w') as f:
            f.write(SYNTHETIC_CSV)

        gaia_tge.csv2h5(csv_fname, self._h5_fname)

        # The published map is not needed by most of these tests
        real_h5 = os.path.join(data_dir(), 'gaia_tge', 'gaia_tge.h5')
        self._real_h5 = real_h5 if os.path.exists(real_h5) else None

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self._tmp_dir)

    def test_csv2h5(self):
        """
        Test that ``csv2h5`` keeps every column of the catalog, with the
        extinction columns in magnitudes.
        """
        table = Table.read(self._h5_fname, path='data')

        self.assertEqual(table.colnames, COLUMNS)
        self.assertEqual(len(table), 54)
        self.assertEqual(table['a0'].unit, units.mag)
        self.assertEqual(table['a0_min'].unit, units.mag)

        np.testing.assert_allclose(
            table['a0'][:5], [0.1, 0.2, 0.3, 0.4, 0.5], rtol=1.e-6)
        np.testing.assert_array_equal(
            table['num_tracers_used'][:5], [10, 20, 30, 40, 50])
        np.testing.assert_array_equal(
            table['optimum_hpx_flag'][:5], [True]*5)
        self.assertFalse(bool(table['optimum_hpx_flag'][5]))

    def test_optimum(self):
        """
        Test loading the optimum HEALPix level at each position. This is the
        code path in which the int8 HEALPix level used to overflow, leaving an
        empty map (issue #70).
        """
        q = gaia_tge.GaiaTGEQuery(map_fname=self._h5_fname)

        self.assertEqual(q._nside, 512)
        self.assertEqual(len(q._pix_val), 12*4**9)

        # At the pole, the finest (level 9) pixel wins
        a0, flags = q(nest_coords(512, 0), return_flags=True)
        self.assertAlmostEqual(float(a0[0]), 0.5)
        self.assertEqual(int(flags['num_tracers_used'][0]), 50)

        # Level-6 pixels cover 64 level-9 pixels each; the second one is not
        # touched by any finer level, so it keeps its own value. (It also shows
        # that the row which is not flagged as optimum was ignored.)
        self.assertAlmostEqual(float(q(nest_coords(512, 64))[0]), 0.2)

        # Pixels that no level covers are left empty
        a0, flags = q(nest_coords(512, [128, 1000]), return_flags=True)
        self.assertTrue(np.all(np.isnan(a0)))
        self.assertTrue(np.all(flags['num_tracers_used'] == -1))
        self.assertFalse(np.any(flags['optimum_hpx_flag']))

    def test_healpix_level(self):
        """
        Test loading a single, specified HEALPix level.
        """
        q = gaia_tge.GaiaTGEQuery(map_fname=self._h5_fname, healpix_level=1)

        self.assertEqual(q._nside, 2)
        self.assertEqual(len(q._pix_val), 48)
        np.testing.assert_allclose(
            q(nest_coords(2, np.arange(48))), 1.0 + 0.01*np.arange(48))

    def test_invalid_healpix_level(self):
        """
        Test the error paths of the ``healpix_level`` argument.
        """
        # A level that is not stored in the map
        with self.assertRaises(ValueError):
            gaia_tge.GaiaTGEQuery(map_fname=self._h5_fname, healpix_level=5)

        # Neither an integer nor "optimum"
        with self.assertRaises(ValueError):
            gaia_tge.GaiaTGEQuery(map_fname=self._h5_fname,
                                  healpix_level='nope')

    def test_shape(self):
        """
        Test that the output shape matches the input shape.
        """
        q = gaia_tge.GaiaTGEQuery(map_fname=self._h5_fname)

        for _ in range(5):
            n_dim = np.random.randint(1, 4)
            shape = tuple(np.random.randint(1, 5, size=(n_dim,)))
            ra = -180. + 360. * np.random.random(shape)
            dec = -90. + 180. * np.random.random(shape)
            c = coords.SkyCoord(ra, dec, frame='icrs', unit='deg')
            self.assertEqual(q(c).shape, shape)

    def test_real_map(self):
        """
        Test the published map, if it has been downloaded.
        """
        if self._real_h5 is None:
            self.skipTest(
                'Gaia TGE map is not in the data directory; run '
                'dustmaps.gaia_tge.fetch() to download it.')

        q = gaia_tge.GaiaTGEQuery()

        self.assertEqual(q._nside, 512)
        self.assertEqual(len(q._pix_val), 12*4**9)

        # Most, but not all, of the sky has a value: the map has no estimate
        # where too few tracer stars were found (e.g. at low Galactic
        # latitude)
        n_finite = np.count_nonzero(np.isfinite(q._pix_val))
        self.assertGreater(n_finite, 0.9 * len(q._pix_val))
        self.assertLess(n_finite, len(q._pix_val))

        a0 = q(coords.SkyCoord(l=0., b=60., unit='deg', frame='galactic'))
        self.assertTrue(np.isfinite(a0))
        self.assertLess(float(a0), 1.)

        # The published catalog has complete levels 6-9, so a single level can
        # also be loaded
        q6 = gaia_tge.GaiaTGEQuery(healpix_level=6)
        self.assertEqual(q6._nside, 64)
        self.assertEqual(len(q6._pix_val), 12*4**6)


if __name__ == '__main__':
    unittest.main()
