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

import contextlib
import io
import os
import tempfile
import time
import warnings

import numpy as np
import astropy.coordinates as coords
import astropy.io.fits as fits
import astropy.wcs as wcs
from scipy.ndimage import map_coordinates

try:
    import ujson as json
except ImportError as error:
    import json

from .. import sfd
from ..std_paths import *

# The base names of the FITS files published for each SFD-style data product on
# Dataverse (doi:10.7910/DVN/EWCNL5). These names are fixed by the data release,
# so they are spelled out here rather than derived from `sfd.SFD_MAPS`.
EXPECTED_MAP_FILES = {
    ('SFD', 'dust'): 'SFD_dust_4096',
    ('SFD', 'i100'): 'SFD_i100_4096',
    ('SFD', 'i60'): 'SFD_i60_4096',
    ('SFD', 'mask'): 'SFD_mask_4096',
    ('SFD', 'temp'): 'SFD_temp',
    ('SFD', 'xmap'): 'SFD_xmap',
    ('FINK', 'Rmap'): 'FINK_Rmap',
    ('Haslam', 'clean'): 'Haslam_clean',
    ('Synch', 'Beta'): 'Synch_Beta',
}


class TestSFD(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        t0 = time.time()

        # Test data comes from NED
        with open(os.path.join(test_dir, 'ned_output.json'), 'r') as f:
            self._test_data = json.load(f)

        # Set up SFD query object
        self._sfd = sfd.SFDQuery()

        # The bit mask is not downloaded by the default `fetch()` call, so it
        # may well be absent.
        mask_files = [
            os.path.join(data_dir(), 'sfd',
                         'SFD_mask_4096_{}.fits'.format(pole))
            for pole in sfd.SFDBase.poles
        ]
        if all(os.path.exists(f) for f in mask_files):
            self._mask = sfd.SFDQuery(component='mask')
        else:
            self._mask = None

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

    def _gal_coords(self, n=256, seed=90210):
        """
        A reproducible set of Galactic coordinates, covering both hemispheres.

        The coordinates are given in the Galactic frame because
        ``ensure_flat_galactic`` then passes them to the query unchanged, which
        makes it possible to compare against a direct lookup in the FITS files.
        """
        rng = np.random.default_rng(seed)
        return coords.SkyCoord(
            rng.uniform(0., 360., n), rng.uniform(-90., 90., n),
            frame='galactic', unit='deg')

    def _require_mask(self):
        """
        Skip the calling test if the SFD bit mask is not available locally.
        """
        if self._mask is None:
            self.skipTest(
                'SFD bit mask is not in the data directory; run '
                "dustmaps.sfd.fetch(component='mask') to download it.")

    def test_map_variant_validation(self):
        """
        Test that an unrecognized ``map_variant`` or ``component`` (or an
        unrecognized combination of the two) raises a ValueError, and that the
        error message names the value that was not recognized.
        """
        for kwargs in (
                {'map_variant': 'nope'},
                {'component': 'nope'},
                {'map_variant': 'FINK', 'component': 'dust'},
                {'map_variant': 'SFD', 'component': 'Rmap'}):
            with self.assertRaises(ValueError) as raised:
                sfd.SFDQuery(**kwargs)
            for value in kwargs.values():
                self.assertIn(value, str(raised.exception))

        # fetch() should validate the arguments before attempting a download
        with self.assertRaises(ValueError):
            sfd.fetch(map_variant='nope')
        with self.assertRaises(ValueError):
            sfd.fetch(component='nope')

    def test_fetch_file_names(self):
        """
        Test that every data product in ``sfd.SFD_MAPS`` is fetched from the
        pair of FITS files (one per Galactic pole) that it is actually
        published as. The download itself is stubbed out, so no network access
        is needed.
        """
        # every product should have a hard-coded expectation, and vice versa
        products = set(
            (map_variant, component)
            for map_variant, components in sfd.SFD_MAPS.items()
            for component in components)
        self.assertEqual(products, set(EXPECTED_MAP_FILES))

        requested = []

        def stub_download(doi, fname, file_requirements=None):
            requested.append(file_requirements['filename'])

        orig_download = sfd.fetch_utils.dataverse_download_doi
        sfd.fetch_utils.dataverse_download_doi = stub_download
        try:
            for (map_variant, component), base_fname \
                    in EXPECTED_MAP_FILES.items():
                del requested[:]
                sfd.fetch(map_variant=map_variant, component=component)
                np.testing.assert_equal(
                    requested,
                    ['{}_{}.fits'.format(base_fname, pole)
                     for pole in sfd.SFDBase.poles])
        finally:
            sfd.fetch_utils.dataverse_download_doi = orig_download

    def test_default_map_variant(self):
        """
        Test that the defaults still select the original Schlegel, Finkbeiner &
        Davis (1998) dust reddening map; adding the other SFD-style data
        products should not have changed the default behavior.
        """
        q = sfd.SFDQuery(map_variant='SFD', component='dust')
        self.assertEqual(q.map_variant, 'SFD')
        self.assertEqual(q.component, 'dust')

        c = self._gal_coords()
        np.testing.assert_equal(self._sfd(c), q(c))

    def test_default_order(self):
        """
        Test that non-mask components default to linear interpolation
        (``order=1``), and that ``order`` is passed through to the
        interpolator.
        """
        c = self._gal_coords()

        e = self._sfd(c)
        self.assertTrue(np.issubdtype(e.dtype, np.floating))
        # an explicit order=1 is the default
        np.testing.assert_equal(e, self._sfd(c, order=1))
        # ... and order=0 (nearest neighbor) gives different values
        self.assertFalse(np.array_equal(e, self._sfd(c, order=0)))

    def test_scalar_coords(self):
        """
        Test that scalar input coordinates give scalar output, rather than a
        one-element array.
        """
        c = coords.SkyCoord(10., 20., frame='galactic', unit='deg')

        e = self._sfd(c)
        self.assertTrue(np.isscalar(e))

        if self._mask is not None:
            m = self._mask(c)
            self.assertTrue(np.isscalar(m))
            self.assertTrue(np.issubdtype(np.asarray(m).dtype, np.integer))

    def test_mask_returns_integers(self):
        """
        Test that the bit mask is returned as an integer array. The FITS files
        store the mask as 8-bit unsigned integers (``BITPIX=8``), but the
        interpolation used to sample it returns floats.
        """
        self._require_mask()

        mask = self._mask(self._gal_coords())

        self.assertTrue(np.issubdtype(mask.dtype, np.integer))
        self.assertEqual(mask.dtype, np.uint8)
        self.assertGreaterEqual(mask.min(), 0)
        self.assertLessEqual(mask.max(), 255)
        # the mask should not be constant over this coordinate set
        self.assertGreater(np.unique(mask).size, 1)

    def test_mask_matches_fits(self):
        """
        Test that the mask equals a direct nearest-neighbor (``order=0``)
        lookup in the FITS files, i.e. that no interpolation has taken place.
        """
        self._require_mask()

        c = self._gal_coords()
        mask = self._mask(c)

        expected = np.empty(len(c.l.deg), dtype=np.uint8)
        for pole in sfd.SFDBase.poles:
            m = (c.b.deg >= 0) if pole == 'ngp' else (c.b.deg < 0)
            if not np.any(m):
                continue
            fname = os.path.join(data_dir(), 'sfd',
                                 'SFD_mask_4096_{}.fits'.format(pole))
            with fits.open(fname) as hdulist:
                w = wcs.WCS(hdulist[0].header)
                x, y = w.wcs_world2pix(c.l.deg[m], c.b.deg[m], 0)
                expected[m] = map_coordinates(
                    hdulist[0].data, [y, x], order=0, mode='nearest')

        np.testing.assert_equal(mask, expected)

        # the number of HCONs is the two lowest bits (values 0-3)
        self.assertTrue(np.all((mask & 0b11) <= 3))

    def test_mask_order_override(self):
        """
        Test that the mask is always queried with nearest-neighbor sampling:
        interpolating a bit field is meaningless, so an explicit non-zero
        ``order`` is ignored (and a warning is issued).
        """
        self._require_mask()

        c = self._gal_coords()

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            default = self._mask(c)
        self.assertEqual(len(caught), 0)    # the default should not warn

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            forced = self._mask(c, order=1)
        self.assertEqual(len(caught), 1)
        self.assertIn('order=0', str(caught[0].message))

        np.testing.assert_equal(default, forced)
        np.testing.assert_equal(default, self._mask(c, order=0))

    def test_missing_data_message(self):
        """
        Test that the message printed when a map is missing names the requested
        data product, and gives a ``fetch()`` command that would actually
        download it.
        """
        with tempfile.TemporaryDirectory() as empty_dir:
            for kwargs, expected in (
                    ({}, 'dustmaps.sfd.fetch()'),
                    ({'component': 'mask'},
                     "dustmaps.sfd.fetch(map_variant='SFD', component='mask')"),
                    ({'map_variant': 'Synch', 'component': 'Beta'},
                     "dustmaps.sfd.fetch(map_variant='Synch', component='Beta')")):
                buf = io.StringIO()
                with contextlib.redirect_stdout(buf):
                    with self.assertRaises(IOError):
                        sfd.SFDQuery(map_dir=empty_dir, **kwargs)
                message = buf.getvalue()
                self.assertIn(expected, message)
                # the requested product should be named, not just "SFD'98"
                for value in kwargs.values():
                    self.assertIn(value, message)

if __name__ == '__main__':
    unittest.main()
