#!/usr/bin/env python
#
# test_bayestar.py
# Test query code for the Green, Schlafly, Finkbeiner et al. (2015) dust map.
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
from unittest import mock

import numpy as np
import astropy.coordinates as coords
import astropy.units as units
import h5py
import os
import re
import shutil
import tempfile
import time

from .. import bayestar
from ..std_paths import *

def parse_argonaut_output(fname, max_samples=None):
    with open(fname, 'r') as f:
        txt = f.read()

    output = []

    meta_fields = ({
        'l':
            (re.compile(r'(?:l = )([-]?[0-9]*[.]?[0-9]*)'), float),
        'b':
            (re.compile(r'(?:b = )([-]?[0-9]*[.]?[0-9]*)'), float),
        'ra':
            (re.compile(r'(?:ra = )([-]?[0-9]*[.]?[0-9]*)'), float),
        'dec':
            (re.compile(r'(?:dec = )([-]?[0-9]*[.]?[0-9]*)'), float),
        'DM_min':
            (re.compile(r'(?:min: )([-]?[0-9]*[.]?[0-9]*)'), float),
        'DM_max':
            (re.compile(r'(?:max: )([-]?[0-9]*[.]?[0-9]*)'), float),
        'converged':
            (re.compile(r'(?:converged: )(True|False)'), lambda x: x == 'True'),
        'n_stars':
            (re.compile(r'(?:stars: )([-]?[0-9]*[.]?[0-9]*)'), int)
    })

    p_dm = re.compile(r'(?:DistanceModulus[\s]*\|)(.*)')
    p_best = re.compile(r'(?:BestFit[\s]*\|)(.*)')
    p_sample = re.compile(r'(?:[0-9]+[\s]*\|)(.*)')

    for block in txt.split('# Line-of-Sight Reddening Results'):
        if not len(block.rstrip()):
            continue

        row = {}

        # Metadata
        for key in meta_fields.keys():
            p, f = meta_fields[key]
            row[key] = f(p.search(block).group(1))

        # Distance slices
        row['DM_bin_edges'] = np.array([
            float(s) for s in p_dm.search(block).group(1).split()
        ])

        # Best fit
        row['best'] = np.array([
            float(s) for s in p_best.search(block).group(1).split()
        ])

        # Samples
        row['samples'] = np.array([
            [float(s) for s in m.group(1).split()]
            for m in re.finditer(p_sample, block)
        ])

        if max_samples != None:
            row['samples'] = row['samples'][:max_samples,:]

        output.append(row)

    return output

class TestBayestar(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        print('Loading Bayestar query object ...')
        t0 = time.time()

        max_samples = 4

        # Test data comes from argonaut.skymaps.info
        fname = os.path.join(test_dir, 'argonaut_output_v1.txt')
        self._test_data = parse_argonaut_output(fname, max_samples=max_samples)

        # Set up Bayestar query object
        self._bayestar = bayestar.BayestarQuery(version='bayestar2015',
                                                max_samples=max_samples)

        t1 = time.time()
        print('Loaded Bayestar test data in {:.5f} s.'.format(t1-t0))

    def _get_equ(self, d, dist=None):
        """
        Get Equatorial (ICRS) coordinates of test data point.
        """
        return coords.SkyCoord(
            d['ra']*units.deg,
            d['dec']*units.deg,
            distance=dist,
            frame='icrs'
        )

    def _get_gal(self, d, dist=None):
        """
        Get Galactic coordinates of test data point.
        """
        return coords.SkyCoord(
            d['l']*units.deg,
            d['b']*units.deg,
            distance=dist,
            frame='galactic'
        )

    def atest_plot_samples(self):
        dm = np.linspace(4., 19., 1001)
        samples = []

        for dm_k in dm:
            d = 10.**(dm_k/5.-2.)
            samples.append(self._interp_ebv(self._test_data[0], d))

        samples = np.array(samples).T
        # print samples

        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        for s in samples:
            ax.plot(dm, s, lw=2., alpha=0.5)

        plt.show()

    def test_equ_med_far_scalar(self):
        """
        Test that median reddening is correct in the far limit, using a single
        location on the sky at a time as input.
        """
        for d in self._test_data:
            c = self._get_gal(d, dist=1.e3*units.kpc)
            # print d['samples']
            ebv_data = np.nanmedian(d['samples'][:,-1])
            ebv_calc = self._bayestar(c, mode='median')
            # print 'ebv_data:', ebv_data
            # print 'ebv_calc:', ebv_calc
            # print ''
            # print r'% residual: {:.6f}'.format((ebv_calc - ebv_data) / (0.001 + 0.001 * ebv_data))
            np.testing.assert_allclose(ebv_data, ebv_calc, atol=0.001, rtol=0.0001)

    def test_equ_med_far_vector(self):
        """
        Test that median reddening is correct in the far limit, using a vector
        of coordinates as input.
        """
        l = [d['l']*units.deg for d in self._test_data]
        b = [d['b']*units.deg for d in self._test_data]
        dist = [1.e3*units.kpc for bb in b]
        c = coords.SkyCoord(l, b, distance=dist, frame='galactic')

        ebv_data = np.array([np.nanmedian(d['samples'][:,-1]) for d in self._test_data])
        ebv_calc = self._bayestar(c, mode='median')

        # print 'vector:'
        # print r'% residual:'
        # for ed,ec in zip(ebv_data, ebv_calc):
        #     print '  {: >8.3f}'.format((ec - ed) / (0.02 + 0.02 * ed))

        np.testing.assert_allclose(ebv_data, ebv_calc, atol=0.001, rtol=0.0001)

    def test_equ_med_vector(self):
        """
        Test that median reddening is correct at arbitary distances, using a
        vector of coordinates as input.
        """
        for reps in range(10):
            l = [d['l']*units.deg for d in self._test_data]
            b = [d['b']*units.deg for d in self._test_data]
            dm = 3. + (25.-3.)*np.random.random(len(self._test_data))
            dist = [d*units.kpc for d in 10.**(dm/5.-2.)]
            dist_unitless = [d for d in 10.**(dm/5.-2.)]
            c = coords.SkyCoord(l, b, distance=dist, frame='galactic')

            ebv_samples = np.array([
                self._interp_ebv(datum, d)
                for datum,d in zip(self._test_data, dist_unitless)
            ])
            ebv_data = np.nanmedian(ebv_samples, axis=1)
            ebv_calc = self._bayestar(c, mode='median')

            # print 'vector arbitrary distance:'
            # print r'% residual:'
            # for ed,ec in zip(ebv_data, ebv_calc):
            #     print '  {: >8.3f}'.format((ec - ed) / (0.02 + 0.02 * ed))

            np.testing.assert_allclose(ebv_data, ebv_calc, atol=0.001, rtol=0.0001)

    def test_equ_med_scalar(self):
        """
        Test that median reddening is correct in at arbitary distances, using
        individual coordinates as input.
        """
        for d in self._test_data:
            l = d['l']*units.deg
            b = d['b']*units.deg

            for reps in range(10):
                dm = 3. + (25.-3.)*np.random.random()
                dist = 10.**(dm/5.-2.)
                c = coords.SkyCoord(l, b, distance=dist*units.kpc, frame='galactic')

                ebv_samples = self._interp_ebv(d, dist)
                ebv_data = np.nanmedian(ebv_samples)
                ebv_calc = self._bayestar(c, mode='median')

                np.testing.assert_allclose(ebv_data, ebv_calc, atol=0.001, rtol=0.0001)

    def test_equ_random_sample_vector(self):
        """
        Test that random sample of reddening at arbitary distance is actually
        from the set of possible reddening samples at that distance. Uses vector
        of coordinates/distances as input.
        """

        # Prepare coordinates (with random distances)
        l = [d['l']*units.deg for d in self._test_data]
        b = [d['b']*units.deg for d in self._test_data]
        dm = 3. + (25.-3.)*np.random.random(len(self._test_data))

        dist = [d*units.kpc for d in 10.**(dm/5.-2.)]
        dist_unitless = [d for d in 10.**(dm/5.-2.)]
        c = coords.SkyCoord(l, b, distance=dist, frame='galactic')

        ebv_data = np.array([
            self._interp_ebv(datum, d)
            for datum,d in zip(self._test_data, dist_unitless)
        ])
        ebv_calc = self._bayestar(c, mode='random_sample')

        d_ebv = np.min(np.abs(ebv_data[:,:] - ebv_calc[:,None]), axis=1)

        # print 'vector arbitrary distance random sample:'
        # print r'% residual:'
        # for ec,de in zip(ebv_calc, d_ebv):
        #     print '  {: >8.3f}  {: >8.5f}'.format(ec, de)

        np.testing.assert_allclose(d_ebv, 0., atol=0.001, rtol=0.0001)

    def test_equ_samples_vector(self):
        """
        Test that full set of samples of reddening at arbitary distance is
        correct. Uses vector of coordinates/distances as input.
        """

        # Prepare coordinates (with random distances)
        l = [d['l']*units.deg for d in self._test_data]
        b = [d['b']*units.deg for d in self._test_data]
        dm = 3. + (25.-3.)*np.random.random(len(self._test_data))

        dist = [d*units.kpc for d in 10.**(dm/5.-2.)]
        dist_unitless = [d for d in 10.**(dm/5.-2.)]
        c = coords.SkyCoord(l, b, distance=dist, frame='galactic')

        ebv_data = np.array([
            self._interp_ebv(datum, d)
            for datum,d in zip(self._test_data, dist_unitless)
        ])
        ebv_calc = self._bayestar(c, mode='samples')

        # print 'vector arbitrary distance random sample:'
        # print r'% residual:'
        # for ed,ec in zip(ebv_data, ebv_calc):
        #     resid = (ed-ec) / (0.001 + 0.0001*ed)
        #     print '  {: >8.3f}  {: >8.3f}  {: >8.5f}  {: >8.5f}'.format(
        #         ed[0], ed[1],
        #         resid[0], resid[1]
        #     )

        np.testing.assert_allclose(ebv_data, ebv_calc, atol=0.001, rtol=0.0001)

    def test_equ_samples_scalar(self):
        """
        Test that full set of samples of reddening at arbitary distance is
        correct. Uses single set of coordinates/distance as input.
        """
        for d in self._test_data:
            # Prepare coordinates (with random distances)
            l = d['l']*units.deg
            b = d['b']*units.deg
            dm = 3. + (25.-3.)*np.random.random()

            dist = 10.**(dm/5.-2.)
            c = coords.SkyCoord(l, b, distance=dist*units.kpc, frame='galactic')

            ebv_data = self._interp_ebv(d, dist)
            ebv_calc = self._bayestar(c, mode='samples')

            np.testing.assert_allclose(ebv_data, ebv_calc, atol=0.001, rtol=0.0001)

    def test_equ_random_sample_scalar(self):
        """
        Test that random sample of reddening at arbitary distance is actually
        from the set of possible reddening samples at that distance. Uses vector
        of coordinates/distances as input. Uses single set of
        coordinates/distance as input.
        """
        for d in self._test_data:
            # Prepare coordinates (with random distances)
            l = d['l']*units.deg
            b = d['b']*units.deg
            dm = 3. + (25.-3.)*np.random.random()

            dist = 10.**(dm/5.-2.)
            c = coords.SkyCoord(l, b, distance=dist*units.kpc, frame='galactic')

            ebv_data = self._interp_ebv(d, dist)
            ebv_calc = self._bayestar(c, mode='random_sample')

            d_ebv = np.min(np.abs(ebv_data[:] - ebv_calc))

            np.testing.assert_allclose(d_ebv, 0., atol=0.001, rtol=0.0001)

    def test_equ_samples_nodist_vector(self):
        """
        Test that full set of samples of reddening vs. distance curves is
        correct. Uses vector of coordinates as input.
        """

        # Prepare coordinates
        l = [d['l']*units.deg for d in self._test_data]
        b = [d['b']*units.deg for d in self._test_data]

        c = coords.SkyCoord(l, b, frame='galactic')

        ebv_data = np.array([d['samples'] for d in self._test_data])
        ebv_calc = self._bayestar(c, mode='samples')

        # print 'vector random sample:'
        # print 'ebv_data.shape = {}'.format(ebv_data.shape)
        # print 'ebv_calc.shape = {}'.format(ebv_calc.shape)
        # print ebv_data[0]
        # print ebv_calc[0]

        np.testing.assert_allclose(ebv_data, ebv_calc, atol=0.001, rtol=0.0001)

    def test_equ_random_sample_nodist_vector(self):
        """
        Test that a random sample of the reddening vs. distance curve is drawn
        from the full set of samples. Uses vector of coordinates as input.
        """

        # Prepare coordinates
        l = [d['l']*units.deg for d in self._test_data]
        b = [d['b']*units.deg for d in self._test_data]

        c = coords.SkyCoord(l, b, frame='galactic')

        ebv_data = np.array([d['samples'] for d in self._test_data])
        ebv_calc = self._bayestar(c, mode='random_sample')

        # print 'vector random sample:'
        # print 'ebv_data.shape = {}'.format(ebv_data.shape)
        # print 'ebv_calc.shape = {}'.format(ebv_calc.shape)
        # print ebv_data[0]
        # print ebv_calc[0]

        d_ebv = np.min(np.abs(ebv_data[:,:,:] - ebv_calc[:,None,:]), axis=1)
        np.testing.assert_allclose(d_ebv, 0., atol=0.001, rtol=0.0001)

    def test_bounds(self):
        """
        Test that out-of-bounds coordinates return NaN reddening, and that
        in-bounds coordinates do not return NaN reddening.
        """

        for mode in (['random_sample', 'random_sample_per_pix',
                      'median', 'samples', 'mean']):
            # Draw random coordinates, both above and below dec = -30 degree line
            n_pix = 1000
            ra = -180. + 360.*np.random.random(n_pix)
            dec = -75. + 90.*np.random.random(n_pix)    # 45 degrees above/below
            c = coords.SkyCoord(ra, dec, frame='icrs', unit='deg')

            ebv_calc = self._bayestar(c, mode=mode)

            nan_below = np.isnan(ebv_calc[dec < -35.])
            nan_above = np.isnan(ebv_calc[dec > -25.])
            pct_nan_above = np.sum(nan_above) / float(nan_above.size)

            # print r'{:s}: {:.5f}% nan above dec=-25 deg.'.format(mode, 100.*pct_nan_above)

            self.assertTrue(np.all(nan_below))
            self.assertTrue(pct_nan_above < 0.05)

    def test_shape(self):
        """
        Test that the output shapes are as expected with input coordinate arrays
        of different shapes.
        """

        for mode in ['random_sample', 'median', 'mean', 'samples']:
            for reps in range(5):
                # Draw random coordinates, with different shapes
                n_dim = np.random.randint(1,4)
                shape = np.random.randint(1,7, size=(n_dim,))

                ra = -180. + 360.*np.random.random(shape)
                dec = -90. + 180. * np.random.random(shape)
                c = coords.SkyCoord(ra, dec, frame='icrs', unit='deg')

                ebv_calc = self._bayestar(c, mode=mode)

                np.testing.assert_equal(ebv_calc.shape[:n_dim], shape)

                if mode == 'samples':
                    self.assertEqual(len(ebv_calc.shape), n_dim+2) # sample, distance
                else:
                    self.assertEqual(len(ebv_calc.shape), n_dim+1) # distance

    def _interp_ebv(self, datum, dist):
        """
        Calculate samples of E(B-V) at an arbitrary distance (in kpc) for one
        test coordinate.
        """
        dm = 5. * (np.log10(dist) + 2.)
        idx_ceil = np.searchsorted(datum['DM_bin_edges'], dm)
        if idx_ceil == 0:
            dist_0 = 10.**(datum['DM_bin_edges'][0]/5. - 2.)
            return dist/dist_0 * datum['samples'][:,0]
        elif idx_ceil == len(datum['DM_bin_edges']):
            return datum['samples'][:,-1]
        else:
            dm_ceil = datum['DM_bin_edges'][idx_ceil]
            dm_floor = datum['DM_bin_edges'][idx_ceil-1]
            a = (dm_ceil - dm) / (dm_ceil - dm_floor)
            return (
                (1.-a) * datum['samples'][:,idx_ceil]
                +    a * datum['samples'][:,idx_ceil-1]
            )



# Format of the ``pixel_info`` dataset of the synthetic map used below.
TEST_PIXEL_INFO_DTYPE = [
    ('nside', 'u4'),
    ('healpix_index', 'u8'),
    ('converged', 'u1'),
    ('DM_reliable_min', 'f4'),
    ('DM_reliable_max', 'f4'),
    ('n_stars', 'u4'),
    ('n_good', 'u4'),
    ('n_dwarfs', 'u4')
]


def make_test_map(fname, n_samples=5, n_distances=6, seed=0):
    """
    Writes a small HDF5 file in the same format as the published Bayestar
    maps, for testing.

    Two nside levels are used. ``nside = 1`` has 12 pixels and ``nside = 2``
    has 48, and pixels ``4p`` to ``4p+3`` of ``nside = 2`` are the four children
    of ``nside = 1`` pixel ``p``. Keeping pixels 0 to 5 and 0 to 23 therefore
    leaves the same half of the sky empty at both levels, so that queries
    outside of the map are exercised.

    Args:
        fname (:obj:`str`): Filename to write the map to.
        n_samples (Optional[:obj:`int`]): Number of samples per pixel.
        n_distances (Optional[:obj:`int`]): Number of distance bins.
        seed (Optional[:obj:`int`]): Seed for the random reddening.

    Returns:
        The ``pixel_info``, ``DM_bin_edges``, ``samples`` and ``best_fit``
        arrays that were written.
    """
    rng = np.random.RandomState(seed)

    nside = np.concatenate([np.full(6, 1, dtype='u4'),
                            np.full(24, 2, dtype='u4')])
    healpix_index = np.concatenate([np.arange(6, dtype='u8'),
                                    np.arange(24, dtype='u8')])
    n_pix = nside.size

    pixel_info = np.empty(n_pix, dtype=TEST_PIXEL_INFO_DTYPE)
    pixel_info['nside'] = nside
    pixel_info['healpix_index'] = healpix_index
    pixel_info['converged'] = rng.randint(0, 2, n_pix)
    pixel_info['DM_reliable_min'] = rng.uniform(2., 6., n_pix)
    pixel_info['DM_reliable_max'] = rng.uniform(10., 18., n_pix)
    pixel_info['n_stars'] = rng.randint(1, 100, n_pix)
    pixel_info['n_good'] = rng.randint(1, 50, n_pix)
    pixel_info['n_dwarfs'] = rng.randint(0, 5, n_pix)

    # Pixels with unreliable distance estimates, to exercise the NaNs
    pixel_info['DM_reliable_min'][3] = np.nan
    pixel_info['DM_reliable_max'][7] = np.nan

    DM_bin_edges = np.linspace(4., 20., n_distances).astype('f4')
    samples = rng.uniform(0., 2., size=(n_pix, n_samples, n_distances))
    samples = samples.astype('f4')
    best_fit = np.median(samples, axis=1).astype('f4')

    with h5py.File(fname, 'w') as f:
        dset = f.create_dataset('pixel_info', data=pixel_info)
        dset.attrs['DM_bin_edges'] = DM_bin_edges
        f.create_dataset('samples', data=samples, chunks=True,
                         compression='gzip')
        f.create_dataset('best_fit', data=best_fit, chunks=True,
                         compression='gzip')

    return pixel_info, DM_bin_edges, samples, best_fit


class TestBayestarRepack(unittest.TestCase):
    """
    Tests of the repacking of the Bayestar maps, and of the memory-mapped
    query path. A small synthetic map is used, so that these tests do not
    require any of the maps to have been downloaded.
    """

    n_samples = 5
    n_distances = 6

    @classmethod
    def setUpClass(cls):
        print('Building a synthetic Bayestar map ...')

        cls._tmpdir = tempfile.mkdtemp()
        cls._orig_fname = os.path.join(cls._tmpdir, 'orig.h5')
        cls._repacked_fname = os.path.join(cls._tmpdir, 'repacked.h5')

        (cls._pixel_info,
         cls._DM_bin_edges,
         cls._samples,
         cls._best_fit) = make_test_map(cls._orig_fname,
                                        n_samples=cls.n_samples,
                                        n_distances=cls.n_distances)

        bayestar.h5_repack(cls._orig_fname, cls._repacked_fname)

        cls._memmap = bayestar.BayestarQuery(cls._repacked_fname, memmap=True)
        cls._in_memory = bayestar.BayestarQuery(cls._repacked_fname,
                                                memmap=False)

        rng = np.random.RandomState(1234)
        cls._n_coords = 200
        l = rng.uniform(0., 360., cls._n_coords)
        b = np.degrees(np.arcsin(rng.uniform(-1., 1., cls._n_coords)))
        cls._coords = coords.SkyCoord(l, b, unit='deg', frame='galactic')
        cls._coords_dist = coords.SkyCoord(
            l, b,
            distance=rng.uniform(0.05, 12., cls._n_coords),
            unit='deg', frame='galactic')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._tmpdir)

    def assert_same(self, coords, mode, pct=None, return_flags=False):
        """
        Checks that a memory-mapped query gives exactly the same result as one
        that has the whole map in memory.
        """
        kwargs = dict(mode=mode, pct=pct, return_flags=return_flags)

        np.random.seed(99)
        result_memmap = self._memmap.query(coords, **kwargs)
        np.random.seed(99)
        result_in_memory = self._in_memory.query(coords, **kwargs)

        if return_flags:
            ret_memmap, flags_memmap = result_memmap
            ret_in_memory, flags_in_memory = result_in_memory
            for name in flags_in_memory.dtype.names:
                np.testing.assert_equal(flags_memmap[name],
                                        flags_in_memory[name])
        else:
            ret_memmap = result_memmap
            ret_in_memory = result_in_memory

        np.testing.assert_allclose(ret_memmap, ret_in_memory,
                                   rtol=0., atol=0., equal_nan=True)

    def test_repacked_layout(self):
        """
        The repacked file stores blocks of pixels together, and keeps the
        pixel information unchanged.
        """
        with h5py.File(self._repacked_fname, 'r') as f:
            self.assertEqual(f['samples'].chunks,
                             (bayestar.CHUNK_PIXELS, self.n_samples,
                              self.n_distances))
            self.assertEqual(f['best_fit'].chunks,
                             (bayestar.CHUNK_PIXELS, self.n_distances))

            # Field by field, because the pixel information contains NaNs
            pixel_info = f['pixel_info'][:]
            self.assertEqual(pixel_info.dtype.names,
                             self._pixel_info.dtype.names)
            for name in self._pixel_info.dtype.names:
                np.testing.assert_allclose(pixel_info[name],
                                           self._pixel_info[name],
                                           equal_nan=True)
            np.testing.assert_equal(f['pixel_info'].attrs['DM_bin_edges'],
                                    self._DM_bin_edges)

    def test_repack_preserves_data(self):
        """Repacking does not change the reddening."""
        with h5py.File(self._repacked_fname, 'r') as f:
            np.testing.assert_equal(f['samples'][:], self._samples)
            np.testing.assert_equal(f['best_fit'][:], self._best_fit)

    def test_h5_is_repacked(self):
        """Only a file in the repacked layout is recognized as repacked."""
        self.assertTrue(bayestar.h5_is_repacked(self._repacked_fname))
        self.assertFalse(bayestar.h5_is_repacked(self._orig_fname))
        self.assertFalse(bayestar.h5_is_repacked(
            os.path.join(self._tmpdir, 'does_not_exist.h5')))

    def test_memmap_is_lazy(self):
        """Memory mapping does not read the reddening until it is queried."""
        self.assertIsInstance(self._memmap._samples, h5py.Dataset)
        self.assertIsInstance(self._memmap._best_fit, h5py.Dataset)
        self.assertIsInstance(self._in_memory._samples, np.ndarray)
        self.assertIsInstance(self._in_memory._best_fit, np.ndarray)

    def test_memmap_matches_in_memory(self):
        """
        Memory-mapped queries agree with in-memory queries, in every mode,
        with and without distances, and with and without flags.
        """
        modes = ['random_sample', 'random_sample_per_pix', 'samples',
                 'median', 'mean', 'best', 'percentile']

        for coords in (self._coords, self._coords_dist):
            for mode in modes:
                pct = 33.3 if mode == 'percentile' else None
                for return_flags in (False, True):
                    self.assert_same(coords, mode, pct=pct,
                                     return_flags=return_flags)

    def test_out_of_bounds(self):
        """
        Coordinates outside of the map give NaN, and do not break the reading
        of the pixels that are inside of it.
        """
        ret = self._memmap.query(self._coords, mode='best')
        missed = np.isnan(ret).all(axis=1)

        # The synthetic map covers only half of the sky
        self.assertTrue(missed.any())
        self.assertTrue((~missed).any())

        # A query in which every coordinate is outside of the map
        outside = self._coords[missed]
        ret = self._memmap.query(outside, mode='best')
        self.assertTrue(np.isnan(ret).all())

        self.assert_same(outside, 'best')
        self.assert_same(outside, 'samples')
        self.assert_same(self._coords, 'samples', return_flags=True)

    def test_max_samples(self):
        """Only ``max_samples`` samples are returned when memory mapping."""
        q = bayestar.BayestarQuery(self._repacked_fname, max_samples=3,
                                   memmap=True)
        few = q.query(self._coords, mode='samples')
        full = self._memmap.query(self._coords, mode='samples')

        self.assertEqual(few.shape[1:], (3, self.n_distances))
        np.testing.assert_allclose(few, full[:,:3,:], equal_nan=True)

    def test_unsorted_pixels_rejected(self):
        """A map whose pixels are not sorted is rejected."""
        fname = os.path.join(self._tmpdir, 'unsorted.h5')
        make_test_map(fname)

        with h5py.File(fname, 'r+') as f:
            pixel_info = f['pixel_info'][:][::-1]
            del f['pixel_info']
            dset = f.create_dataset('pixel_info', data=pixel_info)
            dset.attrs['DM_bin_edges'] = self._DM_bin_edges

        self.assertRaises(ValueError, bayestar.BayestarQuery, fname)

    def test_fetch_repacks_file_on_disk(self):
        """
        A map that is already on disk in the original layout is repacked by
        :obj:`fetch`, rather than being downloaded again.
        """
        tmpdir = tempfile.mkdtemp()

        try:
            os.makedirs(os.path.join(tmpdir, 'bayestar'))
            fname = os.path.join(tmpdir, 'bayestar', 'bayestar2019.h5')
            make_test_map(fname)
            with h5py.File(fname, 'r') as f:
                n_pix = f['samples'].shape[0]
            self.assertFalse(bayestar.h5_is_repacked(fname))

            with mock.patch.object(bayestar, 'data_dir', lambda: tmpdir):
                with mock.patch.dict(bayestar.REPACKED_DSETS, {
                        'bayestar2019': {
                            'samples': (n_pix, self.n_samples,
                                        self.n_distances),
                            'best_fit': (n_pix, self.n_distances),
                            'pixel_info': (n_pix,)
                        }}):
                    with mock.patch.object(
                            bayestar.fetch_utils, 'dataverse_download_doi',
                            side_effect=AssertionError('should not download')):
                        bayestar.fetch('bayestar2019')

            self.assertTrue(bayestar.h5_is_repacked(fname))

            with h5py.File(fname, 'r') as f:
                np.testing.assert_equal(f['samples'][:], self._samples)
                np.testing.assert_equal(f['best_fit'][:], self._best_fit)
        finally:
            shutil.rmtree(tmpdir)


if __name__ == '__main__':
    unittest.main()
