#!/usr/bin/env python
#
# test_decaps.py
# Test query code for the Green, Schlafly, Finkbeiner et al. (2015) dust map.
#
# Copyright (C) 2025 Catherine Zucker, Gregory M. Green
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

from __future__ import print_function, division

import unittest

import numpy as np
import astropy.coordinates as coords
import astropy.units as units
import h5py
import healpy as hp
import os
import re
import shutil
import tempfile
import time
import pickle

from .. import decaps
from ..std_paths import *


def parse_decaps_output(fname):
    with open(fname, "rb") as f:
        output = pickle.load(f)
        
    return output


class TestDECaPS(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        print('Loading DECaPS query object ...')
        t0 = time.time()

        fname = os.path.join(test_dir, 'decaps_test_data.pkl')
        self._test_data = parse_decaps_output(fname)

        # Set up DECaPS query object
        self._decaps = decaps.DECaPSQueryLite()


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


    def test_equ_mean_far_scalar(self):
        """
        Test that mean reddening is correct in the far limit, using a single
        location on the sky at a time as input.
        """
        for d in self._test_data:
            c = self._get_gal(d, dist=1.e3*units.kpc)
            ebv_data = (d['mean'][-1])
            ebv_calc = self._decaps.query(c, mode='mean')
            np.testing.assert_allclose(ebv_data, ebv_calc, atol=0.001, rtol=0.0001)

    def test_equ_mean_far_vector(self):
        """
        Test that mean reddening is correct in the far limit, using a vector
        of coordinates as input.
        """
        l = [d['l']*units.deg for d in self._test_data]
        b = [d['b']*units.deg for d in self._test_data]
        dist = [1.e3*units.kpc for bb in b]
        c = coords.SkyCoord(l, b, distance=dist, frame='galactic')

        ebv_data = np.array([(d['mean'][-1]) for d in self._test_data])
        ebv_calc = self._decaps.query(c, mode='mean')

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
        ebv_calc = self._decaps.query(c, mode='random_sample')

        d_ebv = np.min(np.abs(ebv_data[:,:] - ebv_calc[:,None]), axis=1)


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
        ebv_calc = self._decaps.query(c, mode='samples')

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
            ebv_calc = self._decaps.query(c, mode='samples')

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
            ebv_calc = self._decaps.query(c, mode='random_sample')

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
        ebv_calc = self._decaps.query(c, mode='samples')


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
        ebv_calc = self._decaps.query(c, mode='random_sample')
        
        d_ebv = np.min(np.abs(ebv_data[:,:,:] - ebv_calc[:,None,:]), axis=1)
        np.testing.assert_allclose(d_ebv, 0., atol=0.001, rtol=0.0001)

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

                ebv_calc = self._decaps.query(c, mode=mode)

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



def pixel_centers(nside, healpix_index):
    """
    Returns the Galactic coordinates of the centers of the given HEALPix
    pixels. The HEALPix pixels are taken to be in the Galactic frame, as in
    :obj:`decaps.lb2pix`.
    """
    theta, phi = hp.pix2ang(nside, healpix_index, nest=True)

    l = np.degrees(phi)
    b = 90. - np.degrees(theta)

    return l, b


def make_test_map(fname, nside=16, parent_nside=2, parent=20, n_samples=5,
                  n_distances=6, with_samples=True, seed=0):
    """
    Writes a small HDF5 file in the same format as the published DECaPS maps,
    for testing.

    The map covers a single low-resolution HEALPix pixel: its ``nside``
    children, which occupy a contiguous block of the nested pixel numbering.
    The map therefore has a single ``nside`` and is sorted by
    ``healpix_index``, just as the published map is. The rest of the sky is
    left empty, so that queries outside of the map are also exercised.

    Args:
        fname (:obj:`str`): Filename to write the map to.
        nside (Optional[:obj:`int`]): The HEALPix :obj:`nside` of the map.
        parent_nside (Optional[:obj:`int`]): The :obj:`nside` of the
            low-resolution pixel that the map covers.
        parent (Optional[:obj:`int`]): The index of the low-resolution pixel
            that the map covers.
        n_samples (Optional[:obj:`int`]): Number of samples per pixel.
        n_distances (Optional[:obj:`int`]): Number of distance bins.
        with_samples (Optional[:obj:`bool`]): If ``True`` (the default), the
            file contains the samples, as in ``decaps_mean_and_samples.h5``.
            Otherwise, only the mean map is written, as in ``decaps_mean.h5``.
        seed (Optional[:obj:`int`]): Seed for the random reddening.

    Returns:
        The ``healpix_index``, ``DM_bin_edges``, ``mean`` and ``samples``
        arrays that were written. ``samples`` is :obj:`None` if
        ``with_samples`` is :obj:`False`.
    """
    rng = np.random.RandomState(seed)

    n_children = (nside // parent_nside)**2
    healpix_index = np.arange(parent * n_children, (parent + 1) * n_children)
    n_pix = healpix_index.size

    DM_bin_edges = np.linspace(4., 20., n_distances).astype('f4')
    mean = rng.uniform(0., 2., size=(n_pix, n_distances)).astype('f4')
    samples = rng.uniform(0., 2., size=(n_pix, n_samples, n_distances))
    samples = samples.astype('f4')

    DM_reliable_min = rng.uniform(2., 6., n_pix).astype('f4')
    DM_reliable_max = rng.uniform(10., 18., n_pix).astype('f4')
    converged = rng.randint(0, 2, n_pix)
    infilled = rng.randint(0, 2, n_pix)

    # Pixels with no reliable distance estimate, to exercise the NaNs
    DM_reliable_min[1] = np.nan
    DM_reliable_max[3] = np.nan

    with h5py.File(fname, 'w') as f:
        pixel_info = f.create_group('pixel_info')
        pixel_info.attrs['DM_bin_edges'] = DM_bin_edges
        pixel_info.attrs['n_samples'] = n_samples
        pixel_info.attrs['nside'] = nside
        pixel_info.create_dataset('healpix_index', data=healpix_index)
        pixel_info.create_dataset('DM_reliable_min', data=DM_reliable_min)
        pixel_info.create_dataset('DM_reliable_max', data=DM_reliable_max)
        pixel_info.create_dataset('converged', data=converged)
        pixel_info.create_dataset('infilled', data=infilled)

        f.create_dataset('mean', data=mean, chunks=True, compression='gzip')
        if with_samples:
            f.create_dataset('samples', data=samples, chunks=True,
                             compression='gzip')

    if not with_samples:
        samples = None

    return healpix_index, DM_bin_edges, mean, samples


class TestDECaPSMemmap(unittest.TestCase):
    """
    Tests of the memory-mapped query path of :obj:`DECaPSQuery`. A small
    synthetic map is used, so that these tests do not require any map to have
    been downloaded.
    """

    nside = 16
    n_samples = 5
    n_distances = 6

    @classmethod
    def setUpClass(cls):
        print('Building a synthetic DECaPS map ...')

        cls._tmpdir = tempfile.mkdtemp()
        cls._fname = os.path.join(cls._tmpdir, 'decaps.h5')
        cls._mean_fname = os.path.join(cls._tmpdir, 'decaps_mean.h5')

        (cls._healpix_index,
         cls._DM_bin_edges,
         cls._mean,
         cls._samples) = make_test_map(cls._fname,
                                       nside=cls.nside,
                                       n_samples=cls.n_samples,
                                       n_distances=cls.n_distances)

        make_test_map(cls._mean_fname, nside=cls.nside, with_samples=False,
                      n_distances=cls.n_distances)

        cls._memmap = decaps.DECaPSQuery(cls._fname, memmap=True)
        cls._in_memory = decaps.DECaPSQuery(cls._fname, memmap=False)

        rng = np.random.RandomState(1234)

        # Coordinates inside the map: the centers of its pixels, each queried
        # a few times over, in random order.
        l, b = pixel_centers(cls.nside, cls._healpix_index)
        sel = rng.randint(0, l.size, 200)
        l_in, b_in = l[sel], b[sel]

        # Coordinates outside of the map
        outside_pix = cls._healpix_index[0] + 1000 + np.arange(20)
        l_out, b_out = pixel_centers(cls.nside, outside_pix)

        cls._dist_in = 10.**rng.uniform(-2., 5., l_in.size)
        cls._dist_out = 10.**rng.uniform(-2., 5., l_out.size)

        cls._coords_in = coords.SkyCoord(l_in, b_in, unit='deg',
                                         frame='galactic')
        cls._coords_in_dist = coords.SkyCoord(l_in, b_in,
                                              distance=cls._dist_in,
                                              unit='deg', frame='galactic')

        l = np.concatenate([l_in, l_out])
        b = np.concatenate([b_in, b_out])
        dist = np.concatenate([cls._dist_in, cls._dist_out])

        cls._coords = coords.SkyCoord(l, b, unit='deg', frame='galactic')
        cls._coords_dist = coords.SkyCoord(l, b, distance=dist, unit='deg',
                                           frame='galactic')

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

    def test_memmap_is_lazy(self):
        """Memory mapping does not read the reddening until it is queried."""
        self.assertIsInstance(self._memmap._mean, h5py.Dataset)
        self.assertIsInstance(self._memmap._samples, h5py.Dataset)
        self.assertIsInstance(self._in_memory._mean, np.ndarray)
        self.assertIsInstance(self._in_memory._samples, np.ndarray)

    def test_memmap_matches_in_memory(self):
        """
        Memory-mapped queries agree with in-memory queries, in every mode,
        with and without distances, and with and without flags.
        """
        modes = ['random_sample', 'random_sample_per_pix', 'samples',
                 'median', 'mean', 'percentile']

        for coords in (self._coords, self._coords_dist):
            for mode in modes:
                pct = 33.3 if mode == 'percentile' else None
                for return_flags in (False, True):
                    self.assert_same(coords, mode, pct=pct,
                                     return_flags=return_flags)

    def test_single_coordinate(self):
        """A single coordinate can be queried while memory mapping."""
        l, b = pixel_centers(self.nside, self._healpix_index[:1])
        single = coords.SkyCoord(l, b, unit='deg', frame='galactic')
        single_dist = coords.SkyCoord(l, b, distance=self._dist_in[:1],
                                      unit='deg', frame='galactic')

        for c in (single, single_dist):
            for mode in ('mean', 'samples', 'median'):
                np.testing.assert_allclose(
                    self._memmap.query(c, mode=mode),
                    self._in_memory.query(c, mode=mode),
                    rtol=0., atol=0., equal_nan=True)

    def test_out_of_bounds(self):
        """
        Coordinates outside of the map give NaN, and do not break the reading
        of the pixels that are inside of it.
        """
        ret = self._memmap.query(self._coords, mode='mean')
        missed = np.isnan(ret).all(axis=1)

        # The synthetic map covers only a small part of the sky
        self.assertTrue(missed.any())
        self.assertTrue((~missed).any())

        # A query in which every coordinate is outside of the map
        outside = self._coords[missed]
        ret = self._memmap.query(outside, mode='mean')
        self.assertTrue(np.isnan(ret).all())

        self.assert_same(outside, 'mean')
        self.assert_same(outside, 'samples')
        self.assert_same(self._coords, 'samples', return_flags=True)

    def test_max_samples(self):
        """Only ``max_samples`` samples are returned when memory mapping."""
        q = decaps.DECaPSQuery(self._fname, max_samples=3, memmap=True)
        few = q.query(self._coords, mode='samples')
        full = self._memmap.query(self._coords, mode='samples')

        self.assertEqual(few.shape[1:], (3, self.n_distances))
        np.testing.assert_allclose(few, full[:,:3,:], equal_nan=True)

    def test_mean_only(self):
        """
        A file that contains only the mean map can be queried in ``mean`` mode,
        but no other mode.
        """
        q = decaps.DECaPSQuery(self._mean_fname, mean_only=True, memmap=True)

        np.testing.assert_allclose(
            q.query(self._coords, mode='mean'),
            self._memmap.query(self._coords, mode='mean'),
            rtol=0., atol=0., equal_nan=True)

        for mode in ('random_sample', 'random_sample_per_pix', 'samples',
                     'median', 'percentile'):
            self.assertRaises(ValueError, q.query, self._coords, mode=mode)

    def test_lite_query_agrees(self):
        """
        The memory-mapped :obj:`DECaPSQuery` agrees with
        :obj:`DECaPSQueryLite`, which memory-maps the file in a different way.
        Only the modes that :obj:`DECaPSQueryLite` supports are compared here:
        its ``mean`` mode raises an exception, because the mean map lacks the
        samples axis.
        """
        lite = decaps.DECaPSQueryLite(self._fname, contiguous=True)

        for coords in (self._coords_in, self._coords_in_dist):
            for mode in ('median', 'samples'):
                np.testing.assert_allclose(
                    self._memmap.query(coords, mode=mode),
                    lite.query(coords, mode=mode),
                    rtol=0., atol=0., equal_nan=True)

    def test_unsorted_pixels_rejected(self):
        """A map whose pixels are not sorted is rejected."""
        fname = os.path.join(self._tmpdir, 'unsorted.h5')
        healpix_index = make_test_map(fname, nside=self.nside,
                                      n_distances=self.n_distances)[0]

        with h5py.File(fname, 'r+') as f:
            del f['pixel_info/healpix_index']
            f['pixel_info'].create_dataset('healpix_index',
                                           data=healpix_index[::-1])

        self.assertRaises(ValueError, decaps.DECaPSQuery, fname)


if __name__ == '__main__':
    unittest.main()
