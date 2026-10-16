#!/usr/bin/env python
#
# test_edenhofer2023_repack.py
# Tests the memory-mapped (repacked FITS) reader for the Edenhofer2023 dust
# map.
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

from __future__ import division, print_function

import os
import shutil
import tempfile
import unittest

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

from .. import edenhofer2023

NSIDE = 16
N_RADII = 100
N_SAMPLES = 4
N_PIX = 12 * NSIDE**2
CHECKPOINT_STEP = edenhofer2023.CHECKPOINT_STEP
N_CHECK = int(np.ceil(N_RADII / CHECKPOINT_STEP))


def smooth(shape, seed=0):
    """A positive, smooth field, more like a dust map than noise is."""
    rng = np.random.default_rng(seed)
    field = np.zeros(shape, dtype='f8')
    for axis, n in enumerate(shape):
        t = np.arange(n) / n * 2. * np.pi
        one_axis = [1] * len(shape)
        one_axis[axis] = n
        field += (np.sin(t) + 0.5 * np.cos(3. * t)).reshape(one_axis)
    field /= len(shape)
    return (1e-3 * (2. + field) + 1e-6 * rng.random(shape)).astype('f4')


def make_test_map(fname, with_samples=False, seed=0):
    """
    A small FITS file with the HDU structure of the published maps, written
    with the same float32 precision that they use.
    """
    coo_bounds = np.linspace(0., 2000., N_RADII + 1).astype('f4')
    radii = 0.5 * (coo_bounds[1:] + coo_bounds[:-1]).astype('f4')
    if with_samples:
        name = 'samples'
        density = smooth((N_SAMPLES, N_RADII, N_PIX), seed)
        data0 = smooth((N_SAMPLES, N_PIX), seed + 1)
        data0_unc = smooth((N_SAMPLES, N_PIX), seed + 2)
    else:
        name = 'mean'
        density = smooth((N_RADII, N_PIX), seed)
        data0 = smooth((N_PIX,), seed + 1)
        data0_unc = smooth((N_PIX,), seed + 2)
    std = smooth(density.shape, seed + 3)

    main_header = fits.Header([('NSIDE', NSIDE), ('ORDERING', 'NEST'),
                               ('CUNIT', 'E_per_pc')])
    other_header = fits.Header([('NSIDE', NSIDE), ('ORDERING', 'NEST')])
    hdus = [
        fits.PrimaryHDU(),
        fits.ImageHDU(data=density, name=name, header=main_header),
        fits.ImageHDU(data=std, name='std.', header=other_header.copy()),
        fits.BinTableHDU.from_columns(
            [fits.Column(name='radial pixel centers', format='E',
                         array=radii)]),
        fits.BinTableHDU.from_columns(
            [fits.Column(name='radial pixel boundaries', format='E',
                         array=coo_bounds)]),
        fits.ImageHDU(
            data=data0, name='mean of integrated inner 68.8 pc',
            header=fits.Header([
                ('NSIDE', NSIDE), ('ORDERING', 'NEST'),
                ('CTYPE', 'inner density integrated within 68.8pc')])),
        fits.ImageHDU(
            data=data0_unc, name='std. of integrated inner 68.8 pc',
            header=fits.Header([
                ('NSIDE', NSIDE), ('ORDERING', 'NEST'),
                ('CTYPE', 'std. of inner density integrated within 68.8pc')])),
    ]
    fits.HDUList(hdus).writeto(fname, overwrite=True)


def random_coords(seed=0, n=64, min_r=60., max_r=1900.):
    rng = np.random.default_rng(seed)
    lon = rng.uniform(0., 360., n)
    lat = np.degrees(np.arcsin(rng.uniform(-1., 1., n)))
    dist = rng.uniform(min_r, max_r, n)
    return SkyCoord(l=lon * u.deg, b=lat * u.deg, distance=dist * u.pc,
                    frame='galactic')


class RepackTestCase(unittest.TestCase):
    """Repacking a map, and reading it back the same way as the FITS file."""
    with_samples = False

    @classmethod
    def setUpClass(cls):
        cls._dir = tempfile.mkdtemp(prefix='edenhofer2023_repack_')
        if cls.with_samples:
            cls._fits = os.path.join(cls._dir, 'samples_healpix.fits')
        else:
            cls._fits = os.path.join(cls._dir, 'mean_and_std_healpix.fits')
        cls._repacked = edenhofer2023._repacked_fname(cls._fits)
        make_test_map(cls._fits, with_samples=cls.with_samples)
        edenhofer2023.fits_repack(cls._fits, cls._repacked)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._dir, ignore_errors=True)

    def query(self, fname, **kwargs):
        return edenhofer2023.Edenhofer2023Query(map_fname=fname, **kwargs)

    def mmap(self, **kwargs):
        return self.query(self._repacked, memmap=True, **kwargs)

    def fits(self, **kwargs):
        """The published FITS map, transformed in memory rather than mapped."""
        return self.query(self._fits, memmap=False, **kwargs)

    def test_the_file_is_marked_as_repacked(self):
        self.assertTrue(edenhofer2023.fits_is_repacked(self._repacked))
        self.assertFalse(edenhofer2023.fits_is_repacked(self._fits))

    def test_checkpoints_are_the_reference_cumsum(self):
        """A checkpoint is the value the reference implementation computes."""
        rec = edenhofer2023._get_sphere(self._fits)
        work = np.array(rec.data, dtype='f4')
        dvol = np.diff(rec.coo_bounds)
        np.multiply(dvol[:, np.newaxis], work, out=work)
        work[..., 0, :] += rec.data0
        np.cumsum(work, axis=-2, out=work)
        expected = edenhofer2023._stored(work[..., ::CHECKPOINT_STEP, :])
        with fits.open(self._repacked, 'readonly', memmap=True) as hdul:
            np.testing.assert_array_equal(hdul['CHECKPOINTS'].data, expected)

    def test_the_memory_budget_does_not_change_the_file(self):
        """However wide the ranges it reads, the repacked file is the same."""
        chunked = self._repacked + '.chunked'
        try:
            edenhofer2023.fits_repack(self._fits, chunked, max_memory=1 << 15)
            with open(self._repacked, 'rb') as whole:
                with open(chunked, 'rb') as pieces:
                    self.assertEqual(whole.read(), pieces.read())
        finally:
            if os.path.isfile(chunked):
                os.remove(chunked)

    def test_the_density_is_stored_as_published(self):
        rec = edenhofer2023._get_sphere(self._fits)
        with fits.open(self._repacked, 'readonly', memmap=True) as hdul:
            np.testing.assert_array_equal(
                hdul['DATA'].data, edenhofer2023._stored(rec.data))
            np.testing.assert_array_equal(
                hdul['STD'].data,
                edenhofer2023._stored(rec.data_uncertainty))

    def test_a_pixels_values_are_stored_together(self):
        samples = (N_SAMPLES,) if self.with_samples else ()
        with fits.open(self._repacked, 'readonly', memmap=True) as hdul:
            data, check = hdul['DATA'], hdul['CHECKPOINTS']
            self.assertEqual(data.shape, (N_PIX, N_RADII) + samples)
            self.assertEqual(check.shape, (N_PIX, N_CHECK) + samples)
            # The FITS axes run the other way around, so the innermost axis of
            # the stored array, whose values a query reads together, is NAXIS1
            self.assertEqual(data.header['NAXIS1'], data.shape[-1])
            self.assertEqual(data.header['NAXIS2'], data.shape[-2])
            self.assertEqual(data.header['NAXIS'], len(data.shape))
            if self.with_samples:
                self.assertEqual(data.header['NAXIS3'], N_PIX)

    def test_the_radial_axis_matches_the_published_map(self):
        rec = edenhofer2023._get_sphere(self._fits)
        query = self.mmap()
        np.testing.assert_array_equal(
            query.distances.value, np.asarray(rec.radii))
        np.testing.assert_array_equal(
            query.distance_bounds.value, np.asarray(rec.coo_bounds))

    def test_density_queries_match_the_fits_reader(self):
        ref = self.fits()
        mmap = self.mmap()
        coords = random_coords(seed=0)
        for mode in ('mean', 'std'):
            np.testing.assert_array_equal(
                mmap(coords, mode=mode), ref(coords, mode=mode))

    def test_integrated_queries_match_the_fits_reader(self):
        ref = self.fits(integrated=True)
        mmap = self.mmap(integrated=True)
        coords = random_coords(seed=1)
        # Not merely to float rounding: a checkpoint is the value that the
        # reference implementation accumulates, and the terms between a
        # checkpoint and a queried bin are added one at a time, in the order
        # that the reference accumulates them, so the result is identical
        np.testing.assert_array_equal(mmap(coords), ref(coords))

    def test_integrated_can_be_asked_for_per_query(self):
        coords = random_coords(seed=2)
        mmap = self.mmap()
        ref = self.mmap(integrated=True)
        np.testing.assert_array_equal(mmap(coords, integrated=True),
                                      ref(coords))
        # and the density is still available from the same object
        np.testing.assert_array_equal(
            mmap(coords), mmap(coords, integrated=False))

    def test_integrated_cannot_be_changed_for_a_fits_map(self):
        query = self.fits()
        with self.assertRaises(ValueError):
            query(random_coords(seed=3), integrated=True)

    def test_the_repacked_file_can_be_named_directly(self):
        coords = random_coords(seed=8)
        query = self.query(self._repacked)
        np.testing.assert_array_equal(query(coords, integrated=True),
                                      self.mmap(integrated=True)(coords))

    def test_the_repacked_file_read_into_memory_matches_the_fits_reader(self):
        """It holds the same map, so it can be read into memory as well."""
        coords = random_coords(seed=11)
        in_memory = self.query(self._repacked, memmap=False)
        ref = self.fits()
        for mode in ('mean', 'std'):
            np.testing.assert_array_equal(in_memory(coords, mode=mode),
                                          ref(coords, mode=mode))

    def test_the_repacked_file_read_into_memory_integrates_the_same(self):
        coords = random_coords(seed=12)
        in_memory = self.query(self._repacked, memmap=False, integrated=True)
        np.testing.assert_array_equal(
            in_memory(coords), self.fits(integrated=True)(coords))

    def test_the_repacked_map_is_the_default(self):
        coords = random_coords(seed=9)
        query = edenhofer2023.Edenhofer2023Query(map_fname=self._fits)
        np.testing.assert_array_equal(query(coords), self.mmap()(coords))

    def test_out_of_range_queries_are_nan(self):
        far = SkyCoord(l=[0., 10.] * u.deg, b=[0., 5.] * u.deg,
                       distance=[1e6, 1e6] * u.pc, frame='galactic')
        mmap = self.mmap(integrated=True)
        self.assertTrue(np.all(np.isnan(mmap(far))))
        self.assertTrue(np.all(np.isnan(mmap(far, integrated=False))))


class TestEdenhofer2023Repack(RepackTestCase):
    with_samples = False

    def test_a_map_without_samples_has_no_samples(self):
        self.assertIsNone(self.mmap().n_samples)
        self.assertIsNone(self.fits().n_samples)

    def test_std_is_not_available_from_an_integrated_map_without_samples(self):
        mmap = self.mmap(integrated=True)
        with self.assertRaises(ValueError):
            mmap(random_coords(seed=4), mode='std')


class TestEdenhofer2023RepackSamples(RepackTestCase):
    with_samples = True

    def test_the_map_has_samples(self):
        self.assertEqual(self.mmap().n_samples, N_SAMPLES)
        self.assertEqual(self.fits().n_samples, N_SAMPLES)

    def test_every_mode_matches_the_fits_reader(self):
        ref = self.fits()
        mmap = self.mmap()
        coords = random_coords(seed=5)
        for mode in ('mean', 'std', 'samples'):
            np.testing.assert_array_equal(
                mmap(coords, mode=mode), ref(coords, mode=mode))
        # The same sample is drawn from the same seed by both readers
        np.testing.assert_array_equal(
            mmap(coords, mode='random_sample', seed_or_rng=42),
            ref(coords, mode='random_sample', seed_or_rng=42))

    def test_integrated_queries_match_the_fits_reader(self):
        ref = self.fits(integrated=True)
        mmap = self.mmap(integrated=True)
        coords = random_coords(seed=6)
        for mode in ('mean', 'std', 'samples', 'random_sample'):
            np.testing.assert_array_equal(
                mmap(coords, mode=mode, seed_or_rng=7),
                ref(coords, mode=mode, seed_or_rng=7))

    def test_every_mode_matches_the_fits_reader_read_into_memory(self):
        coords = random_coords(seed=13)
        in_memory = self.query(self._repacked, memmap=False)
        ref = self.fits()
        for mode in ('mean', 'std', 'samples'):
            np.testing.assert_array_equal(in_memory(coords, mode=mode),
                                          ref(coords, mode=mode))


class TestEdenhofer2023RepackMissing(unittest.TestCase):
    """Asking for a memory map that has not been written yet."""

    def test_memmap_needs_a_repacked_file(self):
        d = tempfile.mkdtemp(prefix='edenhofer2023_repack_missing_')
        try:
            fname = os.path.join(d, 'mean_and_std_healpix.fits')
            make_test_map(fname)
            with self.assertRaises(FileNotFoundError) as caught:
                edenhofer2023.Edenhofer2023Query(map_fname=fname, memmap=True)
            self.assertIn('fetch(', str(caught.exception))
            # The message names the function to call, spelling out its
            # arguments, so match the name rather than the call itself
            self.assertIn('edenhofer2023.fetch(', str(caught.exception))
            # The same map is still readable without a memory map
            edenhofer2023.Edenhofer2023Query(map_fname=fname, memmap=False)
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_a_map_that_was_never_downloaded(self):
        d = tempfile.mkdtemp(prefix='edenhofer2023_repack_absent_')
        try:
            fname = os.path.join(d, 'mean_and_std_healpix.fits')
            with self.assertRaises(FileNotFoundError):
                edenhofer2023.Edenhofer2023Query(map_fname=fname)
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_a_map_repacked_before_the_inner_map_was_kept(self):
        """
        Reading such a file into memory would quietly integrate without the
        inner map, so it is refused; the memory map does not need it.
        """
        d = tempfile.mkdtemp(prefix='edenhofer2023_repack_older_')
        try:
            fname = os.path.join(d, 'mean_and_std_healpix.fits')
            make_test_map(fname)
            repacked = edenhofer2023._repacked_fname(fname)
            edenhofer2023.fits_repack(fname, repacked)
            with fits.open(repacked, mode='update') as hdul:
                del hdul[0].header['HASINNER']

            with self.assertRaises(ValueError):
                edenhofer2023.Edenhofer2023Query(map_fname=repacked,
                                                 memmap=False)
            edenhofer2023.Edenhofer2023Query(
                map_fname=repacked, memmap=True)(random_coords(seed=14))
        finally:
            shutil.rmtree(d, ignore_errors=True)


class TestEdenhofer2023RepackOnly(unittest.TestCase):
    """
    Once a file has been repacked, the FITS file that it was made from is no
    longer needed: the repacked file holds the same map, and both ways of
    reading it give what the FITS file gave.
    """

    def query(self, fname, **kwargs):
        return edenhofer2023.Edenhofer2023Query(map_fname=fname, **kwargs)

    def test_the_repacked_file_stands_in_for_the_fits_file(self):
        d = tempfile.mkdtemp(prefix='edenhofer2023_repack_only_')
        try:
            fname = os.path.join(d, 'mean_and_std_healpix.fits')
            make_test_map(fname)
            coords = random_coords(seed=10)
            density = self.query(fname, memmap=False)(coords)
            integrated = self.query(fname, memmap=False,
                                    integrated=True)(coords)

            # A budget small enough that the map is read in many ranges
            edenhofer2023.repack_and_remove(fname, max_memory=1 << 15)
            self.assertFalse(os.path.isfile(fname))
            self.assertFalse(os.path.isfile(
                edenhofer2023._repacked_fname(fname) + '.part'))
            self.assertTrue(edenhofer2023.fits_is_repacked(
                edenhofer2023._repacked_fname(fname)))

            for memmap in (True, False):
                query = self.query(fname, memmap=memmap)
                np.testing.assert_array_equal(query(coords), density)
                if memmap:
                    # Only a mapped map can be integrated per query
                    np.testing.assert_array_equal(
                        query(coords, integrated=True), integrated)
                np.testing.assert_array_equal(
                    self.query(fname, memmap=memmap,
                               integrated=True)(coords), integrated)
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
