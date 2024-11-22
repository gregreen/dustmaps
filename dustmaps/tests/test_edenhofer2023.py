#!/usr/bin/env python
#
# test_edenhofer2023.py
# Tests the query code for the Edenhofer2023 dust map (Edenhofer et al. 2023).
#
# Copyright (C) 2023  Gordian Edenhofer and Gregory M. Green
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

import sys
import time
import unittest
from functools import partial

import numpy as np
from astropy.coordinates import SkyCoord

from .. import edenhofer2023

log = partial(print, file=sys.stderr)


def random_coords(rng, n_dim, min_r=40., max_r=4e+3, n_max_elements=7):
    shp = rng.integers(1, n_max_elements, size=(n_dim,))
    l = rng.uniform(-180., +180., size=shp)
    b = rng.uniform(-90., +90., size=shp)
    dist = rng.uniform(min_r, max_r, size=shp)
    return SkyCoord(
        l=l, b=b, distance=dist, unit=('deg','deg','pc'), frame='galactic'
    )


class TestEdenhofer2023(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        msg = 'Loading all possible combinations of the data for {}...'
        log(msg.format(cls.__name__))

        msg_timed = 'Loaded {} {} in {:.5f} s'
        fmt_timed = partial(msg_timed.format, cls.__name__)

        t0 = time.time()
        cls._query_wo_smpls = edenhofer2023.Edenhofer2023Query(
            load_samples=False
        )
        log(fmt_timed("density w/o samples", time.time() - t0))

        t0 = time.time()
        cls._query_w_smpls = edenhofer2023.Edenhofer2023Query(
            load_samples=True
        )
        log(fmt_timed("density w/ samples", time.time() - t0))

        t0 = time.time()
        cls._query_wo_smpls_int = edenhofer2023.Edenhofer2023Query(
            load_samples=False, integrated=True
        )
        log(fmt_timed("integrated w/o samples", time.time() - t0))

        # Do not test integrated samples b/c it is too memory intensive :(
        # t0 = time.time()
        # cls._query_w_smpls_int = edenhofer2023.Edenhofer2023Query(
        #     load_samples=True, integrated=True
        # )
        # log(fmt_timed("integrated w/ samples", time.time() - t0))

    def test_samples_no_samples_consistency(self):
        for seed in (42, 314, 271828):
            rng = np.random.default_rng(seed)
            for mode in ("mean", "std"):
                n_dim = rng.integers(1, 5)
                coords = random_coords(rng, n_dim, min_r=59.0, max_r=1.3e+3)
                r1 = self._query_wo_smpls(coords, mode=mode)
                assert r1.shape == coords.shape
                r2 = self._query_w_smpls(coords, mode=mode)
                assert r2.shape == coords.shape
                # It makes a difference whether we inerpolate the mean or the
                # samples. Hence, allow for some "significant" tolerance.
                np.testing.assert_allclose(r1, r2, atol=2e-4, rtol=0.1)

                if mode == "mean":
                    r1i = self._query_wo_smpls_int(coords, mode=mode)
                    assert r1i.shape == coords.shape
                    np.testing.assert_array_equal(
                        (r1i > r1)[~np.isnan(r1)], True
                    )
                # The following would be too expensive memory-wise :(
                # r2i = self._query_w_smpls_int(coords, mode=mode)
                # assert r2i.shape == coords.shape
                # np.testing.assert_allclose(r1i, r2i, atol=1e-3, rtol=1e-3)

    def test_samples_shape(self):
        mode = "samples"
        for seed in (42, 314, 271828):
            rng = np.random.default_rng(seed)
            n_dim = rng.integers(1, 5)
            coords = random_coords(rng, n_dim)
            r_density = self._query_w_smpls(coords, mode=mode)
            assert r_density.shape[:-1] == coords.shape
            # The following would be too expensive memory-wise :(
            # r_int = self._query_w_smpls_int(coords, mode=mode)
            # assert r_int.shape[:-1] == coords.shape
            # np.testing.assert_equal(r_int > r_density, True)

    def test_random_sample_shape(self):
        mode = "random_sample"
        for seed in (42, 314, 271828):
            rng = np.random.default_rng(seed)
            n_dim = rng.integers(1, 5)
            coords = random_coords(rng, n_dim)
            r_density = self._query_w_smpls(coords, mode=mode)
            assert r_density.shape == coords.shape

    def test_monotonicty_of_integrated(self):
        mode = "mean"
        for seed in (42, 314, 271828):
            rng = np.random.default_rng(seed)
            n_dim = rng.integers(1, 4)
            coords = random_coords(rng, n_dim, min_r=65.0, max_r=75.0)
            for query in (self._query_wo_smpls_int,):
                r1 = query(coords, mode=mode)
                for _ in range(4):
                    coords = SkyCoord(
                        l=coords.l,
                        b=coords.b,
                        distance=2. * coords.distance,
                        frame="galactic"
                    )
                    r2 = query(coords, mode=mode)
                    np.testing.assert_equal((r2 > r1)[~np.isnan(r1)], True)
                    r1 = r2


if __name__ == '__main__':
    unittest.main()
