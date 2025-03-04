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

import numpy as np
import astropy.coordinates as coords
import astropy.units as units
import os
import re
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



if __name__ == '__main__':
    unittest.main()
