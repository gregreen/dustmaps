#!/usr/bin/env python
#
# test_planck.py
# Test query code for the Planck Collaboration (2013) dust reddening map.
#
# Copyright (C) 2018  Gregory M. Green
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

from .. import planck
from ..std_paths import *


class TestPlanck(unittest.TestCase):
    component = 'extragalactic'
    load_errors = False

    @classmethod
    def setUpClass(self):
        print('Loading Planck {} dust map ...'.format(self.component))
        t0 = time.time()

        # Set up Planck query object
        if self.component == 'GNILC':
            self._planck = planck.PlanckGNILCQuery(load_errors=self.load_errors)
        else:
            self._planck = planck.PlanckQuery(component=self.component)

        t1 = time.time()
        print('Loaded Planck test data in {:.5f} s.'.format(t1-t0))

    def test_shape(self):
        """
        Test that the output shapes are as expected with input coordinate arrays
        of different shapes.
        """

        for reps in range(5):
            # Draw random coordinates, with different shapes
            n_dim = np.random.randint(1,4)
            shape = np.random.randint(1,7, size=(n_dim,))

            ra = np.random.uniform(-180., 180., size=shape) * units.deg
            dec = np.random.uniform(-90., 90., size=shape) * units.deg
            c = coords.SkyCoord(ra, dec, frame='icrs')

            E = self._planck(c)

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
        E0 = self._planck(c)

        for fr in frames:
            cc = c.transform_to(fr)
            E = self._planck(cc)
            np.testing.assert_equal(E, E0)

        u,v,w = np.random.uniform(0., 5., size=(3,100))
        try:
            c = coords.SkyCoord(
                u=u, v=v, w=w,
                unit='kpc',
                representation_type='cartesian',
                frame='galactic'
            )
        except ValueError as err:
            # Astropy version < 3.0
            c = coords.SkyCoord(
                u=u, v=v, w=w,
                unit='kpc',
                representation='cartesian',
                frame='galactic'
            )
        E0 = self._planck(c)

        for fr in frames:
            cc = c.transform_to(fr)
            E = self._planck(cc)
            np.testing.assert_equal(E, E0)



class TestPlanckTau(TestPlanck):
    component = 'tau'


class TestPlanckRadiance(TestPlanck):
    component = 'radiance'


class TestPlanckTemperature(TestPlanck):
    component = 'temperature'


class TestPlanckBeta(TestPlanck):
    component = 'beta'


class TestPlanckTemperature(TestPlanck):
    component = 'temperature'


class TestPlanckBeta(TestPlanck):
    component = 'beta'


class TestPlanckTemperatureErr(TestPlanck):
    component = 'err_temp'


class TestPlanckBetaErr(TestPlanck):
    component = 'err_beta'


class TestPlanckBetaErr(TestPlanck):
    component = 'GNILC'


class TestPlanckBetaErr(TestPlanck):
    component = 'GNILC'
    load_errors = True


if __name__ == '__main__':
    unittest.main()
