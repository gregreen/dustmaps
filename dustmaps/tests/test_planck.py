#!/usr/bin/env python
#
# test_planck.py
# Test query code for the Planck Collaboration (2013) dust reddening map.
#
# Copyright (C) 2018  Gregory M. Green
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
import os
import re
import time

from .. import planck
from ..std_paths import *


class TestPlanck(unittest.TestCase):
    component = 'extragalactic'

    @classmethod
    def setUpClass(self):
        print('Loading Planck {} dust map ...'.format(self.component))
        t0 = time.time()

        # Set up Planck query object
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

            ra = (-180. + 360.*np.random.random(shape)) * units.deg
            dec = (-90. + 180. * np.random.random(shape)) * units.deg
            c = coords.SkyCoord(ra, dec, frame='icrs')

            E = self._planck(c)

            np.testing.assert_equal(E.shape, shape)


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


if __name__ == '__main__':
    unittest.main()
