#!/usr/bin/env python
#
# test_serializers.py
# Test code that serializes and deserializes numpy and Astropy objects.
#
# Copyright (C) 2020  Gregory M. Green
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
import json
import astropy.units as units
from astropy.coordinates import SkyCoord

from .. import json_serializers


class TestSerializers(unittest.TestCase):
    def test_numpy_readable(self):
        """
        Test serializing/deserializing an array using human-readable serializer.
        """
        x = np.random.random(size=(15,7))
        o = json_serializers.serialize_ndarray_readable(x)
        y = json_serializers.deserialize_ndarray(o)
        np.testing.assert_allclose(x, y, atol=1.e-5, rtol=1.e-5)

    def test_numpy_b64(self):
        """
        Test serializing/deserializing an array using b64 serializer.
        """
        x = np.random.random(size=(15,7))
        o = json_serializers.serialize_ndarray_b64(x)
        y = json_serializers.deserialize_ndarray(o)
        np.testing.assert_allclose(x, y, atol=1.e-5, rtol=1.e-5)

    def test_numpy_npy(self):
        """
        Test serializing/deserializing an array using npy serializer.
        """
        x = np.random.random(size=(15,7))
        d = json_serializers.serialize_ndarray_npy(x)
        y = json_serializers.deserialize_ndarray(d)
        np.testing.assert_allclose(x, y, atol=1.e-5, rtol=1.e-5)

    def test_skycoord(self):
        """
        Test serializing/deserializing SkyCoord objects.
        """
        lon = np.random.uniform(0., 360., 23) * units.deg
        lat = np.random.uniform(-90., 90., lon.size) * units.deg
        d = np.random.uniform(0.1, 10., lon.size) * units.kpc

        decoder = json_serializers.MultiJSONDecoder

        for mode in ('b64', 'readable', 'npy'):
            for frame in ('galactic', 'icrs'):
                encoder = json_serializers.get_encoder(ndarray_mode=mode)

                # Without distance
                c = SkyCoord(lon, lat, frame=frame)
                s = json.dumps(c, cls=encoder)
                c_dec = json.loads(s, cls=decoder)
                sep = c.separation(c_dec).to('rad').value
                np.testing.assert_allclose(sep, np.zeros_like(sep), atol=1.e-7, rtol=0.)

                # With distance
                c = SkyCoord(lon, lat, distance=d, frame=frame)
                s = json.dumps(c, cls=encoder)
                c_dec = json.loads(s, cls=decoder)
                sep = c.separation_3d(c_dec).to('kpc').value
                np.testing.assert_allclose(sep, np.zeros_like(sep), atol=1.e-7, rtol=0.)


if __name__ == '__main__':
    unittest.main()

