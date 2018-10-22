#!/usr/bin/env python
#
# unstructured_map.py
# Implements a class for querying dust maps with unstructured pixels. Sky
# coordinates are assigned to the nearest pixel.
#
# Copyright (C) 2016  Gregory M. Green
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

import numpy as np
import astropy.coordinates as coordinates
import astropy.units as units
from scipy.spatial import cKDTree as KDTree

from .map_base import DustMap


class UnstructuredDustMap(DustMap):
    """
    A class for querying dust maps with unstructured pixels. Sky coordinates are
    assigned to the nearest pixel.
    """

    def __init__(self, pix_coords, max_pix_scale, metric_p=2, frame=None):
        """
        Args:
            pix_coords (array-like :obj:`astropy.coordinates.SkyCoord`): The sky
                coordinates of the pixels.
            max_pix_scale (scalar :obj:`astropy.units.Quantity`): Maximum angular
                extent of a pixel. If no pixel is within this distance of a
                query point, NaN will be returned for that query point.
            metric_p (Optional[:obj:`float`]): The metric to use. Defaults to 2, which
                is the Euclidean metric. A value of 1 corresponds to the
                Manhattan metric, while a value approaching infinity yields the
                maximum component metric.
            frame (Optional[:obj:`str`]): The coordinate frame to use internally. Must
                be a frame understood by :obj:`astropy.coordinates.SkyCoord`.
                Defaults to :obj:`None`, meaning that the frame will be inferred
                from :obj:`pix_coords`.
        """
        self._n_pix = pix_coords.shape[0]
        self._metric_p = metric_p

        if frame is None:
            self._frame = pix_coords.frame
        else:
            self._frame = frame

        # Tesselate the space
        self._pix_vec = self._coords2vec(pix_coords)
        self._kd = KDTree(self._pix_vec)

        # Don't query more than this distance from any point
        self._max_pix_scale = max_pix_scale.to('rad').value

    def _coords2vec(self, coords):
        """
        Converts from sky coordinates to unit vectors. Before conversion to unit
        vectors, the coordiantes are transformed to the coordinate system used
        internally by the :obj:`UnstructuredDustMap`, which can be set during
        initialization of the class.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): Input coordinates to
                convert to unit vectors.

        Returns:
            Cartesian unit vectors corresponding to the input coordinates, after
            transforming to the coordinate system used internally by the
            :obj:`UnstructuredDustMap`.
        """

        # c = coords.transform_to(self._frame)
        # vec = np.empty((c.shape[0], 2), dtype='f8')
        # vec[:,0] = coordinates.Longitude(coords.l, wrap_angle=360.*units.deg).deg[:]
        # vec[:,1] = coords.b.deg[:]
        # return np.radians(vec)

        c = coords.transform_to(self._frame).represent_as('cartesian')
        vec_norm = np.sqrt(c.x**2 + c.y**2 + c.z**2)

        vec = np.empty((c.shape[0], 3), dtype=c.x.dtype)
        vec[:,0] = (c.x / vec_norm).value[:]
        vec[:,1] = (c.y / vec_norm).value[:]
        vec[:,2] = (c.z / vec_norm).value[:]

        return vec

    def _coords2idx(self, coords):
        """
        Converts from sky coordinates to pixel indices.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): Sky coordinates.

        Returns:
            Pixel indices of the coordinates, with the same shape as the input
            coordinates. Pixels which are outside the map are given an index
            equal to the number of pixels in the map.
        """

        x = self._coords2vec(coords)
        idx = self._kd.query(x, p=self._metric_p,
                             distance_upper_bound=self._max_pix_scale)
        return idx[1]
