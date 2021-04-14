#!/usr/bin/env python
#
# healpix_map.py
# A set of HEALPix map classes.
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
import six

import numpy as np
import healpy as hp
import astropy.io.fits as fits

from .map_base import DustMap, coord2healpix


class HEALPixQuery(DustMap):
    """
    A class for querying HEALPix maps.
    """

    def __init__(self, pix_val, nest, coord_frame):
        """
        Args:
            pix_val (array): Value of the map in every pixel. The length of the
                array must be of the form `12 * nside**2`, where `nside` is a
                power of two.
            nest (bool): `True` if the map uses nested ordering. `False` if
                ring ordering is used.
            coord_frame (str): The coordinate system that the HEALPix map is in.
                Should be one of the frames supported by `astropy.coordinates`.
        """
        self._nside = hp.pixelfunc.npix2nside(len(pix_val))
        self._pix_val = pix_val
        self._nest = nest
        self._frame = coord_frame
        super(HEALPixQuery, self).__init__()

    def query(self, coords):
        """
        Args:
            coords (`astropy.coordinates.SkyCoord`): The coordinates to query.

        Returns:
            A float array of the value of the map at the given coordinates. The
            shape of the output is the same as the shape of the coordinates
            stored by `coords`.
        """
        pix_idx = coord2healpix(coords, self._frame,
                                self._nside, nest=self._nest)
        return self._pix_val[pix_idx]


class HEALPixFITSQuery(HEALPixQuery):
    """
    A HEALPix map class that is initialized from a FITS file.
    """

    def __init__(self, fname, coord_frame, hdu=0, field=None,
                                           dtype='f8', scale=None):
        """
        Args:
            fname (str, HDUList, TableHDU or BinTableHDU): The filename, HDUList
                or HDU from which the map should be loaded.
            coord_frame (str): The coordinate system in which the HEALPix map is
                defined. Must be a coordinate frame which ``astropy``
                understands.
            hdu (Optional[int or str]): Specifies which HDU to load the map
                from. Defaults to ``0``.
            field (Optional[int or str]): Specifies which field (column) to load
                the map from. Defaults to ``None``, meaning that ``hdu.data[:]``
                is used.
            dtype (Optional[str or type]): The data will be coerced to this
                datatype. Can be any type specification that numpy understands,
                including a structured datatype, if multiple fields are to be
                loaded. Defaults to ``'f8'``, for IEEE754 double precision.
            scale (Optional[:obj:`float`]): Scale factor to be multiplied into
                the data.
        """
        close_file = False

        if isinstance(fname, six.string_types):
            close_file = True
            hdulist = fits.open(fname)
            print(hdulist.info())
            hdu = hdulist[hdu]
        elif isinstance(fname, fits.HDUList):
            hdu = fname[hdu]
        elif (isinstance(fname, fits.TableHDU)
                or isinstance(fname, fits.BinTableHDU)):
            hdu = fname
        else:
            raise TypeError('`fname` must be a `str`, `HDUList`, `TableHDU` or '
                            '`BinTableHDU`.')

        if field is None:
            pix_val = np.array(hdu.data[:].ravel().astype(dtype))
        else:
            pix_val = np.array(hdu.data[field][:].ravel().astype(dtype))

        if scale is not None:
            names = pix_val.dtype.names
            if names is None:
                pix_val *= scale
            else:
                for n in names:
                    pix_val[n] *= scale

        nest = hdu.header.get('ORDERING', 'NESTED').strip() == 'NESTED'

        if close_file:
            hdulist.close()

        super(HEALPixFITSQuery, self).__init__(pix_val, nest, coord_frame)
