#!/usr/bin/env python
#
# healpix_map.py
# A set of HEALPix map classes.
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
import six

import numpy as np
import healpy as hp
import astropy.io.fits as fits

from .map_base import DustMap, coord2healpix


class HEALPixQuery(DustMap):
    """
    A class for querying HEALPix maps.
    """

    def __init__(self, pix_val, nest, coord_frame, flags=None):
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
        self._flags = flags
        if (flags is not None) and (flags.shape[0] != pix_val.shape[0]):
            raise ValueError((
                'The shape of `flags` ({}) must match the shape '
                'of `pix_val` ({}) along the first axis.'
            ).format(flags.shape, pix_val.shape))
        super(HEALPixQuery, self).__init__()

    def query(self, coords, return_flags=False):
        """
        Args:
            coords (`astropy.coordinates.SkyCoord`): The coordinates to query.
            return_flags ([Optional[:obj:`bool`]): If `True`, return flags at
                each pixel. Only possible if flags were provided during
                initialization.

        Returns:
            A float array of the value of the map at the given coordinates. The
            shape of the output is the same as the shape of the coordinates
            stored by `coords`. If `return_flags` is `True`, then a second
            array, containing flags at each pixel, is also returned.
        """
        pix_idx = coord2healpix(coords, self._frame,
                                self._nside, nest=self._nest)
        sel_pix = self._pix_val[pix_idx]

        if return_flags:
            if self._flags is None:
                raise ValueError(
                    '`return_flags` is True, but the class was initialized '
                    'without flags.'
                )
            return sel_pix, self._flags[pix_idx]

        return self._pix_val[pix_idx]


class HEALPixFITSQuery(HEALPixQuery):
    """
    A HEALPix map class that is initialized from a FITS file.
    """

    def __init__(self, fname, coord_frame, hdu=0, field=None,
                                           dtype='f8', scale=None, memmap=True):
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
            memmap (Optional[:obj:`bool`]): If ``True`` (the default) and
                ``fname`` is a filename, the FITS file is memory-mapped, so the
                map data is not loaded into memory. If ``False``, the data is
                read into memory, allowing the file to be closed, modified or
                deleted afterwards. Has no effect when ``fname`` is an
                already-open :obj:`HDUList` or an HDU object.
        """
        self._out_dtype = np.dtype(dtype)
        self._scale = scale
        close_file = False

        if isinstance(fname, six.string_types):
            close_file = True
            hdulist = fits.open(fname, memmap=memmap)
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
            # ``np.asarray`` turns the ``FITS_rec`` into a plain
            # (memory-mapped) ndarray view. Keeping the ``FITS_rec`` breaks
            # fancy-indexing with multi-dimensional pixel indices.
            pix_val = np.asarray(hdu.data[:]).ravel()
        else:
            pix_val = hdu.data[field]

        nest = hdu.header.get('ORDERING', 'NESTED').strip() == 'NESTED'

        if close_file:
            hdulist.close()

        super(HEALPixFITSQuery, self).__init__(pix_val, nest, coord_frame)

    def query(self, coords, return_flags=False):
        result = super(HEALPixFITSQuery, self).query(coords,
                                                     return_flags=return_flags)
        if return_flags:
            sel_pix, flags = result
        else:
            sel_pix = result
        sel_pix = sel_pix.astype(self._out_dtype)
        if self._scale is not None:
            # ``astype`` above already copied the data, so scaling in place is
            # safe. Structured dtypes (e.g. GNILC with error estimates) must be
            # scaled field-by-field, since ufuncs do not support record dtypes.
            names = sel_pix.dtype.names
            if names is None:
                sel_pix *= self._scale
            else:
                for n in names:
                    sel_pix[n] *= self._scale
        return (sel_pix, flags) if return_flags else sel_pix
