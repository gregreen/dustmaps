#!/usr/bin/env python
#
# map2d.py
# A generic interface to a 2D dust reddening map.
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

import astropy.coordinates as coordinates

def ensure_coord_type(f):
    def _wrapper_func(self, coords, **kwargs):
        if not isinstance(coords, coordinates.SkyCoord):
            raise TypeError('`coords` must be an astropy.coordinates.SkyCoord object.')
        return f(self, coords, **kwargs)
    return _wrapper_func

class DustMap(object):
    def __init__(self):
        pass

    @ensure_coord_type
    def __call__(self, coords, **kwargs):
        return self.query(coords, **kwargs)

    def query(self, coords, **kwargs):
        pass

    def query_gal(self, l, b, **kwargs):
        coords = coordinates.SkyCoord(l, b, frame='galactic', unit='deg')
        return self.query(coords, **kwargs)

    def query_equ(self, ra, dec, frame='icrs', **kwargs):
        valid_frames = ['icrs', 'fk4', 'fk5', 'fk4noeterms']
        if frame not in valid_frames:
            raise ValueError('`frame` not understood. Must be one of [{}].'.format(', '.join(valid_frames)))
        coords = coordinates.SkyCoord(ra, dec, frame='icrs', unit='deg')
        return self.query(coords, **kwargs)
