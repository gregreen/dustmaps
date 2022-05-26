#!/usr/bin/env python
#
# gaia_tge.py
# Reads the Gaia TGE dust reddening maps.
#
# Copyright (C) 2022  Gregory M. Green
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

import os
import numpy as np
import healpy as hp
import astropy.io.fits as fits
import astropy.units as units

from .std_paths import *
from .healpix_map import HEALPixQuery
from . import fetch_utils
from . import dustexceptions


class GaiaTGEQuery(HEALPixQuery):
    """
    Queries the Gaia TGE (2022) dust map.
    """

    def __init__(self, map_fname=None, healpix_level='optimum'):
        """
        Args:
            map_fname (Optional[`str`]): Filename of the Gaia TGE map.
                Defaults to ``None``, meaning that the default location is
                used.
            healpix_level (Optional[`int` or `str`]): Which HEALPix
                level to load into the map. If "optimum" (the default), loads
                the optimum HEALPix level available at each location. If an
                `int`, instead loads the specified HEALPix level.
        """

        if map_fname is None:
            map_fname = os.path.join(
                data_dir(),
                'gaia_tge',
                'gaia_tge.fits.zip'
            )

        try:
            with fits.open(map_fname) as hdulist:
                d = hdulist[1].data[:]
        except IOError as error:
            print(dustexceptions.data_missing_message('gaia_tge',
                                                      'Gaia TGE'))
            raise error

        if isinstance(healpix_level, int):
            idx = (d['HealpixLevel'] == healpix_level)
            n_pix = np.count_nonzero(idx)
            if n_pix == 0:
                levels_avail = np.unique(d['HealpixLevel']).tolist()
                raise ValueError(
                    'Requested HEALPix level not stored in map. Available '
                    'levels: {}'.format(levels_avail)
                )
            hpx_sort_idx = np.argsort(d['HealpixId'][idx])
            idx = np.where(idx)[0]
            idx = idx[hpx_sort_idx]
        elif healpix_level == 'optimum':
            idx_opt = (d['OptimumHpxFlag'] == 'true') # Uses str, not bool
            # Upscale to highest HEALPix level
            hpx_level = d['HealpixLevel'][idx_opt]
            hpx_level_max = np.max(hpx_level)
            n_pix = 12 * 4**hpx_level_max
            # Index from original array to use in each pixel of final map
            idx = np.full(n_pix, -1, dtype='i8') # Empty pixel -> index=-1
            # Get the ring-ordered index of the optimal pixels
            idx_opt = np.where(idx_opt)[0]
            hpx_idx = d['HealpixId'][idx_opt]
            # Add pixels of each level to the map
            for level in np.unique(hpx_level):
                nside = 2**level
                idx_lvl = (hpx_level == level)
                # Get the nest-ordered index of optimal pixels at this level
                hpx_idx_ring = hpx_idx[idx_lvl]
                hpx_idx_nest = hp.pixelfunc.ring2nest(nside, hpx_idx_ring)
                # Fill in index (in orig arr) of these pixels
                mult_factor = 4**(hpx_level_max-level)
                hpx_idx_base = hpx_idx_nest*mult_factor
                for offset in range(mult_factor):
                    idx[hpx_idx_base+offset] = idx_opt[idx_lvl]
            # Reorder from nested to ring
            idx = hp.pixelfunc.reorder(idx, n2r=True)
        else:
            raise ValueError(
                '`healpix_level` must be either an integer or "optimum"'
            )

        bad_mask = (idx == -1)

        pix_val = d['A0Tge'][idx]
        pix_val[bad_mask] = np.nan

        dtype = [
            #('A0TgeChiSqr', 'f8'),
            #('A0TgeRange', 'S13'),
            ('A0TgeUncertainty', 'f8'),
            ('NbrTracersUsed', 'i8'),
            ('OptimumHpxFlag', 'bool'),
            #('R0TgeChiSqr', 'f8'),
            #('R0TgeRange', 'f8'),
            #('R0TgeUncertainty', 'f8'),
            #('SolutionId', 'i8'),
            #('SourceIdNearestTracer', 'i8'),
            #('Status', 'i8')
        ]
        flags = np.empty(n_pix, dtype=dtype)
        for key,dt in dtype:
            if key == 'OptimumHpxFlag':
                flags[key] = (d[key][idx] == 'true')
            else:
                flags[key] = d[key][idx]
            flags[key][bad_mask] = {'f8':np.nan, 'i8':-1, 'bool':False}[dt]

        super(GaiaTGEQuery, self).__init__(
            pix_val, False, 'icrs', flags=flags
        )

    def query(self, coords, **kwargs):
        """
        Returns a numpy array containing A0 at the specified
        location(s) on the sky. Optionally, returns a 2nd array containing
        flags at the same location(s).

        Args:
            coords (`astropy.coordinates.SkyCoord`): The coordinates to
                query.
            return_flags (Optional[`bool`]): If `True`, returns a 2nd array
                containing flags at each coordinate. Defaults to `False`.

        Returns:
            A numpy array containing A0 at the specified
            coordinates. The shape of the output is the same as the shape of
            the input coordinate array, ``coords``. If `return_flags` is
            `True`, a 2nd record array containing flags at each coordinate
            is also returned.
        """
        return super(GaiaTGEQuery, self).query(coords, **kwargs)


def fetch():
    """
    Downloads the Gaia Total Galactic Extinction (TGE) dust maps, placing
    it in the default ``dustmaps`` directory.
    """
    props = {
        'url': (
            'https://www.dropbox.com/s/264nhm6mnww04x8/TGE_sim.fits.zip'
            '?dl=1&file_subpath=%2FTGE_sim.fits'
        ),
        'md5': '9bc93c6786951bb0b7a646d080021e4c',
        'fname': 'gaia_tge.fits.zip'
    }
    fname = os.path.join(data_dir(), 'gaia_tge', props['fname'])
    fetch_utils.download_and_verify(props['url'], props['md5'], fname=fname)


def main():
    from astropy.coordinates import SkyCoord
    q = GaiaTGEQuery()
    c = SkyCoord([0., 180., 0.], [0., 0., 90.], frame='galactic', unit='deg')
    print(q(c))


if __name__ == '__main__':
    main()
