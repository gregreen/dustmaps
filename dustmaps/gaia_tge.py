#!/usr/bin/env python
#
# gaia_tge.py
# Reads the Gaia TGE dust reddening maps.
#
# Copyright (C) 2022  Gregory M. Green
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

import os
import numpy as np
import healpy as hp
from astropy.table import Table
import astropy.units as units

from .std_paths import *
from .healpix_map import HEALPixQuery
from . import fetch_utils
from . import dustexceptions


class GaiaTGEQuery(HEALPixQuery):
    """
    Queries the Gaia Total Galactic Extinction (Delchambre 2022) dust map,
    which contains estimates of monochromatic extinction, A0, in mags.
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
                'TotalGalacticExtinctionMap_001.csv.gz'
            )

        try:
            # Cannot use astropy ECSV reader, due to bug in processing
            # null values
            dtype = [
                ('solution_id', 'i8'),
                ('healpix_id', 'i8'),
                ('healpix_level', 'i1'),
                ('a0', 'f4'),
                ('a0_uncertainty', 'f4'),
                ('a0_min', 'f4'),
                ('a0_max', 'f4'),
                ('num_tracers_used', 'i4'),
                ('optimum_hpx_flag', '?'),
                ('status', 'i2')
            ]
            converters = {8: lambda x: x == '"True"'}
            d = np.genfromtxt(
                map_fname, comments='#', delimiter=',',
                encoding='utf-8', converters=converters,
                dtype=dtype
            )[1:]
        except IOError as error:
            print(dustexceptions.data_missing_message('gaia_tge',
                                                      'Gaia TGE'))
            raise error

        if isinstance(healpix_level, int):
            idx = (d['healpix_level'] == healpix_level)
            n_pix = np.count_nonzero(idx)
            if n_pix == 0:
                levels_avail = np.unique(d['healpix_level']).tolist()
                raise ValueError(
                    'Requested HEALPix level not stored in map. Available '
                    'levels: {}'.format(levels_avail)
                )
            hpx_sort_idx = np.argsort(d['healpix_id'][idx])
            idx = np.where(idx)[0]
            idx = idx[hpx_sort_idx]
        elif healpix_level == 'optimum':
            idx_opt = d['optimum_hpx_flag']
            # Upscale to highest HEALPix level
            hpx_level = d['healpix_level'][idx_opt]
            hpx_level_max = np.max(hpx_level)
            n_pix = 12 * 4**hpx_level_max
            # Index from original array to use in each pixel of final map
            idx = np.full(n_pix, -1, dtype='i8') # Empty pixel -> index=-1
            # Get the ring-ordered index of the optimal pixels
            idx_opt = np.where(idx_opt)[0]
            hpx_idx = d['healpix_id'][idx_opt]
            # Add pixels of each level to the map
            for level in np.unique(hpx_level):
                nside = 2**level
                idx_lvl = (hpx_level == level)
                # Get the nest-ordered index of optimal pixels at this level
                hpx_idx_nest = hpx_idx[idx_lvl]
                # Fill in index (in orig arr) of these pixels
                mult_factor = 4**(hpx_level_max-level)
                hpx_idx_base = hpx_idx_nest*mult_factor
                for offset in range(mult_factor):
                    idx[hpx_idx_base+offset] = idx_opt[idx_lvl]
        else:
            raise ValueError(
                '`healpix_level` must be either an integer or "optimum"'
            )

        bad_mask = (idx == -1)

        pix_val = d['a0'][idx]
        pix_val[bad_mask] = np.nan

        dtype = [
            ('a0_uncertainty', 'f4'),
            ('num_tracers_used', 'i4'),
            ('optimum_hpx_flag', 'bool')
        ]
        flags = np.empty(n_pix, dtype=dtype)
        for key,dt in dtype:
            flags[key] = d[key][idx]
            flags[key][bad_mask] = {'f4':np.nan, 'i4':-1, 'bool':False}[dt]

        super(GaiaTGEQuery, self).__init__(
            pix_val, True, 'icrs', flags=flags
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
            'http://cdn.gea.esac.esa.int/Gaia/gdr3/Astrophysical_parameters/'
            'total_galactic_extinction_map/TotalGalacticExtinctionMap_001.csv.gz'
        ),
        'md5': '5f6271869b7e60960a955f08ca11dc37',
        'fname': 'TotalGalacticExtinctionMap_001.csv.gz'
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
