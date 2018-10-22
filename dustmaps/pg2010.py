#!/usr/bin/env python
#
# pg2010.py
# Reads the Peek & Graves (2010) correction to the SFD'98 dust reddening map.
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

import os
import numpy as np

from .std_paths import *
from .map_base import DustMap, WebDustMap, ensure_flat_galactic
from . import fetch_utils
from . import dustexceptions
from .sfd import SFDBase


class PG2010Query(SFDBase):
    """
    Queries the Peek & Graves (2010) correction to the SFD'98 dust reddening map.
    """

    map_name = 'pg2010'
    map_name_long = "Peek & Graves (2010)"
    poles = ['ngp']

    def __init__(self, map_dir=None, component='dust'):
        """
        Args:
            map_dir (Optional[:obj:`str`]): The directory containing the SFD map.
                Defaults to :obj:`None`, which means that :obj:`dustmaps` will look in its
                default data directory.
            component (Optional[:obj:`str`]): :obj:`'dust'` (the default) to load the correction
                to E(B-V), or :obj:`'err'` to load the uncertainty in the correction.
        """
        
        if map_dir is None:
            map_dir = os.path.join(data_dir(), 'pg2010')
        
        if component not in ['dust', 'err']:
            raise ValueError('`component` must be either "dust" or "err"')

        base_fname = os.path.join(map_dir, 'PG_{}_4096'.format(component))
        
        super(PG2010Query, self).__init__(base_fname)
    
    def query(self, coords, order=1):
        """
        Returns the P&G (2010) correction to the SFD'98 E(B-V) at the specified
        location(s) on the sky. If component is 'err', then return the
        uncertainty in the correction.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.
            order (Optional[:obj:`int`]): Interpolation order to use. Defaults to ``1``,
                for linear interpolation.

        Returns:
            A float array containing the P&G (2010) correction (or its
            uncertainty) to SFD'98 at every input coordinate. The shape
            of the output will be the same as the shape of the coordinates
            stored by :obj:`coords`.
        """
        return super(PG2010Query, self).query(coords, order=order)
    

def fetch():
    """
    Downloads the Peek & Graves (2010) dust map, placing it in
    the data directory for :obj:`dustmap`.
    """
    doi = '10.7910/DVN/VBSI4A'
    
    for component in ['dust', 'err']:
        requirements = {'filename': 'PG_{}_4096_ngp.fits'.format(component)}
        local_fname = os.path.join(
            data_dir(),
            'pg2010', 'PG_{}_4096_ngp.fits'.format(component))
        print('Downloading P&G (2010) {} data file to {}'.format(
            component, local_fname))
        fetch_utils.dataverse_download_doi(
            doi,
            local_fname,
            file_requirements=requirements)
