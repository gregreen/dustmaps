#!/usr/bin/env python
#
# sfd.py
# Reads the Schlegel, Finkbeiner & Davis (1998; SFD) dust reddening map.
#
# Copyright (C) 2016-2018  Gregory M. Green, Edward F. Schlafly
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

import astropy.wcs as wcs
import astropy.io.fits as fits
from scipy.ndimage import map_coordinates

from .std_paths import *
from .map_base import DustMap, WebDustMap, ensure_flat_galactic
from . import fetch_utils
from . import dustexceptions


class SFDBase(DustMap):
    """
    Queries maps stored in the same format as Schlegel, Finkbeiner & Davis (1998).
    """

    map_name = ''
    map_name_long = ''
    poles = ['ngp', 'sgp']

    def __init__(self, base_fname):
        """
        Args:
            base_fname (str): The map should be stored in two FITS files, named
                ``base_fname + '_' + X + '.fits'``, where ``X`` is ``'ngp'`` and
                ``'sgp'``.
        """
        self._data = {}

        for pole in self.poles:
            fname = '{}_{}.fits'.format(base_fname, pole)
            try:
                with fits.open(fname) as hdulist:
                    self._data[pole] = [hdulist[0].data, wcs.WCS(hdulist[0].header)]
            except IOError as error:
                print(dustexceptions.data_missing_message(self.map_name,
                                                          self.map_name_long))
                raise error

    @ensure_flat_galactic
    def query(self, coords, order=1):
        """
        Returns the map value at the specified location(s) on the sky.

        Args:
            coords (`astropy.coordinates.SkyCoord`): The coordinates to query.
            order (Optional[int]): Interpolation order to use. Defaults to `1`,
                for linear interpolation.

        Returns:
            A float array containing the map value at every input coordinate.
            The shape of the output will be the same as the shape of the
            coordinates stored by `coords`.
        """
        out = np.full(len(coords.l.deg), np.nan, dtype='f4')

        for pole in self.poles:
            m = (coords.b.deg >= 0) if pole == 'ngp' else (coords.b.deg < 0)

            if np.any(m):
                data, w = self._data[pole]
                x, y = w.wcs_world2pix(coords.l.deg[m], coords.b.deg[m], 0)
                out[m] = map_coordinates(data, [y, x], order=order, mode='nearest')

        return out


class SFDQuery(SFDBase):
    """
    Queries maps stored in the same format as Schlegel, Finkbeiner & Davis (1998) dust reddening map, providing results from the Schlegel, Finkbeiner & Davis (1998) dust map by default.
    """

    map_name = 'sfd'
    map_name_long = "SFD'98"

    def __init__(
        self, 
        map_dir=None,
        whichparent='SFD',
        which='dust'
    ):
        """
        Args:
            map_dir (Optional[str]): The directory containing the SFD map.
                Defaults to `None`, which means that `dustmaps` will look in its
                default data directory.
            whichparent (Optional[:obj:`str`]): The parent name of the SFD style data product to download.
                Should be either ``SFD`` (the default), ``Synch``, ``FINK``, or ``Haslam``.
            which (Optional[:obj:`str`]): The name of the SFD style data product to download.
                Should be either ``dust`` (the default), ``i100``, ``i60``, ``mask``, ``temp``, or ``xmap`` with whichparent ``SFD``. Should be ``Beta`` for whichparent ``Synch``. Should be ``Rmap`` for whichparent ``FINK``. Should be ``clean`` for whichparent ``Haslam``.
        """
        
        if map_dir is None:
            map_dir = os.path.join(data_dir(), 'sfd')

        if (whichparent == 'SFD') and (which == 'dust' or which == 'i100' or which == 'i60' or which == 'mask'):
            base_fname = os.path.join(map_dir, '{}_{}_4096'.format(whichparent,which))
        else:
            base_fname = os.path.join(map_dir, '{}_{}'.format(whichparent,which))
        
        super(SFDQuery, self).__init__(base_fname)
    
    def query(self, coords, order=1):
        """
        Returns value of SFD-style map defined by ``whichparent`` and ``which`` during the intialization of the `SFDQuery` object at the specified location(s) on the sky.
        
        Defaults return E(B-V) at the specified location(s) on the sky. See Table 6 of Schlafly & Finkbeiner (2011) for instructions on how to convert this quantity to extinction in various passbands.
        
        For whichparent ``SFD`` and which ``dust``, the map is the Schlegel, Finkbeiner & Davis (1998) dust reddening map (E(B-V)).
        
        For whichparent ``SFD`` and which ``i100``, the map is the Schlegel, Finkbeiner & Davis (1998) 100 micron intensity map (MJy/sr).
        
        For whichparent ``SFD`` and which ``i60``, the map is the Schlegel, Finkbeiner & Davis (1998) 60 micron intensity map (MJy/sr).
        
        For whichparent ``SFD`` and which ``mask``, the map is the Schlegel, Finkbeiner & Davis (1998) bit mask map.
            Bit 0, 1: The first two bits express (in binary) the number of HCONs (0, 1, 2, or 3)
            Bit 2: Asteroid removed
            Bit 3: Small no-data region replaced
            Bit 4: Source removed (any)
            Bit 5: No source removal
            Bit 6: Large objects - LMC, SMC or M31
            Bit 7: No IRAS data (excluded zone OR Saturn)
            
        For whichparent ``SFD`` and which ``temp``, the map is the Schlegel, Finkbeiner & Davis (1998) dust temperature map (K).
        
        For whichparent ``SFD`` and which ``xmap``, the map is the Schlegel, Finkbeiner & Davis (1998) X-factor map. This map contains a temperature correction factor derived from the 100mu/240mu ratio.  Multiply the 100mu map by this factor to obtain temperature-corrected emission in regions that are expected to have an unusual dust temperature. In some cases (e.g. high Galactic latitude) this factor is poorly constrained and should be used with caution. The mean value for this quantity in "normal" parts of the sky is 1.
        
        For whichparent ``FINK`` and which ``Rmap``, the map is the Finkbeiner-Davis-Schlegel (1999) DIRBE 100/240mu RATIO map. This map is described in "Extrapolation of Galactic Dust Emission at 100 Microns to CMBR Frequencies using FIRAS" by Finkbeiner, Davis, & Schlegel (1999). Please note that this 100/240mu R map differs from the R map described in Schlegel, Finkbeiner, & Davis, ApJ 500, 525 (1998).
        
        # check citation ##FIXME
        For whichparent ``Haslam`` and which ``clean``, the map is the Haslam et al. (1982) 408 MHz all-sky continuum survey, cleaned of bright sources (K).The map has bright point sources removed, and has been Fourier destriped using a method similar to that applied to the IRAS/ISSA data in Schlegel, Finkbeiner, & Davis 1998, Apj, 500, 525. Due to this reprocessing, the effective beam (PSF) of the map has increased from 0.85 deg to 1.0 deg. A CMB monopole (2.73K) has been subtracted from the map. 
        
        # check citation ##FIXME
        For whichparent ``Synch`` and which ``Beta``, the map is the Finkbeiner & Davis (1999) synchrotron spectral index map. This map is based on the 408 MHz Haslam et al. (1982) map, 1.42 GHz Reich & Reich (1986) map, and 2.326 GHz Jonas, Baart, & Nicolson (1998) map.

        Args:
            coords (`astropy.coordinates.SkyCoord`): The coordinates to query.
            order (Optional[int]): Interpolation order to use. Defaults to `1`,
                for linear interpolation.

        Returns:
            A float array containing the values of SFD-style map at every input coordinate. The shape of the output will be the same as the shape of the
            coordinates stored by `coords`.
        """
        return super(SFDQuery, self).query(coords, order=order)
    

class SFDWebQuery(WebDustMap):
    """
    Remote query over the web for the Schlegel, Finkbeiner & Davis (1998) dust
    map.

    This query object does not require a local version of the data, but rather
    an internet connection to contact the web API. The query functions have the
    same inputs and outputs as their counterparts in ``SFDQuery``, but
    are limited in keywords to the SFD dustmap.
    """

    def __init__(self, api_url=None):
        super(SFDWebQuery, self).__init__(
            api_url=api_url,
            map_name='sfd')


def fetch(whichparent='SFD',which='dust'):
    """
    Downloads maps in the format of the Schlegel, Finkbeiner & Davis (1998) dust map, placing them in the data directory for `dustmaps`. By default, it downloads the Schlegel, Finkbeiner & Davis (1998) dust map.
    
    Args:
        whichparent (Optional[:obj:`str`]): The parent name of the SFD style data product to download.
            Should be either ``SFD`` (the default), ``Synch``, ``FINK``, or ``Haslam``.
        which (Optional[:obj:`str`]): The name of the SFD style data product to download.
            Should be either ``dust`` (the default), ``i100``, ``i60``, ``mask``, ``temp``, or ``xmap`` with whichparent ``SFD``. Should be ``Beta`` for whichparent ``Synch``. Should be ``Rmap`` for whichparent ``FINK``. Should be ``clean`` for whichparent ``Haslam``.
    """
    doi = '10.7910/DVN/EWCNL5'

    for pole in ['ngp', 'sgp']:
        if (whichparent == 'SFD') and (which == 'dust' or which == 'i100' or which == 'i60' or which == 'mask'):
            requirements = {'filename': '{}_{}_4096_{}.fits'.format(whichparent,which,pole)}
            local_fname = os.path.join(
            data_dir(),
            'sfd', '{}_{}_4096_{}.fits'.format(whichparent,which,pole))
        else:
            requirements = {'filename': '{}_{}_{}.fits'.format(whichparent,which,pole)}
            local_fname = os.path.join(
            data_dir(),
            'sfd', '{}_{}_{}.fits'.format(whichparent,which,pole))

        print('Downloading {} {} data file to {}'.format(whichparent,which,local_fname))
        fetch_utils.dataverse_download_doi(
            doi,
            local_fname,
            file_requirements=requirements)
