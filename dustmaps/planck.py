#!/usr/bin/env python
#
# planck.py
# Reads the Planck Collaboration dust reddening maps.
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

import os
import numpy as np
import healpy as hp
import astropy.io.fits as fits
import astropy.units as units

from .std_paths import *
from .healpix_map import HEALPixFITSQuery
from . import fetch_utils
from . import dustexceptions


class PlanckQuery(HEALPixFITSQuery):
    """
    Queries the Planck Collaboration (2013) dust map.
    """

    def __init__(self, map_fname=None, component='extragalactic'):
        """
        Args:
            map_fname (Optional[:obj:`str`]): Filename of the Planck map.
                Defaults to ```None``, meaning that the default location is
                used.
            component (Optional[:obj:`str`]): Which measure of reddening to use. There
                are seven valid components. Three denote reddening measures:
                ``'extragalactic'``, ``'tau'`` and ``'radiance'``. Four refer
                to dust properties: ``'temperature'``, ``'beta'``,
                ``'err_temp'`` and ``'err_beta'``. Defaults to
                ``'extragalactic'``.
        """

        if map_fname is None:
            map_fname = os.path.join(
                data_dir(),
                'planck',
                'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
            )

        if component.lower() in ('ebv','extragalactic'):
            field = 'EBV'
            self._scale = 1.
        elif component.lower() in ('tau','tau353','tau_353','optical depth'):
            field = 'TAU353'
            self._scale = 1.49e4
        elif component.lower() in ('radiance','r'):
            field = 'RADIANCE'
            self._scale = 5.4e5
        elif component.lower() in ('temperature','temp','t'):
            field = 'TEMP'
            self._scale = units.Kelvin
        elif component.lower() in ('sigma_temp','sigma_t','err_temp','err_t'):
            field = 'ERR_TEMP'
            self._scale = units.Kelvin
        elif component.lower() in ('beta','b'):
            field = 'BETA'
            self._scale = 1.
        elif component.lower() in ('sigma_beta','sigma_b','err_beta','err_b'):
            field = 'ERR_BETA'
            self._scale = 1.
        else:
            raise ValueError((
                "Invalid `component`: '{}'\n"
                "Valid components for reddening are 'extragalactic', 'tau', "
                "and 'radiance'. Valid components for dust properties are "
                "'temperature', 'err_temp', 'beta' and 'err_beta'."
                ).format(component))

        try:
            with fits.open(map_fname) as hdulist:
                super(PlanckQuery, self).__init__(
                    hdulist, 'galactic',
                    hdu='COMP-MAP',
                    field=field,
                    dtype='f4',
                    scale=self._scale
                )
        except IOError as error:
            print(dustexceptions.data_missing_message('planck',
                                                      'Planck Collaboration'))
            raise error

    def query(self, coords, **kwargs):
        """
        Returns E(B-V) (or a different Planck dust inference, depending on how
        the class was intialized) at the specified location(s) on the sky.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to
                query.

        Returns:
            A float array of the selected Planck component, at the given
            coordinates. The shape of the output is the same as the shape of the
            input coordinate array, ``coords``. If extragalactic E(B-V), tau_353
            or radiance was chosen, then the output has units of magnitudes of
            E(B-V). If the selected Planck component is temperature (or
            temperature error), then an :obj:`astropy.Quantity` is returned,
            with units of Kelvin. If beta (or beta error) was chosen, then the
            output is unitless.
        """
        return super(PlanckQuery, self).query(coords, **kwargs)


class PlanckGNILCQuery(HEALPixFITSQuery):
    """
    Queries the Planck Collaboration (2016) GNILC dust map.
    """
    def __init__(self, map_fname=None, load_errors=False):
        """
        Args:
            map_fname (Optional[:obj:`str`]): Filename of the Planck map.
                Defaults to ```None``, meaning that the default location is
                used.
            load_errors (Optional[:obj:`str`]): If ``True``, then the error
                estimates will be loaded as well, and returned with any query.
                If ``False`` (the default), then queries will only return the
                the reddening estimate, without any error estimate.
        """

        if load_errors:
            self._has_errors = True
            field = None
            dtype = [('EBV','f4'), ('EBV_err','f4')]
        else:
            self._has_errors = False
            field = 'TAU353'
            dtype = 'f4'

        if map_fname is None:
            map_fname = os.path.join(
                data_dir(),
                'planck',
                'COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits'
            )

        try:
            with fits.open(map_fname) as hdulist:
                super(PlanckGNILCQuery, self).__init__(
                    hdulist, 'galactic',
                    hdu=1,
                    field=field,
                    dtype=dtype,
                    scale=1.49e4
                )
        except IOError as error:
            print(dustexceptions.data_missing_message('planck',
                                                      'Planck GNILC'))
            raise error

    def has_errors(self):
        """
        Returns ``True`` if the error estimates have been loaded.
        """
        return self._has_errors

    def query(self, coords, **kwargs):
        """
        Returns E(B-V) at the specified location(s) on the sky.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to
                query.

        Returns:
            If the error estimates have been loaded, then a structured array
            containing ``'EBV'`` and ``'EBV_err'`` is returned. Otherwise,
            returns a float array of E(B-V), at the given coordinates. The
            shape of the output is the same as the shape of the input
            coordinate array, ``coords``.
        """
        return super(PlanckGNILCQuery, self).query(coords, **kwargs)


def fetch(which='2013'):
    """
    Downloads the Planck dust maps, placing them in the default ``dustmaps``
    directory. There are two different Planck dust maps that can be
    downloaded: the Planck Collaboration (2013) map and the "GNILC" (Planck
    Collaboration 2016) map.
    
    Args:
        which (Optional[:obj:`str`]): The name of the dust map to download.
            Should be either ``2013`` (the default) or ``GNILC``.
    """
    planck_maps = {
        '2013': {
            'url': 'http://pla.esac.esa.int/pla/aio/product-action?MAP.MAP_ID=HFI_CompMap_ThermalDustModel_2048_R1.20.fits',
            'md5': '8d804f4e64e709f476a63f0dfed1fd11',
            'fname': 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
        },
        'GNILC': {
            'url': 'http://pla.esac.esa.int/pla/aio/product-action?MAP.MAP_ID=COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits',
            'md5': 'fc385c2ee5e82edf039cbca6e82d6872',
            'fname': 'COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits'
        }
    }
    if which not in planck_maps:
        raise ValueError(
            'Unknown map: "{}". Must be one of {}.'.format(
                which, tuple(planck_maps.keys())
            )
        )
    props = planck_maps[which]
    fname = os.path.join(data_dir(), 'planck', props['fname'])
    fetch_utils.download_and_verify(props['url'], props['md5'], fname=fname)


def main():
    from astropy.coordinates import SkyCoord
    q = PlanckQuery()
    c = SkyCoord([0., 180., 0.], [0., 0., 90.], frame='galactic', unit='deg')
    print(q(c))


if __name__ == '__main__':
    main()
