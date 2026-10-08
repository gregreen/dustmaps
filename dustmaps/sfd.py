#!/usr/bin/env python
#
# sfd.py
# Reads the Schlegel, Finkbeiner & Davis (1998; SFD) dust reddening map.
#
# Copyright (C) 2016-2018  Gregory M. Green, Edward F. Schlafly
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
import warnings
import numpy as np

import astropy.wcs as wcs
import astropy.io.fits as fits
from scipy.ndimage import map_coordinates

from .std_paths import *
from .map_base import DustMap, WebDustMap, ensure_flat_galactic
from . import fetch_utils
from . import dustexceptions

# The SFD-style data products that `dustmaps` is able to download and query.
# Maps ``map_variant`` --> ``component`` --> ``(base_fname, description)``,
# where ``base_fname`` is the name of the FITS files, without the
# ``_ngp.fits``/``_sgp.fits`` suffix, and ``description`` is a short,
# human-readable name for the map.
SFD_MAPS = {
    'SFD': {
        'dust': ('SFD_dust_4096', "SFD'98"),
        'i100': ('SFD_i100_4096', "SFD'98 100 micron intensity"),
        'i60': ('SFD_i60_4096', "SFD'98 60 micron intensity"),
        'mask': ('SFD_mask_4096', "SFD'98 bit mask"),
        'temp': ('SFD_temp', "SFD'98 temperature"),
        'xmap': ('SFD_xmap', "SFD'98 X-factor"),
    },
    'FINK': {
        'Rmap': ('FINK_Rmap', "Finkbeiner, Davis & Schlegel (1999) 100/240 micron R"),
    },
    'Haslam': {
        'clean': ('Haslam_clean', 'Haslam et al. (1982) cleaned 408 MHz'),
    },
    'Synch': {
        'Beta': ('Synch_Beta', 'Finkbeiner & Davis (1999) synchrotron'),
    },
}


def _lookup_map(map_variant, component):
    """
    Look up an SFD-style data product.

    Args:
        map_variant (str): The parent name of the data product (e.g., ``'SFD'``).
        component (str): The specific map within ``map_variant`` (e.g.,
            ``'dust'``).

    Returns:
        A tuple ``(base_fname, description)``. See `SFD_MAPS`.

    Raises:
        ValueError: If ``map_variant`` or ``component`` is not one of the
            data products listed in `SFD_MAPS`.
    """
    if map_variant not in SFD_MAPS:
        raise ValueError(
            "Unrecognized map_variant: '{}'. Valid options are: {}.".format(
                map_variant,
                ', '.join("'{}'".format(v) for v in sorted(SFD_MAPS))))

    if component not in SFD_MAPS[map_variant]:
        raise ValueError(
            "Unrecognized component: '{}' for map_variant '{}'. "
            "Valid options are: {}.".format(
                component,
                map_variant,
                ', '.join("'{}'".format(c) for c in sorted(SFD_MAPS[map_variant]))))

    return SFD_MAPS[map_variant][component]


class SFDBase(DustMap):
    """
    Queries maps stored in the same format as Schlegel, Finkbeiner & Davis
    (1998).
    """

    map_name = ''
    map_name_long = ''
    fetch_args = ''
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
                print(dustexceptions.data_missing_message(
                    self.map_name, self.map_name_long, self.fetch_args))
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
    Queries an SFD-style map, providing results from the Schlegel, Finkbeiner &
    Davis (1998) dust reddening map by default.

    The map that is queried is selected by the ``map_variant`` and
    ``component`` arguments to the constructor.
    """

    map_name = 'sfd'
    map_name_long = "SFD'98"

    def __init__(self, map_dir=None, map_variant='SFD', component='dust'):
        """
        Args:
            map_dir (Optional[str]): The directory containing the SFD map.
                Defaults to `None`, which means that `dustmaps` will look in its
                default data directory.
            map_variant (Optional[str]): The parent name of the SFD-style data
                product to query. Should be one of ``SFD`` (the default),
                ``Synch``, ``FINK`` or ``Haslam``.
            component (Optional[str]): The specific map to query. For
                ``map_variant='SFD'``, should be one of ``dust`` (the default),
                ``i100``, ``i60``, ``mask``, ``temp`` or ``xmap``. For
                ``map_variant='Synch'``, should be ``Beta``. For
                ``map_variant='FINK'``, should be ``Rmap``. For
                ``map_variant='Haslam'``, should be ``clean``.

        Raises:
            ValueError: If ``map_variant`` or ``component`` is not one of the
                data products listed above.
        """
        base_fname, name = _lookup_map(map_variant, component)

        if map_dir is None:
            map_dir = os.path.join(data_dir(), 'sfd')

        self.map_variant = map_variant
        self.component = component
        self.map_name_long = name

        # Used to generate an accurate error message if the data are missing.
        if (map_variant, component) == ('SFD', 'dust'):
            self.fetch_args = ''
        else:
            self.fetch_args = "map_variant='{}', component='{}'".format(
                map_variant, component)

        super(SFDQuery, self).__init__(os.path.join(map_dir, base_fname))

    def query(self, coords, order=None):
        """
        Returns the value of the SFD-style map (selected by ``map_variant`` and
        ``component`` when the `SFDQuery` object was created) at the specified
        location(s) on the sky.

        By default, this returns E(B-V). See Table 6 of Schlafly & Finkbeiner
        (2011) for instructions on how to convert this quantity to extinction in
        various passbands.

        For ``map_variant='SFD'``:

        ``component='dust'``
            The Schlegel, Finkbeiner & Davis (1998) dust reddening map (E(B-V)).

        ``component='i100'``
            The Schlegel, Finkbeiner & Davis (1998) 100 micron intensity map
            (MJy/sr).

        ``component='i60'``
            The Schlegel, Finkbeiner & Davis (1998) 60 micron intensity map
            (MJy/sr).

        ``component='mask'``
            The Schlegel, Finkbeiner & Davis (1998) bit mask map. The bits are:

            - Bit 0, 1: The first two bits express (in binary) the number of
              HCONs (0, 1, 2, or 3)
            - Bit 2: Asteroid removed
            - Bit 3: Small no-data region replaced
            - Bit 4: Source removed (any)
            - Bit 5: No source removal
            - Bit 6: Large objects - LMC, SMC or M31
            - Bit 7: No IRAS data (excluded zone OR Saturn)

        ``component='temp'``
            The Schlegel, Finkbeiner & Davis (1998) dust temperature map (K).

        ``component='xmap'``
            The Schlegel, Finkbeiner & Davis (1998) X-factor map. This map
            contains a temperature correction factor derived from the
            100mu/240mu ratio. Multiply the 100mu map by this factor to obtain
            temperature-corrected emission in regions that are expected to have
            an unusual dust temperature. In some cases (e.g. high Galactic
            latitude) this factor is poorly constrained and should be used with
            caution. The mean value for this quantity in "normal" parts of the
            sky is 1.

        For ``map_variant='FINK'``:

        ``component='Rmap'``
            The Finkbeiner-Davis-Schlegel (1999) DIRBE 100/240mu RATIO map. This
            map is described in "Extrapolation of Galactic Dust Emission at 100
            Microns to CMBR Frequencies using FIRAS" by Finkbeiner, Davis, &
            Schlegel (1999). Please note that this 100/240mu R map differs from
            the R map described in Schlegel, Finkbeiner, & Davis, ApJ 500, 525
            (1998).

        For ``map_variant='Haslam'``:

        ``component='clean'``
            An unpublished version of the Haslam et al. (1982) 408 MHz all-sky
            continuum survey, cleaned of bright sources (K). The map has bright
            point sources removed, and has been Fourier destriped using a method
            similar to that applied to the IRAS/ISSA data in Schlegel,
            Finkbeiner, & Davis 1998, Apj, 500, 525. Due to this reprocessing,
            the effective beam (PSF) of the map has increased from 0.85 deg to
            1.0 deg. A CMB monopole (2.73K) has been subtracted from the map.

        For ``map_variant='Synch'``:

        ``component='Beta'``
            An unpublished work by Finkbeiner & Davis (1999) to derive a
            synchrotron spectral index map. This map is based on the 408 MHz
            Haslam et al. (1982) map, 1.42 GHz Reich & Reich (1986) map, and
            2.326 GHz Jonas, Baart, & Nicolson (1998) map.

        Args:
            coords (`astropy.coordinates.SkyCoord`): The coordinates to query.
            order (Optional[int]): Interpolation order to use. Defaults to
                `None`, which means `0` (nearest neighbor) for the ``mask``
                component, and `1` (linear interpolation) for all other
                components. An `order` other than `0` is not meaningful for the
                ``mask`` component, and is ignored (with a warning).

        Returns:
            An array containing the value of the SFD-style map at every input
            coordinate. The shape of the output will be the same as the shape
            of the coordinates stored by `coords`. The ``mask`` component
            returns integers (the same 8-bit unsigned integers stored in the
            FITS files); all other components return floats.
        """
        if self.component == 'mask':
            # The mask stores a bit field, in which the individual bits are not
            # independent of one another, so interpolating it would produce
            # meaningless values. Force nearest-neighbor (order 0) sampling.
            if order is None:
                order = 0
            elif order != 0:
                warnings.warn(
                    "The SFD bit mask cannot be interpolated; ignoring "
                    "order={} and using order=0 (nearest neighbor).".format(
                        order))
                order = 0
        elif order is None:
            order = 1

        result = super(SFDQuery, self).query(coords, order=order)

        if self.component == 'mask':
            # The mask is stored as an unsigned 8-bit integer (BITPIX=8), but
            # the interpolation above returns floats. Cast the result back to
            # the integer type in which the mask is stored.
            result = result.astype(self._data[self.poles[0]][0].dtype)

        return result


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


def fetch(map_variant='SFD', component='dust'):
    """
    Downloads an SFD-style map, placing it in the data directory for
    `dustmaps`. By default, it downloads the Schlegel, Finkbeiner & Davis (1998)
    dust reddening map.

    Args:
        map_variant (Optional[str]): The parent name of the SFD-style data
            product to download. Should be one of ``SFD`` (the default),
            ``Synch``, ``FINK`` or ``Haslam``.
        component (Optional[str]): The specific map to download. For
            ``map_variant='SFD'``, should be one of ``dust`` (the default),
            ``i100``, ``i60``, ``mask``, ``temp`` or ``xmap``. For
            ``map_variant='Synch'``, should be ``Beta``. For
            ``map_variant='FINK'``, should be ``Rmap``. For
            ``map_variant='Haslam'``, should be ``clean``.

    Raises:
        ValueError: If ``map_variant`` or ``component`` is not one of the data
            products listed above.
    """
    base_fname, _ = _lookup_map(map_variant, component)
    doi = '10.7910/DVN/EWCNL5'

    for pole in SFDBase.poles:
        filename = '{}_{}.fits'.format(base_fname, pole)
        local_fname = os.path.join(data_dir(), 'sfd', filename)

        print('Downloading {} {} data file to {}'.format(
            map_variant, component, local_fname))

        fetch_utils.dataverse_download_doi(
            doi,
            local_fname,
            file_requirements={'filename': filename})
