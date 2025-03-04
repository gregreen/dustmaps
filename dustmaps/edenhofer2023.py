#!/usr/bin/env python
#
# edenhofer2023.py
# Reads the Edenhofer2023 dust map, which is described in
# Edenhofer et al. (2023).
#
# Copyright (C) 2023  Gordian Edenhofer and Gregory M. Green
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

from __future__ import division, print_function

import os
import sys
from collections import namedtuple
from functools import partial

import astropy.units as units
import numpy as np

from .map_base import DustMap, ensure_flat_galactic
from .std_paths import data_dir

_DustSphere = namedtuple(
    "_DustSphere", (
        "data", "nside", "nest", "radii", "coo_bounds", "radii0", "data0",
        "units", "data_uncertainty", "data0_uncertainty"
    )
)

DATA_DIR_SUBDIR = "edenhofer_2023"


def _removeprefix(s, prefix):
    return s[len(prefix):] if s.startswith(prefix) else s[:]


def _removesuffix(s, suffix):
    return s[:-len(suffix)] if s.endswith(suffix) else s[:]


def _get_sphere(filepath):
    from astropy.io import fits

    nside = None
    nest = None
    radii = None
    coo_bounds = None
    radii0 = None
    rec_dust0 = None
    units = None
    rec_uncertainty = None
    rec0_uncertainty = None
    with fits.open(filepath, "readonly") as hdul:
        for hdu in hdul:
            nm = hdu.name.lower()
            if isinstance(hdu, fits.PrimaryHDU) and hdu.data is None:
                continue
            elif nm in ("primary", "mean", "samples", "healpix times distance"):
                dust_density = hdu.data
                nside = hdu.header["NSIDE"]
                nest = hdu.header["ORDERING"].lower().startswith("nest")
                units = hdu.header.get("CUNIT")
                if dust_density.shape[-1] != 12 * nside**2:
                    ve = "invalid shape of dust density {!r}"
                    raise ValueError(ve.format(dust_density.shape))
            elif isinstance(hdu, fits.BinTableHDU):
                if hdu.data.names == ['radial pixel centers']:
                    radii = hdu.data["radial pixel centers"]
                elif hdu.data.names == ['radial pixel boundaries']:
                    coo_bounds = hdu.data["radial pixel boundaries"]
                else:
                    ve = "unrecognized entry in BinTableHDU {!r}"
                    raise ValueError(ve.format(hdu.data.names))
            elif isinstance(hdu, fits.ImageHDU) and (
                nm.startswith("mean of integrated inner") or
                nm.startswith("integrated inner")
            ):
                prfx = "inner density integrated within"
                ctp = hdu.header.get("CTYPE")
                ctp = hdu.header.get("CTYPE1") if ctp is None else ctp
                if ctp.lower().startswith(prfx):
                    radii0 = _removeprefix(ctp.lower(), prfx)
                    if not radii0.endswith("pc"):
                        ve = "unrecognized units {!r}".format(radii0)
                        raise ValueError(ve)
                    radii0 = _removesuffix(radii0, "pc")
                    radii0 = float(radii0.strip())
                    rec_dust0 = hdu.data
                    if nest is None or nside is None:
                        ve = (
                            "main reconstruction needs to come before inner"
                            " in FITS file"
                        )
                        raise ValueError(ve)
                    if nest != hdu.header["ORDERING"].lower(
                    ).startswith("nest"):
                        raise ValueError("ordering mismatch")
                    if rec_dust0.shape[-1] != 12 * nside**2:
                        ve = "incompatible shape of dust density {!r}"
                        raise ValueError(ve.format(dust_density.shape))
            elif isinstance(hdu, fits.ImageHDU
                           ) and nm.startswith("std. of integrated inner"):
                prfx = "std. of inner density integrated within"
                ctp = hdu.header.get("CTYPE")
                ctp = hdu.header.get("CTYPE1") if ctp is None else ctp
                if hdu.header["CTYPE"].lower().startswith(prfx):
                    radii0_unc = _removeprefix(ctp.lower(), prfx)
                    if not radii0_unc.endswith("pc"):
                        ve = "unrecognized units {!r}".format(radii0)
                        raise ValueError(ve)
                    radii0_unc = _removesuffix(radii0_unc, "pc")
                    radii0_unc = float(radii0_unc.strip())
                    rec0_uncertainty = hdu.data
                    if nest is None or nside is None or radii0 is None:
                        ve = (
                            "main reconstruction and inner needs to come before"
                            " inner std. in FITS file"
                        )
                        raise ValueError(ve)
                    if radii0_unc != radii0:
                        raise ValueError("radii mismatch")
                    if nest != hdu.header["ORDERING"].lower(
                    ).startswith("nest"):
                        raise ValueError("ordering mismatch")
                    if rec0_uncertainty.shape[-1] != 12 * nside**2:
                        ve = "incompatible shape of dust density uncertainty {!r}"
                        raise ValueError(ve.format(rec0_uncertainty.shape))
            elif isinstance(hdu, fits.ImageHDU) and nm.lower() == "std.":
                rec_uncertainty = hdu.data
                if nest is None or nside is None:
                    ve = "main reconstruction needs to come before std."
                    raise ValueError(ve)
                if nest != hdu.header["ORDERING"].lower().startswith("nest"):
                    raise ValueError("ordering mismatch")
                if rec_uncertainty.shape[-1] != 12 * nside**2:
                    ve = "incompatible shape of dust density uncertainty {!r}"
                    raise ValueError(ve.format(rec_uncertainty.shape))
            else:
                raise ValueError("unrecognized HDU\n{!r}".format(hdu.header))

    return _DustSphere(
        dust_density,
        nside=nside,
        nest=nest,
        radii=radii,
        coo_bounds=coo_bounds,
        radii0=radii0,
        data0=rec_dust0,
        units=units,
        data_uncertainty=rec_uncertainty,
        data0_uncertainty=rec0_uncertainty,
    )


def _interp_hpxr2lbd(data, radii, nside, nest, lon, lat, dist):
    """Interpolate a 3D map of HEALPix times radii to arbitrary longitude,
    latitude, distance positions."""
    from healpy.pixelfunc import get_interp_weights

    assert lon.shape == lat.shape == dist.shape
    final_shape = data.shape[:-2] + lon.shape
    lon, lat, dist = lon.ravel(), lat.ravel(), dist.ravel()

    idx_pos, wgt_pos = get_interp_weights(
        nside, lon, lat, nest=nest, lonlat=True
    )
    idx_r = np.searchsorted(radii, dist)
    idx_l = idx_r - 1
    mask = (idx_l < 0) | (idx_r >= radii.size)
    idx_l, idx_r = idx_l.clip(0, radii.size - 1), idx_r.clip(0, radii.size - 1)
    wgt_los = np.abs(np.stack((radii[idx_r] - dist, radii[idx_l] - dist)))
    wgt_los /= wgt_los.sum(axis=0)

    assert wgt_pos.ndim == 2
    data = data[...,
                np.stack((idx_l, idx_r))[:, np.newaxis],
                np.stack((idx_pos, ) * 2)]
    # Manual POS and LOS interpolation
    data = (data * wgt_pos).sum(axis=-wgt_pos.ndim)
    data = (data * wgt_los).sum(axis=-wgt_los.ndim)
    return np.where(mask, np.nan, data).reshape(final_shape)


class Edenhofer2023Query(DustMap):
    """
    A class for querying the Edenhofer et al. (2023) 3D dust map.

    The map is in units of E of Zhang, Green, and Rix (2023) per parsec but can
    be translated to an extinction at any given wavelength by using the
    extinction curve published at https://doi.org/10.5281/zenodo.6674521 . For
    further details on the map, see the original paper.

    If you use this map in your work, please cite Edenhofer et al. (2023).

    The data is deposited at https://doi.org/10.5281/zenodo.8187943.
    """
    def __init__(
        self,
        map_fname=None,
        load_samples=False,
        integrated=False,
        flavor='main',
        seed=None,
    ):
        """
        Args:
            map_fname (Optional[str]): Filename of the map. Defaults
                to :obj:`None`, meaning that the default location
                is used.
            load_samples (Optional[bool]): Whether to load the posterior samples
                of the extinction density. The samples give more accurate
                interpolation resoluts and are required for standard deviations
                of integrated extinctions. Defaults to :obj:`False`.
            integrated (Optional[bool]): Whether to return integrated extinction
                density. In this case, the units of the map are E of Zhang,
                Green, and Rix (2023). The pre-processing for efficient access
                to the integrated extinction map can take a couple of minutes
                but subsequent queries will be fast. Defaults to :obj:`False`.
            flavor (Optional[str]): Flavor of the map to use. Must be in
                ('main', 'less_data_but_2kpc').
            seed (Optional[int]): A random seed, used when drawing random
                samples from the map. Set this seed in order to make the
                pseudo-random draws reproducible.
        """
        if not load_samples in (True, False):
            raise TypeError("`load_samples` must be bool")
        if not isinstance(flavor, str):
            raise TypeError("`flavor` must be str")
        if map_fname is None:
            if flavor.lower() == 'main':
                if load_samples is False:
                    fn = 'mean_and_std_healpix.fits'
                else:
                    fn = 'samples_healpix.fits'
            elif flavor.lower() == 'less_data_but_2kpc':
                if load_samples is False:
                    fn = 'validation_with_less_data_but_2kpc_mean_and_std_healpix.fits'
                else:
                    fn = 'validation_with_less_data_but_2kpc_samples_healpix.fits'
            else:
                raise ValueError("unrecognized flavor {!r}".format(flavor))
            map_fname = os.path.join(data_dir(), DATA_DIR_SUBDIR, fn)

            if not os.path.isfile(map_fname):
                from .dustexceptions import data_missing_message

                msg = data_missing_message(
                    "edenhofer2023", "Edenhofer et al. (2023)"
                )
                print(msg, file=sys.stderr)
                err = "{} does not exist".format(repr(map_fname))
                raise FileNotFoundError(err)

        self._flavor = flavor

        self._rec = _get_sphere(map_fname)
        self._has_samples = (self._rec.data.ndim == 3)
        if not self._has_samples and load_samples:
            raise ValueError("failed to load samples")
        if map_fname is None and self._rec.data0 is None:
            ve = "default map should have zeroth integrated layer but it doesn't"
            raise ValueError(ve)

        self._allowed_modes = ("mean",)
        # Replace the data in-place with its log respectively square because the
        # interpolation is done in log- respectively squared-space.
        if integrated is True:
            dvol = np.diff(self._rec.coo_bounds)
            np.multiply(dvol[:, np.newaxis], self._rec.data, out=self._rec.data)
            if self._rec.data0 is not None:
                self._rec.data[..., 0, :] += self._rec.data0
            msg = "Integrating extinction map (this might take a couple of minutes)..."
            print(msg, file=sys.stderr)
            np.cumsum(self._rec.data, axis=-2, out=self._rec.data)
            msg = "Optimizing map for querying (this might take a couple of seconds)..."
            print(msg, file=sys.stderr)
            np.log(self._rec.data, out=self._rec.data)
            if self._has_samples:
                self._allowed_modes += ("std", "samples", "random_sample")
        elif integrated is False:
            msg = "Optimizing map for querying (this might take a couple of seconds)..."
            print(msg, file=sys.stderr)
            np.log(self._rec.data, out=self._rec.data)
            if self._rec.data_uncertainty is not None:
                np.square(
                    self._rec.data_uncertainty, out=self._rec.data_uncertainty
                )
            self._allowed_modes += ("std",)
            if self._has_samples:
                self._allowed_modes += ("samples", "random_sample")
        else:
            te = "`integrated` must be bool; got {}".format(integrated)
            raise TypeError(te)
        self._integrated = integrated

        # If samples are loaded, initialize random number generator
        self._rng = None
        if self._has_samples:
            self._rng = np.random.default_rng(seed)

    @ensure_flat_galactic
    def query(self, coords, mode="mean"):
        """
        Returns the 3D dust extinction density or integrated extinction
        (depending on how the class was initialized) from Edenhofer et al.
        (2023) at the given coordinates. The map is in units of E of Zhang,
        Green, and Rix (2023) per parsec or if integrated simply in units of E.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): Coordinates at which
                to query the extinction. Must be 3D (i.e., include distance
                information).
            mode (str): Which mode to return. Allowed values are 'mean' (for the
                mean), 'std' (for the standard deviation), 'samples' (for all
                posterior samples), or 'random_sample' (for a single sample).
                Defaults to 'mean'.

        Notes:
            To query integrated extinction, set `integrated=True` during
            initizalization.

        Returns:
            Depending on how the class was initialized, either extinction
            density or extinction are returned. The output is either a numpy
            array or float, with the same shape as the input :obj:`coords`.
        """
        if not isinstance(mode, str):
            te = "`mode` must be str; got {}".format(type(mode))
            raise TypeError(te)
        mode = mode.strip().lower()
        if mode not in ("mean", "std", "samples", "random_sample"):
            ve = (
                "`mode` must be 'mean', 'std', 'samples', or 'random_sample'"
                "; got {!r}"
            )
            raise ValueError(ve.format(mode))
        if mode not in self._allowed_modes:
            ve = "`mode={!r}` requires samples but none are available"
            raise ValueError(ve.format(mode))

        interp = partial(
            _interp_hpxr2lbd,
            radii=self._rec.radii,
            nside=self._rec.nside,
            nest=self._rec.nest,
            lon=coords.galactic.l.deg,
            lat=coords.galactic.b.deg,
            dist=coords.galactic.distance.to("pc").value
        )
        if self._has_samples:
            if mode == "random_sample":
                # Same sample for all coordinates
                samp_idx = self._rng.choice(self.n_samples)
                res = interp(self._rec.data[samp_idx])
            else:
                res = interp(self._rec.data)
            res = np.exp(res, out=res)
            if mode == "mean":
                res = res.mean(axis=0)
            elif mode == "std":
                res = res.std(axis=0)
            elif mode == "random_sample":
                pass  # No sample-axis reduction necessary
            else:
                assert mode == "samples"
                # Swap sample and coordinate axes to be consistent with other
                # 3D dust maps. The output shape will be (..., sample).
                res = np.swapaxes(res, 0, -1)
        else:
            if mode == "mean":
                res = interp(self._rec.data)
                res = np.exp(res, out=res)
            else:
                assert mode == "std"
                res = interp(self._rec.data_uncertainty)
                res = np.sqrt(res, out=res)
        return res

    @property
    def distance_bounds(self):
        """
        Returns the distance bin edges that the map uses. The return type is
        :obj:`astropy.units.Quantity`, which stores unit-full quantities.
        """
        return self._rec.coo_bounds * units.pc

    @property
    def distances(self):
        """
        Returns the distance bin centers that the map uses. The return type is
        :obj:`astropy.units.Quantity`, which stores unit-full quantities.
        """
        return self._rec.radii * units.pc

    @property
    def integrated(self):
        """
        Returns ``True`` if the map contains integrated extinction, or
        ``False`` if it contains extinction density.
        """
        return self._integrated

    @property
    def n_samples(self):
        """
        Returns the number of samples stored in the map. If no samples have
        been loaded, then ``None`` is returned.
        """
        if self._has_samples:
            return self._rec.data.shape[0]
        else:
            return None

    @property
    def flavor(self):
        """
        Returns the flavor of the map. This is a string that is either
        ``"main"`` or ``"less_data_but_2kpc"``.
        """
        return self._flavor


def fetch(clobber=False, fetch_samples=False, fetch_2kpc=False):
    """
    Downloads the 3D dust map of Edenhofer et al. (2023).

    Args:
        clobber (Optional[bool]): If ``True``, any existing file will be
            overwritten, even if it appears to match. If ``False`` (the
            default), ``fetch()`` will attempt to determine if the dataset
            already exists. This determination is not 100\% robust against data
            corruption.
        fetch_samples (Optional[bool]): If ``True``, the samples will also be
            downloaded. If ``False`` (the default), only the mean and standard
            deviation will be downloaded taking up about 3.2 GB in size while
            the samples take up 19 GB.
        fetch_2kpc (Optional[bool]): If ``True``, the validation run using
            less data though which extends out to 2kpc in distance will also be
            downloaded.
    """
    from . import fetch_utils

    dest_dir = os.path.join(data_dir(), DATA_DIR_SUBDIR)

    file_spec = [
        ('mean_and_std_healpix.fits', '10c823a5fcf81b47b6e15530bcdf54dc')
    ]
    if fetch_2kpc:
        file_spec += [
            (
                'validation_with_less_data_but_2kpc_mean_and_std_healpix.fits',
                '0826b2f4cdfccad69bdc5e2d85bb4c54'
            )
        ]
    if fetch_samples:
        file_spec += [
            ('samples_healpix.fits', 'aa0aa435e013784fe18a5cb24e379b05')
        ]
        if fetch_2kpc:
            file_spec += [
                (
                    'validation_with_less_data_but_2kpc_samples_healpix.fits',
                    '1ef8ca17a68e724d21c357554951959c'
                )
            ]

    for fn, md5sum in file_spec:
        fname = os.path.join(dest_dir, fn)

        # Download from the server
        url = 'https://zenodo.org/record/8187943/files/{}'.format(fn)
        fetch_utils.download_and_verify(url, md5sum, fname, clobber=clobber)
