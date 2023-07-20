# Copyright(C) 2023 Gordian Edenhofer

from __future__ import division, print_function

import os
from collections import namedtuple
from functools import partial

import astropy.units as units
import numpy as np

from . import fetch_utils
from .map_base import DustMap, ensure_flat_galactic
from .std_paths import data_dir

DustSphere = namedtuple(
    "DustSphere", (
        "data", "nside", "nest", "radii", "coo_bounds", "radii0", "data0",
        "units", "data_uncertainty", "data0_uncertainty"
    )
)

DATA_DIR_SUBDIR = "edenhofer2023"


def _get_sphere(filepath):
    from astropy.io import fits

    radii0 = None
    rec_dust0 = None
    rec_uncertainty = None
    rec0_uncertainty = None
    with fits.open(filepath, "readonly") as hdul:
        for hdu in hdul:
            if isinstance(hdu, fits.PrimaryHDU):
                dust_density = hdu.data
                nside = hdu.header["NSIDE"]
                nest = hdu.header["ORDERING"].lower() == "nest"
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
                    raise ValueError(ve.format(hdu.names))
            elif isinstance(hdu, fits.ImageHDU) and hdu.name.lower(
            ).startswith("mean of integrated inner"):
                prfx = "inner density integrated within"
                ctp = hdu.header.get("CTYPE")
                ctp = hdu.header.get("CTYPE1") if ctp is None else ctp
                if hdu.header["CTYPE"].lower().startswith(prfx):
                    radii0 = ctp.lower().removeprefix(prfx)
                    if not radii0.endswith("pc"):
                        ve = "unrecognized units {!r}".format(radii0)
                        raise ValueError(ve)
                    radii0 = radii0.removesuffix("pc")
                    radii0 = float(radii0.strip())
                    rec_dust0 = hdu.data
                    if nest != (hdu.header["ORDERING"].lower() == "nest"):
                        raise ValueError("ordering mismatch")
                    if rec_dust0.shape[-1] != 12 * nside**2:
                        ve = "incompatible shape of dust density {!r}"
                        raise ValueError(ve.format(dust_density.shape))
            elif isinstance(hdu, fits.ImageHDU) and hdu.name.lower(
            ).startswith("std. of integrated inner"):
                prfx = "std. of inner density integrated within"
                ctp = hdu.header.get("CTYPE")
                ctp = hdu.header.get("CTYPE1") if ctp is None else ctp
                if hdu.header["CTYPE"].lower().startswith(prfx):
                    radii0_unc = ctp.lower().removeprefix(prfx)
                    if not radii0_unc.endswith("pc"):
                        ve = "unrecognized units {!r}".format(radii0)
                        raise ValueError(ve)
                    radii0_unc = radii0_unc.removesuffix("pc")
                    radii0_unc = float(radii0_unc.strip())
                    rec0_uncertainty = hdu.data
                    if radii0_unc != radii0:
                        raise ValueError("radii mismatch")
                    if nest != (hdu.header["ORDERING"].lower() == "nest"):
                        raise ValueError("ordering mismatch")
                    if rec0_uncertainty.shape[-1] != 12 * nside**2:
                        ve = "incompatible shape of dust density uncertainty {!r}"
                        raise ValueError(ve.format(rec0_uncertainty.shape))
            elif isinstance(hdu, fits.ImageHDU) and hdu.name.lower() == "std.":
                rec_uncertainty = hdu.data
                if nest != (hdu.header["ORDERING"].lower() == "nest"):
                    raise ValueError("ordering mismatch")
                if rec_uncertainty.shape[-1] != 12 * nside**2:
                    ve = "incompatible shape of dust density uncertainty {!r}"
                    raise ValueError(ve.format(rec_uncertainty.shape))
            else:
                raise ValueError("unrecognized HDU {!r}".format(hdu))

    return DustSphere(
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
    shp = lon.shape
    lon, lat, dist = lon.ravel(), lat.ravel(), dist.ravel()

    idx_pos, wgt_pos = get_interp_weights(
        nside, lon, lat, nest=nest, lonlat=True
    )
    idx_r = np.searchsorted(radii, dist)
    idx_l = idx_r - 1
    mask = (idx_l < 0) | (idx_r >= radii.size)
    idx_l, idx_r = idx_l.clip(0, radii.size - 1), idx_r.clip(0, radii.size - 1)
    wgt_los = np.abs(np.stack((radii[idx_l] - dist, radii[idx_r] - dist)))
    wgt_los /= wgt_los.sum(axis=0)

    assert wgt_pos.ndim == 2
    data = data[...,
                np.stack((idx_l, idx_r))[:, np.newaxis],
                np.stack((idx_pos, ) * 2)]
    # Manual POS and LOS interpolation
    data = (data * wgt_pos).sum(axis=-wgt_pos.ndim)
    data = (data * wgt_los).sum(axis=-wgt_los.ndim)
    return np.where(mask, np.nan, data).reshape(shp)


class Edenhofer2023Query(DustMap):
    """
    A class for querying the Edenhofer et al. (2023) 3D dust map.

    For details on how to use this map, see the original paper:
    TODO.
    The map is in units of E of Zhang, Green, and Rix (2023) but can be
    translated to an extinction at any given wavelength by using the extinction
    curve published at https://doi.org/10.5281/zenodo.6674521 .

    The data is deposited at Zenodo: TODO.
    """
    def __init__(
        self,
        map_fname=None,
        load_samples=True,
        integrated=False,
        version='pre_release',
    ):
        """
        Args:
            map_fname (Optional[str]): Filename of the map. Defaults
                to :obj:`None`, meaning that the default location
                is used.
            load_samples (Optional[bool]): Whether to load the posterior samples
                of the extinction density. The samples give more accurate
                interpolation resoluts and are required for standard deviations
                of integrated extinctions. Defaults to :obj:`True`.
            integrated (Optional[bool]): Whether to return integrated extinction
                density. The pre-processing for efficient access to the
                integrated extinction map can take a couple of minutes but
                subsequent queries will be fast. Defaults to :obj:`False`.
            version (Optional[str]): Version of the map to use. Must be in
                ('pre_release', 'pre_release_samples').
        """
        if map_fname is None:
            if version.lower() == 'pre_release':
                fn = 'mean_and_std_healpix.fits'
            elif version.lower() == 'pre_release_samples':
                fn = 'sampels_healpix.fits'
            else:
                raise ValueError("unrecognized version {!r}".format(version))
            map_fname = os.path.join(data_dir(), DATA_DIR_SUBDIR, fn)
        self._rec = _get_sphere(map_fname)

        # Replace the data in-place with its log respectively square because the
        # interpolation is done in log- respectively squared-space.
        if integrated is True:
            dvol = np.diff(self._rec.coo_bounds)
            np.multiply(dvol[:, np.newaxis], self._rec.data, out=self._rec.data)
            self._rec.data[..., 0, :] += self._rec.data0
            np.cumsum(self._rec.data, axis=-2, out=self._rec.data)
            np.log(self._rec.data, out=self._rec.data)
        elif integrated is False:
            np.log(self._rec.data, out=self._rec.data)
            np.log(self._rec.data0, out=self._rec.data0)
            if self._rec.data_uncertainty is not None:
                np.square(
                    self._rec.data_uncertainty, out=self._rec.data_uncertainty
                )
                np.square(
                    self._rec.data0_uncertainty,
                    out=self._rec.data0_uncertainty
                )
        else:
            te = "`integrated` must be bool; got {}".format(integrated)
            raise TypeError(te)

        self._integrated = integrated
        self._has_samples = (self._rec.data.ndim == 4)

    @ensure_flat_galactic
    def query(self, coords, component="mean"):
        """
        Returns the 3D dust extinction from Edenhofer et al. (2023) at the given
        coordinates. The map is in units of E of Zhang, Green, and Rix (2023).

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): Coordinates at which
                to query the extinction. Must be 3D (i.e., include distance
                information).
            component (str): Which component to return. Allowable values are
                'mean' (for the mean extinction density), 'std' (for the
                standard deviation of extinction density), and 'samples' (for
                the posterior samples of the extinction density). Defaults to
                'mean'.

        Notes:
            To query integrated extinction, set `integrated=True` during
            initizalization.

        Returns:
            The extinction density, in units of e-foldings / pc, as either a
            numpy array or float, with the same shape as the input
            :obj:`coords`.
        """
        if not isinstance(component, str):
            te = "`component` must be str; got {}".format(type(component))
            raise TypeError(te)
        component = component.strip().lower()
        if not component in ("mean", "std", "samples"):
            ve = "`component` must be 'mean', 'std', or 'samples'; got {!r}"
            raise ValueError(ve.format(component))
        if component == "std" and self._integrated and not self._has_samples:
            raise ValueError("need samples for std. of integrated density")
        if component == "samples" and not self._has_samples:
            raise ValueError("no samples available")

        interp = partial(
            _interp_hpxr2lbd,
            radii=self._rec.radii,
            nside=self._rec.nside,
            nest=self._rec.nest,
            lon=coords.galactic.l.deg,
            lat=coords.galactic.b.deg,
            dist=coords.galactic.distance.to("pc").value
        )
        if component == "samples" or (
            component == "mean" and not self._has_samples
        ):
            res = interp(self._rec.data)
            res = np.exp(res, out=res)
        elif component == "mean" and self._has_samples:
            res = interp(self._rec.data)
            res = np.exp(res, out=res)
            res = res.mean(axis=0)
        elif component == "std" and not self._has_samples:
            res = interp(self._rec.data_uncertainty)
            res = np.sqrt(res, out=res)
        elif component == "std" and self._has_samples:
            res = interp(self._rec.data)
            res = np.sqrt(res, out=res)
            res = res.std(axis=0)
        else:
            raise AssertionError("wait! how?!")
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


def fetch(clobber=False, fetch_samples=False):
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
    """
    dest_dir = os.path.join(data_dir(), DATA_DIR_SUBDIR)

    file_spec = [
        ('mean_and_std_healpix.fits', '8be2ac6343f236dcb3a12bb15de8d9f2')
    ]
    if fetch_samples:
        file_spec += [
            ('samples_healpix.fits', '309e5e5a6c5001ee5c675a6750f7daaa')
        ]

    for fn, md5sum in file_spec:
        fname = os.path.join(dest_dir, fn)

        # Check if the file already exists
        if (not clobber) and fetch_utils.check_md5sum(fname, md5sum):
            print(
                'File "{}" appears to exist already. Call '.format(fn) +
                '`fetch(clobber=True)` to force overwriting of existing ' +
                'file.'
            )
            continue

        # Download from the server
        url = 'https://faun.rc.fas.harvard.edu/gedenhofer/volatile/E+23_v0p1/{}'.format(
            fn
        )
        fetch_utils.download_and_verify(url, md5sum, fname)
