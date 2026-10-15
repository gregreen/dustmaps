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
from astropy.io import fits

from .map_base import DustMap, ensure_flat_galactic
from .std_paths import data_dir
from .rng_utils import get_rng

_DustSphere = namedtuple(
    "_DustSphere", (
        "data", "nside", "nest", "radii", "coo_bounds", "radii0", "data0",
        "units", "data_uncertainty", "data0_uncertainty"
    )
)

# A healpix map that stays in a memory-mapped FITS file, read a few values at
# a time by `.Edenhofer2023Query`. This is what `_get_sphere()` gives for a
# published FITS file, except that the maps themselves are never read into
# memory, and that the cumulative sum over the radial axis is checkpointed.
_SphereFile = namedtuple(
    "_SphereFile", (
        "fname", "data_shape", "nside", "nest", "radii", "coo_bounds",
        "radii0", "units", "n_pix", "n_radii", "n_samples", "n_check",
        "checkpoint_step", "has_uncertainty"
    )
)

DATA_DIR_SUBDIR = "edenhofer_2023"

# A repacked map is written next to the published FITS file that it is made
# from, under this suffix.
REPACKED_SUFFIX = ".repacked.fits"

# Number of radial bins between the checkpoints of the integrated map. A query
# at bin `b` reads the checkpoint at the bin `b // CHECKPOINT_STEP * CHECKPOINT_STEP`
# and sums the density between the two, so this is also how many extra values a
# query reads: at the cost of 1/CHECKPOINT_STEP of the file size, integrating a
# map at query time costs a handful of values instead of a whole radial column.
CHECKPOINT_STEP = 32

# Pixels written at a time while repacking, and coordinates read at a time
# while querying. These only affect how much memory either uses.
REPACK_BLOCK_PIXELS = 4096
QUERY_BLOCK_COORDS = 8192

_FITS_BITPIX = {"f4": -32, "f8": -64}


_InterpWeights = namedtuple(
    "_InterpWeights",
    ("idx_l", "idx_r", "idx_pos", "wgt_pos", "wgt_los", "mask",
     "coord_shape")
)


def _removeprefix(s, prefix):
    return s[len(prefix):] if s.startswith(prefix) else s[:]


def _removesuffix(s, suffix):
    return s[:-len(suffix)] if s.endswith(suffix) else s[:]


def _get_sphere(filepath):
    nside = None
    nest = None
    radii = None
    coo_bounds = None
    radii0 = None
    rec_dust0 = None
    units = None
    rec_uncertainty = None
    rec0_uncertainty = None
    # The maps are left memory-mapped, never read: a query of a repacked map
    # reads blocks of them, and `fits_repack()` blocks of the published ones.
    with fits.open(filepath, "readonly", memmap=True) as hdul:
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


def _integral_block(density, dvol, data0, cp_bins):
    """
    The cumulative sum of a block of a map along the radial axis, at the
    checkpointed radii.

    This is the sequence that `.Edenhofer2023Query` runs over a whole map when
    it is initialized with `integrated=True`, applied to a block of pixels
    instead, in the same order and in the same precision: the density is
    multiplied and summed in the precision it is stored in, so that the
    checkpoints are bit for bit the values that the reference implementation
    would reach at those radii.
    """
    work = np.empty(density.shape, dtype=density.dtype)
    np.multiply(dvol[:, np.newaxis], density, out=work)
    if data0 is not None:
        work[..., 0, :] += data0
    np.cumsum(work, axis=-2, out=work)
    return work[..., cp_bins, :]


def _stored(arr):
    """
    A published map, ``(radii, pixels)`` or ``(samples, radii, pixels)``, in
    the layout that a repacked file stores: ``(pixels, radii)`` or
    ``(pixels, radii, samples)``. The values that a query wants for a pixel,
    and that a FITS image keeps next to each other, are its radial column,
    followed by its samples.
    """
    arr = np.moveaxis(arr, -1, 0)
    if arr.ndim == 3:
        arr = np.moveaxis(arr, 1, 2)
    return arr


def _published(arr):
    """A stored map, back in the layout that the published files use."""
    if arr.ndim == 2:
        return arr.T
    return np.swapaxes(arr, 0, 2)


def _repacked_fname(fits_fname):
    """The name of the repacked form of a published FITS file."""
    return os.path.splitext(fits_fname)[0] + REPACKED_SUFFIX


def _fits_dtype(arr):
    """The FITS name of the dtype that a published array is stored in."""
    itemsize = np.asarray(arr).dtype.itemsize
    if itemsize == 4:
        return "f4"
    if itemsize == 8:
        return "f8"
    raise ValueError("unrecognized dtype {!r}".format(np.asarray(arr).dtype))


def _primary_header(shape, dtype, keywords=()):
    """
    The header of a primary HDU with an image of the given numpy shape.
    `shape` is the shape of the array as it is stored, whose last axis varies
    fastest, so the FITS axes come out reversed, as FITS requires.
    """
    header = fits.Header()
    header["SIMPLE"] = (True, "conforms to FITS standard")
    header["BITPIX"] = (_FITS_BITPIX[dtype], "number of bits per data pixel")
    header["NAXIS"] = (len(shape), "number of data axes")
    for axis, size in enumerate(reversed(shape), start=1):
        header["NAXIS{}".format(axis)] = (size, "length of data axis")
    header["EXTEND"] = (True, "FITS dataset may contain extensions")
    for key, value, comment in keywords:
        header[key] = (value, comment)
    return header


def _image_header(shape, dtype, name=None, keywords=()):
    """The header of an image HDU of the given numpy shape, no data."""
    header = fits.Header()
    header["XTENSION"] = ("IMAGE", "Image extension")
    header["BITPIX"] = (_FITS_BITPIX[dtype], "number of bits per data pixel")
    header["NAXIS"] = (len(shape), "number of data axes")
    for axis, size in enumerate(reversed(shape), start=1):
        header["NAXIS{}".format(axis)] = (size, "length of data axis")
    header["PCOUNT"] = (0, "required keyword; must = 0")
    header["GCOUNT"] = (1, "required keyword; must = 1")
    if name is not None:
        header["EXTNAME"] = (name, "name of this HDU")
    for key, value, comment in keywords:
        header[key] = (value, comment)
    return header


def _write_hdu(fh, header, blocks, dtype):
    """
    Write one image HDU: its header, the blocks of its data in the order they
    are stored, and the padding that ends it.

    `blocks` is an iterator, and the blocks are written as they are given, so
    that a map never has to be in memory as a whole.
    """
    fh.write(header.tostring().encode("ascii"))
    n_bytes = 0
    for block in blocks:
        block = np.ascontiguousarray(block, dtype=">" + dtype)
        n_bytes += block.nbytes
        fh.write(block.tobytes())
    fh.write(b"\0" * (-n_bytes % 2880))


def _published_blocks(arr, n_pix, block_pixels):
    """The stored blocks of a published map, in pixel order."""
    for start in range(0, n_pix, block_pixels):
        stop = min(start + block_pixels, n_pix)
        yield _stored(arr[..., start:stop])


def fits_repack(fits_in, fits_out, checkpoint_step=CHECKPOINT_STEP,
                block_pixels=REPACK_BLOCK_PIXELS):
    """
    Repack an Edenhofer et al. (2023) FITS file into the layout that
    `.Edenhofer2023Query` reads through a memory map.

    The maps are stored pixel-major: the radial column of a pixel, and its
    samples, are contiguous, so that the values a query interpolates between
    lie next to each other in the file, in the pages that the query reads.
    Alongside the density, the cumulative sum over the radial axis is stored at
    every `checkpoint_step` radii. That is what lets a query integrate on the
    fly: to integrate at a radial bin, read the checkpoint below it and add the
    density between the two, a handful of values rather than a whole radial
    column, at a cost of :math:`1 / \\text{checkpoint_step}` of the file size.

    The density is stored as published, neither logged nor integrated, because
    a query takes the log of the two values that each coordinate interpolates
    between, which is cheaper than storing them logged, and no less accurate.
    The checkpoints are computed exactly as `.Edenhofer2023Query` computes them
    for a map held in memory, so that the two agree bit for bit.

    The maps are written a block of pixels at a time, so a file of samples can
    be repacked on a machine that could not hold it; only the checkpoints,
    which are 1/checkpoint_step of the map, are held in memory.

    Args:
        fits_in (:obj:`str`): The published FITS file.
        fits_out (:obj:`str`): Filename to write.
        checkpoint_step (Optional[:obj:`int`]): Radii between checkpoints of
            the integrated map. Defaults to :obj:`CHECKPOINT_STEP`.
        block_pixels (Optional[:obj:`int`]): Pixels held in memory at a time
            while repacking. Defaults to :obj:`REPACK_BLOCK_PIXELS`.
    """
    rec = _get_sphere(fits_in)
    n_pix = 12 * rec.nside**2
    n_radii = int(np.asarray(rec.radii).size)
    n_samples = int(rec.data.shape[0]) if rec.data.ndim == 3 else None
    cp_bins = np.arange(0, n_radii, checkpoint_step)
    dvol = np.diff(rec.coo_bounds)

    if rec.data.shape[-1] != n_pix:
        ve = "invalid shape of map {!r}"
        raise ValueError(ve.format(rec.data.shape))

    samples = () if n_samples is None else (n_samples,)
    data_shape = (n_pix, n_radii) + samples
    check_shape = (n_pix, cp_bins.size) + samples
    units = "" if rec.units is None else str(rec.units)
    keywords = (
        ("REPACKED", True, "repacked for memory-mapped reading"),
        ("NSIDE", int(rec.nside), "HEALPix resolution of the map"),
        ("ORDERING", "NESTED" if rec.nest else "RING", "HEALPix ordering"),
        ("CUNIT", units, "units of the map"),
        ("RADII0", (np.nan if rec.radii0 is None else float(rec.radii0)),
         "radius of the inner map, in pc"),
        ("CPSTEP", int(checkpoint_step), "radii between checkpoints"),
        ("SRCFILE", os.path.basename(fits_in), "the map this was made from"),
        ("SRCSIZE", os.path.getsize(fits_in), "size of that file"),
        ("HASINNER", rec.data0 is not None,
         "the published map carried an inner map"),
    )
    bunit = (("BUNIT", units, "units of the map"),)
    cpunit = (("BUNIT", units + " pc" if units else "", "integrated units"),)

    # The checkpoints are small enough to hold, and holding them means reading
    # the published map once rather than once per HDU
    checkpoints = np.empty(check_shape, dtype="f4")
    for start in range(0, n_pix, block_pixels):
        stop = min(start + block_pixels, n_pix)
        data0 = None
        if rec.data0 is not None:
            data0 = np.asarray(rec.data0)[..., start:stop]
        block = _integral_block(
            rec.data[..., start:stop], dvol, data0, cp_bins
        )
        checkpoints[start:stop] = _stored(block)

    with open(fits_out, "wb") as fh:
        _write_hdu(fh, _primary_header((), "f4", keywords), (), "f4")
        _write_hdu(fh, _image_header(data_shape, "f4", "DATA", bunit),
                   _published_blocks(rec.data, n_pix, block_pixels), "f4")
        if rec.data_uncertainty is not None:
            _write_hdu(
                fh, _image_header(data_shape, "f4", "STD", bunit),
                _published_blocks(rec.data_uncertainty, n_pix, block_pixels),
                "f4"
            )
        if rec.data0 is not None:
            # Kept because a repacked map can also be read into memory and
            # integrated there, which needs the inner map
            inner = np.asarray(rec.data0)
            _write_hdu(fh, _image_header(inner.shape, _fits_dtype(inner),
                                         "DATA0", bunit),
                       (inner,), _fits_dtype(inner))
        _write_hdu(fh, _image_header(check_shape, "f4", "CHECKPOINTS", cpunit),
                   (checkpoints,), "f4")
        # The radial tables are stored as the published file stores them, so
        # that a query selects the same radial bins either way
        for name, arr in (("RADII", rec.radii), ("COO_BOUNDS", rec.coo_bounds)):
            arr = np.asarray(arr)
            _write_hdu(fh, _image_header((arr.size,), _fits_dtype(arr), name),
                       (arr,), _fits_dtype(arr))


def fits_is_repacked(fname):
    """Whether the given file is a map written by :obj:`fits_repack`."""
    try:
        with fits.open(fname, "readonly", memmap=False) as hdul:
            return bool(hdul[0].header.get("REPACKED", False))
    except (OSError, ValueError, KeyError, IndexError):
        return False


def _get_sphere_repacked(fname):
    """
    The description of a repacked map, as :obj:`_get_sphere` gives one for a
    published FITS file. The maps themselves are not read: they stay in the
    file, where a query maps them into memory.
    """
    with fits.open(fname, "readonly", memmap=True) as hdul:
        header = hdul[0].header
        data_shape = tuple(int(s) for s in hdul["DATA"].shape)
        units = str(header.get("CUNIT", ""))
        radii0 = float(header.get("RADII0", np.nan))
        return _SphereFile(
            fname=fname,
            data_shape=data_shape,
            nside=int(header["NSIDE"]),
            nest=str(header.get("ORDERING", "")).lower().startswith("nest"),
            radii=np.asarray(hdul["RADII"].data),
            coo_bounds=np.asarray(hdul["COO_BOUNDS"].data),
            radii0=None if not np.isfinite(radii0) else radii0,
            units=None if units == "" else units,
            n_pix=data_shape[0],
            n_radii=data_shape[1],
            n_samples=int(data_shape[2]) if len(data_shape) > 2 else None,
            n_check=int(hdul["CHECKPOINTS"].shape[1]),
            checkpoint_step=int(header.get("CPSTEP", CHECKPOINT_STEP)),
            has_uncertainty="STD" in hdul,
        )


def _get_sphere_repacked_data(fname):
    """
    The maps of a repacked file, in the layout that the published files use,
    for a query that wants them read into memory rather than mapped.
    :obj:`_get_sphere` gives the same thing for a published FITS file.
    """
    with fits.open(fname, "readonly", memmap=True) as hdul:
        header = hdul[0].header
        data = np.ascontiguousarray(_published(hdul["DATA"].data))
        uncertainty = None
        if "STD" in hdul:
            uncertainty = np.ascontiguousarray(_published(hdul["STD"].data))
        has_inner = header.get("HASINNER", None)
        if has_inner is None or (has_inner and "DATA0" not in hdul):
            # Without the inner map, integrating a map read into memory would
            # quietly leave out what the reference adds to its first bin
            ve = ("{} was repacked before the inner map was kept: repack it"
                  " again.")
            raise ValueError(ve.format(fname))
        units = str(header.get("CUNIT", ""))
        radii0 = float(header.get("RADII0", np.nan))
        return _DustSphere(
            data,
            nside=int(header["NSIDE"]),
            nest=str(header.get("ORDERING", "")).lower().startswith("nest"),
            radii=np.asarray(hdul["RADII"].data),
            coo_bounds=np.asarray(hdul["COO_BOUNDS"].data),
            radii0=None if not np.isfinite(radii0) else radii0,
            data0=(np.array(hdul["DATA0"].data) if has_inner else None),
            units=None if units == "" else units,
            data_uncertainty=uncertainty,
            data0_uncertainty=None,
        )


def repack_and_remove(fits_fname):
    """
    Repack a FITS file into the memory-mapped layout that
    `.Edenhofer2023Query` reads, and remove the FITS file: it holds the same
    map, and the repacked one is all that a query needs.

    The repacked file is written under a temporary name and renamed once it is
    complete, so a repacked file that is there is a finished one, and the FITS
    file is removed only once that has happened.
    """
    repacked_fname = _repacked_fname(fits_fname)
    part_fname = repacked_fname + ".part"
    try:
        fits_repack(fits_fname, part_fname)
        os.replace(part_fname, repacked_fname)
    except BaseException:
        if os.path.isfile(part_fname):
            os.remove(part_fname)
        raise
    os.remove(fits_fname)


def _gather(data, pix, bins, samples=None):
    """
    The stored values at the four pixels and the two radial bins (or
    checkpoints) that each coordinate falls between, in the shape that
    `.Edenhofer2023Query` interpolates: ``(samples, 2, 4, coordinates)``, or
    ``(2, 4, coordinates)`` for a map without samples.

    `data` is a memory-mapped map, laid out as a repacked file stores it, so
    the values for a pixel come from one place in the file. Fancy indexing
    copies what it reads, which is also what keeps a query from ever writing
    to the map.
    """
    if samples is None:
        return data[pix[np.newaxis], bins[:, np.newaxis]]
    return data[
        pix[np.newaxis, np.newaxis],
        bins[np.newaxis, :, np.newaxis],
        samples[:, np.newaxis, np.newaxis, np.newaxis],
    ]


def _gather_integral(hdul, pix, bins, samples, dvol, checkpoint_step, n_radii):
    """
    The cumulative sum of the density at each coordinate's two bracketing bins,
    assembled from the stored checkpoints.

    The radii between a checkpoint and the queried bin are read once, as part
    of the block of radii of one pixel, and their terms are then added to the
    checkpoint one radius at a time, in the order that the reference
    implementation accumulates them, so that the result is the same.
    """
    radii = (bins // checkpoint_step) * checkpoint_step
    offsets = bins - radii
    term_radii = np.minimum(
        radii[:, :, np.newaxis] + np.arange(1, checkpoint_step + 1),
        n_radii - 1
    )

    value = _gather(
        hdul["CHECKPOINTS"].data, pix, bins // checkpoint_step, samples
    )

    # The terms, in the shape of the checkpoints that they are added to: the
    # samples (if any), the two bins, the four pixels, the coordinates, and
    # then the radii between the checkpoint and the queried bin
    data = hdul["DATA"].data
    scale = dvol[term_radii]
    if samples is None:
        terms = data[pix[np.newaxis, :, :, np.newaxis],
                     term_radii[:, np.newaxis, :, :]]
        terms *= scale[:, np.newaxis]
    else:
        terms = data[
            pix[np.newaxis, np.newaxis, :, :, np.newaxis],
            term_radii[np.newaxis, :, np.newaxis, :, :],
            samples[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis],
        ]
        terms *= scale[np.newaxis, :, np.newaxis]

    for term in range(1, checkpoint_step + 1):
        keep = term <= offsets
        if not keep.any():
            break
        if samples is None:
            keep = keep[:, np.newaxis, :]
        else:
            keep = keep[np.newaxis, :, np.newaxis, :]
        value = value + np.where(keep, terms[..., term - 1], np.float32(0.))
    return value


def _interp_weights(radii, nside, nest, lon, lat, dist):
    """
    The pixels, radial bins and weights that each coordinate interpolates
    between, shared by :obj:`_interp_hpxr2lbd` and the memory-mapped query.
    """
    from healpy.pixelfunc import get_interp_weights

    assert lon.shape == lat.shape == dist.shape
    coord_shape = lon.shape
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
    return _InterpWeights(idx_l, idx_r, idx_pos, wgt_pos, wgt_los, mask,
                          coord_shape)


def _interp_reduce(block, weights, final_shape):
    """
    The weighted mean of the four pixels and two radial bins that each
    coordinate falls between.

    `block` holds those values, with the two radial bins and the four pixels on
    its last two axes, and `weights` those of the coordinates that it holds.
    """
    block = (block * weights.wgt_pos).sum(axis=-weights.wgt_pos.ndim)
    block = (block * weights.wgt_los).sum(axis=-weights.wgt_los.ndim)
    return np.where(weights.mask, np.nan, block).reshape(final_shape)


def _slice_weights(weights, start, stop):
    """The weights of the coordinates in `[start, stop)`."""
    return _InterpWeights(
        weights.idx_l[start:stop], weights.idx_r[start:stop],
        weights.idx_pos[:, start:stop], weights.wgt_pos[:, start:stop],
        weights.wgt_los[:, start:stop], weights.mask[start:stop],
        weights.coord_shape
    )


def _interp_hpxr2lbd(data, radii, nside, nest, lon, lat, dist):
    """Interpolate a 3D map of HEALPix times radii to arbitrary longitude,
    latitude, distance positions."""
    lead_shape = data.shape[:-2]
    weights = _interp_weights(radii, nside, nest, lon, lat, dist)
    data = data[...,
                np.stack((weights.idx_l, weights.idx_r))[:, np.newaxis],
                np.stack((weights.idx_pos, ) * 2)]
    # Manual POS and LOS interpolation
    return _interp_reduce(data, weights, lead_shape + weights.coord_shape)


class Edenhofer2023Query(DustMap):
    """
    A class for querying the Edenhofer et al. (2023) 3D dust map.

    The map is in units of E of Zhang, Green, and Rix (2023) per parsec but can
    be translated to an extinction at any given wavelength by using the
    extinction curve published at https://doi.org/10.5281/zenodo.6674521 . For
    further details on the map, see the original paper.

    If you use this map in your work, please cite Edenhofer et al. (2023).

    The data is deposited at https://doi.org/10.5281/zenodo.8187943.

    By default the map is read from a repacked copy of the FITS file, which
    `fetch()` writes and which a query maps into memory a few values at a
    time; see the `memmap` argument of `.Edenhofer2023Query.__init__`.
    """
    def __init__(
        self,
        map_fname=None,
        load_samples=False,
        integrated=False,
        flavor='main',
        memmap=True,
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
                Green, and Rix (2023). A published FITS map is pre-processed
                for this, which can take a couple of minutes but makes
                subsequent queries fast; a repacked one integrates as it is
                queried. Defaults to :obj:`False`.
            flavor (Optional[str]): Flavor of the map to use. Must be in
                ('main', 'less_data_but_2kpc').
            memmap (Optional[bool]): Whether to read the repacked, memory-
                mapped form of the map, which `dustmaps.edenhofer2023.fetch()`
                writes next to the FITS file. Nothing is loaded or transformed,
                so initialization is immediate and the map never occupies
                memory, and a query can return either density or integrated
                extinction. Defaults to :obj:`True`; with :obj:`False`, the map
                is read into memory and transformed in place instead, which
                costs a copy of the whole map in memory but spares every query
                the reads. Either way the map comes from the repacked file,
                which holds the same data as the FITS file that it is made
                from, so that file is not needed once it has been repacked.
        """
        if not load_samples in (True, False):
            raise TypeError("`load_samples` must be bool")
        if not isinstance(flavor, str):
            raise TypeError("`flavor` must be str")
        if not memmap in (True, False):
            raise TypeError("`memmap` must be bool")
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

        self._flavor = flavor

        # A repacked file holds the same map as the FITS file that it is made
        # from, so it is read either way: mapped into memory, or, with
        # `memmap=False`, read into it. The FITS file is only needed until
        # `fetch()` has repacked it.
        repacked = fits_is_repacked(map_fname)
        if not repacked:
            repacked_fname = _repacked_fname(map_fname)
            if fits_is_repacked(repacked_fname):
                map_fname, repacked = repacked_fname, True
            elif memmap and os.path.isfile(map_fname):
                ve = (
                    "{} does not exist. Call `dustmaps.edenhofer2023.fetch()`"
                    " to repack the FITS file into the memory-mapped layout"
                    " that `memmap=True` reads, or query with `memmap=False`,"
                    " which reads the FITS map itself."
                )
                raise FileNotFoundError(ve.format(repacked_fname))

        if not repacked and not os.path.isfile(map_fname):
            # Neither a repacked copy nor the FITS file that it is made from
            from .dustexceptions import data_missing_message

            msg = data_missing_message(
                "edenhofer2023", "Edenhofer et al. (2023)"
            )
            print(msg, file=sys.stderr)
            err = "{} does not exist".format(repr(map_fname))
            raise FileNotFoundError(err)

        # Whether the maps stay in the file and are mapped into memory, rather
        # than read into it
        self._memmap = repacked and memmap

        if self._memmap:
            if integrated not in (True, False):
                te = "`integrated` must be bool; got {}".format(integrated)
                raise TypeError(te)
            self._rec = _get_sphere_repacked(map_fname)
            self._has_samples = self._rec.n_samples is not None
            if not self._has_samples and load_samples:
                raise ValueError("failed to load samples")
            self._integrated = integrated
            self._allowed_modes = self._modes_for(integrated)
            if not self._rec.has_uncertainty and not self._has_samples:
                # Nothing to take a standard deviation of
                self._allowed_modes = tuple(
                    mode for mode in self._allowed_modes if mode != "std"
                )
            # Nothing is transformed here. The maps stay in the file, and a
            # query maps into them and reads the few values that it needs.
            return

        self._rec = (_get_sphere_repacked_data(map_fname) if repacked
                     else _get_sphere(map_fname))
        self._has_samples = (self._rec.data.ndim == 3)
        if not self._has_samples and load_samples:
            raise ValueError("failed to load samples")
        if map_fname is None and self._rec.data0 is None:
            ve = "default map should have zeroth integrated layer but it doesn't"
            raise ValueError(ve)

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
        elif integrated is False:
            msg = "Optimizing map for querying (this might take a couple of seconds)..."
            print(msg, file=sys.stderr)
            np.log(self._rec.data, out=self._rec.data)
            if self._rec.data_uncertainty is not None:
                np.square(
                    self._rec.data_uncertainty, out=self._rec.data_uncertainty
                )
        else:
            te = "`integrated` must be bool; got {}".format(integrated)
            raise TypeError(te)
        self._integrated = integrated
        self._allowed_modes = self._modes_for(integrated)

    def _modes_for(self, integrated):
        """
        The modes that a query may ask for. These only depend on `integrated`
        for a map without samples, which has no standard deviation of an
        integrated quantity to offer (see the reference implementation).
        """
        modes = ("mean",)
        if integrated is False or self._has_samples:
            modes += ("std",)
        if self._has_samples:
            modes += ("samples", "random_sample")
        return modes

    @ensure_flat_galactic
    def query(self, coords, mode="mean", seed_or_rng=None, integrated=None):
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
            seed_or_rng (Optional[:obj:`int` or :obj:`numpy.random.Generator`]):
                A seed for, or an existing generator of, the random numbers
                drawn by the 'random_sample' mode. An integer makes each call
                draw the same sample, so that a result can be repeated; a
                generator is drawn from as it stands, so that successive calls
                give different samples. Defaults to :obj:`None`, meaning a
                fresh, unpredictable generator.
            integrated (Optional[:obj:`bool`]): Whether to return integrated
                extinction. Defaults to :obj:`None`, meaning the setting from
                initialization. A map read through a memory map can be queried
                either way; one read from FITS has already been transformed in
                place, so can only be queried as it was initialized.

        Notes:
            To query integrated extinction, set `integrated=True` at
            initialization, or per query if the map is read through a memory
            map.

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

        if integrated is None:
            integrated = self._integrated
            allowed_modes = self._allowed_modes
        else:
            if integrated not in (True, False):
                te = "`integrated` must be bool or None; got {}".format(
                    integrated)
                raise TypeError(te)
            if integrated != self._integrated and not self._memmap:
                ve = ("`integrated` cannot be changed for a map that has been"
                      " read into memory: it has already been transformed in"
                      " place. Reinitialize the query, or read the map through"
                      " a memory map (`memmap=True`), whose queries integrate"
                      " as they are made.")
                raise ValueError(ve)
            allowed_modes = self._modes_for(integrated)

        if mode not in allowed_modes:
            ve = "`mode={!r}` requires samples but none are available"
            raise ValueError(ve.format(mode))

        if self._memmap:
            # Read only the values that each coordinate interpolates between
            return self._query_repacked(coords, mode, integrated, seed_or_rng)

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
                samp_idx = get_rng(seed_or_rng).choice(self.n_samples)
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

    def _query_repacked(self, coords, mode, integrated, seed_or_rng):
        """
        Read the few values that each coordinate interpolates between out of
        the repacked, memory-mapped file, rather than transforming whole maps
        in memory.

        Only the four pixels and two radial bins that a coordinate falls
        between are ever read, and an integrated value is assembled from the
        checkpoint below the queried bin and the density between the two.
        """
        rec = self._rec
        weights = _interp_weights(
            rec.radii, rec.nside, rec.nest,
            coords.galactic.l.deg, coords.galactic.b.deg,
            coords.galactic.distance.to("pc").value
        )
        n_coords = weights.idx_l.size
        dvol = np.diff(rec.coo_bounds)

        # Which samples the query needs: all of them, unless a single random
        # sample was asked for
        if self._has_samples and mode != "random_sample":
            samples = np.arange(rec.n_samples)
        elif self._has_samples:
            samples = np.array([get_rng(seed_or_rng).choice(rec.n_samples)])
        else:
            samples = None
        lead_shape = (() if not self._has_samples or mode == "random_sample"
                      else (samples.size,))

        value = np.empty(lead_shape + (n_coords,), dtype='f8')
        with fits.open(rec.fname, "readonly", memmap=True) as hdul:
            for start in range(0, n_coords, QUERY_BLOCK_COORDS):
                stop = min(start + QUERY_BLOCK_COORDS, n_coords)
                block = self._read_repacked_block(
                    hdul, weights, samples, start, stop, mode, integrated, dvol
                )
                value[..., start:stop] = _interp_reduce(
                    block, _slice_weights(weights, start, stop),
                    lead_shape + (stop - start,)
                )
        res = value.reshape(lead_shape + weights.coord_shape)

        if self._has_samples:
            res = np.exp(res, out=res)
            if mode == "mean":
                res = res.mean(axis=0)
            elif mode == "std":
                res = res.std(axis=0)
            elif mode == "random_sample":
                pass  # No sample-axis reduction necessary
            else:
                assert mode == "samples"
                res = np.swapaxes(res, 0, -1)
        else:
            if mode == "mean":
                res = np.exp(res, out=res)
            else:
                assert mode == "std"
                res = np.sqrt(res, out=res)
        return res

    def _read_repacked_block(self, hdul, weights, samples, start, stop, mode,
                             integrated, dvol):
        """
        The values whose weighted mean each coordinate in ``[start, stop)``
        interpolates: the two radial bins that bracket it, the four pixels it
        falls between, and each of the samples the query needs.

        These are the values that the reference implementation interpolates,
        so the log of the density, the square of the uncertainty, or the log
        of the cumulative sum: a query takes the log of the two values that it
        interpolates between, rather than of a whole map.
        """
        bins = np.stack((weights.idx_l[start:stop], weights.idx_r[start:stop]))
        pix = weights.idx_pos[:, start:stop]

        if mode == "std" and not self._has_samples:
            block = _gather(hdul["STD"].data, pix, bins)
            if not integrated:
                # The reference squares the uncertainty unless the map is
                # integrated (see `query`)
                np.square(block, out=block)
            return block

        if integrated:
            block = _gather_integral(
                hdul, pix, bins, samples, dvol, self._rec.checkpoint_step,
                self._rec.n_radii
            )
        else:
            block = _gather(hdul["DATA"].data, pix, bins, samples)
        np.log(block, out=block)
        return block

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
        if not self._has_samples:
            return None
        if self._memmap:
            return self._rec.n_samples
        return self._rec.data.shape[0]

    @property
    def flavor(self):
        """
        Returns the flavor of the map. This is a string that is either
        ``"main"`` or ``"less_data_but_2kpc"``.
        """
        return self._flavor


def fetch(clobber=False, fetch_samples=False, fetch_2kpc=False, repack=True):
    """
    Downloads the 3D dust map of Edenhofer et al. (2023).

    Each file is also repacked into the memory-mapped layout that
    `.Edenhofer2023Query` reads when it is initialized with `memmap=True`,
    which is what allows a query to interpolate between four pixels and two
    radial bins, and to integrate on the fly, without transforming a whole map
    in memory. Repacking takes a few minutes and about as much disk space as
    the FITS file.

    Args:
        clobber (Optional[bool]): If ``True``, any existing file will be
            overwritten, even if it appears to match. If ``False`` (the
            default), ``fetch()`` will attempt to determine if the dataset
            already exists. This determination is not 100% robust against data
            corruption.
        fetch_samples (Optional[bool]): If ``True``, the samples will also be
            downloaded. If ``False`` (the default), only the mean and standard
            deviation will be downloaded taking up about 3.2 GB in size while
            the samples take up 19 GB.
        fetch_2kpc (Optional[bool]): If ``True``, the validation run using
            less data though which extends out to 2kpc in distance will also be
            downloaded.
        repack (Optional[bool]): If ``True`` (the default), each downloaded
            file is also repacked into the memory-mapped layout used by
            `.Edenhofer2023Query`, and the FITS file that it was made from is
            then removed: it holds the same map, and only the repacked file is
            needed to query it. A file that has already been repacked is left
            alone, so nothing is downloaded twice. If ``False``, the maps are
            left as FITS files, which take no less disk space but are read a
            whole column at a time.
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

        if repack and not clobber and fits_is_repacked(_repacked_fname(fname)):
            # The repacked file holds the same map, and there is nothing to
            # download
            continue

        # Download from the server
        url = 'https://zenodo.org/record/8187943/files/{}'.format(fn)
        fetch_utils.download_and_verify(url, md5sum, fname, clobber=clobber)

        if repack:
            msg = "Repacking {} (this might take a couple of minutes)..."
            print(msg.format(fn), file=sys.stderr)
            repack_and_remove(fname)
