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

from __future__ import print_function, division

import numpy as np
import healpy as hp
import astropy.coordinates as coordinates
import astropy.units as units

from functools import wraps
import inspect
import requests
import json

from . import json_serializers
from . import dustexceptions

# import time


def coord2healpix(coords, frame, nside, nest=True):
    """
    Calculate HEALPix indices from an astropy SkyCoord. Assume the HEALPix
    system is defined on the coordinate frame ``frame``.

    Args:
        coords (:obj:`astropy.coordinates.SkyCoord`): The input coordinates.
        frame (:obj:`str`): The frame in which the HEALPix system is defined.
        nside (:obj:`int`): The HEALPix nside parameter to use. Must be a power of 2.
        nest (Optional[:obj:`bool`]): ``True`` (the default) if nested HEALPix ordering
            is desired. ``False`` for ring ordering.

    Returns:
        An array of pixel indices (integers), with the same shape as the input
        SkyCoord coordinates (:obj:`coords.shape`).

    Raises:
        :obj:`dustexceptions.CoordFrameError`: If the specified frame is not supported.
    """
    if coords.frame.name != frame:
        c = coords.transform_to(frame)
    else:
        c = coords

    if hasattr(c, 'ra'):
        phi = c.ra.rad
        theta = 0.5*np.pi - c.dec.rad
        return hp.pixelfunc.ang2pix(nside, theta, phi, nest=nest)
    elif hasattr(c, 'l'):
        phi = c.l.rad
        theta = 0.5*np.pi - c.b.rad
        return hp.pixelfunc.ang2pix(nside, theta, phi, nest=nest)
    elif hasattr(c, 'x'):
        x,y,z = [v.to('kpc').value for v in (c.x,c.y,c.z)]
        return hp.pixelfunc.vec2pix(nside, x, y, z, nest=nest)
    elif hasattr(c, 'w'):
        x,y,z = [v.to('kpc').value for v in (c.u,c.v,c.w)]
        return hp.pixelfunc.vec2pix(nside, x, y, z, nest=nest)
    else:
        raise dustexceptions.CoordFrameError(
            'No method to transform from coordinate frame "{}" to HEALPix.'.format(
                frame))


def ensure_coord_type(f):
    """
    A decorator for class methods of the form

    .. code-block:: python

        Class.method(self, coords, **kwargs)

    where ``coords`` is an :obj:`astropy.coordinates.SkyCoord` object.

    The decorator raises a :obj:`TypeError` if the ``coords`` that gets passed to
    ``Class.method`` is not an :obj:`astropy.coordinates.SkyCoord` instance.

    Args:
        f (class method): A function with the signature
            ``(self, coords, **kwargs)``, where ``coords`` is a :obj:`SkyCoord`
            object containing an array.

    Returns:
        A function that raises a :obj:`TypeError` if ``coords`` is not an
        :obj:`astropy.coordinates.SkyCoord` object, but which otherwise behaves
        the same as the decorated function.
    """
    @wraps(f)
    def _wrapper_func(self, coords, **kwargs):
        if not isinstance(coords, coordinates.SkyCoord):
            raise TypeError('`coords` must be an astropy.coordinates.SkyCoord object.')
        return f(self, coords, **kwargs)
    return _wrapper_func


def reshape_coords(coords, shape):
    pos_attr = ['l', 'b', 'ra', 'dec', 'x', 'y', 'z', 'w', 'u', 'v', 'distance']
    pos_kwargs = {}

    for attr in pos_attr:
        if hasattr(coords, pos_attr):
            pos_kwargs[attr] = np.reshape()
            # TODO: finish reshape


def coords_to_shape(gal, shape):
    l = np.reshape(gal.l.deg, shape) * units.deg
    b = np.reshape(gal.b.deg, shape) * units.deg

    has_dist = hasattr(gal.distance, 'kpc')
    d = np.reshape(gal.distance.kpc, shape) * units.kpc if has_dist else None

    return coordinates.SkyCoord(l, b, distance=d, frame='galactic')


def ensure_flat_frame(f, frame=None):
    def _wrapper_func(self, coords, **kwargs):
        if (frame is not None) and (coords.frame.name != frame):
            coords_transf = coords.transform_to(frame)
        else:
            coords_transf = coords

        is_array = not coords.isscalar
        if is_array:
            orig_shape = coords.shape
            shape_flat = (np.prod(orig_shape),)
            coords_transf = coords_to_shape(coords_transf, shape_flat)
        else:
            coords_transf = coords_to_shape(coords_transf, (1,))

        out = f(self, coords_transf, **kwargs)

        if is_array:
            out.shape = orig_shape + out.shape[1:]
        else:
            out = out[0]

        return out

    return _wrapper_func


def gal_to_shape(gal, shape):
    l = np.reshape(gal.l.deg, shape)*units.deg
    b = np.reshape(gal.b.deg, shape)*units.deg

    has_dist = hasattr(gal.distance, 'kpc')
    d = np.reshape(gal.distance.kpc, shape)*units.kpc if has_dist else None

    return coordinates.SkyCoord(l, b, distance=d, frame='galactic')


def ensure_flat_galactic(f):
    """
    A decorator for class methods of the form

    .. code-block:: python

        Class.method(self, coords, **kwargs)

    where ``coords`` is an :obj:`astropy.coordinates.SkyCoord` object.

    The decorator ensures that the ``coords`` that gets passed to
    ``Class.method`` is a flat array of Galactic coordinates. It also reshapes
    the output of ``Class.method`` to have the same shape (possibly scalar) as
    the input ``coords``. If the output of ``Class.method`` is a tuple or list
    (instead of an array), each element in the output is reshaped instead.

    Args:
        f (class method): A function with the signature
            ``(self, coords, **kwargs)``, where ``coords`` is a :obj:`SkyCoord`
            object containing an array.

    Returns:
        A function that takes :obj:`SkyCoord` input with any shape (including
        scalar).
    """

    @wraps(f)
    def _wrapper_func(self, coords, **kwargs):
        # t0 = time.time()

        if coords.frame.name != 'galactic':
            gal = coords.transform_to('galactic')
        else:
            gal = coords

        # t1 = time.time()

        is_array = not coords.isscalar
        if is_array:
            orig_shape = coords.shape
            shape_flat = (np.prod(orig_shape),)
            # print 'Original shape: {}'.format(orig_shape)
            # print 'Flattened shape: {}'.format(shape_flat)
            gal = gal_to_shape(gal, shape_flat)
        else:
            gal = gal_to_shape(gal, (1,))

        # t2 = time.time()

        out = f(self, gal, **kwargs)

        # t3 = time.time()

        if is_array:
            if isinstance(out, list) or isinstance(out, tuple):
                # Apply to each array in output list
                for o in out:
                    o.shape = orig_shape + o.shape[1:]
            else:   # Only one array in output
                out.shape = orig_shape + out.shape[1:]
        else:
            if isinstance(out, list) or isinstance(out, tuple):
                out = list(out)

                # Apply to each array in output list
                for k,o in enumerate(out):
                    out[k] = o[0]
            else:   # Only one array in output
                out = out[0]

        # t4 = time.time()

        # print('')
        # print('time inside ensure_flat_galactic: {:.4f} s'.format(t4-t0))
        # print('{: >7.4f} s : {: >6.4f} s : transform_to("galactic")'.format(t1-t0, t1-t0))
        # print('{: >7.4f} s : {: >6.4f} s : reshape coordinates'.format(t2-t0, t2-t1))
        # print('{: >7.4f} s : {: >6.4f} s : execute query'.format(t3-t0, t3-t2))
        # print('{: >7.4f} s : {: >6.4f} s : reshape output'.format(t4-t0, t4-t3))
        # print('')

        return out

    return _wrapper_func


def ensure_flat_coords(f):
    """
    A decorator for class methods of the form

    .. code-block:: python

        Class.method(self, coords, **kwargs)

    where ``coords`` is an :obj:`astropy.coordinates.SkyCoord` object.

    The decorator ensures that the ``coords`` that gets passed to
    ``Class.method`` is a flat array. It also reshapes
    the output of ``Class.method`` to have the same shape (possibly scalar) as
    the input ``coords``. If the output of ``Class.method`` is a tuple or list
    (instead of an array), each element in the output is reshaped instead.

    Args:
        f (class method): A function with the signature
            ``(self, coords, **kwargs)``, where ``coords`` is a :obj:`SkyCoord`
            object containing an array.

    Returns:
        A function that takes :obj:`SkyCoord` input with any shape (including
        scalar).
    """

    @wraps(f)
    def _wrapper_func(self, coords, **kwargs):
        is_array = not coords.isscalar

        if is_array:
            orig_shape = coords.shape
            shape_flat = (np.prod(orig_shape),)
            coords = coords.reshape(shape_flat)
        else:
            coords = coords.reshape((1,))

        # t2 = time.time()

        out = f(self, coords, **kwargs)

        # t3 = time.time()

        if is_array:
            if isinstance(out, list) or isinstance(out, tuple):
                # Apply to each array in output list
                for o in out:
                    o.shape = orig_shape + o.shape[1:]
            else:   # Only one array in output
                out.shape = orig_shape + out.shape[1:]
        else:
            if isinstance(out, list) or isinstance(out, tuple):
                out = list(out)

                # Apply to each array in output list
                for k,o in enumerate(out):
                    out[k] = o[0]
            else:   # Only one array in output
                out = out[0]

        # t4 = time.time()

        # print('')
        # print('time inside ensure_flat_galactic: {:.4f} s'.format(t4-t0))
        # print('{: >7.4f} s : {: >6.4f} s : transform_to("galactic")'.format(t1-t0, t1-t0))
        # print('{: >7.4f} s : {: >6.4f} s : reshape coordinates'.format(t2-t0, t2-t1))
        # print('{: >7.4f} s : {: >6.4f} s : execute query'.format(t3-t0, t3-t2))
        # print('{: >7.4f} s : {: >6.4f} s : reshape output'.format(t4-t0, t4-t3))
        # print('')

        return out

    return _wrapper_func


def web_api_method(url,
                   encoder=json_serializers.get_encoder(),
                   decoder=json_serializers.MultiJSONDecoder):
    def decorator(f):
        @wraps(f)
        def api_wrapper(self, *args, **kwargs):
            # Collect the arguments
            data = inspect.getcallargs(f, self, *args, **kwargs)
            data.pop('self')
            kw = data.pop('kwargs', {})
            data.update(**kw)

            # Serialize the arguments
            data = json.dumps(data, cls=encoder)

            # POST request to server
            headers = {'content-type': 'application/json'}
            r = requests.post(
                self.base_url.rstrip('/') + '/' + url.lstrip('/'),
                data=data,
                headers=headers)
            try:
                r.raise_for_status()
            except requests.exceptions.HTTPError as err:
                print('Response received from server:')
                print(r.text)
                raise err

            # Deserialize the response
            return json.loads(r.text, cls=decoder)
        return api_wrapper
    return decorator



class DustMap(object):
    """
    Base class for querying dust maps. For each individual dust map, a different
    subclass should be written, implementing the :obj:`query()` function.
    """

    def __init__(self):
        pass

    @ensure_coord_type
    def __call__(self, coords, **kwargs):
        """
        An alias for :obj:`DustMap.query`.
        """
        return self.query(coords, **kwargs)

    def query(self, coords, **kwargs):
        """
        Query the map at a set of coordinates.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates at which to
                query the map.

        Raises:
            :obj:`NotImplementedError`: This function must be defined by derived
                classes.
        """
        raise NotImplementedError(
            '`DustMap.query` must be implemented by subclasses.\n'
            'The `DustMap` base class should not itself be used.')

    def query_gal(self, l, b, d=None, **kwargs):
        """
        Query using Galactic coordinates.

        Args:
            l (:obj:`float`, scalar or array-like): Galactic longitude, in degrees,
                or as an :obj:`astropy.unit.Quantity`.
            b (:obj:`float`, scalar or array-like): Galactic latitude, in degrees,
                or as an :obj:`astropy.unit.Quantity`.
            d (Optional[:obj:`float`, scalar or array-like]): Distance from the Solar
                System, in kpc, or as an :obj:`astropy.unit.Quantity`. Defaults to
                ``None``, meaning no distance is specified.
            **kwargs: Any additional keyword arguments accepted by derived
                classes.

        Returns:
            The results of the query, which must be implemented by derived
            classes.
        """

        if not isinstance(l, units.Quantity):
            l = l * units.deg
        if not isinstance(b, units.Quantity):
            b = b * units.deg

        if d is None:
            coords = coordinates.SkyCoord(l, b, frame='galactic')
        else:
            if not isinstance(d, units.Quantity):
                d = d * units.kpc
            coords = coordinates.SkyCoord(
                l, b,
                distance=d,
                frame='galactic')

        return self.query(coords, **kwargs)

    def query_equ(self, ra, dec, d=None, frame='icrs', **kwargs):
        """
        Query using Equatorial coordinates. By default, the ICRS frame is used,
        although other frames implemented by :obj:`astropy.coordinates` may also be
        specified.

        Args:
            ra (:obj:`float`, scalar or array-like): Galactic longitude, in degrees,
                or as an :obj:`astropy.unit.Quantity`.
            dec (`float`, scalar or array-like): Galactic latitude, in degrees,
                or as an :obj:`astropy.unit.Quantity`.
            d (Optional[:obj:`float`, scalar or array-like]): Distance from the Solar
                System, in kpc, or as an :obj:`astropy.unit.Quantity`. Defaults to
                ``None``, meaning no distance is specified.
            frame (Optional[:obj:`icrs`]): The coordinate system. Can be ``'icrs'`` (the
                default), ``'fk5'``, ``'fk4'`` or ``'fk4noeterms'``.
            **kwargs: Any additional keyword arguments accepted by derived
                classes.

        Returns:
            The results of the query, which must be implemented by derived
            classes.
        """

        valid_frames = ['icrs', 'fk4', 'fk5', 'fk4noeterms']

        if frame not in valid_frames:
            raise ValueError(
                '`frame` not understood. Must be one of {}.'.format(valid_frames))

        if not isinstance(ra, units.Quantity):
            ra = ra * units.deg
        if not isinstance(dec, units.Quantity):
            dec = dec * units.deg

        if d is None:
            coords = coordinates.SkyCoord(ra, dec, frame='icrs')
        else:
            if not isinstance(d, units.Quantity):
                d = d * units.kpc
            coords = coordinates.SkyCoord(
                ra, dec,
                distance=d,
                frame='icrs')

        return self.query(coords, **kwargs)


class WebDustMap(object):
    """
    Base class for querying dust maps through a web API. For each individual
    dust map, a different subclass should be written, specifying the base URL.
    """

    def __init__(self, api_url=None, map_name=''):
        """
        Initialize the :obj:`WebDustMap` object.

        Args:
            api_url (Optional[:obj:`str`]): The base URL for the API. Defaults to
                ``'http://argonaut.skymaps.info/api/v2/'``.
            map_name (Optional[:obj:`str`]): The name of the dust map to query. For
                example, the Green et al. (2015) dust map is hosted at
                ``http://argonaut.skymaps.info/api/v2/bayestar2015``, so the
                correct specifier for that map is ``map_name='bayestar2015'``.
        """
        if api_url is None:
            api_url = 'http://argonaut.skymaps.info/api/v2/'
        self.base_url = api_url.rstrip('/') + '/' + map_name.lstrip('/')

    @ensure_coord_type
    def __call__(self, coords, **kwargs):
        """
        An alias for :obj:`WebDustMap.query()`.
        """
        return self.query(coords, **kwargs)

    @ensure_coord_type
    @web_api_method('/query')
    def query(self, coords, **kwargs):
        """
        A web API version of :obj:`DustMap.query`. See the documentation for the
        corresponding local query object.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates at which to
                query the map.
        """
        pass

    @web_api_method('/query')
    def query_gal(self, l, b, d=None, **kwargs):
        """
        A web API version of :obj:`DustMap.query_gal()`. See the documentation for
        the corresponding local query object. Queries using Galactic
        coordinates.

        Args:
            l (:obj:`float`, scalar or array-like): Galactic longitude, in degrees,
                or as an :obj:`astropy.unit.Quantity`.
            b (:obj:`float`, scalar or array-like): Galactic latitude, in degrees,
                or as an :obj:`astropy.unit.Quantity`.
            d (Optional[:obj:`float`, scalar or array-like]): Distance from the Solar
                System, in kpc, or as an :obj:`astropy.unit.Quantity`. Defaults to
                ``None``, meaning no distance is specified.
            **kwargs: Any additional keyword arguments accepted by derived
                classes.

        Returns:
            The results of the query.
        """
        pass

    @web_api_method('/query')
    def query_equ(self, ra, dec, d=None, frame='icrs', **kwargs):
        """
        A web API version of :obj:`DustMap.query_equ()`. See the documentation for
        the corresponding local query object. Queries using Equatorial
        coordinates. By default, the ICRS frame is used, although other frames
        implemented by :obj:`astropy.coordinates` may also be specified.

        Args:
            ra (:obj:`float`, scalar or array-like): Galactic longitude, in degrees,
                or as an :obj:`astropy.unit.Quantity`.
            dec (:obj:`float`, scalar or array-like): Galactic latitude, in degrees,
                or as an :obj:`astropy.unit.Quantity`.
            d (Optional[:obj:`float`, scalar or array-like]): Distance from the Solar
                System, in kpc, or as an :obj:`astropy.unit.Quantity`. Defaults to
                ``None``, meaning no distance is specified.
            frame (Optional[icrs]): The coordinate system. Can be 'icrs' (the
                default), 'fk5', 'fk4' or 'fk4noeterms'.
            **kwargs: Any additional keyword arguments accepted by derived
                classes.

        Returns:
            The results of the query.
        """
        pass
