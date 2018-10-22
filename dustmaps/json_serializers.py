#!/usr/bin/env python
#
# json_serializers.py
# Contains JSON (de)serializers for:
#     astropy.units.Quantity
#     astropy.coordinates.SkyCoord
#     numpy.ndarray
#     numpy.dtype
#     numpy.floating
#     numpy.integer
#
# Copyright (C) 2016-2017  Gregory M. Green
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


from __future__ import print_function
import six

import json
import base64
import io
import numpy as np
import astropy.units as units
import astropy.coordinates as coords


def deserialize_tuple(d):
    """
    Deserializes a JSONified tuple.

    Args:
        d (:obj:`dict`): A dictionary representation of the tuple.

    Returns:
        A tuple.
    """
    return tuple(d['items'])


def serialize_dtype(o):
    """
    Serializes a :obj:`numpy.dtype`.

    Args:
        o (:obj:`numpy.dtype`): :obj:`dtype` to be serialized.

    Returns:
        A dictionary that can be passed to :obj:`json.dumps`.
    """
    if len(o) == 0:
        return dict(
            _type='np.dtype',
            descr=str(o))
    return dict(
        _type='np.dtype',
        descr=o.descr)
    # res = []
    # for k in range(len(o)):
    #     res.append((o.names[k], str(o[k])))
    # return dict(
    #     _type='np.dtype',
    #     desc=res)


def deserialize_dtype(d):
    """
    Deserializes a JSONified :obj:`numpy.dtype`.

    Args:
        d (:obj:`dict`): A dictionary representation of a :obj:`dtype` object.

    Returns:
        A :obj:`dtype` object.
    """
    if isinstance(d['descr'], six.string_types):
        return np.dtype(d['descr'])
    descr = []
    for col in d['descr']:
        col_descr = []
        for c in col:
            if isinstance(c, six.string_types):
                col_descr.append(str(c))
            elif type(c) is list:
                col_descr.append(tuple(c))
            else:
                col_descr.append(c)
        descr.append(tuple(col_descr))
    return np.dtype(descr)


def serialize_ndarray_b64(o):
    """
    Serializes a :obj:`numpy.ndarray` in a format where the datatype and shape are
    human-readable, but the array data itself is binary64 encoded.

    Args:
        o (:obj:`numpy.ndarray`): :obj:`ndarray` to be serialized.

    Returns:
        A dictionary that can be passed to :obj:`json.dumps`.
    """
    if o.flags['C_CONTIGUOUS']:
        o_data = o.data
    else:
        o_data = np.ascontiguousarray(o).data
    data_b64 = base64.b64encode(o_data)
    return dict(
        _type='np.ndarray',
        data=data_b64.decode('utf-8'),
        dtype=o.dtype,
        shape=o.shape)


def hint_tuples(o):
    """
    Annotates tuples before JSON serialization, so that they can be
    reconstructed during deserialization. Each tuple is converted into a
    dictionary of the form:

        {'_type': 'tuple', 'items': (...)}

    This function acts recursively on lists, so that tuples nested inside a list
    (or doubly nested, triply nested, etc.) will also be annotated.
    """
    if isinstance(o, tuple):
        return dict(_type='tuple', items=o)
    elif isinstance(o, list):
        return [hint_tuples(el) for el in o]
    else:
        return o


def serialize_ndarray_readable(o):
    """
    Serializes a :obj:`numpy.ndarray` in a human-readable format.

    Args:
        o (:obj:`numpy.ndarray`): :obj:`ndarray` to be serialized.

    Returns:
        A dictionary that can be passed to :obj:`json.dumps`.
    """
    return dict(
        _type='np.ndarray',
        dtype=o.dtype,
        value=hint_tuples(o.tolist()))


def serialize_ndarray_npy(o):
    """
    Serializes a :obj:`numpy.ndarray` using numpy's built-in :obj:`save` function.
    This produces totally unreadable (and very un-JSON-like) results (in "npy"
    format), but it's basically guaranteed to work in 100% of cases.

    Args:
        o (:obj:`numpy.ndarray`): :obj:`ndarray` to be serialized.

    Returns:
        A dictionary that can be passed to :obj:`json.dumps`.
    """
    with io.BytesIO() as f:
        np.save(f, o)
        f.seek(0)
        serialized = json.dumps(f.read().decode('latin-1'))
    return dict(
        _type='np.ndarray',
        npy=serialized)


def deserialize_ndarray_npy(d):
    """
    Deserializes a JSONified :obj:`numpy.ndarray` that was created using numpy's
    :obj:`save` function.

    Args:
        d (:obj:`dict`): A dictionary representation of an :obj:`ndarray` object, created
            using :obj:`numpy.save`.

    Returns:
        An :obj:`ndarray` object.
    """
    with io.BytesIO() as f:
        f.write(json.loads(d['npy']).encode('latin-1'))
        f.seek(0)
        return np.load(f)


def deserialize_ndarray(d):
    """
    Deserializes a JSONified :obj:`numpy.ndarray`. Can handle arrays serialized
    using any of the methods in this module: :obj:`"npy"`, :obj:`"b64"`,
    :obj:`"readable"`.

    Args:
        d (`dict`): A dictionary representation of an :obj:`ndarray` object.

    Returns:
        An :obj:`ndarray` object.
    """
    if 'data' in d:
        x = np.fromstring(
            base64.b64decode(d['data']),
            dtype=d['dtype'])
        x.shape = d['shape']
        return x
    elif 'value' in d:
        return np.array(d['value'], dtype=d['dtype'])
    elif 'npy' in d:
        return deserialize_ndarray_npy(d)
    else:
        raise ValueError('Malformed np.ndarray encoding.')


def serialize_quantity(o):
    """
    Serializes an :obj:`astropy.units.Quantity`, for JSONification.

    Args:
        o (:obj:`astropy.units.Quantity`): :obj:`Quantity` to be serialized.

    Returns:
        A dictionary that can be passed to :obj:`json.dumps`.
    """
    return dict(
        _type='astropy.units.Quantity',
        value=o.value,
        unit=o.unit.to_string())


def deserialize_quantity(d):
    """
    Deserializes a JSONified :obj:`astropy.units.Quantity`.

    Args:
        d (:obj:`dict`): A dictionary representation of a :obj:`Quantity` object.

    Returns:
        A :obj:`Quantity` object.
    """
    return units.Quantity(
        d['value'],
        unit=d['unit'])


def serialize_skycoord(o):
    """
    Serializes an :obj:`astropy.coordinates.SkyCoord`, for JSONification.

    Args:
        o (:obj:`astropy.coordinates.SkyCoord`): :obj:`SkyCoord` to be serialized.

    Returns:
        A dictionary that can be passed to :obj:`json.dumps`.
    """
    representation = o.representation.get_name()
    frame = o.frame.name

    r = o.represent_as('spherical')

    d = dict(
        _type='astropy.coordinates.SkyCoord',
        frame=frame,
        representation=representation,
        lon=r.lon,
        lat=r.lat)

    if len(o.distance.unit.to_string()):
        d['distance'] = r.distance

    return d


def deserialize_skycoord(d):
    """
    Deserializes a JSONified :obj:`astropy.coordinates.SkyCoord`.

    Args:
        d (:obj:`dict`): A dictionary representation of a :obj:`SkyCoord` object.

    Returns:
        A :obj:`SkyCoord` object.
    """
    if 'distance' in d:
        args = (d['lon'], d['lat'], d['distance'])
    else:
        args = (d['lon'], d['lat'])

    return coords.SkyCoord(
        *args,
        frame=d['frame'],
        representation='spherical')


def get_encoder(ndarray_mode='b64'):
    """
    Returns a JSON encoder that can handle:
        * :obj:`numpy.ndarray`
        * :obj:`numpy.floating` (converted to :obj:`float`)
        * :obj:`numpy.integer` (converted to :obj:`int`)
        * :obj:`numpy.dtype`
        * :obj:`astropy.units.Quantity`
        * :obj:`astropy.coordinates.SkyCoord`

    Args:
        ndarray_mode (Optional[:obj:`str`]): Which method to use to serialize
            :obj:`numpy.ndarray` objects. Defaults to :obj:`'b64'`, which converts the
            array data to binary64 encoding (non-human-readable), and stores the
            datatype/shape in human-readable formats. Other options are
            :obj:`'readable'`, which produces fully human-readable output, and
            :obj:`'npy'`, which uses numpy's built-in :obj:`save` function and
            produces completely unreadable output. Of all the methods :obj:`'npy'`
            is the most reliable, but also least human-readable. :obj:`'readable'`
            produces the most human-readable output, but is the least reliable
            and loses precision.

    Returns:
        A subclass of :obj:`json.JSONEncoder`.
    """

    # Use specified numpy.ndarray serialization mode
    serialize_fns = {
        'b64': serialize_ndarray_b64,
        'readable': serialize_ndarray_readable,
        'npy': serialize_ndarray_npy}

    if ndarray_mode not in serialize_fns:
        raise ValueError('"ndarray_mode" must be one of {}'.format(
            serialize_fns.keys))

    serialize_ndarray = serialize_fns[ndarray_mode]

    class MultiJSONEncoder(json.JSONEncoder):
        """
        A JSON encoder that can handle:
            * :obj:`numpy.ndarray`
            * :obj:`numpy.floating` (converted to :obj:`float`)
            * :obj:`numpy.integer` (converted to :obj:`int`)
            * :obj:`numpy.dtype`
            * :obj:`astropy.units.Quantity`
            * :obj:`astropy.coordinates.SkyCoord`
        """
        def default(self, o):
            if isinstance(o, coords.SkyCoord):
                return serialize_skycoord(o)
            if isinstance(o, units.Quantity):
                return serialize_quantity(o)
            elif isinstance(o, np.ndarray):
                return serialize_ndarray(o)
            elif isinstance(o, np.dtype):
                return serialize_dtype(o)
            elif isinstance(o, np.floating):
                return float(o)
            elif isinstance(o, np.integer):
                return int(o)
            elif isinstance(o, np.bool_):
                return bool(o)
            elif isinstance(o, np.void):
                try:
                    o = np.array(o)
                except:
                    pass
                else:
                    return o
            return json.JSONEncoder.default(self, o)

    return MultiJSONEncoder


class MultiJSONDecoder(json.JSONDecoder):
    """
    A JSON decoder that can handle:
        * :obj:`numpy.ndarray`
        * :obj:`numpy.dtype`
        * :obj:`astropy.units.Quantity`
        * :obj:`astropy.coordinates.SkyCoord`
    """
    def __init__(self, *args, **kwargs):
        json.JSONDecoder.__init__(
            self,
            object_hook=self.object_hook,
            *args,
            **kwargs)

    def object_hook(self, d):
        if isinstance(d, dict):
            if ('_type' in d):
                if d['_type'] == 'astropy.coordinates.SkyCoord':
                    return deserialize_skycoord(d)
                elif d['_type'] == 'astropy.units.Quantity':
                    return deserialize_quantity(d)
                elif d['_type'] == 'np.ndarray':
                    return deserialize_ndarray(d)
                elif d['_type'] == 'np.dtype':
                    return deserialize_dtype(d)
                elif d['_type'] == 'tuple':
                    return deserialize_tuple(d)
        return d
