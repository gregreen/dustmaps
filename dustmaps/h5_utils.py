#!/usr/bin/env python
#
# h5_utils.py
# Reads individual rows out of memory-mapped (HDF5) maps.
#
# Copyright (C) 2026  Gregory M. Green
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

"""
Reading random rows out of a memory-mapped map.

:obj:`iter_nd_points()` is what the memory-mapped queries of the Bayestar and
DECaPS maps read through. The obvious way to read a set of points is
``dset[rows]``, and h5py hands that to HDF5 as a *point selection*: HDF5 finds
the chunk holding each element separately, walking the chunk index once per
element. For scattered rows that is enormously more expensive than the data
itself. On the repacked DECaPS map (18 GB, 51.4 million pixels, 3.2 million
chunks of 16 rows), reading 4064 rows scattered over the whole sky takes

    dset[rows, :5, :]                     18.3 s
    the same rows as 4063 hyperslabs       0.27 s

for the same rows, the same chunks and the same 78 MB decompressed. The
difference is the per-element index lookup, not I/O or decompression: that
0.27 s is exactly the time needed to stream those 78 MB (the map reads at about
300 MB/s), and none of it is spent locating anything.

Caching cannot rescue the point selection, which is why the read pattern is
changed here rather than the file being opened differently. Measured on the same
4064-row read, averaged over three passes (the third pass is the interesting one,
because the whole 78 MB working set then fits in the cache many times over):

    default settings                       10.1 / 15.8 / 18.1 s
    rdcc_nbytes=1 GB, rdcc_nslots=100003   17.5 / 18.6 / 18.5 s
    rdcc_nbytes=0 (no chunk cache)         16.8 / 17.6 / 18.1 s
    metadata cache raised to 128 MB        18.7 / 21.3 / 20.7 s

(128 MB is the largest metadata cache that HDF5 will accept.) Since the cost is
proportional to the number of scattered rows, no fixed-size cache can fix it.

Two things that do *not* work, and that are worth not retrying:

- Reading one hyperslab from ``min(rows)`` to ``max(rows)``. For 4096 sightlines
  in a 5 deg x 2 deg patch, the queried rows span 602458 rows, and reading that
  span costs 2.7 s (723 MB decompressed) against 0.8 s for the point selection
  and 4.9 MB of wanted data. Runs of rows are read; spans between them are not.
- ``dset.iter_chunks()``. The chunk shape, ``dset.chunks``, is what matters here
  (it says which points can be read together), and iterating millions of chunk
  index tuples in Python to group points by chunk would cost much more than the
  reads it saves. :obj:`iter_nd_points()` groups them with arithmetic on the
  index arrays instead.

Sorting the *coordinates* does not enter into it either: :obj:`iter_nd_points()`
sorts the points by chunk, so the order in which a query is given its coordinates
makes no difference to the reads.

The data is never accumulated. :obj:`iter_nd_points()` hands points over a group
at a time, one row of the map per point and nothing else, and a query reduces
each group into its results as it arrives. The only arrays that grow with the
size of a query are then the coordinates and the results, so a query for an
entire map needs no more memory than a query for a single coordinate. (Gathering
every requested row first, and indexing that array per coordinate, would instead
need 1.2 kB per queried coordinate for the DECaPS samples: 62 GB for the map.)

:obj:`iter_nd_points()` is the same idea for points that are given by an index
expression per axis, as in ``dset[(x, y)]``. It works in chunks, joining runs of
neighbouring chunks along the first per-point axis into a single read, so for a
dense grid it issues as few reads as the row-wise reader does: 62 reads for 65536
points in a patch of the DECaPS map, against 12438 before those runs were joined.
It reads whole chunks rather than exactly the rows that are wanted, so it reads
more data (55952 rows against 12057, for 4096 points in a patch), which on chunks
this small costs a few per cent of the time. It hands over only what was asked
for.

Each axis of a dataset is chosen in one of two ways, and that is what says how
much of the axis is read. A slice is the same range for every point; an integer
array gives each point its own values along that axis. So scattered coordinates
at one distance are read out of the DECaPS samples, of shape
``(n_pixels, n_samples, n_distances)``, with

    iter_nd_points(dset, (pixels, slice(0, 5), bins), in_bounds)

where ``bins`` holds the two distance bins either side of each coordinate's
distance. A read then covers only the chunks that the group's pixels and bins
fall in -- one chunk of the distance axis rather than all of it -- and joins
neighbouring chunks along a per-point axis when the group spans them.
"""

from __future__ import print_function, division

import numpy as np


# The most data that a read may hold in memory at once. The rows of a map are
# passed to a query in groups, and this limits how long those groups may be, so
# that a query for a large part of a map is no more demanding than a query for a
# single coordinate.
DEFAULT_MAX_BYTES = 64 * 1024**2


def iter_nd_points(dset, selection, in_bounds=None, max_gap=None,
                   max_bytes=DEFAULT_MAX_BYTES):
    """
    Iterates over scattered points of a dataset, one chunk at a time.

    The points are given by one index expression per axis, from the first one on,
    and each entry says how its axis is chosen:

    - a slice is the same range for every point, e.g. ``slice(0, 5)`` for five
      samples, or ``slice(None)`` for all of an axis;
    - an integer array gives each point its own values along that axis: one index
      per point, or one row of `k` indices per point, in which case the point
      asks for `k` values there.

    Axes that the selection does not reach are read whole. The points are then
    read a chunk at a time rather than one point at a time, which is what h5py
    does with fancy indexing and what makes it slow, and handed over as they are
    read, trimmed to what was asked for, so that nothing accumulates except the
    caller's own results.

    A read covers only the chunks that the group's points fall in. Along the
    first per-point axis -- the leading one, which for a map is the pixel axis --
    neighbouring chunks are joined into a single read, so that a dense grid costs
    one read per run of chunks rather than one read per chunk. Along every later
    per-point axis a read stops at the chunks that the group asks for, joining
    them when it needs more than one, so a map whose chunks hold a slice of that
    axis is read one chunk deep rather than whole. Along a sliced axis, the range
    that was asked for is read.

    Indices are never read as counting from the end of an axis. A negative index
    is an error, rather than being turned into the last row of the map, so that
    a caller can use one to mark points that are not in the map at all.

    Args:
        dset: The dataset (or array) to read from.
        selection: A tuple of index expressions, one per axis from the first. An
            entry is either a slice, meaning the same range for every point, or
            an integer array holding one index per point (or one row of `k`
            indices per point). At least one entry must be an array, and at most
            one of those arrays may ask for more than one value per point.
        in_bounds (Optional[:obj:`ndarray` of `bool`]): Which points to read.
            Defaults to all of them.
        max_gap (Optional[:obj:`int`]): Read chunks that are within this many
            chunks of one another together. Defaults to 0, which joins only
            neighbouring chunks.
        max_bytes (Optional[:obj:`int`]): The most data to read, and to hand
            over, in one group. Defaults to :obj:`DEFAULT_MAX_BYTES`.

    Yields:
        A tuple `(positions, data)` for each group of points. `positions` gives
        the points that the group is for, as indices into the per-point arrays.
        `data[i]` is what the point `positions[i]` asked for: the axes chosen
        with slices, in order, followed by the values it asked for along the axis
        that asked for more than one.
    """
    selection = tuple(selection)
    ndim = len(dset.shape)

    if len(selection) > ndim:
        raise ValueError(
            'Got {} index expressions for a {}-dimensional dataset.'.format(
                len(selection), ndim))

    # Which axes are chosen per point, and what the points asked for along them
    per_point = []
    values = {}
    n_points = None
    k_axis = None
    n_values = None

    for axis, entry in enumerate(selection):
        if isinstance(entry, slice):
            continue

        array = np.asarray(entry)
        if not np.issubdtype(array.dtype, np.integer):
            raise ValueError(
                'Axis {} must be chosen with a slice or with an integer array, '
                'but got {}.'.format(axis, array.dtype))
        if array.ndim not in (1, 2):
            raise ValueError(
                'The index array for axis {} must hold one index per point, or '
                'one row of them per point, but has shape {}.'.format(
                    axis, array.shape))
        if n_points is None:
            n_points = array.shape[0]
        elif array.shape[0] != n_points:
            raise ValueError(
                'Every per-point array must hold the same number of points, but '
                'axis {} holds {} and axis {} holds {}.'.format(
                    per_point[0], n_points, axis, array.shape[0]))
        if array.ndim == 2:
            if k_axis is not None:
                raise ValueError(
                    'At most one axis may ask for several values per point, but '
                    'axes {} and {} both do.'.format(k_axis, axis))
            k_axis, n_values = axis, int(array.shape[1])

        per_point.append(axis)
        values[axis] = array

    if n_points is None:
        raise ValueError(
            'At least one axis must be chosen with one index per point; every '
            'entry given is a slice.')

    if in_bounds is None:
        positions = np.arange(n_points, dtype='i8')
    else:
        positions = np.nonzero(in_bounds)[0]

    if positions.size == 0:
        return

    for axis in per_point:
        here = values[axis][positions]
        if here.size and (here.min() < 0 or here.max() >= dset.shape[axis]):
            raise IndexError(
                'Some points fall outside the dataset in axis {} (which has '
                'size {}).'.format(axis, dset.shape[axis]))

    if isinstance(dset, np.ndarray):
        # In memory, a gather costs nothing next to a read, so one group will do
        spans = tuple(slice(0, int(dset.shape[axis])) if axis in values
                      else (selection[axis] if axis < len(selection)
                            else slice(None))
                      for axis in range(ndim))
        yield positions, _gather_block(dset[spans], spans, values, positions,
                                       per_point, n_values)
        return

    chunks = getattr(dset, 'chunks', None)
    if chunks is None:
        # A contiguous dataset has no chunks to group by, so invent some, to
        # keep the amount read at once bounded
        chunks = tuple(min(64, int(size)) for size in dset.shape)

    if max_gap is None:
        max_gap = 0

    # The chunk that each point falls in, along every axis that is chosen per
    # point. An axis whose points ask for several values is keyed by the first
    # chunk they reach, since a read spans from that chunk to the last.
    chunk_of_point = {}
    for axis in per_point:
        here = values[axis][positions] // int(chunks[axis])
        chunk_of_point[axis] = here if here.ndim == 1 else here.min(axis=1)

    # One key per point, saying which chunk it falls in, with the leading
    # per-point axis varying fastest, so that runs are joined along it. Every
    # other per-point axis is separated by a gap wide enough that stepping to
    # the next chunk along it always breaks a run, since two chunks that are
    # neighbours along such an axis can only be read in one piece together with
    # the chunks between them.
    inner_axis = per_point[0]
    outer_axes = per_point[1:]
    inner = chunk_of_point[inner_axis].astype('i8')
    inner_chunks = -(-int(dset.shape[inner_axis]) // int(chunks[inner_axis]))

    if not outer_axes:
        key = inner
    else:
        outer_grid = tuple(-(-int(dset.shape[axis]) // int(chunks[axis]))
                           for axis in outer_axes)
        outer = np.ravel_multi_index(
            tuple(chunk_of_point[axis] for axis in outer_axes), outer_grid)
        key = outer * (inner_chunks + max_gap + 2) + inner

    # How much one chunk of the leading per-point axis costs to read, and how
    # much one point's values cost, so that a group can be held inside
    # `max_bytes`
    itemsize = np.dtype(dset.dtype).itemsize
    read_len = []
    for axis in range(ndim):
        if axis not in values:
            entry = selection[axis] if axis < len(selection) else slice(None)
            read_len.append(len(range(*entry.indices(int(dset.shape[axis])))))
        elif axis == k_axis and n_values > int(chunks[axis]):
            # The read spans the chunks that this many values fall in, which is
            # at most one more chunk than they strictly need
            axis_chunks = -(-int(dset.shape[axis]) // int(chunks[axis]))
            read_len.append(min(axis_chunks,
                                -(-n_values // int(chunks[axis])) + 1)
                            * int(chunks[axis]))
        else:
            read_len.append(int(chunks[axis]))

    chunk_bytes = int(chunks[inner_axis]) * itemsize
    point_bytes = itemsize
    for axis in range(ndim):
        if axis == inner_axis:
            continue
        chunk_bytes *= read_len[axis]
        if axis not in values:
            point_bytes *= read_len[axis]
    if n_values is not None:
        point_bytes *= n_values

    # Walk the chunks in order, joining neighbouring chunks along the leading
    # per-point axis
    order = np.argsort(key, kind='stable')
    starts, stops = _group_bounds(
        key[order], max_gap,
        max(1, int(max_bytes // max(1, chunk_bytes))),
        max(1, int(max_bytes // max(1, point_bytes))))

    for start, stop in zip(starts, stops):
        here = order[start:stop]
        points = positions[here]

        # The block to read: a run of chunks along the leading per-point axis,
        # only the chunks the group asks for along the other per-point axes, and
        # the ranges it asked for along the axes that are the same for every
        # point
        spans = []
        for axis in range(ndim):
            if axis not in values:
                spans.append(selection[axis] if axis < len(selection)
                             else slice(None))
                continue

            blocks = values[axis][points] // int(chunks[axis])
            spans.append(slice(int(blocks.min()) * int(chunks[axis]),
                               min(int(blocks.max() + 1) * int(chunks[axis]),
                                   int(dset.shape[axis]))))
        spans = tuple(spans)

        yield points, _gather_block(dset[spans], spans, values, points,
                                    per_point, n_values)


def _gather_block(block, spans, values, points, per_point, n_values):
    """
    The values that each point asked for, out of a block that was read for it.

    `block` holds the axes of the dataset in their own order, each covering the
    span in `spans`. The axes chosen per point are reduced to the points' own
    values, the axes chosen with slices are left as they were read, and the axes
    come out in the order of the dataset -- with the values asked for along a
    per-point axis, if it asked for several, last.
    """
    if len(per_point) == block.ndim:
        ordered = block
    else:
        kept = [axis for axis in range(block.ndim) if axis not in values]
        ordered = np.transpose(block, tuple(per_point) + tuple(kept))

    # Where in the block each point's own values sit
    local = tuple(values[axis][points] - spans[axis].start
                  for axis in per_point)

    if n_values is None:
        return ordered[local]

    # A point that asks for several values pairs each of them with its own
    # position along the other per-point axes, so those are given a trailing axis
    # of one, which broadcasts against the values it does ask for
    local = tuple(here[:, None] if here.ndim == 1 else here
                  for here in local)
    gathered = ordered[local]

    return np.moveaxis(gathered, 1, -1)


def _group_bounds(keys, max_gap, max_span, max_count=None):
    """
    Splits sorted keys into groups, each small enough to read and hand over in
    one piece.

    Args:
        keys (:obj:`ndarray` of `int`): Sorted keys, one per coordinate, saying
            which chunk of the map it falls in. Coordinates may repeat a key.
        max_gap (:obj:`int`): Read coordinates whose keys are within this much of
            one another together.
        max_span (:obj:`int`): The widest spread of keys that one group may
            cover.
        max_count (Optional[:obj:`int`]): The most coordinates that one group may
            hold. Defaults to `max_span`.

    Returns:
        A tuple `(starts, stops)` of index arrays into `keys`, giving the
        coordinates that belong to each group, in order. No group covers a wider
        spread of keys than `max_span`, nor holds more than `max_count`
        coordinates.

    The grouping is done with whole-array operations, rather than by walking the
    groups in Python, which costs a handful of temporary arrays the length of
    the query (about 40 bytes per coordinate) and is about twenty times faster
    when there are very many groups: a million scattered coordinates grouped in
    0.12 s against 2.25 s for a loop over the groups.
    """
    n = keys.size
    if n == 0:
        return np.zeros(0, dtype='i8'), np.zeros(0, dtype='i8')

    if max_count is None:
        max_count = max_span

    # 1. Runs of keys that are close enough to read together. A key that repeats
    #    does not break a run.
    new_run = np.concatenate(([True], np.diff(keys) > (max_gap + 1)))
    run = np.cumsum(new_run) - 1
    run_start = keys[new_run][run]

    # 2. Cut the runs so that no group covers a wider spread than `max_span`: a
    #    key this far past the start of its run belongs to the next piece
    piece_of_key = (keys - run_start) // max_span
    new_piece = np.concatenate(([True], (run[1:] != run[:-1])
                                | (piece_of_key[1:] != piece_of_key[:-1])))
    piece = np.cumsum(new_piece) - 1
    piece_start = np.nonzero(new_piece)[0]

    # 3. Cut again, so that no group holds more than `max_count` coordinates:
    #    the coordinates of a group are contiguous, so count from where the
    #    group starts
    n_within = (np.arange(n) - piece_start[piece]) // max_count

    starts = np.nonzero(np.concatenate(([True], (piece[1:] != piece[:-1])
                                        | (n_within[1:] != n_within[:-1]))))[0]
    stops = np.append(starts[1:], n)

    return starts.astype('i8'), stops.astype('i8')

