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
array per dimension, as in ``dset[(x, y)]``. It works in chunks, joining runs of
neighbouring chunks along the last indexed dimension into a single read, so for a
dense grid it issues as few reads as the row-wise reader does: 62 reads for 65536
points in a patch of the DECaPS map, against 12438 before those runs were joined.
It reads whole chunks rather than exactly the runs of rows that are wanted, so it
reads more data (55952 rows against 12057, for 4096 points in a patch), which on
chunks this small costs a few per cent of the time. Both hand over only what was
asked for.
"""

from __future__ import print_function, division

import numpy as np


# The most data that a read may hold in memory at once. The rows of a map are
# passed to a query in groups, and this limits how long those groups may be, so
# that a query for a large part of a map is no more demanding than a query for a
# single coordinate.
DEFAULT_MAX_BYTES = 64 * 1024**2


def iter_nd_points(dset, indices, in_bounds=None, column_index=(),
                   max_gap=None, max_bytes=DEFAULT_MAX_BYTES):
    """
    Iterates over scattered points of a dataset, one chunk at a time.

    The points are given by one index array per dimension, as in ``dset[(x, y)]``,
    and read a chunk at a time rather than one point at a time, which is what
    h5py does with fancy indexing and what makes it slow. Points are handed over
    as they are read, trimmed to what was asked for, so that nothing accumulates
    except the caller's own results.

    The indexed dimensions are the leading ones. Any dimensions after them are
    read whole, so ``iter_nd_points(dset, (x,))`` on a 3-D dataset hands over
    blocks of shape ``(n, d1, d2)``, just as ``dset[x]`` would.

    Points that fall in neighbouring chunks along the last indexed dimension are
    read together, so that a dense grid costs one read per run of chunks rather
    than one read per chunk. Chunks along any other dimension cannot be joined
    into a single hyperslab without also reading the chunks between them, so
    stepping to the next chunk along one of those is always a break.

    Indices are never read as counting from the end of an axis. A negative index
    is an error, rather than being turned into the last row of the map, so that
    a caller can use one to mark points that are not in the map at all.

    Args:
        dset: The dataset (or array) to read from.
        indices: A tuple of 1-D integer arrays, one for each of the leading
            dimensions that is indexed, all of the same length.
        in_bounds (Optional[:obj:`ndarray` of `bool`]): Which points to read.
            Defaults to all of them.
        column_index (Optional[:obj:`tuple`]): Index expressions for the axes
            after the indexed ones, applied to each group as it is read, e.g.
            `(slice(0, n_samples), slice(None))`.
        max_gap (Optional[:obj:`int`]): Read chunks that are within this many
            chunks of one another together. Defaults to 0, which joins only
            neighbouring chunks.
        max_bytes (Optional[:obj:`int`]): The most data to read, and to hand
            over, in one group. Defaults to :obj:`DEFAULT_MAX_BYTES`.

    Yields:
        A tuple `(positions, data)` for each chunk of points. `positions` gives
        the points that the group is for, as indices into the index arrays, and
        `data` holds their values, such that `data[i]` is the value of the point
        `positions[i]`.
    """
    indices = tuple(np.asarray(idx) for idx in indices)

    if len(indices) == 0:
        raise ValueError('At least one index array is needed.')

    n_points = indices[0].size
    n_indexed = len(indices)

    for dim, idx in enumerate(indices):
        if idx.ndim != 1 or idx.size != n_points:
            raise ValueError(
                'All index arrays must be 1D and of the same length.')
        if not np.issubdtype(idx.dtype, np.integer):
            raise ValueError('Index arrays must be of an integer type, but '
                             'dimension {} is {}.'.format(dim, idx.dtype))
        if n_indexed > len(dset.shape):
            raise ValueError(
                'Got {} index arrays for a {}-dimensional dataset.'.format(
                    n_indexed, len(dset.shape)))

    if in_bounds is None:
        positions = np.arange(n_points, dtype='i8')
    else:
        positions = np.nonzero(in_bounds)[0]

    if positions.size == 0:
        return

    for dim, idx in enumerate(indices):
        values = idx[positions]
        if values.size and (values.min() < 0 or values.max() >= dset.shape[dim]):
            raise IndexError(
                'Some points fall outside the dataset in dimension {} (which '
                'has size {}).'.format(dim, dset.shape[dim]))

    if isinstance(dset, np.ndarray):
        # In memory, fancy indexing is a gather: one group will do
        yield positions, dset[tuple(idx[positions] for idx in indices)
                             + tuple(column_index)]
        return

    chunks = getattr(dset, 'chunks', None)
    if chunks is None:
        # A contiguous dataset has no chunks to group by, so invent some, to
        # keep the amount read at once bounded
        chunks = tuple(min(64, int(size)) for size in dset.shape)

    if max_gap is None:
        max_gap = 0

    grid = tuple(-(-int(dset.shape[dim]) // int(chunks[dim]))
                 for dim in range(n_indexed))
    chunk_of_point = [indices[dim][positions] // chunks[dim]
                      for dim in range(n_indexed)]

    # One key per point, saying which chunk it falls in, with the last indexed
    # dimension varying fastest. Every other dimension is separated by a gap wide
    # enough that stepping to the next chunk along it always breaks a run, since
    # two chunks that are neighbours along such a dimension can only be read in
    # one piece together with the chunks between them.
    inner = chunk_of_point[-1].astype('i8')
    if n_indexed == 1:
        key = inner
    else:
        outer = np.ravel_multi_index(chunk_of_point[:-1], grid[:-1])
        key = outer * (int(grid[-1]) + max_gap + 2) + inner

    # How much one chunk of the read costs, and how much one point's values cost,
    # so that a group can be held inside `max_bytes`
    itemsize = np.dtype(dset.dtype).itemsize
    trailing = [int(dset.shape[dim])
                for dim in range(n_indexed, len(dset.shape))]
    chunk_bytes = int(np.prod([int(chunks[dim]) for dim in range(n_indexed)]
                              + trailing, dtype='i8')) * itemsize
    point_bytes = int(np.prod(trailing, dtype='i8')) * itemsize

    # Walk the chunks in order, joining neighbouring chunks along the last axis
    order = np.argsort(key, kind='stable')
    starts, stops = _group_bounds(
        key[order], max_gap,
        max(1, int(max_bytes // max(1, chunk_bytes))),
        max(1, int(max_bytes // max(1, point_bytes))))

    for start, stop in zip(starts, stops):
        here = order[start:stop]
        points = positions[here]
        chunk = [chunk_of_point[dim][here] for dim in range(n_indexed)]

        # The points of a group share a chunk along every axis but the last, and
        # along that one they cover a run of neighbouring chunks
        slices = tuple(
            slice(int(chunk[dim].min()) * chunks[dim],
                  min(int(chunk[dim].max() + 1) * chunks[dim], dset.shape[dim]))
            for dim in range(n_indexed))
        block = dset[slices + tuple(column_index)]

        # Where in the block each wanted point sits
        local = tuple(indices[dim][points] - slices[dim].start
                      for dim in range(n_indexed))

        yield points, block[local]


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

