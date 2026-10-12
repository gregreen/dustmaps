#!/usr/bin/env python
#
# test_h5_utils.py
# Tests for reading rows out of memory-mapped maps.
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

from __future__ import print_function, division

import os
import shutil
import tempfile
import unittest

import h5py
import numpy as np

from .. import h5_utils


def read_rows_directly(source, row_idx, in_bounds_idx, column_index=()):
    """
    The rows that a query should see, read the obvious way: one fancy-index read
    of the whole map. This is the reference that the streamed groups are checked
    against, and the slow path that :obj:`h5_utils.iter_nd_points` exists to
    avoid, since h5py reads it as a point selection.
    """
    rows = np.unique(row_idx[in_bounds_idx])
    data = source[(rows,) + tuple(column_index)]

    # h5py wants the rows in increasing order, so put them back into the order
    # that the coordinates came in
    return data[np.searchsorted(rows, row_idx[in_bounds_idx])]


def collect(groups):
    """
    Gathers the groups that :obj:`h5_utils.iter_nd_points` yields into a
    dictionary keyed by point, as a query that insisted on holding every row at
    once would have to.
    """
    rows = {}
    for points, data in groups:
        for i, point in enumerate(points):
            rows[int(point)] = data[i]
    return rows


def check_groups(case, source, row_idx, in_bounds_idx, column_index=(),
                 **kwargs):
    """
    Checks that the streamed groups give, for every coordinate that is inside
    the map and for no other, exactly the rows of `source`.
    """
    groups = list(h5_utils.iter_nd_points(
        source, (row_idx,), in_bounds_idx, column_index=column_index, **kwargs))
    got = collect(groups)

    # Every coordinate that falls in the map is handed over, exactly once, and
    # nothing else is
    case.assertEqual(sorted(got), list(np.nonzero(in_bounds_idx)[0]))

    np.testing.assert_allclose(
        np.array([got[i] for i in sorted(got)]),    # noqa: E501
        read_rows_directly(source, row_idx, in_bounds_idx, column_index))

    return groups, got


class SpyDataset(object):
    """
    Stands in for an h5py dataset, recording how it is read.
    """
    def __init__(self, shape, chunks, dtype='f4'):
        self.shape = tuple(shape)
        self.chunks = None if chunks is None else tuple(chunks)
        self.dtype = np.dtype(dtype)
        self.keys = []

    def __getitem__(self, key):
        self.keys.append(key)

        shape = []
        for axis, sub in enumerate(key):
            if isinstance(sub, slice):
                start, stop, step = sub.indices(self.shape[axis])
                shape.append(len(range(start, stop, step)))
            else:
                shape.append(np.size(sub))
        shape += list(self.shape[len(key):])

        return np.zeros(shape, dtype=self.dtype)

    @property
    def spans(self):
        return [(key[0].start, key[0].stop) for key in self.keys]

    @property
    def padding(self):
        return sum(e - s for s, e in self.spans)

    @property
    def chunks_touched(self):
        """
        The chunks that each read touches, as a list of sets, in terms of the
        index of the chunk along the first axis.
        """
        chunk_rows = self.chunks[0]
        return [set(range(s // chunk_rows, (e - 1) // chunk_rows + 1))
                for s, e in self.spans]


class TestGrouping(unittest.TestCase):
    """
    Checks which chunks are read together, through the public iterator.
    """
    def setUp(self):
        self.shape = (2000, 5, 6)

    def _spans(self, rows, chunk_rows=16, **kwargs):
        row_idx = np.asarray(rows, dtype='i8')
        spy = SpyDataset(self.shape, (chunk_rows, 5, 6))
        groups = list(h5_utils.iter_nd_points(
            spy, (row_idx,), np.ones(row_idx.size, dtype=bool), **kwargs))
        return spy.spans, groups

    def test_rows_in_one_chunk_are_read_together(self):
        # Rows 3 to 10 are all inside the first chunk of 16 rows
        spans, _ = self._spans([3, 4, 5, 10])
        self.assertEqual(spans, [(0, 16)])

    def test_neighbouring_chunks_are_read_together(self):
        # Rows 18 and 40 fall in chunks 1 and 2, which are neighbours
        spans, _ = self._spans([18, 40])
        self.assertEqual(spans, [(16, 48)])

    def test_chunks_further_apart_are_read_separately(self):
        # Rows 3 and 70 fall in chunks 0 and 4, which are not neighbours
        spans, _ = self._spans([3, 70])
        self.assertEqual(spans, [(0, 16), (64, 80)])

    def test_max_gap_counts_chunks_between_runs(self):
        # Chunks 0 and 2, with chunk 1 between them: a gap of one chunk joins
        # them, a gap of none does not
        spans, _ = self._spans([3, 35], max_gap=1)
        self.assertEqual(spans, [(0, 48)])

        spans, _ = self._spans([3, 35], max_gap=0)
        self.assertEqual(spans, [(0, 16), (32, 48)])

    def test_a_repeated_row_does_not_break_a_run(self):
        # Several coordinates in the same pixel must not split a run of chunks
        spans, groups = self._spans([0, 0, 0, 20])
        self.assertEqual(spans, [(0, 32)])
        self.assertEqual(sum(c.size for c, _ in groups), 4)

    def test_the_same_chunks_are_read_whatever_the_order(self):
        rows = np.array([3, 70, 71, 900])
        spans, _ = self._spans(rows)
        shuffled, _ = self._spans(rows[::-1])
        self.assertEqual(spans, shuffled)


class TestIterNDPointsInMemory(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.RandomState(1)
        self.data = self.rng.uniform(0., 2., size=(200, 5, 6)).astype('f4')

    def _query(self, size=100, n_out=10):
        row_idx = self.rng.randint(0, 200, size=size)
        row_idx = np.concatenate((row_idx, np.full(n_out, -1)))
        return row_idx.astype('i8'), (row_idx != -1)

    def test_matches_a_direct_read(self):
        for _ in range(5):
            row_idx, in_bounds = self._query()
            check_groups(self, self.data, row_idx, in_bounds,
                         column_index=(slice(0, 3), slice(None)))

    def test_all_out_of_bounds(self):
        row_idx = np.full(10, -1, dtype='i8')
        in_bounds = np.zeros(10, dtype=bool)

        self.assertEqual(list(h5_utils.iter_nd_points(
            self.data, (row_idx,), in_bounds)), [])

    def test_repeated_rows(self):
        row_idx = np.array([3, 3, 3, 7, 3, 7], dtype='i8')
        in_bounds = np.ones(6, dtype=bool)

        _, got = check_groups(self, self.data, row_idx, in_bounds)
        for i in range(6):
            np.testing.assert_allclose(got[i], self.data[row_idx[i]])


class TestIterNDPointsChunked(unittest.TestCase):
    """
    Reads from a real (small) HDF5 file, in the repacked layout of the maps.
    """
    @classmethod
    def setUpClass(cls):
        cls._tmpdir = tempfile.mkdtemp()
        cls.chunk_rows = 16
        cls.row_size = 5 * 6 * np.dtype('f4').itemsize

        rng = np.random.RandomState(2)
        cls.samples = rng.uniform(0., 2., size=(200, 5, 6)).astype('f4')
        cls.best_fit = rng.uniform(0., 2., size=(200, 6)).astype('f4')

        cls._fname = os.path.join(cls._tmpdir, 'test.h5')
        with h5py.File(cls._fname, 'w') as f:
            for name, data in (('samples', cls.samples),
                               ('best_fit', cls.best_fit)):
                f.create_dataset(name, data=data,
                                 chunks=(cls.chunk_rows,) + data.shape[1:],
                                 compression='gzip', compression_opts=4,
                                 shuffle=True)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._tmpdir)

    def setUp(self):
        self.rng = np.random.RandomState(3)
        self.f = h5py.File(self._fname, 'r')

    def tearDown(self):
        self.f.close()

    def _row_idx(self, size=64, scattered=True):
        if scattered:
            rows = self.rng.randint(0, 200, size=size)
            row_idx = np.concatenate((rows, np.full(size // 4, -1)))
        else:
            row_idx = np.arange(0, size // 2)
            row_idx = np.concatenate((row_idx, -np.ones(size // 8,
                                                        dtype='i8')))
        return row_idx.astype('i8')

    def _check(self, row_idx, dataset='samples', **kwargs):
        return check_groups(self, self.f[dataset], row_idx,
                            (row_idx != -1), **kwargs)

    def test_scattered_rows(self):
        for _ in range(5):
            self._check(self._row_idx())

    def test_clustered_rows(self):
        for _ in range(5):
            self._check(self._row_idx(scattered=False))

    def test_row_order_does_not_matter(self):
        row_idx = self._row_idx()
        groups, _ = self._check(row_idx)
        shuffled, _ = self._check(self.rng.permutation(row_idx))

        # The same rows are read, whatever order the coordinates come in
        self.assertEqual(sum(d.shape[0] for _, d in shuffled),
                         sum(d.shape[0] for _, d in groups))

    def test_no_merging(self):
        self._check(np.arange(0, 64, 4), max_gap=0)

    def test_reads_are_trimmed(self):
        # Rows 0, 5, 18 and 40 sit in chunks 0, 0, 1 and 2, which are read
        # together as one hyperslab; only the wanted rows are handed over
        spy = SpyDataset((200, 5, 6), (16, 5, 6))
        row_idx = np.array([0, 5, 18, 40], dtype='i8')
        groups = list(h5_utils.iter_nd_points(
            spy, (row_idx,), np.ones(row_idx.size, dtype=bool)))

        self.assertEqual(spy.spans, [(0, 48)])
        self.assertEqual([data.shape[0] for _, data in groups], [4])

    def test_column_index(self):
        row_idx = self._row_idx()
        groups, got = check_groups(
            self, self.f['samples'], row_idx, (row_idx != -1),
            column_index=(slice(0, 2), slice(None)))

        for _, data in groups:
            self.assertEqual(data.shape[1], 2)

        # The columns that were left out are the ones that were left out
        full = read_rows_directly(self.f['samples'], row_idx, (row_idx != -1))
        np.testing.assert_allclose(
            np.array([got[i] for i in sorted(got)]), full[:, :2, :])

    def test_two_dimensional_dataset(self):
        # The best-fit dataset, which the caller reshapes to add a sample axis
        row_idx = self._row_idx()
        groups, got = self._check(row_idx, dataset='best_fit')

        for _, data in groups:
            self.assertEqual(data.ndim, 2)

        np.testing.assert_allclose(
            np.array([got[i] for i in sorted(got)]),
            self.best_fit[row_idx[row_idx != -1]])

    def test_contiguous_dataset(self):
        # A dataset that is not chunked has no chunks to group by, so the
        # iterator invents them, to keep the amount read at once bounded
        fname = os.path.join(self._tmpdir, 'contiguous.h5')
        with h5py.File(fname, 'w') as f:
            f.create_dataset('samples', data=self.samples)

        with h5py.File(fname, 'r') as f:
            row_idx = self._row_idx()
            check_groups(self, f['samples'], row_idx, (row_idx != -1))

    def test_all_out_of_bounds(self):
        row_idx = np.full(16, -1, dtype='i8')
        in_bounds = np.zeros(16, dtype=bool)

        self.assertEqual(
            list(h5_utils.iter_nd_points(self.f['samples'], (row_idx,),
                                         in_bounds)), [])


class TestReadPattern(unittest.TestCase):
    """
    Checks how the helper reads: whole chunks as hyperslabs, nothing else, and
    never more than `max_bytes` at a time.
    """
    def setUp(self):
        self.shape = (2000, 5, 6)
        self.chunk_rows = 16
        self.row_size = 5 * 6 * np.dtype('f4').itemsize

    def _groups(self, rows, **kwargs):
        row_idx = np.concatenate((np.asarray(rows, dtype='i8'),
                                  np.full(3, -1, dtype='i8')))
        in_bounds = (row_idx != -1)
        spy = SpyDataset(self.shape, (self.chunk_rows, 5, 6))
        groups = list(h5_utils.iter_nd_points(spy, (row_idx,), in_bounds,
                                              **kwargs))
        return spy, groups

    def test_never_reads_by_point_selection(self):
        rows = np.arange(0, 2000, 64)
        spy, _ = self._groups(rows)

        # Each read is a hyperslab over the first axis, aligned to the chunks
        for key in spy.keys:
            self.assertIsInstance(key, tuple)
            self.assertIsInstance(key[0], slice)
            self.assertEqual(key[1:], ())
            self.assertEqual(key[0].start % self.chunk_rows, 0)
            self.assertEqual(key[0].stop % self.chunk_rows, 0)

        # Rows four chunks apart are read on their own
        self.assertEqual(len(spy.keys), rows.size)

    def test_reads_are_in_ascending_order(self):
        rows = np.random.RandomState(4).permutation(np.arange(0, 2000, 64))
        spy, _ = self._groups(rows)

        starts = [s for s, _ in spy.spans]
        self.assertEqual(starts, sorted(starts))

    def test_no_chunk_is_read_twice(self):
        rows = np.sort(np.unique(np.random.RandomState(5).randint(
            0, 2000, size=200)))
        spy, _ = self._groups(rows)

        touched = spy.chunks_touched
        for i in range(len(touched)):
            for j in range(i + 1, len(touched)):
                self.assertFalse(touched[i] & touched[j],
                                 'chunk read by spans {} and {}'.format(i, j))

    def test_holds_at_most_max_bytes(self):
        # One long run of rows is cut into pieces, so that no group reads or
        # hands over more than `max_bytes`
        max_chunks = 25
        rows = np.arange(2000)
        spy, groups = self._groups(
            rows, max_bytes=max_chunks * self.chunk_rows * self.row_size)

        for start, stop in spy.spans:
            self.assertLessEqual(stop - start, max_chunks * self.chunk_rows)
        self.assertEqual(len(spy.spans),
                         2000 // (max_chunks * self.chunk_rows))

        # Every point still comes through
        points = np.concatenate([c for c, _ in groups])
        np.testing.assert_equal(np.sort(points), rows)

    def test_bounds_the_points_that_are_handed_over(self):
        # Many coordinates in one chunk are cut into pieces too, so that no
        # group hands over more than `max_bytes` at once
        max_points = 8
        rows = np.zeros(100, dtype='i8')
        _, groups = self._groups(rows, max_bytes=max_points * self.row_size)

        for points, data in groups:
            self.assertLessEqual(points.size, max_points)
            self.assertEqual(data.shape[0], points.size)

        points = np.concatenate([c for c, _ in groups])
        self.assertEqual(points.size, 100)

    def test_yields_only_the_points_that_were_asked_for(self):
        # The chunks are read whole, but only the wanted rows are handed over
        rows = np.array([0, 5, 18, 40, 41, 900])
        _, groups = self._groups(rows)

        for points, data in groups:
            self.assertEqual(data.shape[0], points.size)

        points = np.concatenate([c for c, _ in groups])
        np.testing.assert_equal(np.sort(points), np.arange(6))

    def test_single_row_chunks_are_never_merged(self):
        # A dataset whose chunks hold one row has nothing to merge
        row_idx = np.array([3, 5, 100, 700], dtype='i8')
        spy = SpyDataset(self.shape, (1, 5, 6))
        list(h5_utils.iter_nd_points(spy, (row_idx,),
                                     np.ones(4, dtype=bool)))
        self.assertEqual(spy.padding, 4)


class TestIterNDPoints(unittest.TestCase):
    """
    The general N-dimensional point iterator.
    """
    @classmethod
    def setUpClass(cls):
        cls._tmpdir = tempfile.mkdtemp()

        rng = np.random.RandomState(7)
        cls.data2d = rng.uniform(0., 2., size=(200, 6)).astype('f4')
        cls.data3d = rng.uniform(0., 2., size=(200, 5, 6)).astype('f4')
        cls.data_grid = rng.uniform(0., 2., size=(64, 32)).astype('f4')

        cls._fname = os.path.join(cls._tmpdir, 'nd.h5')
        with h5py.File(cls._fname, 'w') as f:
            f.create_dataset('d2', data=cls.data2d, chunks=(16, 6),
                             compression='gzip')
            f.create_dataset('d3', data=cls.data3d, chunks=(16, 5, 6),
                             compression='gzip')
            f.create_dataset('grid', data=cls.data_grid, chunks=(16, 8),
                             compression='gzip')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._tmpdir)

    def setUp(self):
        self.rng = np.random.RandomState(8)
        self.f = h5py.File(self._fname, 'r')

    def tearDown(self):
        self.f.close()

    def _collect(self, source, indices, in_bounds=None, **kwargs):
        groups = list(h5_utils.iter_nd_points(source, indices, in_bounds,
                                              **kwargs))
        got = {}
        for positions, data in groups:
            for i, position in enumerate(positions):
                got[int(position)] = data[i]
        return groups, got

    def _check(self, source, indices, reference, in_bounds=None, **kwargs):
        groups, got = self._collect(source, indices, in_bounds, **kwargs)

        if in_bounds is None:
            positions = np.arange(indices[0].size)
        else:
            positions = np.nonzero(in_bounds)[0]

        # Every point that was asked for is handed over, and nothing else
        np.testing.assert_equal(sorted(got), list(positions))
        np.testing.assert_allclose(np.array([got[i] for i in positions]),
                                   reference)

        return groups, got

    def _rows(self, n=50):
        return np.unique(self.rng.randint(0, 200, size=n))

    def test_one_dimension(self):
        rows = self._rows()
        self._check(self.f['d3'], (rows,), self.data3d[rows])

    def test_one_dimension_in_memory(self):
        rows = self._rows()
        self._check(self.data3d, (rows,), self.data3d[rows])

    def test_two_dimensions(self):
        x = self.rng.randint(0, 64, size=60)
        y = self.rng.randint(0, 32, size=60)
        self._check(self.f['grid'], (x, y), self.data_grid[(x, y)])

    def test_two_dimensions_in_memory(self):
        x = self.rng.randint(0, 64, size=60)
        y = self.rng.randint(0, 32, size=60)
        self._check(self.data_grid, (x, y), self.data_grid[(x, y)])

    def test_two_dimensions_with_a_trailing_axis(self):
        x = self.rng.randint(0, 200, size=40)
        y = self.rng.randint(0, 5, size=40)

        groups, _ = self._check(self.f['d3'], (x, y), self.data3d[(x, y)])
        self.assertEqual(groups[0][1].shape[1:], (6,))

    def test_in_bounds_mask(self):
        x = self.rng.randint(0, 200, size=50)
        in_bounds = self.rng.rand(50) < 0.5

        self._check(self.f['d3'], (x,), self.data3d[x[in_bounds]],
                    in_bounds=in_bounds)

    def test_points_outside_the_mask_are_not_read(self):
        spy = SpyDataset((2000, 5, 6), (16, 5, 6))

        list(h5_utils.iter_nd_points(
            spy, (np.array([5, 1500], dtype='i8'),),
            np.array([True, False])))

        # Only the chunk holding row 5 is read, and it is read whole
        self.assertEqual(spy.spans, [(0, 16)])

    def test_empty_selection(self):
        self.assertEqual(list(h5_utils.iter_nd_points(
            self.f['d3'], (np.zeros(0, dtype='i8'),))), [])
        self.assertEqual(list(h5_utils.iter_nd_points(
            self.f['d3'], (np.array([3]),), np.zeros(1, dtype=bool))), [])

    def test_out_of_bounds_raises(self):
        self.assertRaises(IndexError, list, h5_utils.iter_nd_points(
            self.f['d3'], (np.array([200]),)))
        self.assertRaises(IndexError, list, h5_utils.iter_nd_points(
            self.f['d3'], (np.array([-1]),)))
        self.assertRaises(IndexError, list, h5_utils.iter_nd_points(
            self.data3d, (np.array([-1]),)))

    def test_bad_index_arrays_raise(self):
        self.assertRaises(ValueError, list, h5_utils.iter_nd_points(
            self.f['d3'], (np.array([1., 2.]),)))
        self.assertRaises(ValueError, list, h5_utils.iter_nd_points(
            self.f['d3'], (np.zeros(3, dtype='i8'),
                           np.zeros(4, dtype='i8'))))
        self.assertRaises(ValueError, list, h5_utils.iter_nd_points(
            self.f['d3'], (np.zeros((2, 2), dtype='i8'),)))
        self.assertRaises(ValueError, list, h5_utils.iter_nd_points(
            self.f['d2'], (np.zeros(1, dtype='i8'),) * 3))

    def test_reads_whole_chunks(self):
        x = self.rng.randint(0, 64, size=60)
        y = self.rng.randint(0, 32, size=60)

        spy = SpyDataset((64, 32), (16, 8))
        list(h5_utils.iter_nd_points(spy, (x, y)))

        covered = set()
        for key in spy.keys:
            self.assertEqual(len(key), 2)
            self.assertTrue(all(isinstance(sub, slice) for sub in key))

            # Each read is aligned to chunk boundaries
            self.assertEqual(key[0].start % 16, 0)
            self.assertEqual(key[1].start % 8, 0)
            self.assertEqual(key[0].stop % 16, 0)
            self.assertEqual(key[1].stop % 8, 0)

            # And no chunk is read twice
            for row in range(key[0].start // 16, key[0].stop // 16):
                for col in range(key[1].start // 8, key[1].stop // 8):
                    self.assertNotIn((row, col), covered)
                    covered.add((row, col))

        # Every point's chunk was read
        for xi, yi in zip(x, y):
            self.assertIn((xi // 16, yi // 8), covered)

    def test_merges_neighbouring_chunks(self):
        # A dense grid of points costs one read per row of chunks, rather than
        # one read per chunk
        x = np.repeat(np.arange(64), 32)
        y = np.tile(np.arange(32), 64)

        spy = SpyDataset((64, 32), (16, 8))
        groups = list(h5_utils.iter_nd_points(spy, (x, y)))

        self.assertEqual(spy.keys, [(slice(0, 16), slice(0, 32)),
                                    (slice(16, 32), slice(0, 32)),
                                    (slice(32, 48), slice(0, 32)),
                                    (slice(48, 64), slice(0, 32))])
        self.assertEqual(sum(c.size for c, _ in groups), 64 * 32)

    def test_does_not_merge_across_other_axes(self):
        # Chunks that are neighbours along a dimension other than the last can
        # only be read together with the chunks between them, so they are not
        x = np.array([0, 16, 32])
        y = np.array([0, 0, 0])

        spy = SpyDataset((64, 32), (16, 8))
        list(h5_utils.iter_nd_points(spy, (x, y)))

        self.assertEqual(spy.keys, [(slice(0, 16), slice(0, 8)),
                                    (slice(16, 32), slice(0, 8)),
                                    (slice(32, 48), slice(0, 8))])

    def test_respects_max_bytes(self):
        # A budget of two chunks splits a dense row of four
        x = np.repeat(np.arange(64), 32)
        y = np.tile(np.arange(32), 64)

        spy = SpyDataset((64, 32), (16, 8))
        groups = list(h5_utils.iter_nd_points(spy, (x, y), max_bytes=1024))

        for key in spy.keys:
            n_chunks = ((key[0].stop - key[0].start) // 16
                        * ((key[1].stop - key[1].start) // 8))
            self.assertLessEqual(n_chunks * 512, 1024)

        self.assertEqual(sum(c.size for c, _ in groups), 64 * 32)

    def test_column_index(self):
        x = self.rng.randint(0, 200, size=40)

        _, got = self._collect(self.f['d3'], (x,),
                               column_index=(slice(0, 2), slice(None)))

        np.testing.assert_allclose(
            np.array([got[i] for i in sorted(got)]), self.data3d[x][:, :2, :])

    def test_contiguous_dataset(self):
        fname = os.path.join(self._tmpdir, 'contiguous_nd.h5')
        with h5py.File(fname, 'w') as f:
            f.create_dataset('grid', data=self.data_grid)

        with h5py.File(fname, 'r') as f:
            x = self.rng.randint(0, 64, size=60)
            y = self.rng.randint(0, 32, size=60)
            self._check(f['grid'], (x, y), self.data_grid[(x, y)])


if __name__ == '__main__':
    unittest.main()
