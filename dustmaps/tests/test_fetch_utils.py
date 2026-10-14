#!/usr/bin/env python
#
# test_fetch_utils.py
# Tests for the MD5 sums in dustmaps.fetch_utils.
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

import hashlib
import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np

from dustmaps import fetch_utils


def content(size, seed=0):
    """
    Returns a reproducible block of bytes, of the given length.

    The bytes are made with :obj:`hashlib` rather than with numpy, so that the
    sums pinned in the tests below cannot drift with a numpy version: numpy
    guarantees the stream of its legacy generator, but not that of
    :obj:`numpy.random.default_rng`.
    """
    digest = hashlib.sha256(b'dustmaps-test' + str(seed).encode()).digest()
    out = bytearray()
    while len(out) < size:
        digest = hashlib.sha256(digest).digest()
        out += digest

    return bytes(out[:size])


class SparseMd5TestCase(unittest.TestCase):
    """
    Writes the files that the tests checksum into a temporary directory.
    """
    def setUp(self):
        self._tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._tmpdir)

    def write(self, size, seed=0, name='data.bin', sparse=False):
        """
        Writes a file of the given length, full of holes if ``sparse``, and
        returns its filename.
        """
        fname = os.path.join(self._tmpdir, name)
        with open(fname, 'wb') as f:
            if sparse:
                f.truncate(size)
            else:
                f.write(content(size, seed))

        return fname

    def flip_byte(self, fname, pos):
        """Flips the byte at the given offset, in place."""
        with open(fname, 'r+b') as f:
            f.seek(pos)
            byte = f.read(1)
            f.seek(pos)
            f.write(bytes([byte[0] ^ 0xff]))


class TestSparseMd5Offsets(unittest.TestCase):
    """
    The chunks that are sampled have to depend on the length of the file, and
    on nothing else, so that a sum computed here can be compared with one
    computed on another machine, or years ago.
    """
    def test_depends_only_on_the_size(self):
        size = 10 * 2**20
        np.testing.assert_equal(fetch_utils.sparse_md5_offsets(size),
                                fetch_utils.sparse_md5_offsets(size))

    def test_different_sizes_give_different_offsets(self):
        a = fetch_utils.sparse_md5_offsets(10 * 2**20)
        b = fetch_utils.sparse_md5_offsets(10 * 2**20 + 1)
        self.assertFalse(np.array_equal(a, b))

    def test_sizes_two_to_the_32_apart_differ(self):
        """
        :obj:`numpy.random.RandomState` refuses a seed of 2**32 or more, so the
        seed is the length modulo 2**32. The whole length is still used, as the
        range that the offsets are drawn from, so files whose lengths differ by
        a multiple of 2**32 do not share offsets either.
        """
        size = 5 * 2**20
        a = fetch_utils.sparse_md5_offsets(size)
        b = fetch_utils.sparse_md5_offsets(size + 2**32)
        self.assertFalse(np.array_equal(a, b))

    def test_offsets_are_sorted(self):
        offsets = fetch_utils.sparse_md5_offsets(2**30, n_samples=100)
        np.testing.assert_equal(offsets, np.sort(offsets))

    def test_offsets_are_in_bounds(self):
        """
        A chunk starts at least ``chunk_size`` bytes before the end of the
        file, so that every read is a full one.
        """
        size = 3 * 2**20
        offsets = fetch_utils.sparse_md5_offsets(size, n_samples=200)
        self.assertGreaterEqual(offsets.min(), 0)
        self.assertLessEqual(offsets.max(), size - 2**20)
        self.assertEqual(offsets.dtype, np.dtype('i8'))

    def test_the_last_legal_offset_is_available(self):
        """
        ``randint``'s upper bound is exclusive, so it has to be given
        ``size - chunk_size + 1`` for the offset that places the chunk exactly
        at the end of the file to be among those that can be drawn.
        """
        size, n_samples, chunk_size = 10 * 2**20, 4, 2**20
        seen = {}

        class FakeRandomState(object):
            """Records how it was called, and returns offsets of zero."""
            def __init__(self, seed):
                seen['seed'] = seed

            def randint(self, low, high, size=None, dtype=None):
                seen['args'] = (low, high, size, dtype)
                return np.zeros(size, dtype=dtype)

        with mock.patch.object(fetch_utils.np.random, 'RandomState',
                               FakeRandomState):
            fetch_utils.sparse_md5_offsets(size, n_samples, chunk_size)

        self.assertEqual(seen['seed'], size % 2**32)
        self.assertEqual(seen['args'],
                         (0, size - chunk_size + 1, n_samples, 'i8'))


class TestGetSparseMd5sum(SparseMd5TestCase):
    """
    The sums themselves. The three that are pinned below cover the two paths
    through :func:`get_sparse_md5sum` -- a file short enough to read in full,
    and one that is sampled -- and freeze the format of a sum, so that a change
    to it cannot be mistaken for a file that has gone bad.
    """
    def test_a_short_file_read_in_full(self):
        fname = self.write(2**20)
        self.assertEqual(fetch_utils.get_sparse_md5sum(fname),
                         'fd1986c7797380609a91a5ffbe4d44de')

    def test_a_file_that_is_sampled(self):
        fname = self.write(2 * 2**20)
        self.assertEqual(
            fetch_utils.get_sparse_md5sum(fname, n_samples=4,
                                          chunk_size=1 << 16),
            '5b87d3df1aa71df77a46edf4aa5515f6')

    def test_a_huge_file(self):
        """
        A 5 GiB file of holes, so that the offsets come from a seed larger than
        the 2**32 that :obj:`numpy.random.RandomState` takes, and reading them
        costs nothing.
        """
        fname = self.write(5 * 2**30, sparse=True)
        self.assertEqual(fetch_utils.get_sparse_md5sum(fname),
                         'c41a24850047f47d947513caa8dda12e')

    def test_repeatable(self):
        fname = self.write(2 * 2**20)
        self.assertEqual(fetch_utils.get_sparse_md5sum(fname),
                         fetch_utils.get_sparse_md5sum(fname))

    def test_different_contents_differ(self):
        a = self.write(2 * 2**20, seed=0, name='a.bin')
        b = self.write(2 * 2**20, seed=1, name='b.bin')
        self.assertNotEqual(fetch_utils.get_sparse_md5sum(a),
                            fetch_utils.get_sparse_md5sum(b))

    def test_a_truncated_file_has_a_different_sum(self):
        """
        Which is what the sum is for: a download that was interrupted leaves a
        file that is shorter than the one that was expected.
        """
        fname = self.write(2 * 2**20)
        full = fetch_utils.get_sparse_md5sum(fname)

        with open(fname, 'r+b') as f:
            f.truncate(2**20)

        self.assertNotEqual(full, fetch_utils.get_sparse_md5sum(fname))

    def test_files_differing_in_length_have_different_sums(self):
        short = self.write(2**20, name='short.bin')
        long = self.write(2**20 + 1, name='long.bin')
        self.assertNotEqual(fetch_utils.get_sparse_md5sum(short),
                            fetch_utils.get_sparse_md5sum(long))

    def test_the_parameters_are_part_of_the_sum(self):
        """
        A sum computed with other parameters is a different sum, so a mismatch
        cannot hide a comparison between two different samples of a file.
        """
        fname = self.write(2 * 2**20)
        sums = {
            fetch_utils.get_sparse_md5sum(fname),
            fetch_utils.get_sparse_md5sum(fname, n_samples=4,
                                          chunk_size=1 << 16),
            fetch_utils.get_sparse_md5sum(fname, n_samples=8,
                                          chunk_size=1 << 16),
        }
        self.assertEqual(len(sums), 3)

    def test_content_outside_the_sample_is_invisible(self):
        """
        Only the sampled chunks are read, which is the trade that this sum
        makes: corruption elsewhere in the file cannot be detected at all. The
        test flips a byte that is not read, and sees the sum stand still, then
        flips one that is, and sees it move.
        """
        size, n_samples, chunk_size = 2 * 2**20, 4, 1 << 16
        fname = self.write(size)

        def sparse_md5():
            return fetch_utils.get_sparse_md5sum(fname, n_samples=n_samples,
                                                 chunk_size=chunk_size)

        before = sparse_md5()

        offsets = fetch_utils.sparse_md5_offsets(size, n_samples, chunk_size)
        sampled = set()
        for offset in offsets:
            sampled.update(range(int(offset), int(offset) + chunk_size))

        self.flip_byte(fname, next(i for i in range(size) if i not in sampled))
        self.assertEqual(before, sparse_md5())

        self.flip_byte(fname, int(offsets[1]))
        self.assertNotEqual(before, sparse_md5())

    def test_a_missing_file_raises(self):
        fname = os.path.join(self._tmpdir, 'absent.bin')
        with self.assertRaises(FileNotFoundError):
            fetch_utils.get_sparse_md5sum(fname)


class TestCheckSparseMd5sum(SparseMd5TestCase):
    def test_a_file_that_matches(self):
        fname = self.write(2 * 2**20)
        md5sum = fetch_utils.get_sparse_md5sum(fname)
        self.assertTrue(fetch_utils.check_sparse_md5sum(fname, md5sum))

    def test_a_file_that_does_not_match(self):
        fname = self.write(2 * 2**20)
        md5sum = fetch_utils.get_sparse_md5sum(fname)

        with open(fname, 'r+b') as f:
            f.truncate(2**20)      # as an interrupted download would

        self.assertFalse(fetch_utils.check_sparse_md5sum(fname, md5sum))

    def test_a_file_that_is_not_there(self):
        fname = os.path.join(self._tmpdir, 'absent.bin')
        self.assertFalse(fetch_utils.check_sparse_md5sum(fname, '0' * 32))


if __name__ == '__main__':
    unittest.main()
