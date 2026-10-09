#!/usr/bin/env python
"""
Repack a DECaPS 3D dust map file into the chunk layout that random-point queries
want to read.

    python repack_decaps.py decaps_mean_and_samples.h5
    python repack_decaps.py decaps_mean_and_samples.h5 repacked.h5
    python repack_decaps.py --dry-run decaps_mean_and_samples.h5
    python repack_decaps.py --no-verify decaps_mean_and_samples.h5

Only `h5py` and `numpy` are needed. The input file is never modified: the
repacked file is written alongside it, and the last step is left to you.


Why repack?
-----------
As published, `decaps_mean_and_samples.h5` stores its two arrays in chunks of
shape `(200842 pixels, 1 sample, 1 distance bin)`. A chunk therefore holds a
single `(sample, distance)` value for a slab of 200842 pixels. That is efficient
for reading the map from end to end, which is how the file was written, but for a
query at a random coordinate it means touching one ~400 kB chunk for each of the
5 samples x 120 distance bins:

    query                published          repacked
    ------------------   ----------------   ------------------
    1 coordinate         ~0.24 GB read      ~0.03 MB read
    16 random coords     0.84 s, 0.24 GB    0.0006 s, 0.3 MB
    16 scattered coords  13.1 s, 3.86 GB    0.0006 s, 0.3 MB
    256 random coords    reads everything  0.013 s, 5 MB

Repacking so that each chunk holds 16 adjacent pixels with *all* of their
samples and distance bins fixes that, and it also makes the file much smaller,
because the 120 distance bins of the same 16 pixels compress better together than
long runs of a single `(sample, distance)` pair do. Repacking the published
`decaps_mean_and_samples.h5` took it from 33.84 GB to **18.37 GB**, i.e. 46%
smaller, and made a query at a random coordinate about a thousand times faster:

    layout                            size per pixel   whole file
    -------------------------------   --------------   ----------
    published (200842, 1, 1)          0.658 kB         33.84 GB
    16 pixels, deflate 7, shuffle on  0.357 kB         18.37 GB

The saving is almost entirely in the samples (0.491 kB per pixel becomes 0.212);
the mean map, which is already an average, only goes from 0.145 to 0.116.

16 pixels per chunk is the compromise that `dustmaps` uses. A test on a slice of
the map suggested that 64 pixels would save about another 5%, at the cost of
reading five times more data for a query that touches many pixels, so 16 was kept.
The byte shuffling of the published file was kept too (it is worth a few per cent
here), and the deflate level was raised from 3 to 7.


What this script does
---------------------
* opens the input file and reports its layout,
* writes every dataset to a new file, in chunks of 16 pixels x (the rest of the
  dataset's shape), with deflate level 7 and byte shuffling,
* copies every group (including `pixel_info`), every file attribute and every
  dataset attribute unchanged, and adds `repacked = True` and `chunk_pixels = 16`
  to the file attributes,
* optionally re-reads the output and checks that it holds exactly the same data
  as the input, pixel by pixel,
* prints the sizes, the compression achieved and the MD5 of the output file.

Nothing about the values is changed: the pixels stay in the same order, and the
repacked file is a drop-in replacement for the published one.


Time
----
Repacking the 33.8 GB "mean and samples" file took 45 minutes on a 2025 laptop
(about 74 GB of data has to be deflate-compressed, at roughly 28 MB/s), and
verifying 5% of the pixel blocks took a couple of minutes more. Verifying all of
them reads both files from end to end, so allow another 20-30 minutes for it, or
use `--verify-fraction 0.05` to sample the blocks instead, which still touches
the first and last blocks. Progress is printed as the work goes, so the script can
be left to run.
"""
import argparse
import hashlib
import os
import sys
import time

import numpy as np
import h5py


CHUNK_PIXELS = 16
COMPRESSION_OPTS = 7
SHUFFLE = True

# Number of pixels read and written at a time. This only affects how much memory
# the script uses.
BLOCK_PIXELS = 100000


def human(n_bytes):
    """Formats a number of bytes for people to read."""
    for unit, scale in (('GB', 1e9), ('MB', 1e6), ('kB', 1e3)):
        if n_bytes >= scale:
            return '{:.2f} {}'.format(n_bytes / scale, unit)
    return '{} B'.format(n_bytes)


def md5sum(fname, block=1 << 22):
    """MD5 of a file, computed by streaming it."""
    sig = hashlib.md5()
    with open(fname, 'rb') as f:
        for chunk in iter(lambda: f.read(block), b''):
            sig.update(chunk)
    return sig.hexdigest()


def storage_size(dset):
    """Bytes that a dataset occupies in the file, if h5py can tell us."""
    try:
        return dset.id.get_storage_size()
    except AttributeError:
        return None


def describe(fname):
    """
    Prints the layout of an HDF5 file: its attributes, and for each dataset its
    shape, dtype, chunks, filters and size on disk.
    """
    print('  {} ({})'.format(fname, human(os.path.getsize(fname))))
    with h5py.File(fname, 'r') as f:
        for key, value in f.attrs.items():
            print('    attribute {} = {!r}'.format(key, value))
        for name in f:
            obj = f[name]
            if isinstance(obj, h5py.Group):
                print('    {}/ ({} members)'.format(name, len(obj)))
                for key, value in obj.attrs.items():
                    if np.ndim(value) == 0:
                        print('      attribute {} = {!r}'.format(key, value))
                for sub in obj:
                    print('      {:<20} shape {} dtype {}'.format(
                        sub, obj[sub].shape, obj[sub].dtype))
            else:
                comp = obj.compression
                if comp == 'gzip':
                    comp = 'gzip {}'.format(obj.compression_opts)
                print('    {:<20} shape {} dtype {}'.format(
                    name, obj.shape, obj.dtype))
                print('    {:<20} chunks {} filter {}{}'.format(
                    '', obj.chunks, comp or 'none',
                    ', shuffle' if obj.shuffle else ''))
                if storage_size(obj) is not None:
                    print('    {:<20} stored {} ({:.3f} kB per pixel)'.format(
                        '', human(storage_size(obj)),
                        storage_size(obj) / obj.shape[0] / 1e3))


def datasets_of(f):
    """Names of the datasets (as opposed to groups) at the top level."""
    return [name for name in f if not isinstance(f[name], h5py.Group)]


def check_input(fname, chunk_pixels):
    """
    Looks over the input file before doing any work: it has to be readable, and
    its pixels have to be sorted by HEALPix index, which is what queries depend
    on.
    """
    with h5py.File(fname, 'r') as f:
        names = datasets_of(f)
        if not names:
            raise RuntimeError('{} has no datasets at the top level'.format(fname))

        for name in names:
            n_pix = f[name].shape[0]
            if f[name].chunks is not None and f[name].chunks[0] == chunk_pixels:
                print('Note: {} already looks repacked (chunks {}), so this '
                      'will mostly rewrite it.'.format(name, f[name].chunks))
                break

        if 'pixel_info' in f and 'healpix_index' in f['pixel_info']:
            hpi = f['pixel_info/healpix_index']
            n_pix = hpi.shape[0]
            previous = None
            first = True
            for start in range(0, n_pix, 10 * BLOCK_PIXELS):
                stop = min(start + 10 * BLOCK_PIXELS, n_pix)
                block = hpi[start:stop].astype('i8')
                if not np.all(np.diff(block) > 0):
                    raise RuntimeError(
                        'Pixels are not sorted by healpix_index, which queries '
                        'rely on. Was this file truncated?')
                if previous is not None and block[0] <= previous:
                    raise RuntimeError('Pixels are not sorted by healpix_index.')
                previous = block[-1]
            print('  pixels are sorted by healpix_index: ok')


def repack(h5_in, h5_out, chunk_pixels, compression_opts, shuffle):
    """
    Writes `h5_in` out to `h5_out`, with every dataset chunked as
    `(chunk_pixels,) + shape[1:]`, deflate-compressed at `compression_opts`, and
    with the bytes of each value shuffled if `shuffle` is true.

    Groups, file attributes and dataset attributes are copied unchanged, and
    `repacked` and `chunk_pixels` are recorded in the file attributes.
    """
    t_start = time.time()

    with h5py.File(h5_in, 'r') as f_in, h5py.File(h5_out, 'w') as f_out:
        for key in f_in.attrs:
            f_out.attrs[key] = f_in.attrs[key]
        f_out.attrs['repacked'] = True
        f_out.attrs['chunk_pixels'] = chunk_pixels

        for name in f_in:
            obj = f_in[name]

            # Groups, such as `pixel_info`, are small and are copied as they are
            if isinstance(obj, h5py.Group):
                f_in.copy(name, f_out)
                continue

            dset_in = obj
            n_pix = dset_in.shape[0]
            chunks = (chunk_pixels,) + dset_in.shape[1:]

            print('  {}: {} -> chunks {}'.format(name, dset_in.shape, chunks))

            dset_out = f_out.create_dataset(
                name,
                shape=dset_in.shape,
                dtype=dset_in.dtype,
                chunks=chunks,
                compression='gzip',
                compression_opts=compression_opts,
                shuffle=shuffle)

            for key in dset_in.attrs:
                dset_out.attrs[key] = dset_in.attrs[key]

            t0 = time.time()
            for start in range(0, n_pix, BLOCK_PIXELS):
                stop = min(start + BLOCK_PIXELS, n_pix)
                dset_out[start:stop] = dset_in[start:stop]

                done = stop / float(n_pix)
                elapsed = time.time() - t0
                sys.stdout.write(
                    '\r    {:5.1f}%  elapsed {:>5.0f} s  ETA {:>5.0f} s'.
                    format(100. * done, elapsed,
                           elapsed * (1. - done) / max(done, 1e-9)))
                sys.stdout.flush()
            print('\r    {:5.1f}%  written in {:.0f} s{}'.format(
                100., time.time() - t0, ' ' * 30))

    print('  repacked in {:.0f} s'.format(time.time() - t_start))


def verify(h5_in, h5_out, chunk_pixels, fraction, seed):
    """
    Checks that `h5_out` is a repacked but otherwise identical copy of `h5_in`:
    the same datasets with the same shapes, dtypes and attributes, the same file
    attributes, the right chunk shape, and the same values.

    `fraction` of the pixel blocks are compared value by value, drawn at random
    (plus the first and last blocks, which are always compared). A fraction of
    1.0 compares every pixel.

    Returns True if everything matched.
    """
    print('Verifying...')
    ok = True

    with h5py.File(h5_in, 'r') as f_in, h5py.File(h5_out, 'r') as f_out:
        # File attributes: everything from the input, plus the two new ones
        for key, value in f_in.attrs.items():
            if key not in f_out.attrs:
                print('  FAIL: file attribute "{}" is missing'.format(key))
                ok = False
            elif not np.all(f_out.attrs[key] == value):
                print('  FAIL: file attribute "{}" changed'.format(key))
                ok = False
        for key in ('repacked', 'chunk_pixels'):
            if key not in f_out.attrs:
                print('  FAIL: file attribute "{}" is missing'.format(key))
                ok = False
        if f_out.attrs.get('chunk_pixels', None) != chunk_pixels:
            print('  FAIL: chunk_pixels is {}'.format(
                f_out.attrs.get('chunk_pixels', None)))
            ok = False

        # Same groups, with the same members
        groups_in = [n for n in f_in if isinstance(f_in[n], h5py.Group)]
        groups_out = [n for n in f_out if isinstance(f_out[n], h5py.Group)]
        if groups_in != groups_out:
            print('  FAIL: groups changed: {} -> {}'.format(groups_in, groups_out))
            ok = False

        names_in = datasets_of(f_in)
        names_out = datasets_of(f_out)
        if names_in != names_out:
            print('  FAIL: datasets changed: {} -> {}'.format(names_in, names_out))
            return False

        rng = np.random.RandomState(seed)

        for name in names_in:
            a = f_in[name]
            b = f_out[name]

            if a.shape != b.shape or a.dtype != b.dtype:
                print('  FAIL: {}: {} {} -> {} {}'.format(
                    name, a.shape, a.dtype, b.shape, b.dtype))
                ok = False
                continue

            expected = (min(chunk_pixels, a.shape[0]),) + a.shape[1:]
            if b.chunks != expected:
                print('  FAIL: {}: chunks are {}, expected {}'.format(
                    name, b.chunks, expected))
                ok = False

            for key, value in a.attrs.items():
                if key not in b.attrs or not np.all(b.attrs[key] == value):
                    print('  FAIL: {}: attribute "{}" changed'.format(name, key))
                    ok = False

            # Compare values, block by block
            n_pix = a.shape[0]
            n_blocks = -(-n_pix // BLOCK_PIXELS)
            if fraction >= 1.:
                blocks = list(range(n_blocks))
            else:
                chosen = set(rng.randint(0, n_blocks,
                                         max(1, int(fraction * n_blocks))))
                chosen.update((0, n_blocks - 1))
                blocks = sorted(chosen)

            n_bad = 0
            for i, block in enumerate(blocks):
                start = block * BLOCK_PIXELS
                stop = min(start + BLOCK_PIXELS, n_pix)
                if not np.array_equal(a[start:stop], b[start:stop],
                                      equal_nan=True):
                    n_bad += 1
                sys.stdout.write('\r    {:<10} {:5.1f}%  ({} blocks)'.format(
                    name, 100. * (i + 1) / len(blocks), len(blocks)))
                sys.stdout.flush()
            print('\r    {:<10} compared {} of {} pixel blocks, {} differ{}'
                  .format(name, len(blocks), n_blocks, n_bad, ' ' * 20))

            if n_bad:
                ok = False

    print('  {}'.format('OK: the repacked file holds the same data'
                        if ok else 'FAILED: see above'))
    return ok


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Repack a DECaPS dust map file, so that a query for a '
                    'random coordinate reads a few tens of kB instead of '
                    'hundreds of MB, and so that the file is about 30% smaller.',
        epilog='The input file is never modified.')

    parser.add_argument('input', help='the DECaPS file to repack')
    parser.add_argument('output', nargs='?', default=None,
                        help='where to write the repacked file '
                             '(default: <input>.repacked)')
    parser.add_argument('--chunk-pixels', type=int, default=CHUNK_PIXELS,
                        help='pixels per chunk (default: {})'.format(CHUNK_PIXELS))
    parser.add_argument('--compression-opts', type=int, default=COMPRESSION_OPTS,
                        help='deflate level, 0-9 (default: {})'.format(
                            COMPRESSION_OPTS))
    parser.add_argument('--no-shuffle', action='store_true',
                        help='do not shuffle the bytes of each value before '
                             'compressing it (shuffling saves about 3 percent)')
    parser.add_argument('--verify-fraction', type=float, default=1.,
                        help='fraction of the pixel blocks to compare between '
                             'the input and the output, 0-1 (default: 1, all '
                             'of them)')
    parser.add_argument('--no-verify', action='store_true',
                        help='skip verification')
    parser.add_argument('--seed', type=int, default=0,
                        help='seed for choosing which blocks to verify')
    parser.add_argument('--clobber', action='store_true',
                        help='overwrite the output file if it exists')
    parser.add_argument('--dry-run', action='store_true',
                        help='just report the layout of the input file')

    args = parser.parse_args(argv)

    h5_in = args.input
    h5_out = args.output or h5_in + '.repacked'

    if not os.path.isfile(h5_in):
        parser.error('no such file: {}'.format(h5_in))

    if not 0. <= args.verify_fraction <= 1.:
        parser.error('--verify-fraction must be between 0 and 1')

    print('Input file:')
    describe(h5_in)
    print()

    check_input(h5_in, args.chunk_pixels)
    print()

    if args.dry_run:
        print('Dry run: nothing written.')
        return 0

    if os.path.exists(h5_out) and not args.clobber:
        parser.error('{} already exists (use --clobber to overwrite it)'.format(
            h5_out))

    print('Repacking {} -> {}'.format(h5_in, h5_out))
    print('  {}-pixel chunks, deflate {}, shuffle {}'.format(
        args.chunk_pixels, args.compression_opts, not args.no_shuffle))
    repack(h5_in, h5_out,
           chunk_pixels=args.chunk_pixels,
           compression_opts=args.compression_opts,
           shuffle=not args.no_shuffle)

    print()
    print('Output file:')
    describe(h5_out)

    size_in = os.path.getsize(h5_in)
    size_out = os.path.getsize(h5_out)
    print()
    print('  {} -> {} ({:+.1f}%)'.format(
        human(size_in), human(size_out), 100. * (size_out - size_in) / size_in))

    if not args.no_verify:
        print()
        if not verify(h5_in, h5_out, args.chunk_pixels, args.verify_fraction,
                      args.seed):
            print()
            print('The repacked file does NOT match the input. Do not use it.')
            return 1

    print()
    print('MD5 of the repacked file: {}'.format(md5sum(h5_out)))
    print()
    print('Done. Next steps:')
    print('  * check the summary above, and query the file to make sure it '
          'behaves')
    print('  * upload it under the original name ({}), so that it replaces the '
          'published file'.format(os.path.basename(h5_in)))
    print('  * keep this MD5, and the size in bytes ({}), with the upload'.format(
        size_out))
    return 0


if __name__ == '__main__':
    sys.exit(main())
