#!/usr/bin/env python
#
# fetch_utils.py
# Utility funcgtions for fetching data files required by the dustmaps library.
#
# Copyright (C) 2016  Gregory M. Green
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

from . import std_paths

import requests

try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen

import contextlib
import shutil
import hashlib
import json
import os
import os.path
import struct

import numpy as np

from progressbar import ProgressBar, UnknownLength
from progressbar.widgets import DataSize, AdaptiveTransferSpeed, Bar, \
                                AdaptiveETA, Percentage, FormatCustomText
from progressbar.utils import scale_1024


# The URL of the Dataverse to use
dataverse = 'https://dataverse.harvard.edu'

# User-Agent used for HTTP requests. Some servers (e.g., the Harvard Dataverse)
# reject the default ``python-requests`` User-Agent with an HTTP 403 error, so
# we identify ourselves explicitly.
USER_AGENT = 'dustmaps'


def _get(url, **kwargs):
    """
    Wrapper around :func:`requests.get` that identifies the client with a
    :obj:`dustmaps` User-Agent.

    The default ``python-requests`` User-Agent is rejected with an HTTP 403
    error by some servers, including the Harvard Dataverse.
    """
    headers = dict(kwargs.pop('headers', None) or {})
    headers.setdefault('User-Agent', USER_AGENT)
    return requests.get(url, headers=headers, **kwargs)


class Error(Exception):
    pass


class DownloadError(Error):
    """
    An exception that occurs while trying to download a file.
    """


def get_md5sum(fname, chunk_size=1024):
    """
    Returns the MD5 checksum of a file.

    Args:
        fname (str): Filename
        chunk_size (Optional[int]): Size (in Bytes) of the chunks that should be
            read in at once. Increasing chunk size reduces the number of reads
            required, but increases the memory usage. Defaults to 1024.

    Returns:
        The MD5 checksum of the file, which is a string.
    """

    def iter_chunks(f):
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            yield chunk

    sig = hashlib.md5()

    with open(fname, 'rb') as f:
        for chunk in iter_chunks(f):
            sig.update(chunk)

        # data = f.read()
        # return hashlib.md5(data).hexdigest()

    return sig.hexdigest()


def h5_file_exists(fname, size_guess=None, rtol=0.3, atol=1., dsets={}):
    """
    Returns ``True`` if an HDF5 file exists, has the expected file size, and
    contains (at least) the given datasets, with the correct shapes.

    This is for the files that :obj:`dustmaps` builds itself, by repacking a
    downloaded file. Such a file can differ slightly in size depending on the
    versions of the libraries that wrote it (h5py, astropy, ...), so the size
    check is deliberately loose, and the dataset check is the one that matters.

    Args:
        fname (str): Filename to check.
        size_guess (Optional[int]): Expected size (in Bytes) of the file. If
            ``None`` (the default), then filesize is not checked.
        rtol (Optional[float]): Relative tolerance for filesize. Defaults to
            0.3, i.e. 30%.
        atol (Optional[float]): Absolute tolerance (in Bytes) for filesize.
        dsets (Optional[dict]): Dictionary specifying expected datasets. Each
            key is the name of a dataset, while each value is the expected shape
            of the dataset. Defaults to ``{}``, meaning that no datasets are
            checked.

    Returns:
        ``True`` if the file matches by all given criteria.
    """
    # Check if the file exists
    if not os.path.isfile(fname):
        # print('File does not exist.')
        return False

    # Check file size, withe the given tolerances
    if size_guess is not None:
        size = os.path.getsize(fname)
        tol = atol + rtol * size_guess

        if abs(size - size_guess) > tol:
            # print('File size is wrong:')
            # print('  expected: {: >16d}'.format(size_guess))
            # print('     found: {: >16d}'.format(size))
            return False

    # Check the datasets in the file
    if len(dsets):
        import h5py
        try:
            with h5py.File(fname, 'r') as f:
                for key in dsets:
                    # Check that dataset is in file
                    if key not in f:
                        # print('Dataset "{}" not in file.'.format(key))
                        return False
                    # Check that the shape of the dataset is correct
                    if dsets[key] is not None:
                        if f[key].shape != dsets[key]:
                            # print('Dataset "{}" has wrong shape:'.format(key))
                            # print('  expected: {}'.format(dsets[key]))
                            # print('     found: {}'.format(f[key].shape))
                            return False
        except IOError:
            # print('Problem reading file.')
            return False

    return True


class FileTransferProgressBar(ProgressBar):
    def __init__(self, content_length):
        prefixes = ('', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi', 'Yi')
        if content_length is None:
            size_txt = '?'
            content_length = UnknownLength
        else:
            scaled, power = scale_1024(content_length, len(prefixes))
            size_txt = '{:.1f} {:s}B'.format(scaled, prefixes[power])
        widgets = [
            DataSize(),
            FormatCustomText(' of {:s} | '.format(size_txt)),
            AdaptiveTransferSpeed(samples=100),
            FormatCustomText(' '),
            Bar(),
            FormatCustomText(' '),
            Percentage(),
            FormatCustomText(' | '),
            AdaptiveETA(samples=100)]
        super(FileTransferProgressBar, self).__init__(
            max_value=content_length,
            widgets=widgets)


def check_md5sum(fname, md5sum, chunk_size=1 << 16):
    """
    Checks that a file exists, and has the correct MD5 checksum.

    Args:
        fname (str): The filename of the file.
        md5sum (str): The expected MD5 sum.
        chunk_size (Optional[int]): Process in chunks of this size (in Bytes).
            Defaults to 65536 (64 KiB).
    """
    if os.path.isfile(fname):
        md5_existing = get_md5sum(fname, chunk_size=chunk_size)
        return (md5_existing == md5sum)
    return False


# Tag identifying the format of the sparse MD5 sums below, so that a sum
# computed by an older or newer version of `get_sparse_md5sum` cannot be
# mistaken for a file that has been corrupted.
SPARSE_MD5_TAG = b'dustmaps-sparse-md5-v1'


def sparse_md5_offsets(size, n_samples=32, chunk_size=1 << 20):
    """
    Returns the offsets at which :func:`get_sparse_md5sum` samples a file.

    The offsets depend only on ``size``, so that a file's sparse MD5 sum can be
    compared against one computed elsewhere -- on another machine, or in a
    previous version of :obj:`dustmaps`.

    They are drawn from a :obj:`numpy.random.RandomState`, the legacy
    generator, rather than from :obj:`numpy.random.default_rng`: numpy freezes
    the stream of the legacy generator for compatibility, whereas the bit
    generator behind :obj:`default_rng` may change between versions, and these
    sums have to keep working for years. :obj:`RandomState` will not take a seed
    of 2**32 or more, and these files are tens of GB, so only the low 32 bits of
    the file's length seed it; the whole length is still used, as the range the
    offsets are drawn from, and is hashed into the sum as well.

    Args:
        size (int): Length of the file, in Bytes.
        n_samples (Optional[int]): Number of chunks to sample. Defaults to 32.
        chunk_size (Optional[int]): Size of each chunk, in Bytes. Defaults to
            1 MB.

    Returns:
        :obj:`numpy.ndarray`: The offsets of the chunks, in ascending order.
    """
    # The dtype is given explicitly, because the default width of an integer
    # depends on the platform, and these offsets have to come out the same
    # everywhere. The first chunk can start anywhere from 0 to
    # `size - chunk_size`, inclusive.
    rng = np.random.RandomState(size % 2**32)
    offsets = rng.randint(0, size - chunk_size + 1, size=n_samples, dtype='i8')
    offsets.sort()

    return offsets


def get_sparse_md5sum(fname, n_samples=32, chunk_size=1 << 20):
    """
    Returns the MD5 sum of a deterministic sample of the contents of a file.

    What a sparse sum buys is that a file which is not the one that was
    expected -- truncated, incomplete, or of a different version -- is
    vanishingly unlikely to have the same one, which is what neither the file
    size on its own nor the shape of the datasets inside an HDF5 file can tell.
    What it does not buy is a guarantee: only about
    ``n_samples * chunk_size / size`` of the file is read, so corruption
    elsewhere in it is invisible. Use :func:`get_md5sum` when the whole file can
    be read, and this one only to decide whether a file that is already on disk
    has to be fetched again.

    The tag above, the length of the file, and the parameters are hashed along
    with the chunks, so that a sum computed with other parameters, or of a file
    of a different length, cannot be compared against this one.

    Args:
        fname (str): Filename of the file to checksum.
        n_samples (Optional[int]): Number of chunks to read. Defaults to 32.
        chunk_size (Optional[int]): Size of each chunk, in Bytes. Defaults to
            1 MB, i.e. 32 MB read in total. If the file is no longer than
            ``n_samples * chunk_size``, then all of it is read.

    Returns:
        The MD5 sum of the sample, as a string of hex digits.
    """
    size = os.path.getsize(fname)

    hasher = hashlib.md5()
    hasher.update(SPARSE_MD5_TAG)
    hasher.update(struct.pack('<QQQ', size, n_samples, chunk_size))

    if size <= n_samples * chunk_size:
        # Short enough to read in full
        with open(fname, 'rb') as f:
            for chunk in iter(lambda: f.read(1 << 22), b''):
                hasher.update(chunk)
        return hasher.hexdigest()

    with open(fname, 'rb') as f:
        for offset in sparse_md5_offsets(size, n_samples, chunk_size):
            f.seek(int(offset))
            hasher.update(f.read(chunk_size))

    return hasher.hexdigest()


def check_sparse_md5sum(fname, md5sum, **kwargs):
    """
    Checks that a file exists, and has the expected sparse MD5 sum.

    Args:
        fname (str): Filename of the file.
        md5sum (str): The expected sparse MD5 sum, as computed by
            :func:`get_sparse_md5sum` with the same parameters.
        **kwargs: Any other keyword arguments are passed to
            :func:`get_sparse_md5sum`.

    Returns:
        ``True`` if the file exists and matches the sum.
    """
    if not os.path.isfile(fname):
        return False

    return (get_sparse_md5sum(fname, **kwargs) == md5sum)


def download_and_verify(url, md5sum, fname=None,
                        chunk_size=1 << 16, clobber=False,
                        verbose=True):
    """
    Downloads a file and verifies the MD5 sum.

    Args:
        url (str): The URL to download.
        md5sum (str): The expected MD5 sum.
        fname (Optional[str]): The filename to store the downloaded file in.
            If `None`, infer the filename from the URL. Defaults to `None`.
        chunk_size (Optional[int]): Process in chunks of this size (in Bytes).
            Defaults to 65536 (64 KiB).
        clobber (Optional[bool]): If `True`, any existing, identical file will
            be overwritten. If `False`, the MD5 sum of any existing file with
            the destination filename will be checked. If the MD5 sum does not
            match, the existing file will be overwritten. Defaults to `False`.
        verbose (Optional[bool]): If `True` (the default), then a progress bar
            will be shownd during downloads.

    Returns:
        The filename the URL was downloaded to.

    Raises:
        DownloadError: The MD5 sum of the downloaded file does not match
            `md5sum`.
        requests.exceptions.HTTPError: There was a problem connecting to the
            URL.
    """

    # Determine the filename
    if fname is None:
        fname = url.split('/')[-1]

    # Check if the file already exists on disk
    if (not clobber) and os.path.isfile(fname):
        print(('Checking existing file ({}) to see '
               'if MD5 sum matches ...').format(fname))
        md5_existing = get_md5sum(fname, chunk_size=chunk_size)
        if md5_existing == md5sum:
            print('File exists. Not overwriting.')
            return fname
        else:
            print('  -> MD5 sums do not match:')
            print('     * expected: {}'.format(md5sum))
            print('     *    found: {}'.format(md5_existing))

    # Make sure the directory it's going into exists
    dir_name = os.path.dirname(fname)
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    sig = hashlib.md5()

    if verbose:
        print('Downloading {} ...'.format(url))

    if url.startswith('http://') or url.startswith('https://'):
        # Stream the URL as a file, copying to local disk
        with contextlib.closing(_get(url, stream=True)) as r:
            try:
                r.raise_for_status()
            except requests.exceptions.HTTPError as error:
                print('Error connecting to URL: "{}"'.format(url))
                print(r.text)
                raise error

            with open(fname, 'wb') as f:
                content_length = r.headers.get('content-length')
                if content_length is not None:
                    content_length = int(content_length)
                bar = FileTransferProgressBar(content_length)

                for k,chunk in enumerate(r.iter_content(chunk_size=chunk_size)):
                    f.write(chunk)
                    sig.update(chunk)

                    if verbose:
                        bar_val = chunk_size*(k+1)
                        if content_length is not None:
                            bar_val = min(bar_val, content_length)
                        bar.update(bar_val)
    else: # e.g., ftp://
        with contextlib.closing(urlopen(url)) as r:
            content_length = r.headers.get('content-length')
            if content_length is not None:
                content_length = int(content_length)
            bar = FileTransferProgressBar(content_length)

            with open(fname, 'wb') as f:
                k = 0
                while True:
                    chunk = r.read(chunk_size)

                    if not chunk:
                        break

                    f.write(chunk)
                    sig.update(chunk)

                    if verbose:
                        k += 1
                        bar_val = chunk_size*k
                        if content_length is not None:
                            bar_val = min(bar_val, content_length)
                        bar.update(bar_val)


    if sig.hexdigest() != md5sum:
        raise DownloadError('The MD5 sum of the downloaded file is incorrect.\n'
                            + '  download: {}\n'.format(sig.hexdigest())
                            + '  expected: {}\n'.format(md5sum))

    return fname


def download(url, fname=None):
    """
    Downloads a file.

    Args:
        url (str): The URL to download.
        fname (Optional[str]): The filename to store the downloaded file in. If
            `None`, take the filename from the URL. Defaults to `None`.

    Returns:
          The filename the URL was downloaded to.

    Raises:
        requests.exceptions.HTTPError: There was a problem connecting to the
            URL.
    """
    # Determine the filename
    if fname is None:
        fname = url.split('/')[-1]

    # Stream the URL as a file, copying to local disk
    with contextlib.closing(_get(url, stream=True)) as r:
        try:
            r.raise_for_status()
        except requests.exceptions.HTTPError as error:
            print('Error connecting to URL: "{}"'.format(url))
            print(r.text)
            raise error

        with open(fname, 'wb') as f:
            shutil.copyfileobj(r.raw, f)

    return fname


def dataverse_search_doi(doi):
    """
    Fetches metadata pertaining to a Digital Object Identifier (DOI) in the
    Harvard Dataverse.

    Args:
        doi (str): The Digital Object Identifier (DOI) of the entry in the
            Dataverse.

    Raises:
        requests.exceptions.HTTPError: The given DOI does not exist, or there
            was a problem connecting to the Dataverse.
    """

    url = '{}/api/datasets/:persistentId?persistentId=doi:{}'.format(dataverse, doi)
    r = _get(url)

    try:
        r.raise_for_status()
    except requests.exceptions.HTTPError as error:
        print('Error looking up DOI "{}" in the Harvard Dataverse.'.format(doi))
        print(r.text)
        raise error

    return json.loads(r.text)


def dataverse_download_id(file_id, md5sum, **kwargs):
    url = '{}/api/access/datafile/{}'.format(dataverse, file_id)
    download_and_verify(url, md5sum, **kwargs)


def dataverse_download_doi(doi,
                           local_fname=None,
                           file_requirements={},
                           clobber=False):
    """
    Downloads a file from the Dataverse, using a DOI and set of metadata
    parameters to locate the file.

    Args:
        doi (str): Digital Object Identifier (DOI) containing the file.
        local_fname (Optional[str]): Local filename to download the file to. If
            `None`, then use the filename provided by the Dataverse. Defaults to
            `None`.
        file_requirements (Optional[dict]): Select the file containing the
            given metadata entries. If multiple files meet these requirements,
            only the first in downloaded. Defaults to `{}`, corresponding to no
            requirements.

    Raises:
        DownloadError: Either no matching file was found under the given DOI, or
            the MD5 sum of the file was not as expected.
        requests.exceptions.HTTPError: The given DOI does not exist, or there
            was a problem connecting to the Dataverse.

    """
    metadata = dataverse_search_doi(doi)

    def requirements_match(metadata):
        for key in file_requirements.keys():
            if metadata['dataFile'].get(key, None) != file_requirements[key]:
                return False
        return True

    for file_metadata in metadata['data']['latestVersion']['files']:
        if requirements_match(file_metadata):
            file_id = file_metadata['dataFile']['id']
            md5sum = file_metadata['dataFile']['md5']

            if local_fname is None:
                local_fname = file_metadata['dataFile']['filename']

            # Check if the file already exists on disk
            if (not clobber) and os.path.isfile(local_fname):
                print(('Checking existing file ({}) to see '
                       'if MD5 sum matches ...').format(local_fname))
                md5_existing = get_md5sum(local_fname)
                if md5_existing == md5sum:
                    print('File exists. Not overwriting.')
                    return

            print("Downloading data to '{}' ...".format(local_fname))

            dataverse_download_id(file_id, md5sum,
                                  fname=local_fname, clobber=False)

            return

    raise DownloadError(
        'No file found under the given DOI matches the requirements.\n'
        'The metadata found for this DOI was:\n'
        + json.dumps(file_metadata, indent=2, sort_keys=True))


def download_demo():
    doi = '10.5072/FK2/ZSEMG9'
    requirements = {'filename': 'ResizablePng.png'}
    dataverse_download_doi(doi, file_requirements=requirements)
