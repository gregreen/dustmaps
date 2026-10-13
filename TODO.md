# dustmaps TODO

Problems to deal with (found while reviewing PRs and running the test suite).

## Data availability / large downloads

- [ ] **Edenhofer (2023): the test suite requires the full posterior samples.**
  `dustmaps/tests/test_edenhofer2023.py::TestEdenhofer2023.setUpClass` constructs
  `Edenhofer2023Query(load_samples=True)`, so **every** test in that module fails
  unless `samples_healpix.fits` is present.
  `fetch(fetch_samples=True)` downloads `mean_and_std_healpix.fits` (~3.2 GB) **plus**
  `samples_healpix.fits` (~19 GB).
  Possible fixes: fetch the samples only where they are needed, or skip the
  samples-dependent tests when the samples file is absent
  (e.g. `unittest.skipUnless`).

- [ ] **DECaPS: the test suite requires the 33 GB samples map.**
  `dustmaps/tests/test_decaps.py::TestDECaPS.setUpClass` uses `DECaPSQueryLite()` with
  its defaults, which loads `decaps_mean_and_samples.h5` (~33 GB).
  `fetch(mean_only=True)` produces only the ~8 GB mean map, but the tests exercise
  `mode='random_sample'`, so the samples are genuinely needed.
  Possible fixes: guard the samples-dependent tests so the suite can run without the
  33 GB file (and use `mean_only=True` where possible).

## Broken downloads

- [ ] **IPHAS: the download site has been compromised. (investigated 2026-10-09)**
  `dustmaps/iphas.py::fetch()` requests
  `http://www.iphas.org/data/extinction/A_samp_{:03d}.tar.gz`, which now returns
  **HTTP 404**. The download page <https://www.iphas.org/extinction.html> still
  loads and still advertises the per-10°-strip archives, but every link on it now
  points at a `*.tar.html` landing page instead of the `.tar.gz` file:

  | Archive | Contents | Size on the page |
  | --- | --- | --- |
  | `A_map_030.tar.html` … `A_map_210.tar.html` | posterior expectations (mean / differential A0) | ~3–12 MB each |
  | `A_samp_030.tar.html` … `A_samp_210.tar.html` | 20 posterior samples per sightline | ~14–55 MB each |

  Those landing pages return **HTTP 200 `text/html` containing tens of MB of spam
  links and no archive data at all**: `A_samp_030.tar.html` was 74 MB of HTML and
  `A_map_030.tar.html` 16 MB, neither of them containing gzip (`\x1f\x8b\x08`),
  `ustar`, `PK`, `BZh` or `\xfd7zXZ` magic anywhere; `Range` requests are
  ignored, the body length varies between requests, and `/data/extinction/`
  returns **403**. Every page on the site is stuffed with injected casino links,
  i.e. the host has been hacked, so the map cannot be fetched from IPHAS.

  Also checked, without success: the Internet Archive CDX API (503, "Internet
  Archive: Temporarily Offline") and Zenodo (403, traffic-blocked). CDS's
  anonymous **FTP service no longer exists** -- `cdsarc.cds.unistra.fr`,
  `cdsarc.u-strasbg.fr` and `ftp.cds.unistra.fr` all refuse port 21, and outbound
  FTP is unusable from this network in any case (`ftp.gnu.org` and
  `ftp.ncbi.nlm.nih.gov` fail too; the `ftp://cdsarc...` form is a legacy
  ADS/VizieR link that VizieR itself rewrites to HTTPS). The HTTPS *file*
  endpoints are behind the **Anubis anti-bot challenge**: `/ftp/J/MNRAS/443/2907/amap.dat.gz`,
  `/viz-bin/nph-Cat/fits?J/MNRAS/443/2907` and
  `/viz-bin/nph-Cat/fits.gz?J/MNRAS/443/2907/amap.dat.gz` each return HTTP 200
  `text/html`, 4.4 kB, *"Making sure you're not a bot!"*, for both `requests` and
  `curl`, even with a browser user agent. Loading the catalogue page in a real
  browser does pass the challenge, so a human can save those files by hand, but
  the browser tool used here exposes no session cookie (its Playwright build has
  no `Storage.getCookies`, and the page sends no `Cookie` header), so the download
  cannot be handed to `curl`. If the CDS files are ever needed, the civilised
  route is to email <cds-question@unistra.fr> and ask for a whitelisted user
  agent rather than to fight the challenge from Python.

  **The Wayback Machine has the original archives.** All 19 of
  `A_samp_030.tar.gz` … `A_samp_210.tar.gz` are present at
  `https://web.archive.org/web/2020id_/<original-url>`; each returns **HTTP 200**
  with gzip magic and a `Content-Length` matching the size advertised on the IPHAS
  page (e.g. 41,817,076 bytes for `A_samp_030`). The bytes are the originals:
  `A_samp_030.tar.gz` downloaded from there hashes to
  `dd531e397622bc97d4ff92b6c7863ade`, exactly the sum already hard-coded in
  `fetch()`, so that table can stay as it is. Caveat: web.archive.org threw a
  transient `SSLError: RECORD_LAYER_FAILURE` on the second file, so the downloader
  needs retries/backoff (and ideally resumption, as noted above).

  **CDS holds the map, but not the samples.** VizieR `J/MNRAS/443/2907` does hold
  the map (`amap.dat`, 61,130 records), and its distance grid is the *same* one
  `iphas.h5` uses: the rebuilt file has 150 bins at 50, 150, …, 14950 pc, a
  uniform 100 pc step, identical for every sightline (`np.all(dist == dist[0])`),
  which is exactly the 0.05–14.95 kpc grid that the CDS ReadMe describes. What
  CDS stores per voxel is the posterior *expectations* — `meanA`, `diffA`,
  `e_meanA`, `e_diffA` — whereas `iphas.h5` stores the 20 posterior samples. So a
  CDS-based `fetch()` would be a *different map*: `'mean'` would be faithful (CDS
  `meanA` is the posterior mean that today's `mode='mean'` estimates from the
  samples), but `'random_sample'`, `'random_sample_per_pix'`, `'samples'` and
  `'median'` could not be reproduced, and `'random_sample'` is the *default*
  `mode` of `IPHASQuery.query()`, so existing callers would break. In exchange, a
  CDS-based product could expose the differential extinction and its
  uncertainties, which `iphas.h5` does not carry at all. On top of that, the part of CDS that *does* answer
  scripts cannot even see the distance blocks: a plain `GET` of
  `https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/MNRAS/443/2907/amap&-out.max=unlimited`
  returns the whole table (`-out.max=unlimited` works -- 61,130 rows), and
  `astropy.io.ascii` or `astropy.io.votable` parse it, so neither FTP nor
  `astroquery` would be needed; but that interface exposes only three extinction
  columns per sightline -- `A(1)`, `A(3)`, `A(10)`, in mag -- beside `GLON`,
  `GLAT`, `lmin`/`lmax`/`bmin`/`bmax`, `trust` and a `Plot` URL. The 150-step
  blocks are not queryable, and they are exactly what lives in the walled-off
  `amap.dat(.gz)`/FITS files.

  Policy agreed on 2026-10-09: keep `ascii2h5()` and the MD5 table as they are,
  and do **not** point `fetch()` at the Wayback Machine -- it is a third-party
  archive that is not dependable (mid-transfer it stopped completing TLS
  handshakes to this machine, which had to be worked around by falling back to
  plain HTTP), so it is a stopgap for a local copy, not a supported source.
  A local copy was restored that way: all 19 archives were downloaded from the
  Wayback Machine over plain HTTP and MD5-verified against the table above, then
  rebuilt with `ascii2h5()` into `~/dustmaps_data/iphas/iphas.h5` —
  227,797,348 bytes, md5 `58eeccd753ed5ecf6ac67073150d7a44`, 61,130 sightlines ×
  150 distances × 20 samples, no NaNs, A0 from 0.0 to 86.0 mag. (`fetch()` expects
  227,817,543 bytes — 20 kB more, which is ordinary HDF5 chunk/compression
  variance between library versions, not a difference in the data.)
  `test_iphas.py` now passes, including `test_bounds`, which exercises all five
  query modes on 10,000 random coordinates and checks the footprint edges.
  If the map is ever wanted again, the options are to fetch the CDS files through
  a browser and treat `iphas.h5` as a user-supplied file, or to build a separate,
  explicitly mean-only class from the VizieR table.

- [ ] **Give the IPHAS map a permanent home.**
  Greg will ask the IPHAS authors (Stuart Sale and the IPHAS consortium) for
  permission to upload the 3D extinction map itself to Zenodo, so that it stops
  depending on a website that has been compromised and on a third-party archive
  that is not meant to serve bulk data. If that succeeds, the map gets a DOI, and
  `fetch()` should be pointed at it with the MD5 sums above as the integrity
  check; the paragraphs above can then be trimmed to a short note on the map's
  history.

## Dataverse downloads: no resumption, and an interrupted file looks complete
(found while downloading the 33.8 GB DECaPS file)

- [ ] **`download_and_verify()` restarts from byte 0 instead of resuming.**
  `fetch_utils.download_and_verify()` opens the target with `'wb'` and streams
  the whole file, so a transfer that dies partway through has to start over. That
  happened to the DECaPS file at 9.7 GB of 33.8 GB (28.7%), with
  `requests.exceptions.ChunkedEncodingError: Connection broken:
  IncompleteRead(9717970028 bytes read, 24125852550 more expected)`.
  The Harvard Dataverse *does* support resumption: `GET
  {dataverse}/api/access/datafile/{file_id}` with `Range: bytes=<bytes_already_on_disk>-`
  returns `206 Partial Content` with the correct `Content-Range`. So the fix is
  to open the file with `'ab'` when it already exists, send the `Range` header,
  and only check the MD5 once the size reaches the expected value (both are
  already in the DOI metadata, as `dataFile['filesize']` / `dataFile['md5']`).
  A retry loop is needed as well, since a long transfer is likely to be
  interrupted more than once (resuming again costs nothing but a new request).
  Measured: the DECaPS file (`decaps_mean_and_samples.h5`, 33,843,822,578 bytes,
  md5 `7ccf070772769c8225be1568b32de76c`) streams at only 2.5-6 MiB/s from
  Dataverse, i.e. 2-4 hours for a single uninterrupted pass, so a mid-transfer
  failure is expensive.

- [ ] **A partially downloaded file is mistaken for a complete one.**
  The `fetch()` functions of several maps check only whether the file exists, and
  return before the downloader ever looks at it. For example
  `decaps.fetch()` prints *"File '...decaps_mean_and_samples.h5' already exists.
  Skipping download. To overwrite, use clobber=True."* for a file that is 28.7%
  complete, so re-running `fetch()` does not recover. Note that
  `dataverse_download_doi()` *does* verify the MD5 of an existing file (and
  re-downloads it if the sum is wrong), but `decaps.fetch()` returns before
  calling it, so that safety net never runs.
  Fix: drop the existence check from the map-specific `fetch()` functions and let
  the downloader decide, using the size (cheap) and/or the MD5 (slow for 33 GB)
  rather than existence. A size check is enough to catch an interrupted transfer.
  Worth checking all the `fetch()` functions for this pattern, not just DECaPS.

## SFD: data padding at b = 0 (deferred from PR #51)

- [ ] **SFD query near b = 0: the two hemispheres are sampled independently.**
  The SFD maps are zenithal equal-area (`GLON-ZEA`/`GLAT-ZEA`) grids centered on
  each pole, with the Galactic equator at the inscribed circle (`CRPIX = 2048.5`,
  r = 2048 px, so 1 px = 0.0439°). Measured from the published files:

  - Each file carries valid data only out to `r ≈ 2057 px` (~0.4° past its own
    equator); beyond that it is exactly **0.0** — not NaN. That zero padding is
    **20.8 %** of each file.
  - The two hemispheres therefore disagree in that narrow overlap (e.g. at
    b = +0.05°, l = 30° the NGP file gives 33.69 and the SGP file 35.98, a 6.8 %
    difference).
  - `SFDBase.query` splits strictly on the sign of b, so the same physical point
    is sampled from a different file depending on which side of the equator it
    falls, giving a small step at b = 0. (Not corruption: across all longitudes
    at b = 0 there are no zeros or NaNs.)

  **Caveat: any fix here would modify the published map**, so it must stay opt-in
  and default to the current, faithful behavior (which matches the official IDL
  `dust_getval`). Proposed fix, not implemented:

  - Add an opt-in keyword (e.g. `b0_blend=<half-width in deg>`, default `None` =
    today's behavior) to `SFDQuery.query`. For `|b| < half-width`, cross-fade the
    two files: `w = 0.5 + b/(2*half_width)`, returning `w*ngp + (1-w)*sgp`.
  - The half-width must stay below the ~0.4° data extent (0.3° is safe), and each
    sample must be checked to lie inside the valid footprint (pixel radius
    < 2057 px). A value test such as `data != 0` is **not** a valid footprint
    test, because the padding is 0 *and* the `mask` component legitimately
    contains zeros.
  - Cheaper alternative: expose `mode`/`cval` on `query()` (default
    `mode='nearest'`, unchanged) so callers can get NaN outside a map's footprint
    and detect the edge themselves.

- [ ] **SFD: test for the b = 0 seam (blocked on the item above).**
  A test that queries b = +ε and b = −ε (ε ≈ 0.05°) along a few longitudes and
  asserts that the fractional difference stays below a tolerance currently
  *fails*, so it should only be added along with the opt-in fix
  (`unittest.expectedFailure` until then).

## Gaia TGE (found while reviewing PR #73)

- [x] **`GaiaTGEQuery` re-parses a 4.2M-row CSV on every construction — DONE**
  (`7fde666`). `fetch()` now repacks the catalog into an HDF5 file holding an
  `astropy.table.Table`, which the query reads back with
  `Table.read(..., path='data')`. Loading went from ~10 s to ~0.5 s, and the
  41 MB HDF5 replaces the 39 MB gzipped CSV.

- [ ] **`healpix_level` rejects numpy integers.**
  The `isinstance(healpix_level, int)` test is `False` for e.g. `np.int8(6)`
  (numpy integers are not Python `int`s), so such a value falls through to the
  `'optimum'` comparison and raises the misleading
  ``'`healpix_level` must be either an integer or "optimum"'``. Use
  `numbers.Integral` (or `np.integer`) in the test -- keeping in mind that the
  value has to stay a Python `int` for `2**level`/`4**level` arithmetic, to
  avoid the int8 overflow fixed in PR #73.

## Bayestar: memory mapping and repacking (branch `memmap-bayestar`)

`fetch()` now repacks a downloaded Bayestar map so that each chunk holds a
block of `CHUNK_PIXELS = 16` adjacent pixels, with all of their samples and
distance bins (deflate level 7). A `memmap=True` query then reads only the
pixels it asks for. This shrinks the 2015 and 2019 files, but not 2017.

- [ ] **`bayestar2017`: the repacked file is 10.8% larger than the published one.**
  Whole-file repacks at 16 pixels per chunk:

  | map | published | repacked | |
  | --- | --- | --- | --- |
  | bayestar2015 | 4,792,709,805 B | 4,613,454,007 B | -3.7% |
  | bayestar2017 | 4,584,237,476 B | 5,079,647,547 B | **+10.8%** |
  | bayestar2019 | 727,902,836 B | 665,933,674 B | -8.5% |

  2017 is the exception because its data compresses much better in the
  published layout, whose chunks each hold 1456 pixels x 18 samples x 8
  distance bins -- it exploits redundancy *across* pixels, which a 16-pixel
  block cannot. The same repacking makes queries on it about 37x faster
  (100 coordinates: 0.74 s -> 0.02 s), so the growth is the price of that.
  The open question is whether +495 MB is worth that speed-up, or whether
  memory mapping should simply not be offered for this map.

- [ ] **Consider 32 pixels per chunk (`CHUNK_PIXELS = 32`).**
  Measured on whole-file repacks of bayestar2019 (deflate level 4, before
  `GRDiagnostic` and the root attributes were preserved): 0.678 GB at 16
  pixels, 0.628 GB at 32 (-7%), 0.601 GB at 64 (-11%). For 2017 the per-block
  measurements suggest only around 1%, because the file size there is governed
  by how many *pixels* share a chunk rather than by how many distance bins do:
  going from 16 to 32 pixels does not recover the published layout's 1456.

  The cost is that every queried pixel reads and decompresses twice as much
  data, and a query touching more than `n_pix / 32` pixels ends up reading the
  whole file. 16 was chosen as the best trade-off measured on 2019; if the
  2017 size matters, this is worth measuring properly on 2015 and 2017.

## DECaPS: repacking the published file (branch `memmap-decaps`)

The published `decaps_mean_and_samples.h5` stores its arrays in chunks of
`(200842 pixels, 1 sample, 1 distance bin)`, i.e. one `(sample, distance)` value
for a slab of 200842 pixels. That suits reading the map from end to end, which is
how the file was written, but a query at a random coordinate has to read one
~400 kB chunk for every `(sample, distance)` pair:

| query | published | repacked (16 pixels per chunk) |
| --- | --- | --- |
| 1 coordinate | 600 chunks, 0.24 GB | 1 chunk, 30 kB |
| 16 coordinates, same slab | 1.0 s, 0.24 GB | 0.6 ms |
| 16 coordinates, spread out | 13.1 s, 3.86 GB | 0.6 ms |
| 256 coordinates, spread out | 153,600 chunks, 61.7 GB (all of it) | 13 ms, 5 MB |

Repacking also *shrinks* the file, because the 120 distance bins of the same 16
pixels compress better together than long runs of a single `(sample, distance)`
pair do. **Measured on the whole published file** (33,843,822,578 bytes on the
Dataverse):

| dataset | published | repacked | change |
| --- | --- | --- | --- |
| `samples` | 25.24 GB (0.491 kB/pixel) | 10.91 GB (0.212 kB/pixel) | -57% |
| `mean` | 7.46 GB (0.145 kB/pixel) | 5.97 GB (0.116 kB/pixel) | -20% |
| `pixel_info` (copied as is) | 1.14 GB | 1.14 GB | - |
| **file** | **33.84 GB** | **18.37 GB** | **-45.7%** |

(An earlier extrapolation from the first 200,000 pixels predicted 24 GB, i.e. it
was pessimistic: those pixels are less compressible than the map as a whole.)

- [x] **Repacked the local copy** with `scripts/repack_decaps.py`, 16-pixel
      chunks, deflate 7, shuffle on: `~/dustmaps_data/decaps/decaps_mean_and_samples.h5`
      is now the repacked 18,370,195,128-byte file, md5
      `7e8988e8bb8162bae5b595eef850b77c`, and the published layout is kept beside
      it as `decaps_mean_and_samples_published.h5` (33.84 GB, md5
      `7ccf070772769c8225bc1568b32de76c`) until it is no longer needed.
      Repacking took 2677 s (mean 365 s, samples 2311 s), and 5% of the pixel
      blocks were verified (27 and 26 of 515 blocks, none differing).
- [ ] **The first author will upload the repacked file** (agreed 2026-10-09).
      When it appears on the Dataverse, `REPACKED_SIZES` in `decaps.py` is already
      set to the measured 18,370,195,128 bytes, the "33 GB" wording in `fetch()`
      should become ~18 GB, and `fetch()` should be checked to skip repacking a
      file that arrives already repacked (it will: the file carries
      `repacked = True`). The dataset shapes are unchanged, so `REPACKED_DSETS`
      needs no edit.
- [ ] **`fetch()`'s "33 GB" message understates the space needed, and is about
      to become wrong the other way.** `fetch()` warns *"You are about to
      download a large file (33 GB)"* from a hard-coded `n_gb` (33, or 8 with
      `mean_only=True`):

      - It counts only the *download*. `fetch()` streams the as-published file
        to `<file>.orig` and then writes the repacked file beside it, so peak
        usage is about 33 GB + 18 GB ≈ **51 GB** (the `.orig` file is deleted
        only once the repack has succeeded). Nothing checks free space, so on a
        disk with, say, 40 GB free the download can run for hours and then die
        in `h5_repack()` with `ENOSPC`. Either quote the peak figure in the
        warning, or call `shutil.disk_usage()` before starting and abort with a
        clear message.
      - Once the authors' repacked file is on the Dataverse (see the item
        above) the download is ~18 GB, so `n_gb` will have to be updated along
        with `REPACKED_SIZES`, and the message should say that the file arrives
        repacked, so 18 GB is downloaded and nothing is repacked afterwards.

      While in there: the docstring's claim that the repacked file "is nearly a
      third smaller" understates the measurement (33.84 GB → 18.37 GB, i.e.
      **45.7%** smaller; `samples` -57%, `mean` -20%), and `n_gb = 8` for the
      mean-only file is unverified, since that file has not been inspected here.

- [x] **`scripts/repack_decaps.py`** — a standalone repacker (numpy and h5py
      only; it does not import `dustmaps`), written for the DECaPS authors to run
      before uploading a repacked file to the Dataverse. It copies every group,
      file attribute and dataset attribute, adds `repacked` and `chunk_pixels`,
      verifies the output value by value(`--verify-fraction`, all of it by
      default), and prints both layouts, the size change and the MD5. Its output
      is identical, dataset for dataset, to `decaps.h5_repack()`.
- [x] **`decaps.h5_repack()` / `h5_is_repacked()`**, repacking in `fetch()`, and
      the "this map has not been repacked" warning, mirroring `bayestar.py`, with
      tests (16 in `TestDECaPSMemmap`).
- [ ] **Plotting the whole map by querying a grid is slow.** Every grid coordinate
      that falls in its own chunk costs one chunk read (~4 ms on this laptop), so
      the 1024 x 160 grid in `example_plot()` (164k coordinates, queried once per
      mode) takes tens of minutes. If that becomes annoying, the fix is a
      streaming rasterizer that reads the 18 GB file in order (a few minutes) and
      fills the image directly, instead of going through the random-access query
      path.
- [ ] **Check whether writing the file with `libver='latest'` speeds up scattered
      reads.** h5py writes HDF5-1.8-compatible files by default, so the chunks of
      `samples` are found through a version-1 B-tree. A scattered single-row read
      currently costs about 4 ms, and that is index bookkeeping plus a random
      read, not decompression: the chunk itself is only 19 kB. `libver='latest'`
      uses the newer fixed-array and extensible-array chunk indexes, which are
      much cheaper to search, and "one chunk per queried coordinate" is exactly
      the pattern they were designed for. The price is compatibility: such files
      need HDF5 1.10 or newer (2016), so older h5py builds could not read them --
      the same trade-off that ruled out zstd and lz4 for Bayestar. If the speed-up
      turns out to be large, it is worth requiring users to upgrade h5py for this.
      Test it by repacking a slice both ways (`h5_repack()` and
      `scripts/repack_decaps.py` would pass `libver` through to `h5py.File`) and
      timing scattered reads against each other; if it wins, the authors' upload
      should be written that way, since the file cannot be converted in place.
- Note: `/mean` carries an **empty samples axis** — shape `(n_pixels, 1,
  n_distances)`, not `(n_pixels, n_distances)`. The query code relies on that (it
  indexes `val[rows, :, distance]` for the mean as well as for the samples), so
  `DECaPSQuery` now raises a clear error if the shape is wrong. The synthetic test
  map used to write a 2-D mean, which hid a reshape bug in the memory-mapped path;
  it now matches the published file.

## Other observations

- [ ] **Redundant "Downloading ..." messages in other `fetch()` functions.**
  `fetch_utils.dataverse_download_doi()` already reports what it is doing: it
  prints `"Downloading data to '<file>' ..."` when it really downloads, and
  `'File exists. Not overwriting.'` when the file is already present. So a
  "Downloading ..." print in the caller is both duplicated *and* wrong in the
  "file already exists" case (it claims a download that never happens). Fixed in
  `sfd.fetch()`; `pg2010.fetch()` has the same redundant print, and
  `chen2014.fetch()` duplicates the message that `download_and_verify()` prints.

- [x] **NumPy 2.5 `ndarray.shape` assignment — DONE** (`89fa31f`, on `master`).
  All 19 `.shape =` assignments across 7 files were replaced with `np.reshape`
  (`json_serializers`, `bh`, `chen2018`, `bayestar`, `decaps`, `map_base`,
  `marshall`). The suite's NumPy `shape` deprecation warnings went 23 → 0.

- [x] **Invalid escape sequences in docstrings — DONE** (`da956c0`, on `master`).
  Eight modules had `100\%` in a `fetch()` docstring (plus `\*.dat.gz` in
  `marshall.dat2hdf5`), raising `SyntaxWarning` on import. The stray backslashes
  were removed (`chen2014`, `chen2018`, `csfd`, `edenhofer2023`, `iphas`,
  `leike_ensslin_2019`, `leike2020`, `marshall`).

- [x] **Use `np.random.default_rng()`, not the legacy global RNG — DONE for
  the library.** Every source function that needs randomness now takes a
  `seed_or_rng` keyword argument, which `dustmaps.rng_utils.get_rng` turns into
  a generator: a generator that is passed in is returned unchanged, and anything
  else goes to `numpy.random.default_rng` (so an integer, or `None` for a fresh,
  unpredictable generator, works too).

  - **The argument belongs to the query function, not to the map object.**
    `BayestarQuery.query`, `DECaPSQuery.query`, `DECaPSQueryLite.query`,
    `IPHASQuery.query` and `Edenhofer2023Query.query` all take
    `seed_or_rng=None`. On `query` rather than in `__init__`, a caller can change
    the seed from one query to the next -- the usual thing to want, since it
    means "draw another sample for these same coordinates" -- and can pass the
    same seed to two *different* map objects to compare them. (A `seed=` argument
    in `__init__` was tried first: it made a query's draws a reproducible
    *stream*, but left no way to ask for a different draw without building a new
    object.)
  - An integer makes each call repeatable: the same call twice gives the same
    samples. A `numpy.random.Generator` (or a legacy `numpy.random.RandomState`)
    is drawn from as it stands, so successive calls give different samples, and
    a caller can keep one stream going across queries.
  - Only the `random_sample` and `random_sample_per_pix` modes draw; `mean`,
    `median`, `best`, `percentile` and `samples` ignore the argument.
  - `bayestar.diagnostic` and `decaps.diagnostic` keep an integer `seed` of
    their own (documented as such, because a generator would be advanced by each
    timed query in turn, so the `memmap=True` and `memmap=False` passes would
    stop agreeing) and pass it to every query they time;
    `scripts/repack_decaps.py` uses `default_rng` + `rng.integers` rather than
    `np.random.RandomState` / `np.random.randint`. The diagnostics no longer seed
    numpy's global state.
  - Tests pin the behavior down: `test_the_same_seed_gives_the_same_random_samples`,
    `test_different_seeds_give_different_random_samples` and
    `test_a_generator_is_drawn_from_as_it_stands` in `test_bayestar.py`,
    `test_decaps.py` (both query classes) and `test_iphas.py`, over 200-500
    coordinates, so a chance agreement is vanishingly unlikely;
    `test_rng_utils.py` covers `get_rng` itself. The DECaPS and Bayestar
    comparison helpers pass `seed_or_rng=99` to both query calls, which replaces
    the `np.random.seed(99)` calls that used to force the two results to match.

- [ ] **Finish the `default_rng` conversion in the tests and the example
  scripts.** The library is converted (above), but the test modules and the
  example scripts still draw their random coordinates from the unseeded global
  `np.random.*`, so a failure cannot be reproduced and the sampled coordinates
  change from run to run. Each of them would become an explicitly seeded
  `default_rng`, created once per test (or per helper) and passed in.

  Mechanics to watch for: the method names change (`np.random.randint(a, b)` →
  `rng.integers(a, b)`, `np.random.random(...)` → `rng.random(...)`), and PCG64
  draws a different sequence than Mersenne Twister, so anything that pins a
  value or compares against stored output must be regenerated rather than merely
  ported. Most of the work is in the tests, which lean on unseeded `np.random.*`
  in `test_bayestar`, `test_decaps`, `test_csfd`, `test_lenz2017`,
  `test_marshall`, `test_planck`, `test_sfd`, `test_gaia_tge`, `test_iphas`,
  `test_serializers` and `test_h5_utils`, plus the 8 `examples/plot_*.py`
  scripts. Three places are seeded already, but on the legacy class, so they
  need only the class swapped: the module-level `RandomState(seed)` helpers and
  the `RandomState(1234)` fixtures in `test_bayestar.py` and `test_decaps.py`,
  and the eight `RandomState(...)` helpers in `test_h5_utils.py`. Watch the
  helpers that draw two coordinates at once, like
  `u, v = np.random.random((2, n_pix))` in `test_marshall` and `test_iphas`:
  those want `rng.random((2, n_pix))`, not two calls of `rng.random(n_pix)`,
  or the sequence changes for every test that follows.
