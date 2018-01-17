Available Dust Maps
===================


Two-Dimensional Dust Maps
-------------------------


SFD
~~~

A two-dimensional map of dust reddening across the entire sky. The "SFD" dust
map is based on far-infrared emission of dust. The authors model the temperature
and optical depth of the dust, and then calibrate a relationship between the
dust's far-infrared optical depth and optical reddening. This calibration was
updated by
`Schlafly & Finkbeiner (2011) <http://adsabs.harvard.edu/abs/2011ApJ...737..103S>`_.

In order to convert SFD values of E(B-V) to extinction, one should use the
conversions provided in
`Table 6 of Schlafly & Finkbeiner (2011) <http://iopscience.iop.org/0004-637X/737/2/103/article#apj398709t6>`_.

**Reference**: `Schlegel, Finkbeiner & Davis (1998) <http://adsabs.harvard.edu/abs/1998ApJ...500..525S>`_

**Recalibration**: `Schlafly & Finkbeiner (2011) <http://adsabs.harvard.edu/abs/2011ApJ...737..103S>`_


Planck
~~~~~~

Two-dimensional maps of dust reddening across the entire sky. The
`Planck Collaboration (2013) <http://adsabs.harvard.edu/abs/2014A%26A...571A..11P>`_
fits a modified blackbody dust emission model to the Planck and IRAS
far-infrared maps, and provides three different conversions to dust reddening.

The three maps provided by
`Planck Collaboration (2013) <http://adsabs.harvard.edu/abs/2014A%26A...571A..11P>`_
are based on:

1. τ\ :sub:`353`\ : dust optical depth at 353 GHz.
2. ℛ: thermal dust radiance.
3. A recommende extragalactic reddening estimate, based on thermal dust
   radiance, but with point sources removed.

**Reference**: `Planck Collaboration (2013) <http://adsabs.harvard.edu/abs/2014A%26A...571A..11P>`_

**Website**: `Planck Explanatory Supplement <https://wiki.cosmos.esa.int/planckpla/index.php/CMB_and_astrophysical_component_maps#The_.5Bmath.5DE.28B-V.29.5B.2Fmath.5D_map_for_extra-galactic_studies>`_


Burstein & Heiles
~~~~~~~~~~~~~~~~~

Primarily of historical interest, the
`Burstein & Heiles (1982) <http://adsabs.harvard.edu/abs/1982AJ.....87.1165B>`_
dust reddening maps are derived from HI column density and galaxy counts.

**Reference**: `Burstein & Heiles (1982) <http://adsabs.harvard.edu/abs/1982AJ.....87.1165B>`_


Three-Dimensional Dust Maps
---------------------------


Bayestar
~~~~~~~~

A three-dimensional map of Milky Way dust reddening, covering the three quarters
of the sky north of a declination of -30°. The map is probabilistic. containing
samples of the reddening along each line of sight. The "Bayestar" dust map is
inferred from stellar photometry of 800 million stars observed by Pan-STARRS 1,
and 2MASS photometry for a quarter of the stars.

There are two versions of Bayestar, called Bayestar17 and Bayestar15 here. By
default, :code:`dustmaps` will use the latest version, Bayestar17, although the
earlier version of the map can be selected by providing the keyword argument
:code:`version='bayestar2015'` in routines such as
:code:`dustmaps.bayestar.fetch`, :code:`dustmaps.bayestar.BayestarQuery` and
:code:`dustmaps.bayestar.BayestarWebQuery`.

Bayestar17 reports reddening in an arbitrary unit that can be converted to
extinction in different bands using the coefficients given in Table 1 of
`Green, Schlafly, Finkbeiner et al. (2018) <http://adsabs.harvard.edu/abs/2018arXiv180103555G>`_.

Bayestar15 reports reddening in the same units as those used by SFD. Therefore,
in order to convert Bayestar15 reddenings to extinction in different bands, one
should use the conversions provided in
`Table 6 of Schlafly & Finkbeiner (2011) <http://iopscience.iop.org/0004-637X/737/2/103/article#apj398709t6>`_.

**References**: `Green, Schlafly, Finkbeiner et al. (2018) <http://adsabs.harvard.edu/abs/2018arXiv180103555G>`_,
`Green, Schlafly, Finkbeiner et al. (2015) <http://adsabs.harvard.edu/abs/2015arXiv150701005G>`_

**Website**: `argonaut.skymaps.info <http://argonaut.skymaps.info>`_


Chen et al. (2014)
~~~~~~~~~~~~~~~~~~

A three-dimensional map of dust extinction in the Galactic anticenter. The map
covers about 6000 deg\ :sup:`2`\ , from 140° < ℓ < 240° and -60° < b < 40°, and
is based on stellar photometry from the Xuyi Schmidt Telescope Photometric
Survey of the Galactic Anticentre (XSTPS-GAC), 2MASS and *WISE*. The map has an
angular resolution of 3 to 9 arcminutes, and reports *r*-band extinction, along
with Gaussian error estimates.

**Reference**: `Chen et al. (2014) <http://adsabs.harvard.edu/abs/2014MNRAS.443.1192C>`_

**Website**: `http://lamost973.pku.edu.cn <http://lamost973.pku.edu.cn/site/Photometric-Extinctions-and-Distances/>`_


IPHAS
~~~~~

A three-dimensional map of Milky Way dust extinction, covering a 10°-thick strip
of the Galactic plane, between 30° < ℓ < 120°. The map is probabilistic,
containing samples of the cumulative extinction along each line of sight. The
map is based on IPHAS imaging of stars. The map returns A\ :sub:`0`\ , the
monochromatic extinction.

**Reference**: `Sale et al. (2014) <http://adsabs.harvard.edu/abs/2014MNRAS.443.2907S>`_

**Website**: `www.iphas.org/extinction <http://www.iphas.org/extinction/>`_


Marshall et al. (2006)
~~~~~~~~~~~~~~~~~~~~~~

A three-dimensional map of Milky Way dust extinction, covering a 10°-thick strip
of the Galactic plane, between -100° < ℓ < 100°. The map is contains 2MASS
K\ :sub:`s`\ -band extinctions with a Gaussian uncertainty estimates. The map is
based on a comparison of 2MASS colors of stars with expectations from the
Besançon model of the Galaxy.

**Reference**: `Marshall et al. (2006) <http://adsabs.harvard.edu/abs/2006A%26A...453..635M>`_

**Website**: `http://cds.u-strasbg.fr/ <http://cdsarc.u-strasbg.fr/viz-bin/qcat?J/A+A/453/635>`_
