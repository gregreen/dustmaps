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

* **Reference**: `Schlegel, Finkbeiner & Davis (1998) <http://adsabs.harvard.edu/abs/1998ApJ...500..525S>`_
* **Recalibration**: `Schlafly & Finkbeiner (2011) <http://adsabs.harvard.edu/abs/2011ApJ...737..103S>`_


CSFD (Chiang 2023)
~~~~~~~~~~~~~~~~~~

"Corrected SFD," a 2D dust map based on a reanalysis of SFD, using
tomographically constrained templates from WISE galaxy density fields to remove
extragalactic contamination from the cosmic infrared background (CIB).

* **Reference**: `Chiang (2023) <https://ui.adsabs.harvard.edu/abs/2023arXiv230603926C/abstract>`_
* **Website**: `Project description <https://idv.sinica.edu.tw/ykchiang/CSFD.html>`_
* **Data**: `Zenodo <https://doi.org/10.5281/zenodo.8207159>`_


Gaia Total Galactic Extinction (2022)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A two-dimensional map of A0, the monochromatic extinction at 541.4 nm. The map
is based on extinction estimates for giants beyond 300 pc. The individual
exitnction estimates estimates were obtained by fitting Gaia BP/RP spectra,
parallaxes and G-band apparent magnitudes.

The map comes in multiple HEALPix levels (6 to 9). By default, an "optimum"
map is loaded, with an adaptive HEALPix level, based on the local number
of stars (at least 3 stars are required per pixel).

* **Reference**: Delchambre et al. (2022).


Lenz, Hensley & Doré (2017)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A two-dimensional map of dust reddening, covering 40% of the sky with a 16.1'
resolution. This map is derived from emission from low-velocity
(l.o.s. velocity < 90 km/s) HI, which is found to correlate much more strongly
with E(B-V) than emission from high-velocity HI. The underlying data comes from
the HI4PI Survey. This map reports E(B-V) in magnitudes.

* **Reference**: `Lenz, Hensley & Doré (2017) <https://arxiv.org/abs/1706.00011>`_.
* **See also**: `GitHub page <https://github.com/daniellenz/ebv_tools>`_.


Planck Collaboration (2013)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two-dimensional maps of dust reddening across the entire sky. The
`Planck Collaboration (2013) <http://adsabs.harvard.edu/abs/2014A%26A...571A..11P>`_
fits a modified blackbody dust emission model to the Planck and IRAS
far-infrared maps, and provides three different conversions to dust reddening.

The three maps provided by
`Planck Collaboration (2013) <http://adsabs.harvard.edu/abs/2014A%26A...571A..11P>`_
are based on:

#. τ\ :sub:`353`\ : dust optical depth at 353 GHz.
#. ℛ: thermal dust radiance.
#. A recommended extragalactic reddening estimate, based on thermal dust
   radiance, but with point sources removed.

* **Reference**: `Planck Collaboration (2013) <http://adsabs.harvard.edu/abs/2014A%26A...571A..11P>`_
* **Website**: `Planck Explanatory Supplement <https://wiki.cosmos.esa.int/planckpla/index.php/CMB_and_astrophysical_component_maps#The_.5Bmath.5DE.28B-V.29.5B.2Fmath.5D_map_for_extra-galactic_studies>`_


Planck Collaboration (2016; "GNILC")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two-dimensional maps of dust reddening across the entire sky, using the
generalized needlet internal linear combination (GNILC) method to separate
out Galactic dust emission from CIB anisotropies.

This map contains both reddening estimates and estimated uncertainties.

* **Reference**: `Planck Collaboration (2016) <https://ui.adsabs.harvard.edu/abs/2016A%26A...596A.109P/abstract>`_
* **Website**: `Planck Explanatory Supplement <https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/Foreground_maps#GNILC_thermal_dust_and_CIB_products>`_


Peek & Graves (2010)
~~~~~~~~~~~~~~~~~~~~

A correction to the SFD'98 dust map, based on color excess measurements of "standard crayons" -- spectroscopically selected passively evolving galaxies. The maps have an angular resolution of 4.5°, and have a 1σ uncertainty of 1.5 mmag in E(B-V). Subtract this map from SFD'98 to obtain the corrected E(B-V) reddening.

* **Reference**: `Peek & Graves (2010) <http://adsabs.harvard.edu/abs/2010ApJ...719..415P>`_


Burstein & Heiles
~~~~~~~~~~~~~~~~~

Primarily of historical interest, the
`Burstein & Heiles (1982) <http://adsabs.harvard.edu/abs/1982AJ.....87.1165B>`_
dust reddening maps are derived from HI column density and galaxy counts.

* **Reference**: `Burstein & Heiles (1982) <http://adsabs.harvard.edu/abs/1982AJ.....87.1165B>`_


Three-Dimensional Dust Maps
---------------------------


Bayestar
~~~~~~~~

A three-dimensional map of Milky Way dust reddening, covering the three quarters
of the sky north of a declination of -30°. The map is probabilistic. containing
samples of the reddening along each line of sight. The "Bayestar" dust map is
inferred from stellar photometry of 800 million stars observed by Pan-STARRS 1,
and 2MASS photometry for a quarter of the stars. The latest version of Bayestar
also makes use of *Gaia* DR2 parallaxes.

There are three versions of Bayestar, called *Bayestar19*, *Bayestar17* and
*Bayestar15* here. By default, :code:`dustmaps` will use the latest version,
Bayestar19, although the earlier versions of the map can be selected by providing
the keyword argument :code:`version='bayestar2017'` or :code:`version='bayestar2015'`
in routines such as :code:`dustmaps.bayestar.fetch`,
:code:`dustmaps.bayestar.BayestarQuery` and :code:`dustmaps.bayestar.BayestarWebQuery`.
If you want to make sure that your code will always use the same version of the
map, even as new versions of Bayestar are released, then set the :code:`version`
keyword explicitly.

The units of reddening used by each map are slightly different:

#. Bayestar19 reports reddening in an arbitrary unit that can be converted to
   extinction in different bands using the coefficients given in Table 1 of
   Green, Schlafly, Finkbeiner et al. (2019).
#. Bayestar17 reports reddening in an arbitrary unit that can be converted to
   extinction in different bands using the coefficients given in Table 1 of
   `Green, Schlafly, Finkbeiner et al. (2018) <http://adsabs.harvard.edu/abs/2018arXiv180103555G>`_.
#. Bayestar15 reports reddening in the same units as those used by SFD. Therefore,
   in order to convert Bayestar15 reddenings to extinction in different bands, one
   should use the conversions provided in
   `Table 6 of Schlafly & Finkbeiner (2011) <http://iopscience.iop.org/0004-637X/737/2/103/article#apj398709t6>`_.

* **References**: `Green, Schlafly, Finkbeiner et al. (2019) <https://ui.adsabs.harvard.edu/abs/2019ApJ...887...93G>`_,
  `Green, Schlafly, Finkbeiner et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018MNRAS.478..651G>`_
  and `Green, Schlafly, Finkbeiner et al. (2015) <https://ui.adsabs.harvard.edu/abs/2015ApJ...810...25G>`_.
* **Website**: `argonaut.skymaps.info <http://argonaut.skymaps.info>`_

DECaPS
~~~~~~~~

A three-dimensional map of Milky Way dust reddening, covering the southern Galactic plane (239° < l < 6°) in a 20°-thick strip. 
The map has an angular resolution of 1' (on par with Herschel), a typical maximum reliable distance of ~ 10 kpc, and is reliable up to A(V)~12 mag in extinction. 
The map is probabilistic. containing samples of the reddening along each line of sight. The "DECaPS" dust map is
inferred from stellar photometry of 709 million stars observed by the DECaPS2 survey, in combination with 2MASS, VVV, and unWISE. 
The DECaPS map also makes use of *Gaia* DR3 parallaxes where available. When combined with Bayestar, DECaPS completes the coverage
of the Galactic plane within a 20°-thick strip, enabling reddening estimates over the entire disk. 

Unlike the Bayestar map, whose reddening is reported in an arbitrary unit, the DECaPS reddening is reported in units of E(B-V) in mags

* **References**: Zucker, Saydjari, & Speagle et al. (2025)
* **Data**: `Dataverse <https://doi.org/10.7910/DVN/J9JCKO>`_


Chen et al. (2014)
~~~~~~~~~~~~~~~~~~

A three-dimensional map of dust extinction in the Galactic anticenter. The map
covers about 6000 deg\ :sup:`2`\ , from 140° < ℓ < 240° and -60° < b < 40°, and
is based on stellar photometry from the Xuyi Schmidt Telescope Photometric
Survey of the Galactic Anticentre (XSTPS-GAC), 6MASS and *WISE*. The map has an
angular resolution of 3 to 9 arcminutes, and reports *r*-band extinction, along
with Gaussian error estimates.

* **Reference**: `Chen et al. (2014) <http://adsabs.harvard.edu/abs/2014MNRAS.443.1192C>`_
* **Website**: `http://lamost973.pku.edu.cn <http://lamost973.pku.edu.cn/site/Photometric-Extinctions-and-Distances/>`_


Edenhofer et al. (2023)
~~~~~~~~~~~~~~~~~~~~~~~

A three-dimensional map of Milky Way dust extinction, with a Gaussian process
prior on the logarithm of the dust extinction density. The prior is implemented
on a spherical grid. The map starts at 69 pc and extends out to 1.25 kpc in
distance from the Sun. It has an angular resolution of 14' and a maximum
distance voxalization of 0.4 pc at 69 pc and a minimum distance voxalization of
7 pc at 1250 pc. The map is based on the stellar distance and extinction
estimates of Zhang, Green & Rix (2023), and therefore reports extinctions in
their units. Accompanying the main reconstruciton is an additional map that uses
less data but extends out to 2 kpc from the Sun.

* **Reference**: `Edenhofer et al. (2023) <https://ui.adsabs.harvard.edu/abs/2023arXiv230801295E/abstract>`_
* **Data**: `Zenodo <https://doi.org/10.5281/zenodo.8187943>`_


IPHAS
~~~~~

A three-dimensional map of Milky Way dust extinction, covering a 10°-thick strip
of the Galactic plane, between 30° < ℓ < 120°. The map is probabilistic,
containing samples of the cumulative extinction along each line of sight. The
map is based on IPHAS imaging of stars. The map returns A\ :sub:`0`\ , the
monochromatic extinction.

* **Reference**: `Sale et al. (2014) <http://adsabs.harvard.edu/abs/2014MNRAS.443.2907S>`_
* **Website**: `www.iphas.org/extinction <http://www.iphas.org/extinction/>`_


Leike & Enßlin (2019)
~~~~~~~~~~~~~~~~~~~~~~

A three-dimensional map of Milky Way dust extinction, incorporating a Gaussian
process prior on the log of the dust extinction density. The map is based on
the Gaia DR2 catalog parallaxes and G-band extinctions, and spans a (600 pc)³
box centered on the Sun.

* **Reference**: `Leike & Enßlin (2019) <https://ui.adsabs.harvard.edu/abs/2019arXiv190105971L/abstract>`_
* **Data**: `Zenodo <https://doi.org/10.5281/zenodo.2542807>`_


Leike, Glatzle & Enßlin (2020)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A three-dimensional map of Milky Way dust extinction, incorporating a Gaussian
process prior on the log of the dust extinction density, similar to Leike &
Enßlin (2019). The map is based on data from Gaia, 2MASS, Pan-STARRS 1 and
ALLWISE, and is calculated on a Cartesian grid spanning a
(740 pc)×(740 pc)×(540 pc) box (in Galactic *x*, *y* and *z*, respectively)
centered on the Sun.

* **References**: `Leike, Glatzle & Enßlin (2020) <https://ui.adsabs.harvard.edu/abs/2020A%26A...639A.138L/abstract>`_
* **Data**: `Zenodo <https://doi.org/10.5281/zenodo.3993082>`_


Marshall et al. (2006)
~~~~~~~~~~~~~~~~~~~~~~

A three-dimensional map of Milky Way dust extinction, covering a 20°-thick strip
of the Galactic plane, between -100° < ℓ < 100°. The map is contains 2MASS
K\ :sub:`s`\ -band extinctions with a Gaussian uncertainty estimates. The map is
based on a comparison of 2MASS colors of stars with expectations from the
Besançon model of the Galaxy.

* **Reference**: `Marshall et al. (2006) <http://adsabs.harvard.edu/abs/2006A%26A...453..635M>`_
* **Website**: `http://cds.u-strasbg.fr/ <http://cdsarc.u-strasbg.fr/viz-bin/qcat?J/A+A/453/635>`_
