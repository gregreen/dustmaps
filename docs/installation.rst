Installation
============

There are two ways to install :code:`dustmaps`.


1. Using :code:`pip`
--------------------

From the commandline, run

.. code-block :: bash

    pip install dustmaps

You may have to use :code:`sudo`.

Next, we'll configure the package and download the dust maps we'll want to use.
Start up a python interpreter and type:

.. code-block :: python
    
    from dustmaps.config import config
    config['data_dir'] = '/path/to/store/maps/in'
    
    import dustmaps.sfd
    dustmaps.sfd.fetch()
    
    import dustmaps.csfd
    dustmaps.csfd.fetch()
    
    import dustmaps.planck
    dustmaps.planck.fetch()
    
    import dustmaps.planck
    dustmaps.planck.fetch(which='GNILC')
    
    import dustmaps.bayestar
    dustmaps.bayestar.fetch()
    
    import dustmaps.iphas
    dustmaps.iphas.fetch()
    
    import dustmaps.marshall
    dustmaps.marshall.fetch()
    
    import dustmaps.chen2014
    dustmaps.chen2014.fetch()
    
    import dustmaps.lenz2017
    dustmaps.lenz2017.fetch()
    
    import dustmaps.pg2010
    dustmaps.pg2010.fetch()
    
    import dustmaps.leike_ensslin_2019
    dustmaps.leike_ensslin_2019.fetch()
    
    import dustmaps.leike2020
    dustmaps.leike2020.fetch()
    
    import dustmaps.edenhofer2023
    dustmaps.edenhofer2023.fetch()
    
    import dustmaps.gaia_tge
    dustmaps.gaia_tge.fetch()
        
    import dustmaps.decaps
    dustmaps.decaps.fetch()

All the dust maps should now be in the path you gave to
:code:`config['data_dir']`. Note that these dust maps can be very large - some
are several Gigabytes! Only download those you think you'll need.

Note that there are two versions of the Bayestar dust map. By default,
:code:`dustmaps.bayestar.fetch()` will download Bayestar19 (Green et al. 2019).
In order to download earlier version of the map (Green et al. 2015, 2018), you can
provide the keyword argument :code:`version='bayestar2017'` (Green et al. 2018) or
:code:`version='bayestar2015'` (Green et al. 2015).


2. Using :code:`setup.py`
-------------------------

An alternative way to download :code:`dustmaps`, if you don't want to use
:code:`pip`, is to download or clone the respository from
https://github.com/gregreen/dustmaps.


In this case, you will have to manually make sure that the dependencies are
satisfied:

* :code:`numpy`
* :code:`scipy`
* :code:`astropy`
* :code:`h5py`
* :code:`healpy`
* :code:`requests`
* :code:`six`
* :code:`progressbar2`
* :code:`tqdm`


These packages can typically be installed using the Python package manager,
:code:`pip`.

Once these dependencies are installed, run the following command from the root
directory of the :code:`dustmaps` package:

.. code-block :: bash
    
    python setup.py install --large-data-dir=/path/to/store/maps/in

Then, fetch the maps you'd like to use. Depending on which dust maps you choose
to download, this step can take up several Gigabytes of disk space. Be careful
to only download those you think you'll need:

.. code-block :: bash
    
    python setup.py fetch --map-name=sfd
    python setup.py fetch --map-name=csfd
    python setup.py fetch --map-name=planck
    python setup.py fetch --map-name=planckGNILC
    python setup.py fetch --map-name=bayestar
    python setup.py fetch --map-name=iphas
    python setup.py fetch --map-name=marshall
    python setup.py fetch --map-name=chen2014
    python setup.py fetch --map-name=lenz2017
    python setup.py fetch --map-name=leikeensslin2019
    python setup.py fetch --map-name=leike2020
    python setup.py fetch --map-name=edenhofer2023
    python setup.py fetch --map-name=gaia_tge
    python setup.py fetch --map-name=decaps


That's it!

Note that the above code will download the latest version of the Bayestar dust
map (the 2019 version). If you want to download the 2015 and 2017 versions, you
can enter the commands

.. code-block :: bash
    
    python setup.py fetch --map-name=bayestar2015
    python setup.py fetch --map-name=bayestar2017

3. Custom configuration file location (Optional)
------------------------------------------------

By default, a configuration file is stored in :code:`~/.dustmapsrc`. This 
file might look like the following::

    {"data_dir": "/path/to/store/maps/in"}

If you would like :code:`dustmaps` to use a different configuration file, 
then you can set the environmental variable :code:`DUSTMAPS_CONFIG_FNAME`. 
For example, in a :code:`bash` terminal,

.. code-block :: bash

    export DUSTMAPS_CONFIG_FNAME=/path/to/custom/config/file.json
    python script_using_dustmaps.py

The paths listed in the configuration file can also include environmental
variables, which will be expanded when :code:`dustmaps` is loaded. For example,
the configuration file could contain the following::

    {"data_dir": "/path/with/${VARIABLE}/included"}

If the environmental variable :code:`VARIABLE` is set to :code:`"foo"`,
for example, then :code:`dustmaps` will expand :code:`data_dir` to
:code:`"/path/with/foo/included"`.
