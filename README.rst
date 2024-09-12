========================
SARvey - survey with SAR
========================

Open-source InSAR time series analysis software developed within the project SAR4Infra.
**SARvey** aims to analyze InSAR displacement time series for engineering applications.



Documentation
-------------
The documentation with installation instructions, processing steps, and examples with a demo dataset can be found at:
https://luhipi.github.io/sarvey/docs/



Status
------

.. image:: https://github.com/luhipi/sarvey/actions/workflows/ci.yml/badge.svg
        :target: https://github.com/luhipi/sarvey/actions
        :alt: Pipelines
.. image:: https://img.shields.io/static/v1?label=Documentation&message=GitHub%20Pages&color=blue
        :target: https://luhipi.github.io/sarvey/docs/
        :alt: Documentation
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.12544131.svg
        :target: https://doi.org/10.5281/zenodo.12544131
        :alt: DOI


See also the latest coverage_ report and the pytest_ HTML report.


License
-------

**SARvey** is distributed under the GNU General Public License, version 3 (GPLv3).

The following exceptions applies:

This package uses PyMaxFlow. The core of PyMaxflows library is the C++ implementation by Vladimir Kolmogorov. It is also licensed under the GPL, but it REQUIRES that you cite [BOYKOV04] in any resulting publication if you use this code for research purposes.
This requirement extends to **SARvey**.

Please check out the details of the license `here <LICENSE>`_.

How to cite
-----------

If you use **SARvey** in your research, please cite the following.

1. The paper describing the methodology:

   Piter, A., Haghshenas Haghighi, M., Motagh, M.(2024). An in-depth study on Sentinel-1 InSAR for transport infrastructure monitoring. PFG - Journal of Photogrammetry, Remote Sensing and Geoinformation Science. (paper currently under review).

2. The software itself. Please specify the version you use:

   Piter, A., Haghshenas Haghighi, M., FERN.Lab, & Motagh, M. (2024). SARvey - survey with SAR [version]. Zenodo. https://doi.org/10.5281/zenodo.12544131

3. If you use the PUMA method for unwrapping in your research, please cite the following publication as indicated in the license:

   An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Vision. Yuri Boykov and Vladimir Kolmogorov. In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), September 2004. `Link to paper <https://ieeexplore.ieee.org/document/1316848>`_.


Processing overview
-------------------


.. image:: https://seafile.projekt.uni-hannover.de/f/006f702937cd4e618bcb/?dl=1
   :width: 600
   :align: center
   :alt: SARvey workflow

Processing workflow for using the SARvey software to derive displacement time series.


SARvey is a command-line-based software. The major steps for running SARvey are the following:

* **Installation**

  SARvey is a cross-platform python-based software and can be installed on Linux and MacOS. On Windows, SARvey is tested on Windows Subsystem for Linux (WSL_) version 2.
  Details of installation can be found in `installation instruction`_.


* **Preprocessing**

  The software requires a coregistered stack of SLC and the related geometry information in the MiaplPy_  data format.
  The coregistered stack of SLC can be created using an InSAR processor. Currently MiaplPy_ only supports ISCE_. Support for GAMMA and SNAP_ is planned for future.
  After creating the coregistered stack of SLC, run the "load_data" step from Miaplpy_ to create the "inputs" directory which contains "slcStack.h5" and "geometryRadar.h5".
  Details are explained in the preparation_ section


* **Time series analysis**

  Time series analysis is performed using `sarvey`. It consists of 5 steps (steps 0 to 4). The details of each step are explained in `processing steps`_. The processing parameters are handled in a json config file. Visualization and export are handled by `sarvey_plot` and `sarvey_export` packages. Below are the major steps:

  * Go to your working directory:

    .. code-block:: bash

         cd path/to/working_dir/

  * Create a default config file using **"-g"** flag:

    .. code-block:: bash

         sarvey -f config.json 0 4 -g

  * Modify **config.json** to change path to "inputs" directory. Modify other parameters as desired.

  * Run all processing steps (steps 0 to 4):

    .. code-block:: bash

         sarvey -f config.json 0 4

    Different processing steps are explained in `processing`_ section.

  * Plot the resulting displacement time series:

    .. code-block:: bash

         sarvey_plot outputs/p2_coh80_ts.h5 -t

  * Export the results as Shapefiles_:

    .. code-block:: bash

         sarvey_export outputs/p2_coh80_ts.h5 -o outputs/shp/p2_coh80.shp


Feature overview
----------------

**SARvey** has three main components for processing, visualization, and exporting data.

* `sarvey` performs time series analysis.
* `sarvey_plot` plots the outputs.
* `sarvey_export` exports InSAR time series results from to GIS data formats. The GIS data format can be visualized for example in QGIS_.

It also has two components that facilitate transport infrastructure monitoring.

* `sarvey_mask` creates mask from Shapefiles, e.g. for transport infrastructures.
* `sarvey_osm` downloads transport infrastructure information from OSM_ and store as Shapefiles.

You can run each component in the command line with "-h" argument for more information about the usage. For example:

  .. code-block:: bash

       sarvey -h



**SARvey** supports two processing schemes:

* `Two-step unwrapping`_ with atmospheric correction (default).

* `One-step unwrapping`_ for a small area.

History / Changelog
-------------------

You can find the protocol of recent changes in the **SARvey** package
`history`_.

We follow the principle of semantic versioning.
The version number is structured as follows: MAJOR.MINOR.PATCH.
You can find a description of the versioning scheme `here <https://semver.org/>`__.

Credits
-------

This software was developed within the project SAR4Infra (2020-2024) with funds of the German Federal Ministry for Digital and Transport.
The project consortium consists of
the `Institute of Photogrammetry and GeoInformation`_ at Leibniz University Hannover,
`FERN.Lab`_ (innovation and technology transfer lab of the GFZ German Research Centre for Geosciences, Potsdam),
`Landesamt fuer Vermessung und Geoinformation Schleswig-Holstein`_,
and `Landesbetrieb Strassenbau und Verkehr Schleswig-Holstein`_.
The scientific and methodological development was carried out by Andreas Piter (piter@ipi.uni-hannover.de), supervised by Mahmud H. Haghighi (mahmud@ipi.uni-hannover.de) and Mahdi Motagh (motagh@gfz-potsdam.de).
The `FERN.Lab`_ (fernlab@gfz-potsdam.de) contributed to the development, documentation, continuous integration, and testing of the package.


This package was created with Cookiecutter_ and the `fernlab/cookiecutter-pypackage`_ project template.


.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`fernlab/cookiecutter-pypackage`: https://github.com/fernlab/cookiecutter-pypackage
.. _coverage: https://luhipi.github.io/sarvey/coverage/
.. _pytest: https://luhipi.github.io/sarvey/test_reports/report.html
.. _processing: https://luhipi.github.io/sarvey/docs/processing.html
.. _`processing steps`: https://luhipi.github.io/sarvey/docs/processing.html#processing-steps-for-two-step-unwrapping-workflow
.. _preparation: https://luhipi.github.io/sarvey/docs/preparation.html
.. _`Two-step unwrapping`: https://luhipi.github.io/sarvey/docs/processing.html#processing-steps-for-two-step-unwrapping-workflow
.. _`One-step unwrapping`: https://luhipi.github.io/sarvey/docs/processing.html#processing-steps-for-one-step-unwrapping-workflow
.. _`installation instruction`: https://luhipi.github.io/sarvey/docs/installation.html
.. _`history`: https://luhipi.github.io/sarvey/docs/history.html
.. _MiaplPy: https://github.com/insarlab/MiaplPy
.. _ISCE: https://github.com/isce-framework/isce2
.. _SNAP: https://step.esa.int/main/toolboxes/snap
.. _Shapefiles: https://doc.arcgis.com/en/arcgis-online/reference/shapefiles.htm
.. _QGIS: https://qgis.org/en/site/
.. _`PS Time Series Viewer`: https://plugins.qgis.org/plugins/pstimeseries/
.. _OSM: https://www.openstreetmap.org/
.. _WSL: https://learn.microsoft.com/en-us/windows/wsl/
.. _FERN.Lab: https://fernlab.gfz-potsdam.de/
.. _`Institute of Photogrammetry and GeoInformation`: https://www.ipi.uni-hannover.de/en/
.. _`Landesamt fuer Vermessung und Geoinformation Schleswig-Holstein`: https://www.schleswig-holstein.de/DE/landesregierung/ministerien-behoerden/LVERMGEOSH/lvermgeosh_node.html
.. _`Landesbetrieb Strassenbau und Verkehr Schleswig-Holstein`: https://www.schleswig-holstein.de/DE/Landesregierung/LBVSH/lbvsh_node.html
