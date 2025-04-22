=======
History
=======

Future major version
---------------------

* Export data to GIS format in WGS84 coordinates.
* Change file name coordinates_UTM.h5 to coordinates_map.h5.
* Use Transverse Mercator instead of UTM as map coordinates.

Future minor version (release soon)
-----------------------------------

* Update CI docker builder
* Update runner to test installation
* Update documentation with new instruction for installation including pip
* Fix numerical problems when computing grid size.

1.2.0 (2025-02-19)
------------------

* Create the background map and coordinates file each run of step 0.
* Visualize time series of neighbouring points in sarvey_plot -t.
* Ensure that specified grid size is bigger than study area.
* Update runner.
* Visualize amplitude images and interferograms interactively with sarvey_plot -i.

1.1.0 (2024-11-06)
------------------

* Use Scientific colour maps from Crameri.

1.0.0 (2024-08-12) Strawberry Pie
---------------------------------

* First release version on github.
* Change name of files for second-order points from coh* to p2_coh*.
* Check existence of intermediate results before continuing processing.
* Improve parameter names in config.
* Combine all general settings into one section in config.
* Allow adding user comments in config.json file.
* Improve documentation.
* Adapt CI from gitlab to github.
* Mask mean amplitude to avoid zero division warning in log10.
* Set logging level to debug for log file.
