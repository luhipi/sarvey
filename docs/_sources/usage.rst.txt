.. _usage:

=====
Usage
=====

.. image:: https://seafile.projekt.uni-hannover.de/f/39209355cabc4607bf0a/?dl=1
   :alt: SARvey workflow
   :width: 600px
   :align: center

Processing workflow for using the SARvey software to derive displacement time series. The minimal required processing
steps and datasets are depicted in grey. All other steps are optional.



Command-line tools
------------------

The following command-line tools are available and can be run directly in the terminal.

`sarvey`
    A tool to derive displacements from the SLC stack with Multi-Temporal InSAR (MTI).
    A detailed description of the processing steps is given `here <processing.html>`_.

`sarvey_plot`
    A tool to plot the results from `sarvey` processing.

`sarvey_export`
    A tool to export the results from `sarvey` processing to shapefile or geopackage.

`sarvey_mask`
    A tool to create a mask from shapefile containing the area of interest, which can be used in `sarvey` processing.
    The tool reads from an input file, which is a shapefile or geopackage containing the geographic data.
    It supports both 'LineString' and 'Polygon' geometries.
    The tool first gets the spatial extent of the geographic data and searches the location of the polygon/line nodes in the image coordinates of the radar image.
    A buffer around the polygon/line is created specified by a width in pixel.
    The buffer is then used to create the mask.

    Here is an example of how to use the `sarvey_mask` tool:

    .. code-block:: bash

        sarvey_mask --input_file my_shapefile.shp --geom_file ./inputs/geometryRadar.h5 --out_file_name my_mask.h5 --width 5



`sarvey_osm`
    A tool to download OpenStreetMap data for the area of interest specified by the spatial extend of the SLC stack.
    The tool first gets the spatial extent of the SAR image from the geometry file.
    It then uses this spatial extent to download the OpenStreetMap data for the corresponding area.
    The download of railway tracks, highways and bridges is supported.
    After downloading the data, the tool saves it to a shapefile.

    After downloading the OpenStreetMap data with `sarvey_osm`, you can use the `sarvey_mask` tool to create a mask from the shapefile.

    Here is an example of how to use the `sarvey_osm` tool:

    .. code-block:: bash

      sarvey_osm --geom ./geometryRadar.h5 --railway                       # download railway
      sarvey_osm --geom ./geometryRadar.h5 --highway                       # download highway
      sarvey_osm --geom ./geometryRadar.h5 --railway --bridge              # download railway bridge
      sarvey_osm --geom ./geometryRadar.h5 --railway -o mask_railway.shp   # specify output path


Usage of the Python API
-----------------------

To use SARvey in a project:

    .. code-block:: python

        import sarvey

