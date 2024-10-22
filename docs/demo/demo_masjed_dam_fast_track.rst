.. _demo_masjed_dam_fast_track:

Fast Track Guide for Masjed Soleyman Dam
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are an advanced user, you can proceed with this fast track tutorial. If you prefer a more detailed, step-by-step guide, please refer to the :ref:`detailed guide <demo_masjed_dam_detailed_guide>` for this example.

.. note::

    These instructions are based on SARvey version 1.0.0 (Strawberry Pie). Newer versions may differ slightly.


Download the Data
"""""""""""""""""

In this tutorial, a processed stack of data is provided. If you wish to generate data for other areas, please refer to the :ref:`preparation` section.

.. code-block:: bash

    wget https://zenodo.org/records/12189041/files/SARvey_input_data_Masjed_Soleyman_dam_S1_dsc_2015_2018.zip
    unzip SARvey_input_data_Masjed_Soleyman_dam_S1_dsc_2015_2018.zip
    cd SARvey_input_data_Masjed_Soleyman_dam_S1_dsc_2015_2018


Activate SARvey environment
"""""""""""""""""""""""""""

.. code-block:: bash

    conda activate sarvey


Create a Config File
""""""""""""""""""""

.. code-block:: bash

    sarvey -f config.json 0 0 -g

Specify parameters in the config file. Set a reasonable value for **num_cores**.

Run **SARvey**
""""""""""""""

You can run each step individually or a range of steps by specifying the first and last step.

.. code-block:: bash

    sarvey -f config.json 0 4

Check Outputs
"""""""""""""

First, check the output snapshots in the `outputs/pics` directory. You can also use **`sarvey_plot`** to plot various products to assess the quality of the results and decide how to adjust parameters.  Modify the parameters in the config file and rerun the corresponding steps of `sarvey` to improve the results. For instance, changing **`coherence_p2`** from 0.8 to 0.7 and rerunning step 4 can increase the density of the final set of points. However, be cautious that reducing the value too much may include noisy points of low quality in the analysis, potentially leading to poor final results. You can check the details of all parameters using the -p flag in `sarvey` and decide how to tune them. For more explanations, please refer to :ref:`processing`



Plot Time Series Results
""""""""""""""""""""""""

The final products, including the time series, are stored in the coh\*\*_ts.h5 file. The file is named based on the coherence_p2 parameter you used. Plot the time series using the following command:

.. code-block:: bash

    sarvey_plot outputs/p2_coh80_ts.h5 -t

You can visualize velocity and DEM error estimation of second-order points. You can also visualize amplitude, DEM, or temporal coherence as the background. Right-click on any point to see its time series.

.. description of time series options to be added.




Export to GIS Format
""""""""""""""""""""

Export the data to Shapefiles using the following command.


.. code-block:: bash

    sarvey_export outputs/p2_coh80_ts.h5 -o outputs/shp/p2_coh80_ts.shp

You can visualize the data in any GIS software. If you use QGIS, you can use the `PS Time Series Viewer <https://plugins.qgis.org/plugins/pstimeseries/>`_ plugin to draw the time series.



Validate Your Results
"""""""""""""""""""""

You can download a copy of the final SARvey products from `this link <https://doi.org/10.5281/zenodo.12189041>`_. Use these files to compare your results and ensure everything worked correctly.

