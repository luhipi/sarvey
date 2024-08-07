.. _demo_masjed_dam_detailed_guide:

Detailed Guide for Masjed Soleyman Dam
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This tutorial provides a comprehensive guide to SARvey processing. If you are an advanced user, you can proceed directly to the :ref:`fast track for advanced users <demo_masjed_dam_fast_track>`.

.. note::

    This instruction is based on SARvey version 1.0.0 (Strawberry Pie). Newer versions may differ slightly.

Step 1: Before Running SARvey
"""""""""""""""""""""""""""""

Step 1.1: Download the Data
"""""""""""""""""""""""""""

Download the data by running the following commnad in the console:

.. code-block:: bash

    wget https://zenodo.org/records/12189041/files/SARvey_input_data_Masjed_Soleyman_dam_S1_dsc_2015_2018.zip


Unzip the downloaded file and change the directory.

.. code-block:: bash

    unzip SARvey_input_data_Masjed_Soleyman_dam_S1_dsc_2015_2018.zip
    cd SARvey_input_data_Masjed_Soleyman_dam_S1_dsc_2015_2018


Check the downloaded data using `info.py` and `view.py`. For example:

.. code-block:: bash

    info.py inputs/slcStack.h5

.. code-block:: bash

    view.py inputs/geometryRadar.h5


Step 1.2: Activate SARvey and Change Directory
"""""""""""""""""""""""""""""""""""""""""""""""

If you have not installed SARvey, refer to the `installation instructions <installation.html>`_. Activate the SARvey environment:

.. code-block:: bash

    conda activate sarvey

Ensure SARvey can be called from the console.

.. code-block:: bash

    sarvey -h

If you see the following command, it indicates that SARvey cannot be called. Ensure it is installed correctly and the conda environment is activated.

.. code-block:: none

    command not found: sarvey

Step 1.3: Create a Config File
""""""""""""""""""""""""""""""

Create a config file, which is a JSON file containing the parameters for `sarvey`. The config file can be created using the following command:

.. code-block:: bash

    sarvey -f config.json 0 0 -g

Note: The above command only generates a configuration file. Although step 0 is specified, it will not be executed.

Step 1.4: Modify the config.json File
"""""""""""""""""""""""""""""""""""""

1.4.1. Open the config.json file and check the parameters. The first parameters to specify in the config file are **input_path** and **output_path**. For this example dataset, the `slcStack.h5` and `geometryRadar.h5` files are in the `inputs/` directory, which is the default value in the config file. Therefore, you do not need to change it. The **output_path** should be `outputs/` for this example.

.. code-block:: json

    {
        "general": {
            "input_path": "inputs/",
            "output_path": "outputs/"
        }
        // other parameters
    }

1.4.2. Specify the **num_cores**. You can check the number of cores on your computer using the following commands.

In Linux, run:

.. code-block:: bash

    nproc --all

In MacOS, run:

.. code-block:: bash

    sysctl -n hw.ncpu

It is a good practice to specify a number lower than the number of available cores in the config file.

.. code-block:: json

    {
    // other parameters
        "general": {
        "num_cores": 5,
        // other parameters
        },
    //other parameters
    }



Step 2: Running SARvey
""""""""""""""""""""""

SARvey consists of five steps as detailed in :ref:`processing`. You can run all steps by specifying starting step `0` and ending step `4`. In this tutorial, however, we will run the steps separately as follows.

When running `sarvey`, if it finishes normally, you will see a message like the following in the command line:

.. code-block:: none

    2024-06-19 11:05:10,305 - INFO - MTI finished normally.

.. note::
    If you encounter an error, first read all the prompts in the console and carefully track all error and warning messages. If the issue is not clear from the console messages, check the log files stored in the directory specified in the config file. If the error persists and you need assistance, sharing the corresponding log file will help.


Step 2.0: Run Step 0 of SARvey: Preparation
'''''''''''''''''''''''''''''''''''''''''''

The first step creates an interferogram network and calculates the temporal coherence for all pixels. Run the following command:

.. code-block:: bash

    sarvey -f config.json 0 0

In the command line, you will see a list of parameters used by SARvey to run step 0. All parameters that have been changed from the default are indicated:

.. code-block:: none

    ...
    2024-06-19 11:04:28,137 - INFO - Parameter value default
    2024-06-19 11:04:28,137 - INFO - _________ _____ _______
    2024-06-19 11:04:28,138 - INFO - num_cores 5 <--- 50
    2024-06-19 11:04:28,138 - INFO - num_patches 1 1
    2024-06-19 11:04:28,138 - INFO - apply_temporal_unwrapping True True
    2024-06-19 11:04:28,138 - INFO - spatial_unwrapping_method puma puma
    2024-06-19 11:04:28,138 - INFO -
    2024-06-19 11:04:28,138 - INFO - ---------------------------------------------------------------------------------
    2024-06-19 11:04:28,138 - INFO - STEP 0: PREPARATION
    2024-06-19 11:04:28,138 - INFO - ---------------------------------------------------------------------------------
    2024-06-19 11:04:28,138 - INFO - Parameter value default
    2024-06-19 11:04:28,139 - INFO - _________ _____ _______
    2024-06-19 11:04:28,139 - INFO - start_date None None
    2024-06-19 11:04:28,139 - INFO - end_date None None
    2024-06-19 11:04:28,139 - INFO - ifg_network_type sb <--- delaunay
    2024-06-19 11:04:28,139 - INFO - num_ifgs 3 3
    2024-06-19 11:04:28,139 - INFO - max_tbase 100 100
    2024-06-19 11:04:28,139 - INFO - filter_window_size 9 9
    ...

After running this step, a `sbas` directory is created. Inside this directory, you can find the following files:

.. code-block:: none

    outputs/
    ├── temporal_coherence.h5
    ├── ifg_stack.h5
    ├── ifg_network.h5
    ├── coordinates_utm.h5
    ├── config.json
    ├── background_map.h5
    └── pic/
        ├── step_0_temporal_phase_coherence.png
        ├── step_0_interferogram_network.png
        └── step_0_amplitude_image.png


Check the PNG files inside the `outputs/pic` directory and ensure the amplitude image, interferogram network, and temporal coherence look fine. If you are not satisfied with the interferogram network, you can modify the corresponding parameters in the `config.json` file and run step 0 again.

Use the following command to plot the interferograms:

.. code-block:: bash

    sarvey_plot outputs/ifg_stack.h5 -i

This command creates the interferograms as PNG files in the following directory:

.. code-block:: none

    outputs/
    └── pic/
        └── ifgs/
            ├── 0_ifg.png
            ├── 1_ifg.png
            └── ...

Check the interferograms one by one and ensure they look reasonable. In various interferograms, there are fringes associated with deformation approximately at ranges 100-200, azimuth 40-60.


Step 2.1: Run Step 1 of SARvey
''''''''''''''''''''''''''''''

.. code-block:: bash

    sarvey -f config.json 1 1

Outputs of this step are:

.. code-block:: none

    outputs/
    ├── point_network.h5
    ├── p1_ifg_wr.h5
    ├── point_network_parameter.h5
    └── pic/
        ├── selected_pixels_temp_coh_0.8.png
        ├── step_1_mask_p1.png
        ├── step_1_arc_coherence.png
        ├── step_1_arc_coherence_reduced.png
        ├── step_1_rmse_vel_0th_iter.png
        └── step_1_rmse_dem_error_0th_iter.png


Step 2.2: Run Step 2 of SARvey
''''''''''''''''''''''''''''''

.. code-block:: bash

    sarvey -f config.json 2 2


Outputs of this step are:

.. code-block:: none

    outputs/
    ├── p1_ifg_unw.h5
    ├── p1_ts.h5
    └── pic/
        ├── step_2_estimation_dem_error.png
        └── step_2_estimation_velocity.png

Step 2.3: Run Step 3 of SARvey
''''''''''''''''''''''''''''''

.. code-block:: bash

    sarvey -f config.json 3 3


Outputs of this step are:

.. code-block:: none

    outputs/
    ├── coh80_ifg_wr.h5
    ├── coh80_aps.h5
    ├── p1_aps.h5
    ├── p1_ts_filt.h5
    └── pic/
        ├── step_3_temporal_autocorrelation.png
        ├── step_3_stable_points.png
        ├── selected_pixels_temp_coh_0.8.png
        └── step_3_mask_coh80.png


Step 2.4: Run Step 4 of SARvey
''''''''''''''''''''''''''''''

.. code-block:: bash

    sarvey -f config.json 4 4

.. outputs directory structure to be added


The results of step 4 of SARvey, including the time series, are stored in the `coh80_ts.h5` file. The file is named based on the `coherence_p2` parameter in the config.json file.


Step 3: Plot Time Series Results
""""""""""""""""""""""""""""""""

Check the instruction on how to use the `sarvey_plot`.

.. code-block:: bash

    sarvey_plot -h


Plot the time series using the following command. Flag `-t` indicates that you want to plot the time series.

.. code-block:: bash

    sarvey_plot outputs/coh80_ts.h5 -t


You can visualize velocity and DEM error estimation of second-order points. You can also visualize amplitude, DEM, or temporal coherence as the background. Right-click on any point to see its time series. As you will see in the plot, the density of measurement points on the dam is relatively low. In the next section, you will learn how to modify the config file to increase the density of points.


Step 4: Modify Config File and Rerun SARvey
"""""""""""""""""""""""""""""""""""""""""""

Modify the config.json file and change **coherence_p2** from 0.8 to 0.7.

Run steps 3 and 4 using the following command:

.. code-block:: bash

    sarvey -f config.json 3 4


A new file `coh70_ts.h5` is created. You can now visualize this file that has a higher point density.

.. code-block:: bash

    sarvey_plot outputs/coh70_ts.h5 -t


.. note::
    Be cautious that reducing the value of **coherence_p2** too much may include noisy points of low quality in the analysis, potentially leading to poor final results.

    You should carefully read the :ref:`processing` documentation to understand the meaning of each parameter and carefully choose reasonable values. You should also check the details of all parameters using the -p flag in `sarvey` and decide how to tune them.

.. code-block:: bash

    sarvey -f config.json 0 0 -p


Step 5: Export to GIS Format
""""""""""""""""""""""""""""

Export the data to Shapefiles using the following command:

.. code-block:: bash

    sarvey_export outputs/coh70_ts.h5 -o outputs/shp/coh70_ts.shp

You can open the exported data in any GIS software. If you use QGIS, you can use the `PS Time Series Viewer <https://plugins.qgis.org/plugins/pstimeseries/>`_ plugin to draw the time series.


Step 6: Validate Your Results
"""""""""""""""""""""""""""""

You can download a copy of the final SARvey products from `this link <https://doi.org/10.5281/zenodo.12189041>`_. Use these files to compare your results and ensure everything worked correctly.

