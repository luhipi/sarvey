.. _preparation:

===========
Preparation
===========

SARvey requires a coregistered stack of SLC and the related geometry information in the MiaplPy_ data format.
The coregistered stack of SLC can be created using an InSAR processor, such as ISCE, GAMMA, or SNAP.
Currently MiaplPy only supports ISCE_. Support for GAMMA and SNAP_ is planned for future.
After creating the coregistered stack of SLC, run the “load_data” step from MiaplPy to create the “inputs” directory which contains “slcStack.h5” and “geometryRadar.h5”.



Preprocessing
-------------

ISCE
^^^^
... ISCE brief processing to be added

The ISCE products should have the following directory structure that is later in `Loading Data into MiaplPy`_ step.

::

    ISCE_processed_data
    ├─ reference
    │   ├─ IW*.xml
    │   └─ ...
    ├─ merged
    │   ├─ SLC
    │   │   ├─ YYYYMMDD
    │   │   │   ├─ YYYYMMDD.slc.full
    │   │   │   └─ ...
    │   │   ├─ YYYYMMDD
    │   │   ├─ YYYYMMDD
    │   ├─ geom_reference
    │   │   ├─ hgt.rdr.full
    │   │   ├─ lat.rdr.full
    │   │   ├─ lon.rdr.full
    │   │   ├─ los.rdr.full
    │   │   └─ ...
    └─ baselines
        └─ YYYYMMDD_YYYYMMDD
            └─ YYYYMMDD_YYYYMMDD.txt


GAMMA
^^^^^
Support is in progress.


SNAP
^^^^
Support is planned for future.


Loading Data to MiaplPy Format
------------------------------

**Loading Data into MiaplPy**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run the `load_data` step of MiaplPy to convert the preprocessed stack of SLC to `slcStack.h5` and `geometryRadar.h5`.
Refer to MiaplPy_ instruction on how to prepare the stack of coregistered SLC and modify the template file.

.. code-block:: bash

   miaplpyApp miaplpy_template_file.txt --dostep load_data

The output includes the following directory structure that is later used as input in SARvey processing:

::

    inputs
      ├── slcStack.h5
      └── geometryRadar.h5



**Check the data**
^^^^^^^^^^^^^^^^^^

Use `info.py` from MintPy_ to check the files' information.

.. code-block:: bash

   info.py inputs/slcStack.h5
   info.py inputs/geometryRadar.h5


Use `view.py` from MintPy_ to visualize the files and make sure they look fine.

.. code-block:: bash

   view.py inputs/slcStack.h5
   view.py inputs/geometryRadar.h5



Optional Steps
--------------


**Phase Linking**
^^^^^^^^^^^^^^^^^


This step is optional. You can run it if you wish to perform distributed scatterers (DS) analysis.
**Caution:** This step is computationally heavy and might be time-consuming for large datasets.

.. code-block:: bash

   miaplpyApp miaplpy_template_file.txt --dostep phase_linking
    miaplpyApp miaplpy_template_file.txt --dostep concatenate_patches

The output includes the following directory structure that is later used as additional input in SARvey processing if the config file is modified to inclued DS analysis.

::

    MiaplPy working directory
    ├─ inverted
    │   ├── phase_series.h5
    │   ├── ...
    ├── maskPS.h5
    └── ...



Subset Data
^^^^^^^^^^^

Data loaded into MiaplPy can be subset using Mintpy_'s subset function.
This is particularly useful if you have a dataset in MiaplPy format and want to crop a small area of it.
Both slcStack.h5 and geometryRadar.h5 should be subset with the same range and azimuth coordinate ranges.
Also the Phase Linking results (phase_series.h5 and maskPS.h5) should be subset if it has been created.
Please refer to Mintpy_ for more instruction to subset.
Run `subset.py -h` for information about parameters.
The following example crops the data between 500 and 800 in range and 100 and 1000 in azimuth coordinates.


.. code-block:: bash

   subset.py -h

   subset.py inputs/slcStack.h5 -x 500 800 -y 100 1000 -o inputs_crop/slcStack.h5
   subset.py inputs/geometryRadar.h5 -x 500 800 -y 100 1000 -o inputs_crop/geometryRadar.h5

   subset.py inverted/phase_series.h5 -x 500 800 -y 100 1000 -o inverted_crop/phase_series.h5
   subset.py maskPS.h5 -x 500 800 -y 100 1000 -o inverted_crop/maskPS.h5


`Check the data`_ after subsetting it and make sure all products look correct.



Create Manual Mask
^^^^^^^^^^^^^^^^^^
A mask can be created manually using MintPy's `generate_mask.py` tool.
This is particularly useful if you want to limit the MTInSAR processing to certain areas.
Run `generate_mask.py -h` for information about parameters.
The following example allows to draw a polygon on top of the DEM to create a mask.

.. code-block:: bash

   generate_mask.py -h

   generate_mask.py inputs/geometryRadar.h5 height -o mask.h5 --roipoly         # draw polygon on top of the DEM

Alternatively, a mask can be drawn on top of the temporal coherence map, in case step 0 (preparation) of `sarvey` has been executed already.

.. code-block:: bash

   generate_mask.py results_dir/temporal_coherence.h5 -o mask.h5 --roipoly         # draw polygon on top of the temporal coherence image

Follow the instructions in the terminal:

    Select points in the figure by enclosing them within a polygon.
    Press the 'esc' key to start a new polygon.
    Try hold to left key to move a single vertex.
    After complete the selection, close the figure/window to continue.



.. _MiaplPy: https://github.com/insarlab/MiaplPy
.. _ISCE: https://github.com/isce-framework/isce2
.. _SNAP: https://step.esa.int/main/toolboxes/snap
.. _MintPy: https://github.com/insarlab/MintPy

