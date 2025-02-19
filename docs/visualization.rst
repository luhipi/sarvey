.. _visualization:

Visualization
=============


Amplitude image and interferograms
----------------------------------

Even before starting the SARvey processing, you might want to check the amplitude images and interferograms to get to know the study area and the displacement pattern.
This can be done using the command line script `sarvey_plot` which requires the path to the SLC stack as input.
It does not require any processing results or configuration file from SARvey.

.. code-block:: bash

    sarvey_plot inputs/slcStack -i


In the baseline plot (perpendicular and temporal baseline), you can click on a dot corresponding to an image to plot the respective amplitude image.
By switching to the interferogram plot, you can visualize the interferograms from two acquisitions by clicking with LEFT and RIGHT mouse buttons on the two images.



