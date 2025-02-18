.. _tcs:

===========================================================================
Multitemporal InSAR processing workflow for Temporarily Coherent Scatterers
===========================================================================


General information
-------------------
TCS phase unwrapping is only implemented with spatial unwrapping and not for spatiotemporal unwrapping.
In its current implementation, the TCS phase unwrapping cannot be combined with Phase Linking results.
It would require that phase linking is applied on DS pixels only for their coherent lifetime.


Selection of TCS
----------------
CCS are prioritized over TCS, i.e. the whole time span is used if a pixel is both TCS and CCS.
TCS act as additional pixels in the processing.
The file "selection_method.h5" contains a mask for each pixel selection method indicating which pixels are selected by which method.



Processing steps
----------------

Determine which pixel is coherent in which interferogram based on the information about the coherent images.
Accordingly, a spatial network is created for each interferogram seperately with only the points which are coherent in that interferogram.
All interferograms are unwrapped in space.
The network of interferograms is inverted to retrieve the displacement time series.
Notably, each point has a different network of interferograms corresponding to the interferograms in which it is coherent.

