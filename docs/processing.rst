.. _processing:

=======================================
Multitemporal InSAR processing workflow
=======================================

The `sarvey` command line interface executes the multitemporal InSAR processing workflow.
The workflow is described in the paper

    Piter, A., Haghshenas Haghighi, M., Motagh, M.(2024). An in-depth study on Sentinel-1 InSAR for transport infrastructure monitoring. PFG - Journal of Photogrammetry, Remote Sensing and Geoinformation Science. (paper currently under review).

All processing steps are described in detail in the following sections.
Two processing strategies are provided with either one- or two-step unwrapping.
The workflow should be decided based on the characteristics of the displacement (spatial extend, magnitude, temporal behaviour).
The parameters of each step are handled via the `configuration file`_ for which the parameters are named within the description of each step.


Configuration file
------------------
The configuration file is a JSON file containing all the parameters required to run `sarvey`.
This file can be generated using the `sarvey` command with the **"-g"** flag, where you can specify your desired filename.


.. code-block:: bash

     sarvey -f config.json 0 0 -g

Note: The above command only generates a configuration file. Although step 0 is specified, it will not be executed.

The configuration file has various sections, as detailed below:


* General


 This section includes top-level parameters such as the number of cores and the unwrapping method.
 It specifies the paths to the input and output data. The paths can be either absolute or relative.
 Further, it defines the logging level displayed in the command line and the directory path where log files will be stored.



* phase_linking


 This section specifies the Phase Linking parameters. By default, `"use_phase_linking_results": false`.
 If you wish to perform DS analysis, change it to `true`. Note: If `"use_phase_linking_results": true`, you must complete the corresponding step of MiaplPy as described in `preparation <preparation.rst/#Phase Linking>`_. In the configuration file, set `inverted_path` to the path of the inverted directory of MiaplPy data.



* preparation


 This section includes network parameters, such as the start and end dates, network type, and `filter_window_size`, which specifies the window size used to estimate the temporal coherence for each pixel.


* consistency_check


 This section contains parameters related to the first order points.

* unwrapping


 This section will specify parameters related to the unwrapping process.

* filtering


 This section defines the parameters for atmospheric estimation and filtering. Atmospheric filtering is enabled by default. To skip it, set `"apply_aps_filtering": false`.


* densification


 This section includes the settings for second order points.





Processing steps for two-step unwrapping workflow
-------------------------------------------------

Step 0: Preparation
^^^^^^^^^^^^^^^^^^^

- Loading the resampled SLC data:
    The resampled SLC (Single Look Complex) data is read from the inputs/slcStack.h5
    This data is complex-valued and contains both amplitude and phase information.
    The data is subsetted to the specified time span (via **preparation:start_date** and **preparation:end_date** in the config file).
    A description of how to prepare the data and make a spatial subset of the data is described in `data preparation in MiaplPy <preparation.rst>`_.

- Designing the interferogram network:
    From the stack of SLC images, the interferogram network is designed.
    The network of interferograms is designed based on the temporal and perpendicular baselines of the SLC images.
    Different networks can be created (via **preparation:ifg_network_type** in the config file) and should be chosen based on the characteristics of the displacement (spatial extend, magnitude, temporal behaviour).
    Currently five types of networks are supported:

    a) small baseline network ('sb') (Berardino et al. 2002),
    b) small temporal baseline network ('stb') (only consecutive images are used to form interferograms)
    c) small temporal baselines + yearly interferograms ('stb_yearly')
    d) delaunay network ('delaunay')
    e) star network ('star', single-reference network) (Ferretti et al. 2001)


- Generating a stack of interferograms:
    The stack of interferograms is generated based on the specified interferogram network.

- Estimating the temporal coherence:
    The phase noise of each pixel is approximated by the estimation of the temporal phase coherence (Zhao and Mallorqui 2019).
    Thereby, a low-pass filter with a certain window size is used (**preparation:filter_window_size**).
    The temporal coherence is used to select the first- and second-order points in the later steps (**consistency_check:coherence_p1** and **filtering:coherence_p2**).

- Output of this step
    - background_map.h5
    - ifg_stack.h5
    - coordinates_utm.h5
    - ifg_network.h5
    - temporal_coherence.h5


Step 1: Consistency Check
^^^^^^^^^^^^^^^^^^^^^^^^^


- Selecting candidates for first order points:
    Candidates for the first-order points are selected based on the temporal coherence threshold (**consistency_check:coherence_p1**).
    However, not all points with a coherence above the threshold are selected, but only those which have the highest coherence within a grid cell of size **consistency_check:grid_size** (in [m]).
    A mask file can be specified (**consistency_check:mask_p1_file**) to limit the first-order points to the given area of interest.

- Creating a spatial network:
    After selecting the candidates for first order points, the method creates a spatial network to connect the first-order points.
    For each arc in the network, the double difference phase time series is calculated.
    A delaunay network ensures the connectivity in the spatial network and k-nearest neighbors (**consistency_check:num_nearest_neighbours**) can be used to increase the redundancy in the network.
    Arcs with a distance above a threshold (**consistency_check:max_arc_length**) are removed from the network to reduce the impact of the atmospheric effects.

- Temporal unwrapping:
    All arcs in the spatial network are temporally unwrapped based on a phase model consisting of DEM error difference and velocity difference between the two points of the arc.
    The temporal coherence derived from the model fit is maximized by searching within a search space of given bounds (**consistency_check:velocity_bound** and **consistency_check:dem_error_bound**).
    Within the bounds, the search space is discretized (**consistency_check:num_optimization_samples**).
    The final parameters for each arc are derived from a gradient descent refinement of the discrete search space result.

- Performing a consistency check on the data:
    During the atmospheric filtering in step 3, only high quality first-order points are supposed to be used.
    Therefore, outliers among the candidates are removed with a consistency check.
    The consistency check is based on the estimated temporal coherence of the temporal unwrapping of each arc.
    A point is assumed to be an outlier, if it is connected by many arcs having a low temporal coherence from temporal unwrapping.
    Arcs with a temporal coherence below a threshold are removed (**consistency_check:arc_unwrapping_coherence**).
    Similarly, points with mean coherence of all connected arcs are removed (specified by the same parameter **consistency_check:arc_unwrapping_coherence**).
    Moreover, points which are connected by a number of arcs less than a threshold (**consistency_check:min_num_arc**) are removed.
    Afterwards, the consistency within the spatial network is checked.
    For this purpose, the parameters (DEM error difference and velocity difference) of all arcs are integrated in the spatial network relative to an arbitrary reference point using least squares.
    The residuals of the integration are used to identify outliers.

- Output of this step
    - point_network.h5
    - point_network_parameter.h5
    - p1_ifg_wr.h5

Step 2: Unwrapping
^^^^^^^^^^^^^^^^^^

Two unwrapping options (**general:apply_temporal_unwrapping**, also applies to step 4) are implemented and should be chosen based on the characteristics of the displacement (spatial extend, magnitude, temporal behaviour).

- Output of this step
    - p1_ifg_unw.h5
    - p1_ifg_ts.h5

Option 1) Unwrapping in time and space
""""""""""""""""""""""""""""""""""""""

- Integrating parameters from arcs to points:
    The temporal unwrapping results of the spatial network from consistency check in step 1 are used in this step.
    The parameters of the arcs are integrated relative to an arbitrary reference point from the arcs to the points using least squares.

- Removing phase contributions (mean velocity and DEM error):
    After integrating the parameters, the phase contributions are removed from the wrapped interferometric phase of the first-order points.

- Spatial unwrapping of the residuals:
    The residuals in each interferogram are unwrapped in space using a sparse point network unwrapping method (**general:spatial_unwrapping_method**) (Bioucas-Dias and Valadao 2007, Boykov and Kolmogorov 2004).
    The spatial neighbourhood for unwrapping is defined by the arcs of the spatial network.
    There are two options (**unwrapping:use_arcs_from_temporal_unwrapping**).
    Either the spatial network from consistency check (step 2) can be used for unwrapping, i.e. the spatial network after removing arcs with a low temporal coherence from temporal unwrapping.
    Or, the spatial network is re-created with a delaunay network.

- Restore phase contributions to the spatially unwrapped residual phase:
    Finally, the phase contributions are added back to the spatially unwrapped residual phase of each point.

- Adjust reference:
    All restored unwrapped interferograms are referenced to the peak of velocity histogram derived from all points.

- Inverting the interferogram network:
    The interferogram network is inverted for each point to retrieve the displacement time series relative to the first acquisition.

Option 2) Unwrapping in space
"""""""""""""""""""""""""""""

- Spatial unwrapping:
    The interferograms are unwrapped independently in space with a sparse point network unwrapping method (**general:spatial_unwrapping_method**) (Bioucas-Dias and Valadao 2007, Boykov and Kolmogorov 2004).
    The spatial neighbourhood for unwrapping is defined by the arcs of the spatial network.
    There are two options (**unwrapping:use_arcs_from_temporal_unwrapping**).
    Either the spatial network from consistency check (step 2) can be used for unwrapping, i.e. the spatial network after removing arcs with a low temporal coherence from temporal unwrapping.
    Or, the spatial network is re-created with a delaunay network.

- Adjust reference:
    All unwrapped interferograms are referenced to the peak of velocity histogram derived from all points.

- Inverting the interferogram network:
    The interferogram network is inverted for each point to retrieve the displacement time series relative to the first acquisition.

Step 3: Filtering
^^^^^^^^^^^^^^^^^

In this step, the atmospheric phase screen (APS) is estimated from the displacement time series of the first-order points.
Afterwards, the APS is interpolated to the location of the second-order points.
The filtering can be skipped by setting **filtering:apply_aps_filtering** to True.
However, the step 3 has to be executed as the second-order points are selected during this step.

- Selecting pixels with no or linear displacement:
    Among the first-order points, the points with no or merely linear displacement are selected (**filtering:use_moving_points**).
    It is assumed that for these points, the phase consists only of atmospheric effect and noise after removing the mean velocity and DEM error.
    Points with a non-linear displacement behaviour are removed by a threshold on the temporal autocorrelation of the displacement time series (**filtering:max_temporal_autocorrelation**) (Crosetto et al. 2018).
    A regular grid (**filtering:grid_size** in [m]) is applied to select the first-order points with the lowest temporal autocorrelation to reduce the computational complexity during filtering.

- Selecting second-order points:
    Second-order points are selected based on a temporal coherence threshold (**filtering:coherence_p2**) on the temporal phase coherence computed during step 0.
    A mask file can be specified (**filtering:mask_p2_file**) to limit the second-order points to the given area of interest.
    Second-order points can also be selected based on the results of phase-linking (set **phase_linking:use_phase_linking_results** to True) implemented in MiaplPy (Mirzaee et al. 2023).
    More information on Miaplpy and phase-linking can be found `here <preparation>`_.
    The number of siblings (**phase_linking:num_siblings**) used during phase-linking within MiaplPy processing needs to be specified to identify the distributed scatterers (DS) among the pixels selected by MiaplPy.
    A mask file can be specified (**phase_linking:mask_phase_linking_file**) to limit the phase-linking to the given area of interest.
    MiaplPy also provides a selection of persistent scatterers (PS) which can be included as second-order points (set **phase_linking:use_ps** to True) and also specify the path to the maskPS.h5 (**phase_linking:mask_ps_file**) which is also an output of MiaplPy.
    In case the second-order points are selected among the results from MiaplPy, the filtered interferometric phase (MiaplPy result) is used for the respective points.
    The DS pixels from MiaplPy and the pixels selected with the temporal phase coherence from step 0 are both selected with the same coherence threshold (**filtering:coherence_p2**).

- Estimating the atmospheric phase screen (APS):
    The estimation of the APS takes place in time-domain and not interferogram-domain to reduce the computational time.
    The phase contributions are removed from the first-order points which were selected for atmospheric filtering.
    Their residual time series contains atmospheric phase contributions and noise.
    As the APS is assumed to be spatially correlated, the residuals of all points are spatially filtered (**filtering:interpolation_method**) independently for each time step.
    After filtering, the estimated APS is interpolated to the location of the second-order points.

- Output of this step
    - p1_ts_filt.h5
    - p1_aps.h5
    - p2_cohXX_aps.h5
    - p2_cohXX_ifg_wr.h5

The placeholder XX depends on the threshold for the temporal coherence used for selecting the second-order points.
For example, a threshold of 0.8 would result in p2_coh80_aps.h5 and p2_coh80_ifg_wr.h5.

Step 4: Densification
^^^^^^^^^^^^^^^^^^^^^

Two unwrapping options (**general:apply_temporal_unwrapping**, also applies to step 2) are implemented and should be chosen based on the characteristics of the displacement (spatial extend, magnitude, temporal behaviour).

- Output of this step
    - p2_cohXX_ifg_unw.h5
    - p2_cohXX_ts.h5

The placeholder XX depends on the threshold for the temporal coherence used for selecting the second-order points during filtering in step 3.
For example, a threshold of 0.8 would result in p2_coh80_ifg_unw.h5 and p2_coh80_ts.h5.

Option 1: Unwrapping in time and space
""""""""""""""""""""""""""""""""""""""

- Removing APS from interferograms
    The wrapped interferograms are corrected for the interpolated APS for both the first and second order points.

- Densify network:
    The parameters (DEM error and velocity) of each second-order point are estimated independently from the other second-order points.
    The parameters are estimated by temporal unwrapping with respect to the closest first-order points (**densification:num_connections_to_p1**, **densification:max_distance_to_p1**) with a phase model consisting of DEM error and velocity (**densification:velocity_bound** and **densification:dem_error_bound**, **densification:num_optimization_samples**).
    The densification is similar to the approach described by Van Leijen (2014), but jointly maximizes the temporal coherence to find the parameters that fit best to all arcs connecting the second-order point to the first-order points.

- Remove outliers:
    Second-order points which could not be temporally unwrapped with respect to the closest first-order points are removed.
    For this purpose, a threshold on the joint temporal coherence considering the residuals of all arcs connecting the respective second-order point to the closest first-order points is applied (**densification:arc_unwrapping_coherence**).
    First-order points receive a joint temporal coherence value of 1.0 to avoid them being removed from the final set of points.

- Removing phase contributions (mean velocity and DEM error):
    After estimating the parameters of the second-order points, the phase contributions are removed from the wrapped interferometric phase of the first-order points.

- Spatial unwrapping of the residuals:
    The residuals in each interferogram are unwrapped in space using a sparse point network unwrapping method (**general:spatial_unwrapping_method**) (Bioucas-Dias and Valadao 2007, Boykov and Kolmogorov 2004).
    The spatial neighbourhood for unwrapping is defined by spatial network including both first- and second-order points.
    It is created with a delaunay network.

- Restore phase contributions to the spatially unwrapped residual phase:
    Finally, the phase contributions are added back to the spatially unwrapped residual phase of each point.

- Adjust reference:
    All restored unwrapped interferograms are referenced to the peak of velocity histogram derived from all points.

- Inverting the interferogram network:
    The interferogram network is inverted for each point to retrieve the displacement time series relative to the first acquisition.

Option 2: Unwrapping in space
"""""""""""""""""""""""""""""

- Removing APS from interferograms
    The wrapped interferograms are corrected for the interpolated APS for both the first and second order points.

Afterwards, the processing is the same as in the spatial unwrapping during step 2.


Handling big datasets
---------------------
The processing of large datasets can be computationally expensive and time-consuming.
Especially the estimation of the temporal phase coherence in step 0 is a bottleneck, also in terms of memory consumption.
Therefore, it is recommended to set **general:num_cores** for parallel processing.
By setting **general:num_patches** the data is split into spatial patches and processed subsequently to fit into memory.


Processing steps for one-step unwrapping workflow
-------------------------------------------------
The one-step unwrapping workflow is an alternative to the two-step unwrapping workflow.
The steps are similar to the workflow described above, but is only executed until step 2.
This workflow is meant for processing small areas where the atmospheric filtering is not required as the reference point will be selected close to the area of interest.
The idea behind the one-step unwrapping workflow is to apply the consistency check based on the temporal unwrapping (step 1) to all pixels, without differentiating between first and second order points.
This can yield better unwrapping results compared to the two-step unwrapping in case DEM error and/or velocity highly vary in space.
For this purpose, the pixels are selected without gridding (set **preparation:grid_size** to Zero, i.e. all pixels above the specified coherence threshold are selected as final points.
Since the densification step is not performed, you should reduce the coherence threshold (**consistency_check:coherence_p1**) to select the desired number of points.


Literature
----------

* Piter, A., Haghshenas Haghighi, M., Motagh, M.(2024). An in-depth study on Sentinel-1 InSAR for transport infrastructure monitoring. PFG - Journal of Photogrammetry, Remote Sensing and Geoinformation Science. (paper currently under review).

* Zhao F, Mallorqui JJ (2019). A Temporal Phase Coherence Estimation Algorithm and Its Application on DInSAR Pixel Selection. IEEE Transactions on Geoscience and Remote Sensing 57(11):8350–8361, DOI 10.1109/TGRS.2019.2920536

* Ferretti A, Prati C, Rocca F (2001). Permanent scatterers in SAR interferometry. IEEE Transactions on Geoscience and Remote Sensing 39(1):8–20

* Berardino P, Fornaro G, Lanari R, Sansosti E (2002). A new algorithm for surface deformation monitoring based on small baseline differential SAR interferograms. IEEE Transactions on Geoscience and Remote Sensing 40(11):2375–2383

* Bioucas-Dias JM, Valadao G (2007). Phase Unwrapping via Graph Cuts. IEEE Transactions on Image Processing 16(3):698–709, DOI 10.1109/TIP.2006.888351

* Mirzaee S, Amelung F, Fattahi H (2023). Non-linear phase linking using joined distributed and persistent scatterers. Computers & Geosciences 171:105291, DOI 10.1016/j.cageo.2022.105291

* Crosetto M, Devanthéry N, Monserrat O, Barra A, Cuevas-González M, Mróz M, Botey-Bassols J, Vázquez-Suné E, Crippa B (2018). A persistent scatterer interferometry procedure based on stable areas to filter the atmospheric component. Remote Sensing 10(11):1780

* Van Leijen FJ (2014). Persistent scatterer interferometry based on geodetic estimation theory. PhD thesis

* Boykov Y, Kolmogorov V (2004) An experimental comparison of min-cut/max- flow algorithms for energy minimization in vision. IEEE Transactions on Pattern Analysis and Machine Intelligence 26(9):1124–1137, DOI 10.1109/TPAMI.2004.60
