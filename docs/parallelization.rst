Parallel processing (user guide)
================================

This page explains how to tune parallel processing in SARvey from a user perspective.

Main settings
-------------

``general:num_cores``
    Controls process-level parallelism in the CPU-heavy processing stages.

``general:num_patches``
    Splits the image stack into spatial patches to reduce memory pressure.
    This is mostly a memory and I/O control, not a direct multiplier for CPU parallelism.

How SARvey controls native thread libraries
-------------------------------------------

At startup, SARvey applies an *unset-only* policy for native thread environment variables.
If a variable is not set, SARvey sets it to the value configured in ``general:num_cores``. If a variable is already set in your shell,
job scheduler, or container environment, SARvey keeps the existing value.

Variables handled:

* ``OMP_NUM_THREADS``
* ``OPENBLAS_NUM_THREADS``
* ``MKL_NUM_THREADS``
* ``NUMEXPR_NUM_THREADS``

Why this matters
----------------

Some SARvey stages use Python multiprocessing while numerical libraries may also use internal
native threads. If both layers are aggressive, CPU oversubscription can happen and runtime may get
worse even with high CPU utilization.

Recommended tuning workflow
---------------------------

1. Start conservative:

   * set ``general:num_cores`` to about half of physical cores
   * keep ``general:num_patches = 1`` unless memory is insufficient

2. Increase ``general:num_cores`` stepwise and compare runtime for the same step and dataset.

3. Increase ``general:num_patches`` only when memory is a bottleneck.

4. On shared systems, prefer explicit environment settings in your job script if your site requires them.

Troubleshooting high CPU load
-----------------------------

If CPU usage is too high:

1. Lower ``general:num_cores`` first.
2. Keep ``general:num_patches`` low unless needed for memory.
3. Check whether your environment defines high values for BLAS/OpenMP variables.
4. Re-run one processing step and compare wall-clock time, not only CPU percent.

Quick notes
-----------

* Different processing steps scale differently with core count.
* Higher CPU usage does not always mean faster processing.
* Patching can increase I/O overhead; use it primarily for memory control.
