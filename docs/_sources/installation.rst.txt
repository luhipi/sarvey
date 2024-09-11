.. _installation:

============
Installation
============

SARvey is a cross-platform python-based software and can be installed on
  * `Linux`_
  * `MacOS ARM (Apple Silicon M2)`_
  * `Windows using WSL`_


Linux
-----

On Linux, SARvey can be installed `Using Mamba (recommended)`_ or `Using Anaconda or Miniconda`_.

Using Mamba (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^

Using mamba_ (latest version recommended), **SARvey** is installed as follows:


1. Clone the SARvey source code and install SARvey and all dependencies from the environment_sarvey.yml file:

   .. code-block:: bash

    git clone git@gitlab.projekt.uni-hannover.de:ipi-sar4infra/timeseries.git
    cd timeseries


2. Create virtual environment for **SARvey** (optional but recommended):

   .. code-block:: bash

    pip install conda-merge
    wget https://raw.githubusercontent.com/insarlab/MiaplPy/main/conda-env.yml
    conda-merge conda-env.yml tests/CI_docker/context/environment_sarvey.yml > env.yml
    mamba env create -n sarvey -f env.yml
    rm env.yml conda-env.yml
    mamba activate sarvey
    pip install git+https://github.com/insarlab/MiaplPy.git
    pip install .


This is the preferred method to install **SARvey**, as it always installs the most recent stable release and
automatically resolves all the dependencies.


Using Anaconda or Miniconda
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using conda_ (latest version recommended), **SARvey** is installed as follows:


1. Then clone the **SARvey** source code and install **SARvey** and all dependencies from the environment_sarvey.yml file:

   .. code-block:: bash

    git clone git@gitlab.projekt.uni-hannover.de:ipi-sar4infra/timeseries.git
    cd timeseries


1. Create virtual environment for **SARvey** (optional but recommended):

   .. code-block:: bash

    pip install conda-merge
    wget https://raw.githubusercontent.com/insarlab/MiaplPy/main/conda-env.yml
    conda-merge conda-env.yml tests/CI_docker/context/environment_sarvey.yml > env.yml
    conda env create -n sarvey -f env.yml
    rm env.yml conda-env.yml
    conda activate sarvey
    pip install git+https://github.com/insarlab/MiaplPy.git
    pip install .


MacOS ARM (Apple Silicon M2)
----------------------------

This guide provides instructions for installing SARvey on MacOS ARM M2 using conda_.
If you do not have Conda, install `Conda for Mac`_.
Using conda_ (latest version recommended), SARvey is installed as follows:

0. **Create a directory for the SARvey package and navigate to it in the terminal. You can choose any other directory if you prefer.**

    .. code-block:: bash

        mkdir -p ~/software/sarvey

1. **Install MiaplPy before installing SARvey in the same environment where you want to install SARvey.**

    .. code-block:: bash

        cd ~/software/sarvey
        git clone https://github.com/insarlab/MiaplPy.git
        cd MiaplPy

    1.1 Open `conda-env.yml` in an editor of your choice and comment out the line `isce2`. Alternatively, you can run the following command:.

    .. code-block:: bash

        sed -i '' '/isce2/s/^/# /' conda-env.yml

    1.2 Install the package using Conda.

    .. code-block:: bash

        conda env update --name sarvey --file conda-env.yml
        conda activate sarvey
        python -m pip install .

2. **Install SARvey**

   2.1 Download the source code of the SARvey package.

    .. code-block:: bash

        cd ~/software/sarvey
        git clone git@gitlab.projekt.uni-hannover.de:ipi-sar4infra/timeseries.git
        cd timeseries

   2.2 Open `tests/CI_docker/context/environment_sarvey.yml` in an editor of your choice and comment out the lines `isce2` and `gcc_linux-64`. Alternatively, you can run the following commands.

    .. code-block:: bash

         sed -i '' '/isce2/s/^/# /' tests/CI_docker/context/environment_sarvey.yml
         sed -i '' '/gcc_linux-64/s/^/# /' tests/CI_docker/context/environment_sarvey.yml

    Note: As of the time of creation of this document, `isce2` for MacOS ARM64 is not available in Conda repositories. Therefore, it is skipped, but it should not cause any problems for running SARvey. Also, `gcc_linux-64` is not required on ARM64.

   2.3 Install Timeseries using the same environment that you used to install MiaplPy.

    .. code-block:: bash

        conda env update --name sarvey -f tests/CI_docker/context/environment_sarvey.yml
        conda activate sarvey
        pip install .

3. **Set up the PATH for MiaplPy and SARvey.**

   3.1 Run the following commands to set up the path in `~/source_sarvey.sh`.

    .. code-block:: bash

        echo 'export miaplpy_path=~/software/sarvey/MiaplPy/src/' > ~/source_sarvey.sh
        echo 'export PYTHONPATH=${PYTHONPATH:+$PYTHONPATH:}$miaplpy_path' >> ~/source_sarvey.sh
        echo 'export timeseries_path=~/software/sarvey/timeseries' >> ~/source_sarvey.sh
        echo 'export PATH=${PATH}:$timeseries_path:$timeseries_path/sarvey' >> ~/source_sarvey.sh
        echo 'export PYTHONPATH=${PYTHONPATH:+$PYTHONPATH:}:$timeseries_path' >> ~/source_sarvey.sh

4. **Test the installation**

   4.1. Open a new terminal and activate the software.

    .. code-block:: bash

        conda activate sarvey
        source ~/source_sarvey.sh

   4.2. Run the following commands. If the help messages of SARvey and MiaplPy are shown, the installation is correctly done.

    .. code-block:: bash

        sarvey -h


Windows using WSL
-----------------

On Windows, SARvey is tested on Windows Subsystem for Linux (WSL_) version 2. Please follow the `Linux`_ installation.



.. note::

    Timeseries has been tested with Python 3.6+., i.e., should be fully compatible to all Python versions from 3.6 onwards.


.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/
.. _conda: https://conda.io/docs
.. _mamba: https://github.com/mamba-org/mamba
.. _Conda for Mac: https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html
.. _WSL: https://learn.microsoft.com/en-us/windows/wsl/

