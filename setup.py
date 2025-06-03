#!/usr/bin/env python

# SARvey - A multitemporal InSAR time series tool for the derivation of displacements.
#
# Copyright (C) 2021-2025 Andreas Piter (IPI Hannover, piter@ipi.uni-hannover.de)
#
# This software was developed together with FERN.Lab (fernlab@gfz-potsdam.de) in the context
# of the SAR4Infra project with funds of the German Federal Ministry for Digital and
# Transport and contributions from Landesamt fuer Vermessung und Geoinformation
# Schleswig-Holstein and Landesbetrieb Strassenbau und Verkehr Schleswig-Holstein.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Important: This package uses PyMaxFlow. The core of PyMaxflows library is the C++
# implementation by Vladimir Kolmogorov. It is also licensed under the GPL, but it REQUIRES that you
# cite [BOYKOV04] (see LICENSE) in any resulting publication if you use this code for research purposes.
# This requirement extends to SARvey.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

"""The setup script."""


from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

version = {}
with open("sarvey/version.py") as version_file:
    exec(version_file.read(), version)

req = [
    "cython", "numpy<=1.26", "pyproj", "matplotlib", "numba", "scipy",
    "mintpy", "h5py", "overpy", "gstools", "shapely", "pandas", "geopandas", "pymaxflow",
    "pillow", "importlib_resources", "kamui", "json5", "cmcrameri", 'pydantic<=1.10.10',
    "miaplpy @ git+https://github.com/insarlab/MiaplPy.git"


    "kamui @ git+https://github.com/mahmud1/kamui.git@numpy"
]

req_setup = []

req_test = ['pytest>=3', 'pytest-cov', 'pytest-reporter-html1', 'urlchecker']

req_doc = [
    'sphinx>=4.1.1',
    'sphinx-argparse',
    'sphinx-autodoc-typehints',
    'sphinx_rtd_theme'
]

req_lint = ['flake8', 'pycodestyle', 'pydocstyle']

req_dev = ['twine'] + req_setup + req_test + req_doc + req_lint

extra_req = ["gdal"]

setup(
    author="Andreas Piter",
    author_email='piter@ipi.uni-hannover.de',
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: Release Candidate',
        'Intended Audience :: Researchers',
        'None',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10'
    ],
    description="InSAR time series analysis software for SAR4Infra project",
    entry_points={
        'console_scripts': [
            'sarvey=sarvey.sarvey_mti:main',
            'sarvey_plot=sarvey.sarvey_plot:main',
            'sarvey_export=sarvey.sarvey_export:main',
            'sarvey_mask=sarvey.sarvey_mask:main',
            'sarvey_osm=sarvey.sarvey_osm:main',
        ],
    },
    extras_require={
        "doc": req_doc,
        "test": req_test,
        "lint": req_lint,
        "dev": req_dev
    },
    install_requires=req,
    license="GPLv3",
    include_package_data=True,
    keywords='sarvey',
    long_description=readme,
    name='sarvey',
    packages=find_packages(include=['sarvey', 'sarvey.*']),
    setup_requires=req_setup,
    test_suite='tests',
    tests_require=req_test,
    url='https://github.com/luhipi/sarvey',
    version=version['__version__'],
    zip_safe=False,
)
