.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/luhipi/sarvey/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the Github issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the Github issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

SARvey could always use more documentation, whether as part of the
official SARvey docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/luhipi/sarvey/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Commit Changes
--------------

How to
~~~~~~

1. Fork the `sarvey` repo on Github.
2. Clone your fork locally::

    $ git clone https://github.com/luhipi/sarvey.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv sarvey
    $ cd sarvey/
    $ python setup.py develop

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the
   tests, including testing other Python versions with tox::

    $ make pytest
    $ make lint
    $ make urlcheck
    $ tox

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to Github::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a merge request through the Github website.

Sign your commits
~~~~~~~~~~~~~~~~~

Please note that our license terms only allow signed commits.
A guideline how to sign your work can be found here: https://git-scm.com/book/en/v2/Git-Tools-Signing-Your-Work

If you are using the PyCharm IDE, the `Commit changes` dialog has an option called `Sign-off commit` to
automatically sign your work.


License header
~~~~~~~~~~~~~~

If you commit new Python files, please note that they have to contain the following license header:

.. code:: bash

    # SARvey - A multitemporal InSAR time series tool for the derivation of displacements.
    #
    # Copyright (C) 2021-2024 Andreas Piter (IPI Hannover, piter@ipi.uni-hannover.de)
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


Merge Request Guidelines
------------------------

Before you submit a pull request, check that it meets these guidelines:

1. The merge request should include tests.
2. If the merge request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 3.6, 3.7, 3.8 and 3.9. Check
   https://github.com/luhipi/sarvey/pulls
   and make sure that the tests pass for all supported Python versions.

Tips
----

To run a subset of tests::

$ pytest tests.test_processing

