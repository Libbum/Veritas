Veritas Quick Start Guide
=========================

Veritas is an efficient tool for continuum Vlasov-Maxwell simulations.
Implemented in modern c++14, using state of the art numerical schemes for maximal accuracy and an adaptive mesh to minimise memory and runtime requirements.

Dependencies
------------

Veritas requires a *LAPACK* package on your machine.
We highly recommend using the *Intel MKL* library as other performance optimisations have been integrated if *MKL* is enabled.

A c++14 aware compiler is also needed.
This means older complier versions will have issues with the codebase.
We have extensively tested *GCC* and *Intel* compilers.

Configuration and compilation is aided by *CMake*.
Building Veritas is possible on Linux, OSX and Windows; however most testing has been performed on Linux (and to a lesser extent OSX).


Obtaining the Source Code
-------------------------

You may clone the latest build from the Veritas git repository via the command line ::

   $ git clone https://github.com/Libbum/Veritas.git .

or if you have your ssh keys configured with github::

   $ git clone git@github.com:Libbum/Veritas.git .

If you prefer a stable version of the code, release archives can be found `here <https://github.com/Libbum/Veritas/releases>`_, which can be extracted and used in a similar manner as the cloned repository data.

Installation
------------

With the dependencies above configured correctly, and a copy of the Vertias source code either cloned from the git repository or extracted from an archive file, move to the build directory ::

   $ cd veritas/build

The *CMake* build script assumes the environment variables ``CC`` and ``CXX`` indicate the *c* and *c++* compilers you wish to use, which may not be the ones identified in your current path or loaded modules list.
If these values are not set, or you wish to change them, simply call::

   $ export CC=/path/to/cc && export CXX=/path/to/cxx

once configured, run the *CMake* configuration via ::

   $ cmake ..

It's also possible to call cmake with compiler arguments if you do not wish to override your environment.
Using *GCC* for example::

   $ cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CC_COMPILER=gcc ..

Another configuration option one may wish to use is the ``USE_MKL``.
By default this option is on and will search for a configured *MKL* library on your system.
If one is not found, the configuration will revert to a standard *LAPACK* implementation.
However, you may manually control this configuration via::

   $ cmake -DUSE_MKL=off ..

Once the ``CMake`` configuration is successful, you may now make and install Veritas. ::

   $ make
   $ make install

This will place a ``veritas`` binary in the folder ../bin relative to the build directory, as well as set up some output folders for simulations.

If subsequent changes are made to the Vertas source after the initial *CMake* configuration, there is a ``compile.sh`` bash script in the bin directory which can be run to easily recompile when needed.

Usage
------

TODO: Running a Job

TODO: Customising ``veritas.cpp``

Contribute
----------

- Issue Tracker: `github.com/Libbum/Veritas/issues <https://github.com/Libbum/Veritas/issues>`_
- Source Code: `github.com/Libbum/Veritas <https://github.com/Libbum/Veritas>`_

Support
-------

Bugs can be submitted through the `tracker <https://github.com/Libbum/Veritas/issues>`_ at any time.
If you are having other problems, please let us know.
You can contact us directly via the `contact form <http://ft.nephy.chalmers.se/veritas/#three>`_ on the Veritas web page.

License
-------

The project is licensed under the `MIT <https://github.com/Libbum/Veritas/blob/master/LICENSEhttps://github.com/Libbum/Veritas/blob/master/LICENSE>`_ license.
