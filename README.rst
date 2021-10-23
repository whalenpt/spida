=========
SPIDA
=========

.. _kissfft: https://github.com/mborgerding/kissfft
.. _indirect target: kissfft
.. _internal target:

Spectral integration and differentiation algorithms (SPIDA). This project implements several
high-level wrappers for spectral transforms: 

#. Fast-Fourier Transforms kissfft (FFTs) 
#. `Fast Discrete Cosine Transforms <https://www.nayuki.io/page/fast-discrete-cosine-transform-algorithms>` (FDCTs) 
#. Hankel-Transforms  

SPIDA also contains adaptive-step Runge-Kutta numerical propagators for stiff partial-differential-equations (PDEs).
These are the same methods as implemented for the python package `rkstiff <https://github.com/whalenpt/rkstiff>`.
Two types of numerical methods are available, integrating factor methods (IF), and exponential time-differencing (ETD).
In particular, the solvers provided are

#. ETD35 (5th order ETD with 3rd order embedding)
#. ETD34 (4th order ETD with 3rd order embedding) 
#. IF34 (4th order IF with 3rd order embedding)
#. IF45DP (5th order IF with 4th order embedding)

In general, one should prefer ETD35 as it often has the best speed and stability for diagonal systems.
No non-diagonal methods are implemented as of yet (see the python package for that functionality).

This project implements the
functionality of the python rkstiff package in C++ and also includes several high-level wrappers
for third-party Fast-Fourier Transforms (FFTs), among other transforms. The aim of the project
is to provide fast solutions of partial-differential-equations (PDEs) using pseudo-spectral methods.

Dependencies
------------

Third-party packages use MIT licenses or similarly permissive licenses

Libraries built with
* kissfft
* nayukidct
* pwutils
* boost
 
Testing done with
* googletest

Compile and build
* CMake

Usage
-----

Check out the demos. These can be built by configuring CMake with::
    cmake -S . -B build -DCMAKE_DEMOS=ON

License
-------
This project is licensed under the MIT License - see the `LICENSE <./LICENSE>` file for details.

Contact
-------
Patrick Whalen - whalenpt@gmail.com



