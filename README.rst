=========
SPIDA
=========

Spectral integration and differentiation algorithms (SPIDA). This project implements several
high-level wrappers for spectral transforms:

#. Fast-Fourier Transforms (FFTs) 
#. Fast Discrete Cosine Transforms (FDCTs) 
#. Hankel-Transforms  

SPIDA also contains adaptive-step Runge-Kutta numerical propagators for stiff partial-differential-equations (PDEs).
These are the same methods as implemented for the python package `rkstiff <https://github.com/whalenpt/rkstiff>`_.
Two types of numerical methods are available, integrating factor methods (IF), and exponential time-differencing (ETD).
In particular, the solvers provided are

* ETD35 (5th order ETD with 3rd order embedding)
* ETD34 (4th order ETD with 3rd order embedding) 
* IF34 (4th order IF with 3rd order embedding)
* IF45DP (5th order IF with 4th order embedding)

In general, one should prefer ETD35 as it often has the best speed and stability for diagonal systems.
Non-diagonal methods are not implemented as of yet (see the python package for this functionality).
A detailed discussion of these suolvers is provided in the journal article |article|_.

 .. _article: https://www.sciencedirect.com/science/article/pii/S0021999114006743

 .. |article| replace:: *Exponential time-differencing with embedded Runge–Kutta adaptive step control*

Dependencies
------------

Libraries built with

* `KISS FFT <https://github.com/mborgerding/kissfft>`_
* `Nayuki-Lee DCT <https://www.nayuki.io/page/fast-discrete-cosine-transform-algorithms>`_ 
* `pwutils <https://github.com/whalenpt/pwutils>`_
* `Boost <https://www.boost.org>`_
 
Testing done with

* `GoogleTest <https://github.com/google/googletest>`_

Compile and build

* `CMake <https://cmake.org>`_

Usage
-----

Check out the demos. These can be built by configuring CMake with

.. code-block:: none

    cmake -S . -B build -DCMAKE_DEMOS=ON

License
-------
This project is licensed under the MIT License - see the `LICENSE <./LICENSE>_` file for details.
Third-party package dependencies use MIT or similarly permissive licenses

Contact
-------
Patrick Whalen - whalenpt@gmail.com



