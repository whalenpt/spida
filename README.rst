=========
SPIDA
=========

.. |build| image:: https://github.com/whalenpt/spida/actions/workflows/cmake.yml/badge.svg
    :target: https://github.com/whalenpt/spida/actions/workflows/cmake.yml

.. |coverage| image:: https://codecov.io/gh/whalenpt/spida/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/whalenpt/spida

.. |license| image:: https://img.shields.io/badge/License-MIT-yellow.svg
    :target: https://opensource.org/licenses/MIT


|build| |coverage| |license|


Spectral integration and differentiation algorithms (SPIDA). This project implements several
high-level wrappers for spectral transforms:

#. Fast-Fourier Transforms (FFTs)
#. Fast Discrete Cosine Transforms (FDCTs)
#. Hankel-Transforms  

SPIDA also contains adaptive-step Runge-Kutta numerical propagators for stiff partial-differential-equations (PDEs).
These are the same methods as implemented for the python package rkstiff_.

Two types of numerical methods are available: integrating factor methods (IF), and exponential time-differencing (ETD).
In particular, the solvers provided are:

* ETD35 (5th order ETD with 3rd order embedding)
* ETD34 (4th order ETD with 3rd order embedding)
* IF34 (4th order IF with 3rd order embedding)
* IF45DP (5th order IF with 4th order embedding)

In general, one should prefer ETD35 as it often has the best speed and stability for diagonal systems.
Non-diagonal methods are not implemented as of yet (see rkstiff_ for this functionality).

A detailed discussion of these solvers is provided in the journal article |article|_.

.. _article: https://www.sciencedirect.com/science/article/pii/S0021999114006743

.. |article| replace:: *Exponential time-differencing with embedded Runge–Kutta adaptive step control*


Dependencies
------------

Libraries built with:

* `KISS FFT <https://github.com/mborgerding/kissfft>`_
* `Nayuki-Lee DCT <https://www.nayuki.io/page/fast-discrete-cosine-transform-algorithms>`_
* `pwutils <https://github.com/whalenpt/pwutils>`_
* `Boost <https://www.boost.org>`_ (system Boost preferred, fallback included)
* `json11 <https://github.com/dropbox/json11>`_

Testing done with:

* `GoogleTest <https://github.com/google/googletest>`_

Build system:

* CMake (CMake Presets supported)


Hankel transform usage
----------------------

Hankel transforms of order zero computed on a Bessel root grid:

.. code-block:: c

    #include <spida/grid/besselR.h>
    #include <spida/shape/shapeR.h>
    #include <spida/transform/hankelR.h>

    int main()
    {
        int nr = 200;

        spida::BesselRootGridR gridR(nr, 12*w0);
        double A0 = 1.0;
        double w0 = 1.0;

        spida::GaussR shapeR(gridR, A0, w0);
        spida::HankelTransformR transform(gridR);

        std::vector<double> u = shapeR.shapeRV();
        std::vector<double> v(nr);

        transform.R_To_SR(u, v);

        return 0;
    }


Installation
------------

If not installed, first install CMake_.

To build SPIDA using **CMake Presets**:

.. code-block:: none

    git clone https://github.com/whalenpt/spida.git
    cd spida


Build configurations:

**Release build (recommended)**

.. code-block:: none

    cmake --preset release
    cmake --build --preset release --parallel


**Debug build**

.. code-block:: none

    cmake --preset debug
    cmake --build --preset debug --parallel


**Build with demos enabled**

.. code-block:: none

    cmake --preset demos
    cmake --build --preset demos --parallel


Install:

.. code-block:: none

    cmake --install build/release


Demos
-----

Demos are enabled via the `demos` preset:

.. code-block:: none

    cmake --preset demos
    cmake --build --preset demos --parallel


Testing
-------

Testing is enabled via the `release` preset:

.. code-block:: none

    cmake --preset release
    cmake --build --preset release --parallel
    ctest --preset release --output-on-failure

License
-------

This project is licensed under the MIT License - see the LICENSE_ file for details.

Third-party dependencies use MIT or similarly permissive licenses.

Contact
-------

Patrick Whalen - whalenpt@gmail.com
