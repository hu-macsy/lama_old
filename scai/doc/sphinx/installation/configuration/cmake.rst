.. _cmake:

CMake
^^^^^

|CMake| is an open-source, cross-platform family of tools designed to build, test and package software.

.. |CMake| raw:: html

  <a href="https://cmake.org/" target="_blank">CMake</a>

LAMA fully relies on CMake as build system. It should find all other LAMA dependent software packages, so the installation of LAMA takes places by only configuring cmake and starting the build process.

You generally run CMake by calling it with the path to the root source directory and it generates by default on Unix systems (Linux and Mac) **Makefiles** and on Windows **Visual Studio projects**. CMake also supports other generators (Eclipse CDT, Ninja), that can be specified by ``-G"[Generator-Name]"``. Please refer to CMakes webpage for further informations about that.

CMake also comes with a curses interface ``ccmake`` or a ``cmake-gui`` utility, i.e. a GUI where project configuration settings
can be specified interactively. We absolutely recommend its use for the LAMA configuration as it is very easy to adapt the
settings in case of problems. These GUI utilities are not always part of the ``cmake`` package and might be installed separately.

