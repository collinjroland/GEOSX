=================================
Getting Started
=================================

----------
 Welcome
----------

Welcome to GEOSX. In this section, we offer a step-by-step process to install and run GEOSX. The procedure is as follows:

  - 1: Check prerequisites and install external dependencies.
  - 2: Clone GEOSX and the Third Party Libraries (TPL) from GitHub.
  - 3: Compile TPL and GEOSX, in that order.
  - 4: Run simple GEOSX tests to ensure that everything works as intended.

Exact instructions may vary based on individual configurations, but the essential steps should be identical.


Prerequisites
=================================

First, we make sure that all external dependencies are installed. The following packages are required. The versions specified are currently supported, they are not necessarily minimal.


 * A C++ compiler (gcc >=7, clang, Apple LLVM >= 10.0.1),
 * A Fortran compiler (GNU Fortran >= 9.1.0) for some components of the third party libraries (TPL),
 * Git (>= 2.22.0), to clone and track development versions,
 * Git Lfs (>= 2.7.2), for large files Git support,
 * CMake (>= 3.14.5), to create cross-platform project files,
 * MPI such as open-mpi (>= 4.0.1), for parallelization.

The installation of these packages depends on your configuration.

On Linux, the `apt <https://wiki.debian.org/Apt>`__ package manager syntax is shown as an example.

.. code-block:: sh

  sudo apt install git git-lfs gcc g++ gfortran python cmake zlib1g-dev libblas-dev liblapack-dev libopenmpi-dev


For Mac OS, all packages can be installed using `Homebrew <https://docs.brew.sh/Installation>`__. If Homebrew is not installed on your system, and if you are getting ready to install it, it is important to install the Xcode command line tools (as stated on Homebrew's website, but this step can easily be missed).

To install the Xcode command line tools:

.. code-block:: sh

  xcode-select --install



To verify that that the Xcode command line tools are installed and their version:

.. code-block:: sh

  xcode-select -v

We are currently using xcode-select version 2354.

.. code-block:: sh

  brew install git git-lfs gcc python cmake libomp open-mpi



Accessing GitHub
--------------------

GEOSX resides in a git repository hosted at https://github.com/GEOSX/GEOSX. To download the code, it is recommended to setup and use ssh keys as discussed
`here <https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/>`__.

You can test that your SSH configuration works properly `here <https://help.github.com/en/articles/testing-your-ssh-connection>`__.

Alternatively, it is possible to use a less secured https tokens, as described `here <https://help.github.com/en/articles/git-automation-with-oauth-tokens>`__.

Cloning the Third Party Libraries and GEOSX
==================================================================

Before starting, let us create a directory to host the various clones required for an effective development workflow.

1. Setup working directory

.. code-block:: sh

  mkdir geosx
  cd geosx


There are currently two separate repositories that should be downloaded.
The first is the main repository, which may be cloned and initialized by the following steps:

2. Clone the main repository

If you have successfully setup the SSH authentication:

.. code-block:: sh

   git clone git@github.com:GEOSX/GEOSX.git


Otherwise, if you use the https protocol:

.. code-block:: sh

   git clone https://github.com/GEOSX/GEOSX.git

Then:

.. code-block:: sh

  cd GEOSX
  git lfs install
  git submodule init
  git submodule update
  cd ..



3. Clone the third-party libraries

.. code-block:: sh

   git clone git@github.com:GEOSX/thirdPartyLibs.git
   cd thirdPartyLibs
   git lfs install
   git pull
   git submodule init
   git submodule update
   cd ..

Note that git-lfs may not funct-on properly (or may be very slow) if version of git and git-lfs are not current.
If you are using an older version of git/git-lfs you may need to add "git lfs pull" after "git pull" in the above procedures.

Compiling the Code
=================================

GEOSX compilations are typically driven by a hostconfig file, which reside in GEOSX/host-configs.
If your platform does not have a host-config in the repository, you are encouraged to maintain one.
If you are running on an LC system, there is already a hostconfig and copy of the thirdPartyLibs installed, and you can skip step 4.

The first step in compiling GEOSX is to run cmake and generate the makefiles.
Starting with the third-party libraries, the config-build.script will run cmake for you.
Note that the 'make' step should be run serially, as the indiviudal package builds are run in parallel by default.

4. Configure and make the third party libraries

.. code-block:: sh

   cd thirdPartyLibs
   python scripts/config-build.py -hc ../GEOSX/host-configs/your-platform.cmake -bt Release
   cd build-your-platform-release
   make -j1

The next step is to compile the main code.
Again, the config-build sets up cmake for you.

5. Configure and make the main code

.. code-block:: sh

   cd ../../GEOSX
   python scripts/config-build.py -hc host-configs/your-platform.cmake -bt Release
   cd build-your-platform-release
   make -j4


Running the Code
=================================

GEOSX executables read in a XML input file. A simple example XML is located
`here <https://github.com/GEOSX/GEOSX/blob/develop/src/components/core/tests/PhysicsSolvers/LaplaceFEM.xml/>`__.
To execute a serial run enter the following command from a working directory:

.. code-block:: sh

    path-to-geosx-bin/geosx -i ./GEOSX/src/coreComponents/physicsSolvers/SimpleSolvers/integratedTests/10x10x10_LaplaceFEM.xml
