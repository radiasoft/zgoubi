[This document is formatted with GitHub-Flavored Markdown. ]:#
[For better viewing, including hyperlinks, read it online at ]:#
[https://github.com/radiasoft/zgoubi/blob/master/README.md ]:#

zgoubi
======

Building in Docker
------------------

To download and build zgoubi and execute the tests on a system with
Docker installed and the Docker service running, execute the following
in a command-line shell such as bash:
```bash
git clone https://github.com/radiasoft/zgoubi
cd zgoubi
./travis.sh
```
This will use [cmake] to build the master branch of the zgoubi repository.
Then it will use the companion package `ctest` to run the zgoubi test
suite.  Please submit an [issue] describing any problems you
encounter.

Launching an Interactive Docker Session
---------------------------------------

For an interactive, command-line interface to the docker image,
execute the following docker command:

```bash
docker run -it --rm -u vagrant -v "$PWD":/home/vagrant/src/radiasoft/zgoubi "${1:-radiasoft/beamsim}" bash
```

Building from Source
--------------------

To build Zgoubi from source, there are three prerequisites:
  + an up-to-date version of CMake (v 3.10+)
  + an up-to-date Fortran compiler (we recommend gfortran v.8.2+)
  + a BLAS/LAPACK library built with that compiler

If you are the sort of person planning to build Zgoubi from source,
you probably already have CMake installed. But if not, head on over
to https://cmake.org/ to download the latest release for your system.

To get Fortran with the latest Fortran 2018 features, we recommend
installing via [OpenCoarrays](http://www.opencoarrays.org/).
On a typical Linux system, the following steps should do the job
of installing gfortran together with support for coarrays and other
Fortran 2018 features.
Before you start, make sure you have a working internet connection.
```
sudo su -l        # you're logged in as root now, so watch your steps!!
umask 022         # ensure correct permissions for installation
mkdir sw & cd sw  # create and enter a directory for building software

# now get the latest OpenCoarrays
git clone https://github.com/sourceryinstitute/OpenCoarrays.git
cd OpenCoarrays

# build gcc/gfortran (c.15 minutes)
# if necessary, modify the install path---the argument to `-i`---
# to suit your needs (and modify the exported PATH variables also)
./install.sh -j -p gcc -I 8.2.0 -i /opt/gcc/8.2.0 --disable-bootstrap
export PATH=/opt/gcc/8.2.0/bin:$PATH
export LD_LIBRARY_PATH=/opt/gcc/8.2.0/lib64:/opt/gcc/8.2.0/lib:$LD_LIBRARY_PATH

# build mpich (needed for coarray support)
# if necessary, modify the install path---the argument to `-i`---
# to suit your needs (and modify the exported PATH variable also)
./install.sh -j -p mpich -I 3.3 -i /opt/mpich/3.3/gcc/8.2.0
export PATH=/opt/mpich/3.3/gcc/8.2.0/bin:$PATH

# build opencoarrays
# if necessary, modify the install path---the argument to `-i`---
# to suit your needs (and modify the exported PATH variable also)
./install.sh -j -p opencoarrays -i /opt/opencoarrays/master/gcc/8.2.0
export PATH=/opt/opencoarrays/master/gcc/8.2.0/bin:$PATH

# build BLAS/LAPACK with the new compiler
cd ..
git clone https://github.com/reference-LAPACK/lapack.git
cd lapack
mkdir build && cd build
CC=gcc FC=gfortran cmake -DBUILD_SHARED_LIBS=TRUE ..
make -j
make install

# log out from root
exit
```

Users who wants access to the new compilers and libraries should add
the following to their `.bashrc`:
```
source /opt/opencoarrays/master/gcc/8.2.0/setup.sh
export LD_LIBRARY_PATH=/usr/local/lib64:$LD_LIBRARY_PATH
```

You are now ready to build Zgoubi! Assuming you already have the
source code downloaded, you follow the standard cmake build steps:
```
cd <path-to-Zgoubi>
rm -rf build && mkdir build && cd build
cmake ..
make -j

$ to run the tests, execute
ctest
```

That's it!

Examples
--------

These examples are maintained. Efforts are made to ensure that they
are executable with the latest version of zgoubi on the
author's `master` branch. In particular, some are part of the
tools/compare procedure which tests their repeatability.

The example sub-folders found in `exemple` have the similar following
content:

- one or more zgoubi input data files for a given example (a file with
  suffix `.dat`), and the corresponding result files (`.res` suffix),
  and sometimes in addition `*.eps` plot files or other forms of
  computation outputs,
- or just one or more result files (suffix `.res`)

Zgoubi needs a zgoubi.dat as input file. That zgoubi.dat is obtained
by copy-pasting `thisExampleDataFile.dat` or `thisExampleDataFile.res`
file in that folder, to zgoubi.dat. (that works with a .res file since
the first part of a `.res` zgoubi result file is a copy of the
zgoubi.dat file that it stems from).  Note that Zgoubi, instead, also
accepts the command `zgoubi -inFile thisExampleDataFile.dat/.res`.

Some of these examples are discussed in Part C of the [users' guide].

The folder `exemples/tools` is dedicated to the following:

- running `compare` will run some of the examples, and compare their
outcomes (logged in compare.out) with reference ones logged in
compare.out_reference.  Absence of (major) differences between
compare.out_reference and compare.out is a sign of a good installation
of zgoubi.

[CMake]: https://www.cmake.org
[issue]: https://github.com/radiasoft/zgoubi/issues/new
[users' guide]: https://www.bnl.gov/isd/documents/79375.pdf
