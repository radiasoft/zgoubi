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

Examples
--------

These examples are maintained. Efforts are made to ensure that they
are executable with the latest version of zgoubi on the
author's `master` branch. In particular, some are part of the
tools/compare procedure which tests their repeatability.

The example sub-folders found in `exemple` have the similar following
content:

- one or more zgoubi input data files for a given example (a file with
  suffix `.dat`), and the corresponding result files (*.res suffix),
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
