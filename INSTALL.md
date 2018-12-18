Building with cmake
-------------------

First obtain and setup the sirepo envirionment following the
instructions at https://github.com/radiasoft/sirepo/wiki/Development.

Once up and running within the development environment (e.g., after
`vagrant ssh`) follow these steps:

1. Clone zgoubi:

```bash
git clone https://github.com/radiasoft/zgoubi
```

2. Build, test and install zgoubi

```bash
cd zgoubi
rm -rf build
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX="$(pyenv prefix)" ..
make -j $(NPROC)
ctest --output-on-failure
make install
```

The [travis.sh](./travis.sh) script shows an analogous process we use
to build and test zgoubi on
[Travis-CI.org](https://travis-ci.org/radiasoft/zgoubi/branches) via
the radiasoft/beamsim docker image.

Building with make
------------------

Type `make -f <my_Makefile>`, with `<my_Makefile>` one of
`Makefile_zgoubi_gfortran`, `Makefile_zgoubi_ifort`.

`Makefile_zgoubi_gfortran` will build zgoubi, zpop and zgoubi users'
guide.

`Makefile_zgoubi_ifort` will build zgoubi and zgoubi users' guide, it
will not build `zpop`. The latter will require `make -f
Makefile_zpop_ifort`.
