Building with cmake
-------------------

1. Download and install [vagrant](https://www.vagrantup.com).

2. Download and install [VirtalBox](https://www.virtualbox.org).

3. Obtain and launch the RadiaSoft development environment:

```bash
curl radia.run | bash -s vagrant-sirepo-dev
mkdir vagrantdir
vagrant init
vagrant up
vagrant ssh
```

4. Clone zgoubi:

```bash
git clone https://github.com/radiasoft/zgoubi
```

5. Build, test and install zgoubi

```bash
cd zgoubi
rm -rf build
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX="$(pyenv prefix)" ..
make -j
ctest
make install
```
 
Building with make
------------------

Type 'make -f my_Makefile', 
with 'my_Makefile' one of Makefile_zgoubi_gfortran, Makefile_zgoubi_ifort. 

Makefile_zgoubi_gfortran will build zgoubi, zpop and zgoubi users' guide. 

Makefile_zgoubi_ifort will build zgoubi and zgoubi users' guide, it will not build zpop. The latter will require 'make -f Makefile_zpop_ifort'. 


