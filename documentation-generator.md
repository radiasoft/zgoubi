---
project: zgoubi
summary: charged-particle motion in magnetic electric fields
src_dir: toolbox
src_dir: zgoubi
src_dir: zpop
src_dir: common
src_dir: exemples
output_dir: doc
preprocess: true
display: public
         protected
         private
source: true
graph: true
sort: alpha
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
project_github: https://github.com/radiasoft/zgoubi 
project_download: https://github.com/radiasoft/zgoubi/releases
license: gpl
---

[This document is a FORD project file, formatted with Pythonic Markdown                                      ]:#
[See https://github.com/cmacmackin/ford/wiki/Project-File-Options for more info on writing FORD project files]:#

[Descriptions of some of the above commands:]:#
[source: display source code corresponding to item being documented]:#
[graph: generate call graphs, module dependency graphs, derive type composition/inheritance graphs ]:#
[sort: different sorting schemes for the modules or procedures or programs or derived types (alpha = alphabetical see wiki).]:#
[extra_mods: documentation for intrinsic modules]:#

--------------------

Brief description
-----------------

@warning
This archive uses advanced features of Fortran 2018.

### ToDo

X marks indicated documented all files that do not require modification for ford to processes:

 - [X] [zgoubi](./zgoubi)
 - [X] [zpop](./zpop)
 - [ ] [common](./common)
 - [ ] [exemples](./exemples)
 - [ ] [modules](./modules)
 - [ ] [toolbox](./toolbox)


Compilers
---------

This archive has been tested with the GNU Fortran compiler ([gfortran](https://gcc.gnu.org)) 8.1.0 with the parallel
runtime library from [OpenCoarrays](http://www.opencoarrays.org).
