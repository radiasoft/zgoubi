#!/usr/bin/env bash
set -euo pipefail
docker run -i --rm -u vagrant -v "$PWD":/home/vagrant/src/radiasoft/zgoubi "$1" bash -l <<'EOF'
#!/bin/bash
set -xeuo pipefail
cd ~/src/radiasoft/zgoubi
# Needs to match what tests assume
build_d=build
if [[ -d $build_d ]] ; then
    echo "$PWD/$build_d should not exist"
    exit 1
fi
mkdir "$build_d"
cd "$build_d"
cmake -Wdev -DCMAKE_INSTALL_PREFIX=$HOME/.local ..
CTEST_OUTPUT_ON_FAILURE=1 make -j $(nproc) test
# do not install: "make test" creates a new zgoubi that
# has a non-default PARIZ.H. Probably don't need to test
# "install" here.
EOF
