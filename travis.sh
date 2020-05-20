#!/usr/bin/env bash

# Perform CI testing in the desired docker image passed as the first argument
# default to radiasoft/beamsim-part1 if no first argument is passed

set -o pipefail
set -o nounset
set -o errexit

docker run -i --rm -u vagrant -v "$PWD":/home/vagrant/src/radiasoft/zgoubi "${1:-radiasoft/beamsim}" bash <<-'EOF'
    #!/bin/bash
    source ~/.bashrc
    set -veuo pipefail
    cd ~/src/radiasoft/zgoubi
    if [[ -d "${BUILD_DIR:-cmake-build}" ]] ; then
        echo \
            "Warning: Using an old/dirty build directory! Please run 'rm -r ${BUILD_DIR:-cmake-build}' and try again if script fails." >&2
    else
        mkdir "${BUILD_DIR:-cmake-build}"
    fi
    cd "${BUILD_DIR:-cmake-build}"
    cmake -Wdev -DCMAKE_INSTALL_PREFIX=$(pyenv prefix) ..
    make -j $(nproc)
    ctest --output-on-failure
    make install
EOF
