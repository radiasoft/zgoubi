#!/usr/bin/env bash

# Perform CI testing in the desired docker image passed as the first argument
# default to radiasoft/beamsim-part1 if no first argument is passed

set -o pipefail
set -o nounset
set -o errexit

docker run -i --rm -u vagrant -v "$PWD":/home/vagrant/src/radiasoft/zgoubi "${1:-radiasoft/beamsim-part1}" bash  <<-'EOF'
    #!/usr/bin/env bash
    source ~/.bashrc
    set -euo pipefail
    cd ~/src/radiasoft/zgoubi
    if ! mkdir build > /dev/null 2>&1 ; then
        echo \
            "Warning: Using an old/dirty build directory! Please run 'rm -r build' and try again if script fails." >&2
    fi
    cd build
    cmake -Wdev -DCMAKE_INSTALL_PREFIX=$(pyenv prefix) ..
    make -j $(nproc)
    ctest --output-on-failure
    make install
EOF
