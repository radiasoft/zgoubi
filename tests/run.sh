#!/usr/bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

run_cleanup() {
    if [[ -f $run_pariz_h_orig ]]; then
        mv -f "$run_pariz_h_orig" "$run_pariz_h" || true
    fi
}

run_main() {
    # globals
    local run_src_d=$PWD
    local run_build_d=$run_src_d/../build
    local run_out_d=$run_build_d/tests
    local run_warm_snake_gz=$run_build_d/tests/warmSnake.map.gz
    local run_pariz_h=$run_src_d/../zgoubi/PARIZ.H
    local run_pariz_h_orig=$run_pariz_h.zit-orig
    local run_zgoubi=$run_build_d/zgoubi

    mkdir -p "$run_out_d"
    run_pariz_h_setup
    run_warm_snake
    local t s e
    for t in [1-9]*; do
        local s=$run_src_d/$t
        cp "$s/PARIZ.H" "$run_pariz_h"
        cd "$run_build_d"
        make ${MAKEFLAGS:+-}$MAKEFLAGS >& /dev/null
        cd "$run_out_d"
        rm -rf "$t"
        mkdir "$t"
        cd "$t"
        gunzip -d < "$run_warm_snake_gz" > warmSnake.map
        cp "$s/zgoubi.dat" .
        e=$s/zgoubi.res
        "$run_zgoubi"
        (( header=$(wc -l < "$e") - 10 ))
        diff <(head -n "$header" zgoubi.res) <(head -n "$header" "$e")
        cd "$run_out_d"
    done
    echo ZGOUBI_PASSED
}

run_pariz_h_setup() {
    if [[ -f "$run_pariz_h_orig" ]]; then
        echo "$run_pariz_h_orig exists and shouldn't"
        return 1
    fi
    trap run_cleanup EXIT
    mv "$run_pariz_h" "$run_pariz_h_orig"
}

run_warm_snake() {
    local e=$(md5sum "$run_warm_snake_gz" 2>&1 || true)
    if [[ ${e[0]} == 5e1206f09fcc3ac85ce9f4346117f074 ]]; then
        return
    fi
    curl -L -S -s -o "$run_warm_snake_gz" \
         'https://sourceforge.net/p/zgoubi/code/HEAD/tree/trunk/exemples/usersGuide/FIT-and-REBELOTE/warmSnake.map.gz?format=raw'
}

run_main "$@"
