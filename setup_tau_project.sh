#!/usr/bin/env bash

# Script to setup a zgoubi TAU Commander project

set -o nounset
set -o errexit
set -o pipefail
set -o errtrace

err_report() {
    local error_code
    error_code=${?}
    echo "Error in $(basename "$0") in function ${1} on line ${2}" >&2
    exit ${error_code}
}
trap 'err_report "${FUNCNAME:-.}" ${LINENO}' ERR

check_for_tau () {
    if TAU=$(command -v tau) ; then
	export TAU
	echo "TAU Commander found at $TAU"
    else
	if [ -d "${HOME}/taucmdr/bin" ]; then
	    export TAU="${HOME}/taucmdr/bin/tau"
	else
	    echo "TAU Commander could not be located! Try installing it with ./install-update-tau.sh" >&2
	    exit 1
	fi
    fi
}

check_for_project () {
    if [ -d "./.tau" ]; then
	"$TAU" dashboard
	echo "TAU Commander project appears to be setup. Please see the dashboard above for info."
	return 0
    else
	echo "TAU Commander appears to not have been initialized yet."
	return 2
    fi
}

setup_project () {
    echo "Initializing the TAU Commander project with the following command:"
    init_cmd=("$TAU" initialize --linkage dynamic --compilers GNU --application-name zgoubi)
    echo "    ${init_cmd[*]}"
    "${init_cmd[@]}"

    echo "Removing UN-wanted measurement configurations..."
    for measurement in profile trace source-inst ; do
	"$TAU" measurement delete "$measurement"
    done

    select_file="$(pwd)/zgoubi-select.tau"
    if [ -f "$select_file" ]; then
	echo "Setting up compiler instrumentation measurements..."
	"$TAU" measurement edit compiler-inst --select-file "$select_file"
    fi
    select_no_io="$(pwd)/zgoubi-no-io-select.tau"
    if [ -f "$select_no_io" ]; then
	echo "Adding a measurement to more aggressively ignore IO heavy procedures..."
	"$TAU" measurement copy compiler-inst comp-inst-no-io --select-file "$select_no_io"
    fi
    "$TAU" select baseline
    "$TAU" dashboard
    echo "Finished setting up the TAU Commander project for zgoubi. Please see the dashboard above."
}

main () {
    check_for_tau

    if ! check_for_project ; then
	setup_project
    fi
}

main "$@"
