#!/usr/bin/env bash
#
# Script for installing or updating TAU Commander
# taucommander.com
# https://github.com/ParaToolsInc/taucmdr

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

init_defaults () {
    export TAU_REPO='https://github.com/ParaToolsInc/taucmdr'

    export TAU_BRANCH=unstable
    export CODE_DIR="${HOME}/src"
    if [[ -z "${USE_MINICONDA:-}" ]] ; then
       export USE_MINICONDA=false
    fi
    export INSTALLDIR="${HOME}/taucmdr"
}

print_usage () {
    usage="
$(basename "$0") [<options>] -- install TAU Commander

Options:
	-b|--branch <branch|tag> Install from the user specified branch or tag. Default: unstable
	-c|code-dir <path>       Directory in which to clone TAU Commander.     Default: ~/src/
	-h|--help                Print this help message and exit
	-m|--use-miniconda       Use vendored Miniconda instead of current env. Default: false
	-p|--prefix <path>       Installation prefix for TAU Commander          Default: ~/taucmdr
"
    echo "$usage"
}

error_if_arg_missing () {
    if [[ -z "${2:-}" ]] ; then
	echo "Command line flag, ${1:-'<key>'} requires an argument!" >&2
	print_usage >&2
	exit 1
    fi
}

parse_opts () {
    while [[ $# -gt 0 ]]
    do
	key="$1"

	case $key in
	    -b|--taucmdr-branch)
		error_if_arg_missing "$key" "${2}"
		export TAU_BRANCH="${2}"
		shift
		shift
		;;
	    -c|--code-dir)
		error_if_arg_missing "$key" "${2}"
		export CODE_DIR="${2}"
		shift
		shift
		;;
    	    -h|--help)
		print_usage
		exit 0
		;;
	    -m|--use-miniconda)
		export USE_MINICONDA=true
		shift
		;;
	    -p|--prefix)
		error_if_arg_missing "$key" "${2}"
		export INSTALLDIR="${2}"
		shift
		shift
		;;
	    *)
		print_usage >&2
		exit 1
		;;
	esac
    done
}

clone_if_missing () {
    if [ ! -d "${CODE_DIR}/taucmdr" ] ; then
	git clone --branch="${TAU_BRANCH}" "${TAU_REPO}" "${CODE_DIR}/taucmdr" || exit 1
    fi
}

update_tau_src () {
    git checkout "${TAU_BRANCH}"
    git pull origin
}

do_install () {
    make clean || true
    make INSTALLDIR="$INSTALLDIR" USE_MINICONDA="$USE_MINICONDA" install
}

main () {
    init_defaults

    parse_opts "$@"

    clone_if_missing || exit 1

    (cd "${CODE_DIR}/taucmdr" || exit 1
     update_tau_src || exit 1
     do_install || exit 1
    )
}

main "$@"
