#!/usr/bin/env bash
#
# Portions of this script derive from BASH3 Boilerplate and are distributed under
# the following license:
#
# The MIT License (MIT)
#
# Copyright (c) 2014 Kevin van Zonneveld
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#
#  - https://github.com/kvz/bash3boilerplate
#  - http://kvz.io/blog/2013/02/26/introducing-bash3boilerplate/
#
# Version: 2.0.0
#
# Authors:
#
#  - Kevin van Zonneveld (http://kvz.io)
#  - Izaak Beekman (https://izaakbeekman.com/)
#  - Alexander Rathai (Alexander.Rathai@gmail.com)
#  - Dr. Damian Rouson (http://www.sourceryinstitute.org/) (documentation)
#
# Licensed under MIT
# Copyright (c) 2013 Kevin van Zonneveld (http://kvz.io)

# The invocation of bootstrap.sh below performs the following tasks:
# (1) Import several bash3boilerplate helper functions & default settings.
# (2) Set several variables describing the current file and its usage page.
# (3) Parse the usage information (default usage file name: current file's name with -usage appended).
# (4) Parse the command line using the usage information.

### Start of boilerplate -- do not edit this block #######################
export ZGOUBI_SRC_DIR="${ZGOUBI_SRC_DIR:-${PWD%/}}"
if [[ ! -f "${ZGOUBI_SRC_DIR}/zgoubi/zgoubi_main.f" ]]; then
  echo "Please run this script inside the top-level zgoubi source directory or "
  echo "set ZGOUBI_SRC_DIR to the zgoubi source directory path."
  exit 1
fi
export BASH3_BOILERPLATE="${BASH3_BOILERPLATE:-${ZGOUBI_SRC_DIR}/scripts/bash3-boilerplate}"
if [[ ! -f "${BASH3_BOILERPLATE:-}/bootstrap.sh" ]]; then
  echo "Please set BASH3_BOILERPLATE to the bash3boilerplate directory path."
  exit 2
else
    source "${BASH3_BOILERPLATE}/bootstrap.sh" "$@"
fi
### End of boilerplate -- start user edits below #########################

# Set expected value of present flags that take no arguments
export __flag_present=1

# Set up a function to call when receiving an EXIT signal to do some cleanup. Remove if
# not needed. Other signals can be trapped too, like SIGINT and SIGTERM.
function cleanup_before_exit () {
  info "Cleaning up. Done"
}
trap cleanup_before_exit EXIT # The signal is specified here. Could be SIGINT, SIGTERM etc.


### Validation (decide what's required for running your script and error out)
#####################################################################

[ -z "${LOG_LEVEL:-}" ] && emergency "Cannot continue without LOG_LEVEL. "

if [[ "${arg_v}" == "${__flag_present}" ]]; then
   print_debug_only=7
   if [ "$(( LOG_LEVEL < print_debug_only ))" -ne 0 ]; then
     debug "Supressing info and debug messages: -v or --version  present."
     suppress_info_debug_messages
   fi
fi

### Print bootstrapped magic variables to STDERR when LOG_LEVEL
### is at the default value (6) or above.
#####################################################################
{
info "__file: ${__file}"
info "__dir: ${__dir}"
info "__base: ${__base}"
info "__os: ${__os}"
info "__usage: ${__usage}"
info "LOG_LEVEL: ${LOG_LEVEL}"

info  "-c (--with-c):           ${arg_c}"
info  "-d (--debug):            ${arg_d}"
info  "-e (--verbose):          ${arg_e}"
info  "-f (--with-fortran):     ${arg_f}"
info  "-h (--help):             ${arg_h}"
info  "-i (--install-prefix):   ${arg_i}"
info  "-j (--num-threads):      ${arg_j}"
info  "-m (--with-cmake):       ${arg_m}"
info  "-n (--no-color):         ${arg_n}"
info  "-v (--version):          ${arg_v}"
}
# This file is organized into three sections:
# 1. Command-line argument and environment variable processing.
# 2. Function definitions.
# 3. Main body.
# The script depends on several external programs, including a second script that
# builds prerequisite software.  Building prerequisites requires network access
# unless tar balls of the prerequisites are present.

# TODO:
# 1. Collapse the body of the main conditional branches in the find_or_install function
#    into one new function.
# 2. Verify that script-installed packages meet the minimum version number.
# 3. Add a script_transfer function to collapse the stack_pop x; stack_push z y
#    pattern into one statement
# 4. Consider adding mpich and cmake to the dependency stack before passing them to
#    find_or_install to make the blocks inside find_or_install more uniform.
#    Alternatively, check the dependency stacks for the package before entering the
#    main conditional blocks in find_or_install.
#


# __________ Process command-line arguments and environment variables _____________

this_script="$(basename "$0")"
export this_script

export install_prefix="${arg_i%/}"
info "install_path=\"${install_prefix}\""

export num_threads="${arg_j}"
info "num_threads=\"${arg_j}\""

export build_path="${ZGOUBI_SRC_DIR}"/build
info "build_path=\"${ZGOUBI_SRC_DIR}\"/build"

# ________________________________ Start of the Main Body ___________________________________

# Print version information and exit
if [[ "${arg_v}" == "${__flag_present}" ]]; then
  zgoubi_version=$(sed -n '/[0-9]\{1,\}\(\.[0-9]\{1,\}\)\{1,\}/{s/^\([^.]*\)\([0-9]\{1,\}\(\.[0-9]\{1,\}\)\{1,\}\)\(.*\)/\2/p;q;}' "${ZGOUBI_SRC_DIR}/VERSION")
  if [[ "${arg_v}" == "${__flag_present}" ]]; then
    echo "zgoubi ${zgoubi_version}"
  fi
  exit 0
fi

if [[ $(uname) == "Darwin" ]]; then
  # If macOS , install Xcode command line tools (CLT) if not already installed
  source "${ZGOUBI_SRC_DIR}/scripts/install-xcode-clt.sh"
  maybe_install_xcodeCLT
  info "Untested operating system detected: macOS.  Please report any issues at https://github.com/radiasoft/zgoubi."
elif grep -s -q Microsoft /proc/version  ; then
  # Detect and report Windows Subsystem for Linux 
  info "Untested operating system detected: Windows Subsystem for Linux.  Please report any issues at https://github.com/radiasoft/zgoubi/issues."
fi

# Build and install zgoubi and zpop executable programs
installation_record=install-zgoubi.log
source "${ZGOUBI_SRC_DIR}/scripts/set_SUDO_if_needed_to_write_to_directory.sh"
set_SUDO_if_needed_to_write_to_directory "${install_prefix}"

cd ${ZGOUBI_SRC_DIR}
if [[ -d build ]]; then
  rm -rf build
fi
mkdir build
cd build 

if ! type cmake >& /dev/null; then
  emergency "Invoking cmake failed. Please insure cmake is installed and in your PATH.\n"
fi

export CMAKE=${arg_m:-cmake}
FC=${arg_f:-gfortran} CC=${arg_c:-gcc} $CMAKE .. -DCMAKE_INSTALL_PREFIX="$install_prefix" 
make -j ${num_threads:-}
${SUDO:-} make install
if [[ -x "${install_prefix}/bin/zgoubi" && -x "${install_prefix}/bin/zpop" ]]; then
  info "zgoubi and zpop are installed in ${install_prefix}"
else
  emergency "zgoubi and/or zpop failed to install in ${install_prefix}. Please report an issue at https://github.com/radiasoft/zgoubi/issues."
fi
# ____________________________________ End of Main Body ____________________________________________
