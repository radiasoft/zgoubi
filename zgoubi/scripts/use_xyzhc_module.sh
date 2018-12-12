#!/usr/bin/env bash
#
#  use_xyzhc_module.sh
#
# --- Replace XYZHC.H header file inclusion with use association

# Exit on error. (Append ||true if an non-zero exit code is acceptable.)
set -o errexit

# Return an error upon the use of an unset variable
set -o nounset

# Remember & return the highest exit code in a chain of pipes.
# This way you can catch the error in case mysqldump fails in `mysqldump |gzip`
set -o pipefail

sed -e 's/INCLUDE\ \"XYZHC.H\"/USE xyzhc/' *.f
sed -e "s/INCLUDE\ \'XYZHC.H\'/USE xyzhc/" *.f
