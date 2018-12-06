#!/usr/bin/env bash
#
# check-common.sh
#
# -- This script 
#    1. Extracts the DIMENSION statements from all *.f files in the present
#       working directory, 
#    2. Extracts the array names that use the passed argument ($1) as an array dimension,
#    3. Reports whether or not any such array appears in a COMMON block inside ../include.

# Exit on error. (Append ||true if an non-zero exit code is acceptable.)
set -o errexit

# Return an error upon the use of an unset variable
set -o nounset


# Remember & return the highest exit code in a chain of pipes.
# This way you can catch the error in case mysqldump fails in `mysqldump |gzip`
set -o pipefail

: "${1?'Usage: ./scripts/check-common.sh [dimension_name] (example: ./scripts/check-common.sh MMAP)'}"

export dimension="$1"

grep "$dimension" -i *.f | grep DIMENSION  |
{
  while IFS='' read -r line || [[ -n "${line}" ]]; do

    export static_array_list="${line#*DIMENSION}"                 # grab text after DIMENSION
    export static_array_list=${static_array_list//[[:space:]]/}   # delete all white space

    echo "static_array_list=${static_array_list}"

    while [[ ! -z "$static_array_list" ]]; do

      export array_name="${static_array_list%%(*}"                # grab array name (text before first opening parenthesis)
      export remaining_text="${static_array_list#*${array_name}}" # grab text after array name
      export dimension_list="${remaining_text%%)*}"               # grab array dimensions (text before first closing parenthesis)

      if [[  "$dimension_list" == *"$dimension"* ]]; then
        export dynamic_array_list+=("$array_name")
      else
        echo "$array_name does not have $dimension dimension"
      fi

      if [[  "$static_array_list" == *"),"* ]]; then
        export static_array_list="${static_array_list#*),}"       # shorten list to text after "),"
      else
        export static_array_list=""
      fi
    done
  done

  echo "dynamic_array_list=${dynamic_array_list[@]}"
  
  for dynamic_array in "${dynamic_array_list[@]}" ; do
    dynamic_array_in_common=$(grep COMMON -i -r ../include | grep "${dynamic_array}") || true
    if [[ -z "$dynamic_array_in_common" ]]; then
      echo "Array not found in a COMMON statement: ${dynamic_array}"
    else
      echo "Array found in a COMMON statement: ${dynamic_array}"
    fi
  done
}
