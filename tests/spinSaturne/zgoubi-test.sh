#!/usr/bin/env bash
# Execute a zgoubi input file for testing

## usage message
usage='usage: zgoubi-test.sh [-ehHqv] testName'
usageF=$usage$'\n'
usageF=$usageF$'Execute a given test for zgoubi.\n'
#usageF=$usageF$'  -e "ext ..." A quoted list of file extensions.\n'
usageF=$usageF$'  -h   Print a brief help message and exit.\n'
usageF=$usageF$'  -H   Print this help message and exit.\n'
usageF=$usageF$'  -q   Suppress  output from zgoubi.\n'
usageF=$usageF$'  -v   Print all differences (verbose output).\n'

## default flags
quiet=false  # quiet zgoubi
verb=false   # verbose diffs

## initialize the exit status
xstat=0

## parse optional arguments
while getopts "hHqv" opt; do
  case $opt in
  h) echo "$usage"; exit;;
  H) echo "$usageF"; exit;;
  q) quiet=true;;
  v) verb=true;;
  \?) echo "$usage"; exit;;
  esac
done
## strip off the optional argumnents
shift $(( OPTIND - 1 ))

## we require exactly one non-optional argument
if (( $# != 1 )); then
  echo "zgoubi-test.sh: You must supply exactly one test name."
  echo "$usage"
  exit
fi
testName=$1
echo " === Test ${testName} for zgoubi ==="

#cd "@ZGOUBI_TEST_BIN_DIR@" # taken from the corresponding CMake variable
#                           # == argument WORKING_DIRECTORY of add_test()

## execute zgoubi
# zgoubi executable
#zgoubi_exe="@CMAKE_BINARY_DIR@"/zgoubi
zgoubi_exe=/Users/dabell/Projects/zgoubi/zgoubi-code/zgoubi/zgoubi
#zgoubi_exe=/Users/dabell/Applications/bin/zgoubi
# leading part of testName.res.expected contains the required input
zgoubi_cmd="$zgoubi_exe -in ${testName}.res.expected"
# run zgoubi command
if [ "$quiet" = true ] ; then
  zgoubi_res=$(eval $zgoubi_cmd)
else
  eval $zgoubi_cmd
fi

## compare zgoubi.res with the expected result
# we compare all but the last ten lines, which contain date and timing information
footer=10
(( header=$(wc -l < ${testName}.res.expected) - footer ))
result=$(diff <(head -n $header ${testName}.res.expected)  <(head -n $header zgoubi.res))
stat=$?
if [ "$stat" -ne 0 ] ; then
  (( ++xstat ))
  echo "File zgoubi.res differs from ${testName}.res.expected."
  if [ "$verb" = true ] ; then
    echo "$result"
  fi
fi

## compare other files: ${testName}.fe
exts="SPNPRT.Out"
for fe in $exts; do
  result=$(diff ${testName}.${fe}.expected zgoubi.${fe})
  stat=$?
  if [ "$stat" -ne 0 ] ; then
    (( ++xstat ))
    echo "File zgoubi.${fe} differs from ${testName}.${fe}.expected."
    if [ "$verb" = true ] ; then
      echo "$result"
    fi
  fi
done

if [ "$xstat" -eq 0 ]
then
   echo "Test ${testName} passed."
else
   echo "Test ${testName} failed!"
fi

