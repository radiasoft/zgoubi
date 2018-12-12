Automated source analysis and processing
========================================
This script 
  1. Extracts the DIMENSION statements from all *.f files in the present
     working directory, 
  2. Extracts the array names that use the passed argument ($1) as an array dimension,
  3. Reports whether or not any such array appears in a COMMON block inside ../include.

Usage
----- 
```bash
./scripts/check-common-equivalence.sh [dimension_name] 
```
Example: `cd zgoubi && ./scripts/check-common.sh MMAP`
