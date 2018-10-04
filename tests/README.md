[This document is formatted with GitHub-Flavored Markdown. ]:#
[For better viewing, including hyperlinks, read it online at ]:#
[https://github.com/radiasoft/zgoubi/blob/build-test-infrastructure/tests/spinSaturne/README.md ]:#

zgoubi Test Suite
=================

To add a new test,

1. Create a new directory inside the tests directory.
2. Add the new directory name to the `test_directories` list in tests/CMakeLists.txt.
3. Add input files following the exmaple of either the warmSnake or spinSaturne test.
4. Create a CMakeLists.txt file in the new directory by one of two methods:
   a. Copy CMakeLists.txt from warmSnake if the test requires recompilation and uses
      a field map.
   b. Copy CMakeLists.txt from spinSaturne otherwise.
