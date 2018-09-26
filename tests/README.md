[This document is formatted with GitHub-Flavored Markdown. ]:#
[For better viewing, including hyperlinks, read it online at ]:#
[https://github.com/radiasoft/zgoubi/blob/build-test-infrastructure/tests/spinSaturne/README.md ]:#

zgoubi Test Suite
=================

To add a new test,

1. Create a new subdirectory inside the tests subdirectory.
2. Add the name of the new folder to the `test_directories` list in tests/CMakeLists.txt.
3. Create a new CMakeLists.txt file in the new test subdirectory (copy an existing CMakeLists.txt from another test).
4. Possibly edit the bash script (hopefully this step can be eliminated by adding the ability to give the bash
   script arguments specifying paths to files that the script is to use ndiff to compare).
