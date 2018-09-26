spinSaturne Test
================

Description
-----------
This test, spinSaturne, was derived by François Méot from the file
  `xing/Gg7-NuZ_fromVAXSaturne/grorud_2012-Aug.res`.

For more information, see Part C, Section 5, of the [Zgoubi Users' Guide]:
  "Multiturn Spin Tracking in SATURNE 3 GeV Synchrotron"

Files
-----
* Test execution script: `zgoubi-test.sh`
* Principal expected result: `spinSaturne.res.expected`
  - compare with zgoubi.res
  - includes the required input (equivalent to a zgoubi.dat)
* Other expected results: `spinSaturne.SPNPRT.Out.expected`
* gnuplot command script:  `spinSaturne_SzT-gnuplot.cmd`
  - called by a Zgoubi SYSTEM command
 
Expected Output Figure
----------------------
![spinsaturne_szt-save](https://user-images.githubusercontent.com/13108868/46051414-3ba1a980-c0ee-11e8-97ac-ba6a8791a6d4.png)

[Zgoubi Users' Guide]: https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=2ahUKEwjBuJOuvtfdAhUnGDQIHTimDHoQFjAAegQIBRAC&url=https%3A%2F%2Fwww.bnl.gov%2Fisd%2Fdocuments%2F79375.pdf&usg=AOvVaw2b0VWKFd8Hvwrey9kXPXKq
