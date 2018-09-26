[This document is formatted with GitHub-Flavored Markdown. ]:#
[For better viewing, including hyperlinks, read it online at ]:#
[https://github.com/radiasoft/zgoubi/blob/master/tests/spinESRF_GeV/README.md ]:#

Spin ESRF 18 GeV Test
=====================

Description
-----------

This test, spinESRF_18GeV, was derived from earlier simulations of the
ESRF, a 6 GeV light source in Grenoble, France. Here we run at 18 GeV
so as to yield a much shorter damping time for this test. We track just
ten electrons around the 0.8km ring. Synchrotron radiation is ON (though
only in dipole fields), and spin tracking is ON. The 6D motion of the
electrons converges, on average, to the equilibrium emittance of this
ring (zero vertically). The transverse damping time is about 50 turns.
(See the [Graphs] below files for some context.)

NB: To limit the output file size, the `zgoubi` input includes the lines

```
   'OPTIONS'
  1 1
  WRITE OFF
```

which inhibit (most of the) output printed to `zgoubi.res`.


Files
-----

* Test execution script: `zgoubi-test.sh`
* Principal expected result : `spinESRF_18GeV.res.expected`
  - includes the required input (equivalent to a zgoubi.dat)
  - compare with zgoubi.res
* Other expected results: `spinESRF_18GeV.SPNPRT.Out.expected`



Graphs
------

![turn-dp](https://user-images.githubusercontent.com/13108868/46050816-047dc900-c0eb-11e8-93da-c767123cc611.png)

![turn-sz](https://user-images.githubusercontent.com/13108868/46050817-047dc900-c0eb-11e8-898f-318d6da41102.png)

![turn-y](https://user-images.githubusercontent.com/13108868/46050818-047dc900-c0eb-11e8-86a9-fcc27bc672db.png)

![turn-z](https://user-images.githubusercontent.com/13108868/46050819-047dc900-c0eb-11e8-944a-478c16c52b2b.png)


[Graphs]: #graphs
