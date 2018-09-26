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

* Test execution script: 
 - `zgoubi-test.sh`

* Principal expected result (to be compared with zgoubi.res):
 - `spinESRF_18GeV.res.expected` (includes the required input equivalent to a zgoubi.dat)

* Other expected results:
 -  `spinESRF_18GeV.SPNPRT.Out.expected`


Graphs
------

![turn-dp.pdf](https://github.com/radiasoft/zgoubi/files/2417504/turn-dp.pdf)

![turn-SZ.pdf](https://github.com/radiasoft/zgoubi/files/2417505/turn-SZ.pdf)

![turn-Y.pdf](https://github.com/radiasoft/zgoubi/files/2417506/turn-Y.pdf)

![turn-Z.pdf](https://github.com/radiasoft/zgoubi/files/2417507/turn-Z.pdf)

[Graphs]: #graphs
