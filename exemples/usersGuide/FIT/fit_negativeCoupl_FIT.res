MATCHING QUADRUPOLE TRIPLET, INCLUDING COUPLED DRIFTS
 'OBJET'                                                                                                      1
2501.73            RIGIDITY (kG.cm)
5                  11 PARTICLES GENERATED FOR USE OF MATRIX
.01  .01  .01   .01  0. .0001
0.  0.  0.   0.  0.  1.
 'ESL   '                                                                                                     2
200.
 'QUADRUPO'    3                                                                                              3
0
40.  15.  -6.
0.  0.
6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723
0.   0.
6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723
5.
1 0. 0. 0.
 'ESL'                                                                                                        4
30.
 'QUADRUPO'     5                                                                                             5
0
40.  15.  3.
0.  0.
6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723
0.   0.
6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723
5.
1 0. 0. 0.
 'ESL'                                                                                                        6
30.
 'QUADRUPO'      7                                                                                            7
0
40.  15.  -8.
0.  0.
6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723
0.   0.
6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723
5.
1 0. 0. 0.
 'ESL'                                                                                                        8
200.
 'MATRIX'                                                                                                     9
1  0
 'FIT'             ! Vary first and last straight sections                                                   10
2  save            ! # of variables; variables will save in zgoubi.FITVALS.out
2  1 -8.001  .5    ! Coupled: sum of drifts stays constant
5  12  0.   2.     ! Parmtr #12 of elements #3 is the field value
2  1.E-10 200      ! # of constraints; penalty; max. nmbr of iterations
1  1 2 8  17. 1.   ! Cnstrnt #1 : R12=17. after last drift (lmnt #8)
1  3 4 8 -.8  1.   ! Cnstrnt #2 : R34=-.8 after last drift
 'MATRIX'                                                                                                    11
1  0
 'END'                                                                                                       12

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       2501.730 kG*cm

                                         CALCUL  DES  TRAJECTOIRES

                              OBJET  (5)  FORME  DE     11 POINTS 



                                Y (cm)         T (mrd)       Z (cm)        P (mrd)       S (cm)        dp/p 
               Sampling :          0.10E-01      0.10E-01      0.10E-01      0.10E-01       0.0          0.1000E-03
  Refrnce trajectry #      1 :      0.0           0.0           0.0           0.0           0.0           1.000    

************************************************************************************************************************************
      2  Keyword, label(s) :  ESL                               


                              Drift,  length =   200.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00   2.0000000E+02   0.00000E+00

 Cumulative length of optical axis =    2.00000000     m  at  ;  Time  (for ref. rigidity & particle) =   1.068463E-08 s 

************************************************************************************************************************************
      3  Keyword, label(s) :  QUADRUPO    3                     


      -----  QUADRUPOLE  : 
                Length  of  element  =    40.000000      cm
                Bore  radius      RO =    15.000      cm
               B-QUADRUPOLE  = -6.0000000E+00 kG   (i.e.,  -6.0000000E+00 * SCAL)

               Entrance/exit field models are sharp edge
               FINTE, FINTS, gap :    0.0000E+00   0.0000E+00   7.5000E+00

                    Integration step :   5.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.010     0.000     0.000     0.000           40.000     0.011     0.000     0.000     0.000            2
  A    1  1.0000    -0.010     0.000     0.000     0.000           40.000    -0.011     0.000     0.000     0.000            3
  A    1  1.0000     0.000     0.010     0.000     0.000           40.000     0.003     0.000     0.000     0.000            4
  A    1  1.0000     0.000    -0.010     0.000     0.000           40.000    -0.003     0.000     0.000     0.000            5
  A    1  1.0000     0.000     0.000     0.010     0.000           40.000     0.000     0.000     0.009     0.000            6
  A    1  1.0000     0.000     0.000    -0.010     0.000           40.000     0.000     0.000    -0.009     0.000            7
  A    1  1.0000     0.000     0.000     0.000     0.010           40.000     0.000     0.000     0.002     0.000            8
  A    1  1.0000     0.000     0.000     0.000    -0.010           40.000     0.000     0.000    -0.002     0.000            9
  A    1  1.0001     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           10
  A    1  0.9999     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           11

 Cumulative length of optical axis =    2.40000000     m ;  Time  (for ref. rigidity & particle) =   1.282155E-08 s 

************************************************************************************************************************************
      4  Keyword, label(s) :  ESL                               


                              Drift,  length =    30.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00   2.7000000E+02   1.33426E-03

 Cumulative length of optical axis =    2.70000000     m  at  ;  Time  (for ref. rigidity & particle) =   1.442425E-08 s 

************************************************************************************************************************************
      5  Keyword, label(s) :  QUADRUPO    5                     


      -----  QUADRUPOLE  : 
                Length  of  element  =    40.000000      cm
                Bore  radius      RO =    15.000      cm
               B-QUADRUPOLE  =  3.0000000E+00 kG   (i.e.,   3.0000000E+00 * SCAL)

               Entrance/exit field models are sharp edge
               FINTE, FINTS, gap :    0.0000E+00   0.0000E+00   7.5000E+00

                    Integration step :   5.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.010     0.000     0.000     0.000           40.000     0.015     0.000     0.000     0.000            2
  A    1  1.0000    -0.010     0.000     0.000     0.000           40.000    -0.015     0.000     0.000     0.000            3
  A    1  1.0000     0.000     0.010     0.000     0.000           40.000     0.004     0.000     0.000     0.000            4
  A    1  1.0000     0.000    -0.010     0.000     0.000           40.000    -0.004     0.000     0.000     0.000            5
  A    1  1.0000     0.000     0.000     0.010     0.000           40.000     0.000     0.000     0.005     0.000            6
  A    1  1.0000     0.000     0.000    -0.010     0.000           40.000     0.000     0.000    -0.005     0.000            7
  A    1  1.0000     0.000     0.000     0.000     0.010           40.000     0.000     0.000     0.002     0.000            8
  A    1  1.0000     0.000     0.000     0.000    -0.010           40.000     0.000     0.000    -0.002     0.000            9
  A    1  1.0001     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           10
  A    1  0.9999     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           11

 Cumulative length of optical axis =    3.10000000     m ;  Time  (for ref. rigidity & particle) =   1.656117E-08 s 

************************************************************************************************************************************
      6  Keyword, label(s) :  ESL                               


                              Drift,  length =    30.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00   3.4000000E+02   2.66851E-03

 Cumulative length of optical axis =    3.40000000     m  at  ;  Time  (for ref. rigidity & particle) =   1.816387E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  QUADRUPO    7                     


      -----  QUADRUPOLE  : 
                Length  of  element  =    40.000000      cm
                Bore  radius      RO =    15.000      cm
               B-QUADRUPOLE  = -8.0000000E+00 kG   (i.e.,  -8.0000000E+00 * SCAL)

               Entrance/exit field models are sharp edge
               FINTE, FINTS, gap :    0.0000E+00   0.0000E+00   7.5000E+00

                    Integration step :   5.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.010     0.000     0.000     0.000           40.000     0.019     0.000     0.000     0.000            2
  A    1  1.0000    -0.010     0.000     0.000     0.000           40.000    -0.019     0.000     0.000     0.000            3
  A    1  1.0000     0.000     0.010     0.000     0.000           40.000     0.006     0.000     0.000     0.000            4
  A    1  1.0000     0.000    -0.010     0.000     0.000           40.000    -0.006     0.000     0.000     0.000            5
  A    1  1.0000     0.000     0.000     0.010     0.000           40.000     0.000     0.000     0.001     0.000            6
  A    1  1.0000     0.000     0.000    -0.010     0.000           40.000     0.000     0.000    -0.001     0.000            7
  A    1  1.0000     0.000     0.000     0.000     0.010           40.000     0.000     0.000     0.002     0.000            8
  A    1  1.0000     0.000     0.000     0.000    -0.010           40.000     0.000     0.000    -0.002     0.000            9
  A    1  1.0001     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           10
  A    1  0.9999     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           11

 Cumulative length of optical axis =    3.80000000     m ;  Time  (for ref. rigidity & particle) =   2.030079E-08 s 

************************************************************************************************************************************
      8  Keyword, label(s) :  ESL                               


                              Drift,  length =   200.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00   5.8000000E+02   4.00277E-03

 Cumulative length of optical axis =    5.80000000     m  at  ;  Time  (for ref. rigidity & particle) =   3.098542E-08 s 

************************************************************************************************************************************
      9  Keyword, label(s) :  MATRIX                            


  Reference, before change of frame (part #     1)  : 
   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   5.80000000E+02   4.00276914E-03

           Frame for MATRIX calculation moved by :
            XC =    0.000 cm , YC =    0.000 cm ,   A =  0.00000 deg  ( = 0.000000 rad )


  Reference, after change of frame (part #     1)  : 
   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   5.80000000E+02   4.00276914E-03

  Reference particle (#     1), path length :   580.00000     cm  relative momentum :    1.00000    


                  TRANSFER  MATRIX  ORDRE  1  (MKSA units)

           5.25705         16.9606         0.00000         0.00000         0.00000         0.00000    
           1.66150         5.55066         0.00000         0.00000         0.00000         0.00000    
           0.00000         0.00000        -1.15003        -1.03949         0.00000         0.00000    
           0.00000         0.00000       -0.643559        -1.45123         0.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         1.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         0.00000         1.00000    

          DetY-1 =      -0.0000160143,    DetZ-1 =      -0.0000161106

          R12=0 at   -3.056     m,        R34=0 at  -0.7163     m

      First order symplectic conditions (expected values = 0) :
        -1.6014E-05   -1.6111E-05     0.000         0.000         0.000         0.000    

************************************************************************************************************************************
     10  Keyword, label(s) :  FIT2                              

     FIT procedure launched. Method is 2

           variable #            1       IR =            2 ,   ok.
           variable #            1       IP =            1 ,   ok.

          VARIABLE  ELEMENT    2,  PRMTR #  1 :
               COUPLED  WITH  ELEMENT    8,  PRMTR #  1

           variable #            1       XC.=            8 ,   ok.
           variable #            1       .XC=            1 ,   ok.
           variable #            2       IR =            5 ,   ok.
           variable #            2       IP =           12 ,   ok.
           constraint #            1       IR =            8 ,   ok.
           constraint #            2       IR =            8 ,   ok.

                    FIT  variables  in  good  order,  FIT  will proceed. 

                    Final FIT status will be saved in zgoubi.FITVALS.out                                                              


 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
 LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL         MAXIMUM      STEP     NAME       LBL1     LBL2
    2    1      1    100.        200.        261.2965621       300.      2.892E-03 ESL        *          *         
    8    1      1    100.        200.        138.7034379       300.      2.892E-03
    5    2     12   -3.00        3.00        2.475482591       9.00      3.852E-05 QUADRUPO   5          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   1   1   2     8      1.7000000E+01     1.0000E+00    1.6999996E+01     Infinity   ESL        *          *          *   0 : 
   1   3   4     8     -8.0000000E-01     1.0000E+00   -8.0000182E-01     Infinity   ESL        *          *          *   0 : 
 Fit reached penalty value   2.3202E-11

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       2501.730 kG*cm

                                         CALCUL  DES  TRAJECTOIRES

                              OBJET  (5)  FORME  DE     11 POINTS 



                                Y (cm)         T (mrd)       Z (cm)        P (mrd)       S (cm)        dp/p 
               Sampling :          0.10E-01      0.10E-01      0.10E-01      0.10E-01       0.0          0.1000E-03
  Refrnce trajectry #      1 :      0.0           0.0           0.0           0.0           0.0           1.000    

************************************************************************************************************************************
      2  Keyword, label(s) :  ESL                               


                              Drift,  length =   261.29656  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00   2.6129656E+02   0.00000E+00

 Cumulative length of optical axis =    8.41296562     m  at  ;  Time  (for ref. rigidity & particle) =   4.494470E-08 s 

************************************************************************************************************************************
      3  Keyword, label(s) :  QUADRUPO    3                     


      -----  QUADRUPOLE  : 
                Length  of  element  =    40.000000      cm
                Bore  radius      RO =    15.000      cm
               B-QUADRUPOLE  = -6.0000000E+00 kG   (i.e.,  -6.0000000E+00 * SCAL)

               Entrance/exit field models are sharp edge
               FINTE, FINTS, gap :    0.0000E+00   0.0000E+00   7.5000E+00

                    Integration step :   5.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.010     0.000     0.000     0.000           40.000     0.011     0.000     0.000     0.000            2
  A    1  1.0000    -0.010     0.000     0.000     0.000           40.000    -0.011     0.000     0.000     0.000            3
  A    1  1.0000     0.000     0.010     0.000     0.000           40.000     0.003     0.000     0.000     0.000            4
  A    1  1.0000     0.000    -0.010     0.000     0.000           40.000    -0.003     0.000     0.000     0.000            5
  A    1  1.0000     0.000     0.000     0.010     0.000           40.000     0.000     0.000     0.009     0.000            6
  A    1  1.0000     0.000     0.000    -0.010     0.000           40.000     0.000     0.000    -0.009     0.000            7
  A    1  1.0000     0.000     0.000     0.000     0.010           40.000     0.000     0.000     0.003     0.000            8
  A    1  1.0000     0.000     0.000     0.000    -0.010           40.000     0.000     0.000    -0.003     0.000            9
  A    1  1.0001     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           10
  A    1  0.9999     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           11

 Cumulative length of optical axis =    8.81296562     m ;  Time  (for ref. rigidity & particle) =   4.708162E-08 s 

************************************************************************************************************************************
      4  Keyword, label(s) :  ESL                               


                              Drift,  length =    30.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00   3.3129656E+02   1.33426E-03

 Cumulative length of optical axis =    9.11296562     m  at  ;  Time  (for ref. rigidity & particle) =   4.868432E-08 s 

************************************************************************************************************************************
      5  Keyword, label(s) :  QUADRUPO    5                     


      -----  QUADRUPOLE  : 
                Length  of  element  =    40.000000      cm
                Bore  radius      RO =    15.000      cm
               B-QUADRUPOLE  =  2.4754826E+00 kG   (i.e.,   2.4754826E+00 * SCAL)

               Entrance/exit field models are sharp edge
               FINTE, FINTS, gap :    0.0000E+00   0.0000E+00   7.5000E+00

                    Integration step :   5.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.010     0.000     0.000     0.000           40.000     0.015     0.000     0.000     0.000            2
  A    1  1.0000    -0.010     0.000     0.000     0.000           40.000    -0.015     0.000     0.000     0.000            3
  A    1  1.0000     0.000     0.010     0.000     0.000           40.000     0.005     0.000     0.000     0.000            4
  A    1  1.0000     0.000    -0.010     0.000     0.000           40.000    -0.005     0.000     0.000     0.000            5
  A    1  1.0000     0.000     0.000     0.010     0.000           40.000     0.000     0.000     0.005     0.000            6
  A    1  1.0000     0.000     0.000    -0.010     0.000           40.000     0.000     0.000    -0.005     0.000            7
  A    1  1.0000     0.000     0.000     0.000     0.010           40.000     0.000     0.000     0.002     0.000            8
  A    1  1.0000     0.000     0.000     0.000    -0.010           40.000     0.000     0.000    -0.002     0.000            9
  A    1  1.0001     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           10
  A    1  0.9999     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           11

 Cumulative length of optical axis =    9.51296562     m ;  Time  (for ref. rigidity & particle) =   5.082124E-08 s 

************************************************************************************************************************************
      6  Keyword, label(s) :  ESL                               


                              Drift,  length =    30.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00   4.0129656E+02   2.66851E-03

 Cumulative length of optical axis =    9.81296562     m  at  ;  Time  (for ref. rigidity & particle) =   5.242394E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  QUADRUPO    7                     


      -----  QUADRUPOLE  : 
                Length  of  element  =    40.000000      cm
                Bore  radius      RO =    15.000      cm
               B-QUADRUPOLE  = -8.0000000E+00 kG   (i.e.,  -8.0000000E+00 * SCAL)

               Entrance/exit field models are sharp edge
               FINTE, FINTS, gap :    0.0000E+00   0.0000E+00   7.5000E+00

                    Integration step :   5.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.010     0.000     0.000     0.000           40.000     0.020     0.000     0.000     0.000            2
  A    1  1.0000    -0.010     0.000     0.000     0.000           40.000    -0.020     0.000     0.000     0.000            3
  A    1  1.0000     0.000     0.010     0.000     0.000           40.000     0.007     0.000     0.000     0.000            4
  A    1  1.0000     0.000    -0.010     0.000     0.000           40.000    -0.007     0.000     0.000     0.000            5
  A    1  1.0000     0.000     0.000     0.010     0.000           40.000     0.000     0.000     0.001     0.000            6
  A    1  1.0000     0.000     0.000    -0.010     0.000           40.000     0.000     0.000    -0.001     0.000            7
  A    1  1.0000     0.000     0.000     0.000     0.010           40.000     0.000     0.000     0.002     0.000            8
  A    1  1.0000     0.000     0.000     0.000    -0.010           40.000     0.000     0.000    -0.002     0.000            9
  A    1  1.0001     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           10
  A    1  0.9999     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           11

 Cumulative length of optical axis =    10.2129656     m ;  Time  (for ref. rigidity & particle) =   5.456086E-08 s 

************************************************************************************************************************************
      8  Keyword, label(s) :  ESL                               


                              Drift,  length =   138.70344  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00   5.8000000E+02   4.00277E-03

 Cumulative length of optical axis =    11.6000000     m  at  ;  Time  (for ref. rigidity & particle) =   6.197084E-08 s 

************************************************************************************************************************************
      9  Keyword, label(s) :  MATRIX                            


  Reference, before change of frame (part #     1)  : 
   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   5.80000000E+02   4.00276914E-03

           Frame for MATRIX calculation moved by :
            XC =    0.000 cm , YC =    0.000 cm ,   A =  0.00000 deg  ( = 0.000000 rad )


  Reference, after change of frame (part #     1)  : 
   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   5.80000000E+02   4.00276914E-03

  Reference particle (#     1), path length :   580.00000     cm  relative momentum :    1.00000    


                  TRANSFER  MATRIX  ORDRE  1  (MKSA units)

           4.49573         17.0000         0.00000         0.00000         0.00000         0.00000    
           1.78958         6.98947         0.00000         0.00000         0.00000         0.00000    
           0.00000         0.00000       -0.801681       -0.800002         0.00000         0.00000    
           0.00000         0.00000       -0.657618        -1.90360         0.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         1.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         0.00000         1.00000    

          DetY-1 =      -0.0000155969,    DetZ-1 =      -0.0000156851

          R12=0 at   -2.432     m,        R34=0 at  -0.4203     m

      First order symplectic conditions (expected values = 0) :
        -1.5597E-05   -1.5685E-05     0.000         0.000         0.000         0.000    

************************************************************************************************************************************

     10   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
     11  Keyword, label(s) :  MATRIX                            


  Reference, before change of frame (part #     1)  : 
   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   5.80000000E+02   4.00276914E-03

           Frame for MATRIX calculation moved by :
            XC =    0.000 cm , YC =    0.000 cm ,   A =  0.00000 deg  ( = 0.000000 rad )


  Reference, after change of frame (part #     1)  : 
   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   5.80000000E+02   4.00276914E-03

  Reference particle (#     1), path length :   580.00000     cm  relative momentum :    1.00000    


                  TRANSFER  MATRIX  ORDRE  1  (MKSA units)

           4.49573         17.0000         0.00000         0.00000         0.00000         0.00000    
           1.78958         6.98947         0.00000         0.00000         0.00000         0.00000    
           0.00000         0.00000       -0.801681       -0.800002         0.00000         0.00000    
           0.00000         0.00000       -0.657618        -1.90360         0.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         1.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         0.00000         1.00000    

          DetY-1 =      -0.0000155969,    DetZ-1 =      -0.0000156851

          R12=0 at   -2.432     m,        R34=0 at  -0.4203     m

      First order symplectic conditions (expected values = 0) :
        -1.5597E-05   -1.5685E-05     0.000         0.000         0.000         0.000    

************************************************************************************************************************************
     12  Keyword, label(s) :  END                               


                            11 particles have been launched
                     Made  it  to  the  end :     11

************************************************************************************************************************************

           MAIN PROGRAM : Execution ended upon key  END       

************************************************************************************************************************************

                    Saved new version of zgoubi.dat with variables updated.

                    Updated version of zgoubi.dat saved in  zgoubi.FIT.out.dat
   
            ZGOUBI RUN COMPLETED. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  08-02-0015,  at  09:00:49 
  JOB  ENDED  ON    08-02-0015,  AT  09:00:49 

   CPU time, total :    4.000200000000000E-002
