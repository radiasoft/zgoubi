MATCHING A SYMMETRIC QUADRUPOLE TRIPLET, INCLUDING COUPLED QUADS                                              
 'OBJET'                                                                                                      1
2501.73            RIGIDITY (kG.cm)                                                                           
5                  11 PARTICLES GENERATED FOR USE OF MATRIX                                                   
.01  .01  .01   .01  0. .0001                                                                                 
0.  0.  0.   0.  0.  1.                                                                                       
 'ESL   '                                                                                                     2
200.                                                                                                          
 'QUADRUPO'    3                                                                                              3
0                                                                                                             
40.  15.  -7.                                                                                                 
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
40.  15.  -7.                                                                                                 
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
 'FIT'             ! Vary b in quads for FIT of R12 and R34                                                  10
2  save  nofinal   ! # of variables; variables will save in zgoubi.FITVALS.out                                
3  12 7.12  2.     ! Symmetric triplet => quads #1 and #3 are coupled                                         
5  12  0.   2.     ! Parmtr #12 of elements #3, 5 and 7 is the field value                                    
2  1.E-10          ! # of constraints; penalty                                                                
1  1 2 8 16.6 1.   ! Cnstrnt #1 : R12=16.6 after last drift (lmnt #8)                                         
1  3 4 8 -.88 1.   ! Cnstrnt #2 : R34=-.88 after last drift                                                   
 'MATRIX'                                                                                                    11
1  0                                                                                                          
                                                                                                              
 'REBELOTE'                                                                                                  12
3 0.1 0                                                                                                       
                                                                                                              
 'FAISCEAU'                                                                                                  13
 'END'                                                                                                       14

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
               B-QUADRUPOLE  = -7.0000000E+00 kG 

               Entrance/exit field models are sharp edge
               FINTE, FINTS, gap :    0.0000E+00   0.0000E+00   7.5000E+00

                    Integration step :   5.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.010     0.000     0.000     0.000           40.000     0.012     0.000     0.000     0.000            2
  A    1  1.0000    -0.010     0.000     0.000     0.000           40.000    -0.012     0.000     0.000     0.000            3
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
               B-QUADRUPOLE  =  3.0000000E+00 kG 

               Entrance/exit field models are sharp edge
               FINTE, FINTS, gap :    0.0000E+00   0.0000E+00   7.5000E+00

                    Integration step :   5.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.010     0.000     0.000     0.000           40.000     0.016     0.000     0.000     0.000            2
  A    1  1.0000    -0.010     0.000     0.000     0.000           40.000    -0.016     0.000     0.000     0.000            3
  A    1  1.0000     0.000     0.010     0.000     0.000           40.000     0.004     0.000     0.000     0.000            4
  A    1  1.0000     0.000    -0.010     0.000     0.000           40.000    -0.004     0.000     0.000     0.000            5
  A    1  1.0000     0.000     0.000     0.010     0.000           40.000     0.000     0.000     0.004     0.000            6
  A    1  1.0000     0.000     0.000    -0.010     0.000           40.000     0.000     0.000    -0.004     0.000            7
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
               B-QUADRUPOLE  = -7.0000000E+00 kG 

               Entrance/exit field models are sharp edge
               FINTE, FINTS, gap :    0.0000E+00   0.0000E+00   7.5000E+00

                    Integration step :   5.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.010     0.000     0.000     0.000           40.000     0.021     0.000     0.000     0.000            2
  A    1  1.0000    -0.010     0.000     0.000     0.000           40.000    -0.021     0.000     0.000     0.000            3
  A    1  1.0000     0.000     0.010     0.000     0.000           40.000     0.006     0.000     0.000     0.000            4
  A    1  1.0000     0.000    -0.010     0.000     0.000           40.000    -0.006     0.000     0.000     0.000            5
  A    1  1.0000     0.000     0.000     0.010     0.000           40.000     0.000     0.000     0.000     0.000            6
  A    1  1.0000     0.000     0.000    -0.010     0.000           40.000     0.000     0.000     0.000     0.000            7
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

           5.43427         17.0254         0.00000         0.00000         0.00000         0.00000    
           1.67580         5.43425         0.00000         0.00000         0.00000         0.00000    
           0.00000         0.00000        -1.27003       -0.974288         0.00000         0.00000    
           0.00000         0.00000       -0.629171        -1.27004         0.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         1.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         0.00000         1.00000    

          DetY-1 =      -0.0000157172,    DetZ-1 =      -0.0000158156

          R12=0 at   -3.133     m,        R34=0 at  -0.7671     m

      First order symplectic conditions (expected values = 0) :
        -1.5717E-05   -1.5816E-05     0.000         0.000         0.000         0.000    

************************************************************************************************************************************
     10  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1

           variable #            1       IR =            3 ,   ok.
           variable #            1       IP =           12 ,   ok.

          VARIABLE  ELEMENT    3,  PRMTR # 12 :
               COUPLED  WITH  ELEMENT    7,  PRMTR #120

           variable #            1       XC.=            7 ,   ok.
           variable #            1       .XC=          120 ,   ok.
           variable #            2       IR =            5 ,   ok.
           variable #            2       IP =           12 ,   ok.
           constraint #            1       IR =            8 ,   ok.
           constraint #            2       IR =            8 ,   ok.

                    FIT  variables  in  good  order,  FIT  will proceed. 

                    Final FIT status will be saved in zgoubi.FITVALS.out                                                              


 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
 LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL         MAXIMUM      STEP     NAME       LBL1     LBL2
    3    1     12   -21.0       -6.85       -6.847901235       7.00      3.841E-04 QUADRUPO   3          *         
    7    1    120   -6.85       -6.85       -6.847901235       7.00      3.841E-04
    5    2     12   -3.00        3.13        3.130864198       9.00      1.646E-04 QUADRUPO   5          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   1   1   2     8      1.6600000E+01     1.0000E+00    1.6605335E+01   1.4265E-01   ESL        *          *          *   0 : 
   1   3   4     8     -8.8000000E-01     1.0000E+00   -8.6692125E-01   8.5735E-01   ESL        *          *          *   0 : 
 Fit reached penalty value   1.9952E-04

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

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

           5.25981         16.6053         0.00000         0.00000         0.00000         0.00000    
           1.61287         5.28198         0.00000         0.00000         0.00000         0.00000    
           0.00000         0.00000        -1.22982       -0.866921         0.00000         0.00000    
           0.00000         0.00000       -0.623829        -1.25286         0.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         1.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         0.00000         1.00000    

          DetY-1 =      -0.0000155409,    DetZ-1 =      -0.0000156233

          R12=0 at   -3.144     m,        R34=0 at  -0.6920     m

      First order symplectic conditions (expected values = 0) :
        -1.5541E-05   -1.5623E-05     0.000         0.000         0.000         0.000    

************************************************************************************************************************************
     12  Keyword, label(s) :  REBELOTE                          


                                -----  REBELOTE  -----

     End of pass #        1 through the optical structure 

                     Total of         11 particles have been launched

     Multiple pass, 
          from element #     1 : OBJET     /label1=          /label2=           to REBELOTE /label1=          /label2=          
          ending at pass #       4 at element #    12 : REBELOTE  /label1=          /label2=          


     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
 LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL         MAXIMUM      STEP     NAME       LBL1     LBL2
    3    1     12   -21.0       -6.92       -6.919725652       7.00      1.280E-04 QUADRUPO   3          *         
    7    1    120   -6.92       -6.92       -6.919725652       7.00      1.280E-04
    5    2     12   -3.00        3.19        3.188148148       9.00      5.487E-05 QUADRUPO   5          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   1   1   2     8      1.6600000E+01     1.0000E+00    1.6601341E+01   4.9776E-02   ESL        *          *          *   0 : 
   1   3   4     8     -8.8000000E-01     1.0000E+00   -8.7414303E-01   9.5022E-01   ESL        *          *          *   0 : 
 Fit reached penalty value   3.6101E-05

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

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

           5.26566         16.6013         0.00000         0.00000         0.00000         0.00000    
           1.61364         5.27733         0.00000         0.00000         0.00000         0.00000    
           0.00000         0.00000        -1.23675       -0.874143         0.00000         0.00000    
           0.00000         0.00000       -0.623067        -1.24894         0.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         1.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         0.00000         1.00000    

          DetY-1 =      -0.0000157408,    DetZ-1 =      -0.0000158229

          R12=0 at   -3.146     m,        R34=0 at  -0.6999     m

      First order symplectic conditions (expected values = 0) :
        -1.5741E-05   -1.5823E-05     0.000         0.000         0.000         0.000    

************************************************************************************************************************************
     12  Keyword, label(s) :  REBELOTE                          


                                -----  REBELOTE  -----

     End of pass #        2 through the optical structure 

                     Total of         22 particles have been launched

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
 LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL         MAXIMUM      STEP     NAME       LBL1     LBL2
    3    1     12   -21.0       -6.95       -6.950452675       7.00      4.268E-05 QUADRUPO   3          *         
    7    1    120   -6.95       -6.95       -6.950452675       7.00      4.268E-05
    5    2     12   -3.00        3.21        3.212016461       9.00      1.829E-05 QUADRUPO   5          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   1   1   2     8      1.6600000E+01     1.0000E+00    1.6600578E+01   5.1780E-02   ESL        *          *          *   0 : 
   1   3   4     8     -8.8000000E-01     1.0000E+00   -8.7752784E-01   9.4822E-01   ESL        *          *          *   0 : 
 Fit reached penalty value   6.4453E-06

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

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

           5.26851         16.6006         0.00000         0.00000         0.00000         0.00000    
           1.61411         5.27570         0.00000         0.00000         0.00000         0.00000    
           0.00000         0.00000        -1.23982       -0.877528         0.00000         0.00000    
           0.00000         0.00000       -0.622767        -1.24735         0.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         1.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         0.00000         1.00000    

          DetY-1 =      -0.0000158264,    DetZ-1 =      -0.0000159085

          R12=0 at   -3.147     m,        R34=0 at  -0.7035     m

      First order symplectic conditions (expected values = 0) :
        -1.5826E-05   -1.5909E-05     0.000         0.000         0.000         0.000    

************************************************************************************************************************************
     12  Keyword, label(s) :  REBELOTE                          


                                -----  REBELOTE  -----

     End of pass #        3 through the optical structure 

                     Total of         33 particles have been launched


      Next  pass  is  #     4 and  last  pass  through  the  optical  structure


     WRITE statements to zgoubi.res are re-established from now on.

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

 Cumulative length of optical axis =    2.00000000     m  at  ;  Time  (for ref. rigidity & particle) =   4.167004E-08 s 

************************************************************************************************************************************
      3  Keyword, label(s) :  QUADRUPO    3                     


      -----  QUADRUPOLE  : 
                Length  of  element  =    40.000000      cm
                Bore  radius      RO =    15.000      cm
               B-QUADRUPOLE  = -6.9504527E+00 kG 

               Entrance/exit field models are sharp edge
               FINTE, FINTS, gap :    0.0000E+00   0.0000E+00   7.5000E+00

                    Integration step :   5.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.010     0.000     0.000     0.000           40.000     0.012     0.000     0.000     0.000            2
  A    1  1.0000    -0.010     0.000     0.000     0.000           40.000    -0.012     0.000     0.000     0.000            3
  A    1  1.0000     0.000     0.010     0.000     0.000           40.000     0.003     0.000     0.000     0.000            4
  A    1  1.0000     0.000    -0.010     0.000     0.000           40.000    -0.003     0.000     0.000     0.000            5
  A    1  1.0000     0.000     0.000     0.010     0.000           40.000     0.000     0.000     0.009     0.000            6
  A    1  1.0000     0.000     0.000    -0.010     0.000           40.000     0.000     0.000    -0.009     0.000            7
  A    1  1.0000     0.000     0.000     0.000     0.010           40.000     0.000     0.000     0.002     0.000            8
  A    1  1.0000     0.000     0.000     0.000    -0.010           40.000     0.000     0.000    -0.002     0.000            9
  A    1  1.0001     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           10
  A    1  0.9999     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           11

 Cumulative length of optical axis =    2.40000000     m ;  Time  (for ref. rigidity & particle) =   4.380697E-08 s 

************************************************************************************************************************************
      4  Keyword, label(s) :  ESL                               


                              Drift,  length =    30.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00   2.7000000E+02   1.33426E-03

 Cumulative length of optical axis =    2.70000000     m  at  ;  Time  (for ref. rigidity & particle) =   4.540966E-08 s 

************************************************************************************************************************************
      5  Keyword, label(s) :  QUADRUPO    5                     


      -----  QUADRUPOLE  : 
                Length  of  element  =    40.000000      cm
                Bore  radius      RO =    15.000      cm
               B-QUADRUPOLE  =  3.2120165E+00 kG 

               Entrance/exit field models are sharp edge
               FINTE, FINTS, gap :    0.0000E+00   0.0000E+00   7.5000E+00

                    Integration step :   5.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.010     0.000     0.000     0.000           40.000     0.016     0.000     0.000     0.000            2
  A    1  1.0000    -0.010     0.000     0.000     0.000           40.000    -0.016     0.000     0.000     0.000            3
  A    1  1.0000     0.000     0.010     0.000     0.000           40.000     0.004     0.000     0.000     0.000            4
  A    1  1.0000     0.000    -0.010     0.000     0.000           40.000    -0.004     0.000     0.000     0.000            5
  A    1  1.0000     0.000     0.000     0.010     0.000           40.000     0.000     0.000     0.004     0.000            6
  A    1  1.0000     0.000     0.000    -0.010     0.000           40.000     0.000     0.000    -0.004     0.000            7
  A    1  1.0000     0.000     0.000     0.000     0.010           40.000     0.000     0.000     0.002     0.000            8
  A    1  1.0000     0.000     0.000     0.000    -0.010           40.000     0.000     0.000    -0.002     0.000            9
  A    1  1.0001     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           10
  A    1  0.9999     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           11

 Cumulative length of optical axis =    3.10000000     m ;  Time  (for ref. rigidity & particle) =   4.754659E-08 s 

************************************************************************************************************************************
      6  Keyword, label(s) :  ESL                               


                              Drift,  length =    30.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00   3.4000000E+02   2.66851E-03

 Cumulative length of optical axis =    3.40000000     m  at  ;  Time  (for ref. rigidity & particle) =   4.914928E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  QUADRUPO    7                     


      -----  QUADRUPOLE  : 
                Length  of  element  =    40.000000      cm
                Bore  radius      RO =    15.000      cm
               B-QUADRUPOLE  = -7.0000000E+00 kG 

               Entrance/exit field models are sharp edge
               FINTE, FINTS, gap :    0.0000E+00   0.0000E+00   7.5000E+00

                    Integration step :   5.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.010     0.000     0.000     0.000           40.000     0.020     0.000     0.000     0.000            2
  A    1  1.0000    -0.010     0.000     0.000     0.000           40.000    -0.020     0.000     0.000     0.000            3
  A    1  1.0000     0.000     0.010     0.000     0.000           40.000     0.006     0.000     0.000     0.000            4
  A    1  1.0000     0.000    -0.010     0.000     0.000           40.000    -0.006     0.000     0.000     0.000            5
  A    1  1.0000     0.000     0.000     0.010     0.000           40.000     0.000     0.000     0.000     0.000            6
  A    1  1.0000     0.000     0.000    -0.010     0.000           40.000     0.000     0.000     0.000     0.000            7
  A    1  1.0000     0.000     0.000     0.000     0.010           40.000     0.000     0.000     0.002     0.000            8
  A    1  1.0000     0.000     0.000     0.000    -0.010           40.000     0.000     0.000    -0.002     0.000            9
  A    1  1.0001     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           10
  A    1  0.9999     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000           11

 Cumulative length of optical axis =    3.80000000     m ;  Time  (for ref. rigidity & particle) =   5.128621E-08 s 

************************************************************************************************************************************
      8  Keyword, label(s) :  ESL                               


                              Drift,  length =   200.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00   5.8000000E+02   4.00277E-03

 Cumulative length of optical axis =    5.80000000     m  at  ;  Time  (for ref. rigidity & particle) =   6.197084E-08 s 

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

           5.26851         16.6006         0.00000         0.00000         0.00000         0.00000    
           1.61411         5.27570         0.00000         0.00000         0.00000         0.00000    
           0.00000         0.00000        -1.23982       -0.877528         0.00000         0.00000    
           0.00000         0.00000       -0.622767        -1.24735         0.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         1.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         0.00000         1.00000    

          DetY-1 =      -0.0000158264,    DetZ-1 =      -0.0000159085

          R12=0 at   -3.147     m,        R34=0 at  -0.7035     m

      First order symplectic conditions (expected values = 0) :
        -1.5826E-05   -1.5909E-05     0.000         0.000         0.000         0.000    

************************************************************************************************************************************
     10  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1


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
               B-QUADRUPOLE  = -6.9504527E+00 kG 

               Entrance/exit field models are sharp edge
               FINTE, FINTS, gap :    0.0000E+00   0.0000E+00   7.5000E+00

                    Integration step :   5.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.010     0.000     0.000     0.000           40.000     0.012     0.000     0.000     0.000            2
  A    1  1.0000    -0.010     0.000     0.000     0.000           40.000    -0.012     0.000     0.000     0.000            3
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
               B-QUADRUPOLE  =  3.2120165E+00 kG 

               Entrance/exit field models are sharp edge
               FINTE, FINTS, gap :    0.0000E+00   0.0000E+00   7.5000E+00

                    Integration step :   5.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.010     0.000     0.000     0.000           40.000     0.016     0.000     0.000     0.000            2
  A    1  1.0000    -0.010     0.000     0.000     0.000           40.000    -0.016     0.000     0.000     0.000            3
  A    1  1.0000     0.000     0.010     0.000     0.000           40.000     0.004     0.000     0.000     0.000            4
  A    1  1.0000     0.000    -0.010     0.000     0.000           40.000    -0.004     0.000     0.000     0.000            5
  A    1  1.0000     0.000     0.000     0.010     0.000           40.000     0.000     0.000     0.004     0.000            6
  A    1  1.0000     0.000     0.000    -0.010     0.000           40.000     0.000     0.000    -0.004     0.000            7
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
               B-QUADRUPOLE  = -7.0000000E+00 kG 

               Entrance/exit field models are sharp edge
               FINTE, FINTS, gap :    0.0000E+00   0.0000E+00   7.5000E+00

                    Integration step :   5.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000           40.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.010     0.000     0.000     0.000           40.000     0.020     0.000     0.000     0.000            2
  A    1  1.0000    -0.010     0.000     0.000     0.000           40.000    -0.020     0.000     0.000     0.000            3
  A    1  1.0000     0.000     0.010     0.000     0.000           40.000     0.006     0.000     0.000     0.000            4
  A    1  1.0000     0.000    -0.010     0.000     0.000           40.000    -0.006     0.000     0.000     0.000            5
  A    1  1.0000     0.000     0.000     0.010     0.000           40.000     0.000     0.000     0.000     0.000            6
  A    1  1.0000     0.000     0.000    -0.010     0.000           40.000     0.000     0.000     0.000     0.000            7
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

           5.26851         16.6006         0.00000         0.00000         0.00000         0.00000    
           1.61411         5.27570         0.00000         0.00000         0.00000         0.00000    
           0.00000         0.00000        -1.23982       -0.877528         0.00000         0.00000    
           0.00000         0.00000       -0.622767        -1.24735         0.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         1.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         0.00000         1.00000    

          DetY-1 =      -0.0000158264,    DetZ-1 =      -0.0000159085

          R12=0 at   -3.147     m,        R34=0 at  -0.7035     m

      First order symplectic conditions (expected values = 0) :
        -1.5826E-05   -1.5909E-05     0.000         0.000         0.000         0.000    

************************************************************************************************************************************
     10  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1

           variable #            1       IR =            3 ,   ok.
           variable #            1       IP =           12 ,   ok.

          VARIABLE  ELEMENT    3,  PRMTR # 12 :
               COUPLED  WITH  ELEMENT    7,  PRMTR #120

           variable #            1       XC.=            7 ,   ok.
           variable #            1       .XC=          120 ,   ok.
           variable #            2       IR =            5 ,   ok.
           variable #            2       IP =           12 ,   ok.
           constraint #            1       IR =            8 ,   ok.
           constraint #            2       IR =            8 ,   ok.

                    FIT  variables  in  good  order,  FIT  will proceed. 

                    Final FIT status will be saved in zgoubi.FITVALS.out                                                              


 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
 LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL         MAXIMUM      STEP     NAME       LBL1     LBL2
    3    1     12   -21.0       -6.96       -6.963255601       7.00      1.423E-05 QUADRUPO   3          *         
    7    1    120   -6.96       -6.96       -6.963255601       7.00      1.423E-05
    5    2     12   -3.00        3.22        3.221947874       9.00      6.097E-06 QUADRUPO   5          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   1   1   2     8      1.6600000E+01     1.0000E+00    1.6600260E+01   5.8036E-02   ESL        *          *          *   0 : 
   1   3   4     8     -8.8000000E-01     1.0000E+00   -8.7895379E-01   9.4196E-01   ESL        *          *          *   0 : 
 Fit reached penalty value   1.1620E-06

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

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

           5.26969         16.6003         0.00000         0.00000         0.00000         0.00000    
           1.61430         5.27502         0.00000         0.00000         0.00000         0.00000    
           0.00000         0.00000        -1.24110       -0.878954         0.00000         0.00000    
           0.00000         0.00000       -0.622644        -1.24669         0.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         1.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         0.00000         1.00000    

          DetY-1 =      -0.0000158622,    DetZ-1 =      -0.0000159443

          R12=0 at   -3.147     m,        R34=0 at  -0.7050     m

      First order symplectic conditions (expected values = 0) :
        -1.5862E-05   -1.5944E-05     0.000         0.000         0.000         0.000    

************************************************************************************************************************************
     12  Keyword, label(s) :  REBELOTE                          


                         ****  End  of  'REBELOTE'  procedure  ****

      There  has  been          4  passes  through  the  optical  structure 

                     Total of         44 particles have been launched

************************************************************************************************************************************
     13  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                 11 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

O  1   1.0000     0.000     0.000     0.000     0.000       0.0000    1.0000    0.000    0.000    0.000    0.000   5.800000E+02     1
A  1   1.0000     0.010     0.000     0.000     0.000       0.0000    1.0000    0.053    0.161    0.000    0.000   5.800000E+02     2
B  1   1.0000    -0.010     0.000     0.000     0.000       0.0000    1.0000   -0.053   -0.161    0.000    0.000   5.800000E+02     3
C  1   1.0000     0.000     0.010     0.000     0.000       0.0000    1.0000    0.017    0.053    0.000    0.000   5.800000E+02     4
D  1   1.0000     0.000    -0.010     0.000     0.000       0.0000    1.0000   -0.017   -0.053    0.000    0.000   5.800000E+02     5
E  1   1.0000     0.000     0.000     0.010     0.000       0.0000    1.0000    0.000    0.000   -0.012   -0.062   5.800000E+02     6
F  1   1.0000     0.000     0.000    -0.010     0.000       0.0000    1.0000    0.000    0.000    0.012    0.062   5.800000E+02     7
G  1   1.0000     0.000     0.000     0.000     0.010       0.0000    1.0000    0.000    0.000   -0.001   -0.012   5.800000E+02     8
H  1   1.0000     0.000     0.000     0.000    -0.010       0.0000    1.0000    0.000    0.000    0.001    0.012   5.800000E+02     9
I  1   1.0001     0.000     0.000     0.000     0.000       0.0000    1.0001    0.000    0.000    0.000    0.000   5.800000E+02    10
J  1   0.9999     0.000     0.000     0.000     0.000       0.0000    0.9999    0.000    0.000    0.000    0.000   5.800000E+02    11


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

   2.2848E-09   -93.83        305.3        0.000        0.000          11      11    1.000    B-Dim     1      5
   2.2848E-09   -7.837        15.48        0.000        0.000          11      11    1.000    B-Dim     2      5
   1.5462E-12    0.000       1.2031E-10   4.0028E-03    750.0          11      11    1.000    B-Dim     3      5


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  1.818153E-10  1.818153E-10
           AlpX, BetX =  9.382660E+01  3.052583E+02
           AlpZ, BetZ =  7.837328E+00  1.548074E+01

  5.550063E-08  1.705911E-08  0.000000E+00  0.000000E+00

  1.705911E-08  5.244019E-09  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  2.814635E-09  1.424946E-09

  0.000000E+00  0.000000E+00  1.424946E-09  7.331424E-10

************************************************************************************************************************************
     14  Keyword, label(s) :  END                               


************************************************************************************************************************************

           MAIN PROGRAM : Execution ended upon key  END       

************************************************************************************************************************************
   
            Zgoubi run completed. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  29-11-0014,  at  04:47:54 
  Job  ended  on    29-11-0014,  at  04:47:55 

   CPU time, total :    0.664041000000000     
