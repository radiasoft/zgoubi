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
 'FIT2'            ! Vary b in quads for FIT of R12 and R34                                                  10
2  save            ! # of variables; variables will save in zgoubi.FITVALS.out                                
3  12 7.12  2.     ! Symmetric triplet => quads #1 and #3 are coupled                                         
5  12  0.   2.     ! Parmtr #12 of elements #3, 5 and 7 is the field value                                    
2  1.E-10          ! # of constraints; penalty                                                                
1  1 2 8 16.6 1.   ! Cnstrnt #1 : R12=16.6 after last drift (lmnt #8)                                         
1  3 4 8 -.88 1.   ! Cnstrnt #2 : R34=-.88 after last drift                                                   
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
     10  Keyword, label(s) :  FIT2                              

     FIT procedure launched. Method is 2

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
    3    1     12   -21.0       -7.00       -6.972765137       7.00      1.707E-04 QUADRUPO   3          *         
    7    1    120   -6.97       -7.00       -6.972765137       7.00      1.707E-04
    5    2     12   -3.00        3.00        3.229344585       9.00      1.266E-04 QUADRUPO   5          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   1   1   2     8      1.6600000E+01     1.0000E+00    1.6599981E+01     Infinity   ESL        *          *          *   0 : 
   1   3   4     8     -8.8000000E-01     1.0000E+00   -8.8000964E-01     Infinity   ESL        *          *          *   0 : 
 Fit reached penalty value   8.4374E-11

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


                              Drift,  length =   200.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00   2.0000000E+02   0.00000E+00

 Cumulative length of optical axis =    7.80000000     m  at  ;  Time  (for ref. rigidity & particle) =   4.167004E-08 s 

************************************************************************************************************************************
      3  Keyword, label(s) :  QUADRUPO    3                     


      -----  QUADRUPOLE  : 
                Length  of  element  =    40.000000      cm
                Bore  radius      RO =    15.000      cm
               B-QUADRUPOLE  = -6.9727651E+00 kG 

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

 Cumulative length of optical axis =    8.20000000     m ;  Time  (for ref. rigidity & particle) =   4.380697E-08 s 

************************************************************************************************************************************
      4  Keyword, label(s) :  ESL                               


                              Drift,  length =    30.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00   2.7000000E+02   1.33426E-03

 Cumulative length of optical axis =    8.50000000     m  at  ;  Time  (for ref. rigidity & particle) =   4.540966E-08 s 

************************************************************************************************************************************
      5  Keyword, label(s) :  QUADRUPO    5                     


      -----  QUADRUPOLE  : 
                Length  of  element  =    40.000000      cm
                Bore  radius      RO =    15.000      cm
               B-QUADRUPOLE  =  3.2293446E+00 kG 

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

 Cumulative length of optical axis =    8.90000000     m ;  Time  (for ref. rigidity & particle) =   4.754659E-08 s 

************************************************************************************************************************************
      6  Keyword, label(s) :  ESL                               


                              Drift,  length =    30.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00   3.4000000E+02   2.66851E-03

 Cumulative length of optical axis =    9.20000000     m  at  ;  Time  (for ref. rigidity & particle) =   4.914928E-08 s 

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

 Cumulative length of optical axis =    9.60000000     m ;  Time  (for ref. rigidity & particle) =   5.128621E-08 s 

************************************************************************************************************************************
      8  Keyword, label(s) :  ESL                               


                              Drift,  length =   200.00000  cm

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

           5.27056         16.6000         0.00000         0.00000         0.00000         0.00000    
           1.61443         5.27450         0.00000         0.00000         0.00000         0.00000    
           0.00000         0.00000        -1.24205       -0.880010         0.00000         0.00000    
           0.00000         0.00000       -0.622553        -1.24620         0.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         1.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         0.00000         1.00000    

          DetY-1 =      -0.0000158889,    DetZ-1 =      -0.0000159710

          R12=0 at   -3.147     m,        R34=0 at  -0.7062     m

      First order symplectic conditions (expected values = 0) :
        -1.5889E-05   -1.5971E-05     0.000         0.000         0.000         0.000    

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

           5.27056         16.6000         0.00000         0.00000         0.00000         0.00000    
           1.61443         5.27450         0.00000         0.00000         0.00000         0.00000    
           0.00000         0.00000        -1.24205       -0.880010         0.00000         0.00000    
           0.00000         0.00000       -0.622553        -1.24620         0.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         1.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         0.00000         1.00000    

          DetY-1 =      -0.0000158889,    DetZ-1 =      -0.0000159710

          R12=0 at   -3.147     m,        R34=0 at  -0.7062     m

      First order symplectic conditions (expected values = 0) :
        -1.5889E-05   -1.5971E-05     0.000         0.000         0.000         0.000    

************************************************************************************************************************************
     12  Keyword, label(s) :  END                               


                            11 particles have been launched
                     Made  it  to  the  end :     11

************************************************************************************************************************************

           MAIN PROGRAM : Execution ended upon key  END       

************************************************************************************************************************************
   
            Zgoubi run completed. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  28-11-0014,  at  16:21:04 
  Job  ended  on    28-11-0014,  at  16:21:04 

   CPU time, total :    3.600100000000000E-002
