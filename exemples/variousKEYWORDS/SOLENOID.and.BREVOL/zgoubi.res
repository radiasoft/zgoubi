NuFactory. Solenoid pion capture. In-flight decay into muon
 'OBJET'                                                                                1                     1
833.910238
5
0.1  0.1 0.1 0.1 0. .001 'o'
0.  0. 0. 0. 0. 1.
 'DRIFT'                                                                                2                     2
-20.
 'BREVOL'                                                                                3                    3
0 0
0.095933608187 1
TEST CARTE SOLENO
1521                                  IX
sol30m.map
0  0. 0. 0.                            PAS DE DROITE COUPURE
2
1.
1 0 0 0
 'DRIFT'                                                                                4                     4
-20.
 'FAISCEAU'                                                                              5                    5
 'MATRIX'                                                                                6                    6
1  0
 'RESET'                                                                                 7                    7
 'OBJET'                                                                                8                     8
833.910238
5
0.1  0.1 0.1 0.1 0. .001 'o'
0.  0. 0. 0. 0. 1.
 'SOLENOID'                                                                            9                      9
0
3000.  1.  16.
20.  20.
1.
1  0. 0. 0.
 'FAISCEAU'                                                                             10                   10
 'MATRIX'                                                                               11                   11
1  0
 'END'                                                                                 12                    12

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET       1                     

                          MAGNETIC  RIGIDITY =        833.910 kG*cm

                                         CALCUL  DES  TRAJECTOIRES

                              OBJET  (5)  FORME  DE     11 POINTS 



                                Y (cm)         T (mrd)       Z (cm)        P (mrd)       S (cm)        dp/p 
               Sampling :          0.10          0.10          0.10          0.10           0.0          0.1000E-02
  Refrnce trajectry #      1 :      0.0           0.0           0.0           0.0           0.0           1.000    

************************************************************************************************************************************
      2  Keyword, label(s) :  DRIFT       2                     


                              Drift,  length =   -20.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  -2.0000000E+01   0.00000E+00

 Cumulative length of optical axis =  -0.200000000     m  at  ;  Time  (for ref. rigidity & particle) =  -2.591144E-09 s 

************************************************************************************************************************************
      3  Keyword, label(s) :  BREVOL      3                     


           Title in this field map problem :TEST CARTE SOLENO                                                                                                       
  A total of           1  field map(s) to be loaded :
          sol30m.map                                                                      


     Min/max fields seen in map   :   9.980853E-03                 /   1.600000E+01
       @  X,  Y, Z :  3.020E+03   0.00       0.00                  /  1.474E+03   0.00       0.00    

     Given normalisation coeffs on field, x, y, z :   9.593361E-02   1.000000E+00   0.000000E+00   0.000000E+00
       min/max normalised fields (kG) :  9.574992E-04                      1.534938E+00
       @  X (cm),  Y (cm), Z (cm) :  3.020E+03   0.00       0.00   /  1.474E+03   0.00       0.00    

     Length of element,  XL =  3.040000E+03 cm 
                                               from  XI =  -2.000000E+01 cm 
                                               to    XF =   3.020000E+03 cm 

     Nbr of nodes in X =1521;  nbr of nodes in Y =    1
     X-size of mesh =  2.000000E+00 cm ; Y-size =  0.000000E+00 cm

                     OPTION  DE  CALCUL  : 2

                    Integration step :   1.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000         3020.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.100     0.000     0.000     0.000         3020.000     0.069    -0.000    -0.038     0.000            2
  A    1  1.0000    -0.100     0.000     0.000     0.000         3020.000    -0.069     0.000     0.038    -0.000            3
  A    1  1.0000     0.000     0.100     0.000     0.000         3020.000     0.006     0.000    -0.003    -0.000            4
  A    1  1.0000     0.000    -0.100     0.000     0.000         3020.000    -0.006    -0.000     0.003     0.000            5
  A    1  1.0000     0.000     0.000     0.100     0.000         3020.000     0.038    -0.000     0.069    -0.000            6
  A    1  1.0000     0.000     0.000    -0.100     0.000         3020.000    -0.038     0.000    -0.069     0.000            7
  A    1  1.0000     0.000     0.000     0.000     0.100         3020.000     0.003     0.000     0.006     0.000            8
  A    1  1.0000     0.000     0.000     0.000    -0.100         3020.000    -0.003    -0.000    -0.006    -0.000            9
  A    1  1.0010     0.000     0.000     0.000     0.000         3020.000     0.000     0.000     0.000     0.000           10
  A    1  0.9990     0.000     0.000     0.000     0.000         3020.000     0.000     0.000     0.000     0.000           11

 Cumulative length of optical axis =    30.2000000     m ;  Time  (for ref. rigidity & particle) =   3.912628E-07 s 

************************************************************************************************************************************
      4  Keyword, label(s) :  DRIFT       4                     


                              Drift,  length =   -20.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00   3.0000000E+03   1.01403E-01

 Cumulative length of optical axis =    30.0000000     m  at  ;  Time  (for ref. rigidity & particle) =   3.886716E-07 s 

************************************************************************************************************************************
      5  Keyword, label(s) :  FAISCEAU    5                     

0                                             TRACE DU FAISCEAU

                                                 11 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

O  1   1.0000     0.000     0.000     0.000     0.000       0.0000    1.0000    0.000    0.000    0.000    0.000   3.000000E+03     1
A  1   1.0000     0.100     0.000     0.000     0.000       0.0000    1.0000    0.077   -0.394   -0.043    0.218   3.000001E+03     2
B  1   1.0000    -0.100     0.000     0.000     0.000       0.0000    1.0000   -0.077    0.394    0.043   -0.218   3.000001E+03     3
C  1   1.0000     0.000     0.100     0.000     0.000       0.0000    1.0000    0.004    0.077   -0.002   -0.043   3.000000E+03     4
D  1   1.0000     0.000    -0.100     0.000     0.000       0.0000    1.0000   -0.004   -0.077    0.002    0.043   3.000000E+03     5
E  1   1.0000     0.000     0.000     0.100     0.000       0.0000    1.0000    0.043   -0.218    0.077   -0.394   3.000001E+03     6
F  1   1.0000     0.000     0.000    -0.100     0.000       0.0000    1.0000   -0.043    0.218   -0.077    0.394   3.000001E+03     7
G  1   1.0000     0.000     0.000     0.000     0.100       0.0000    1.0000    0.002    0.043    0.004    0.077   3.000000E+03     8
H  1   1.0000     0.000     0.000     0.000    -0.100       0.0000    1.0000   -0.002   -0.043   -0.004   -0.077   3.000000E+03     9
I  1   1.0010     0.000     0.000     0.000     0.000       0.0000    1.0010    0.000    0.000    0.000    0.000   3.000000E+03    10
J  1   0.9990     0.000     0.000     0.000     0.000       0.0000    0.9990    0.000    0.000    0.000    0.000   3.000000E+03    11


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

   2.2848E-07    3.916        7.749        0.000        0.000          11      11    1.000    B-Dim     1      1
   2.2848E-07    3.916        7.749        0.000        0.000          11      11    1.000    B-Dim     2      1
   2.9468E-08  -1.8631E-14   2.0636E-07   0.1014        250.0          11      11    1.000    B-Dim     3      1


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  1.818150E-08  1.818150E-08
           AlpX, BetX = -3.915979E+00  7.749441E+00
           AlpZ, BetZ = -3.915979E+00  7.749441E+00

  1.408964E-07 -7.119836E-08 -2.047013E-20  4.357924E-11

 -7.119836E-08  3.832442E-08 -4.358089E-11  8.444700E-16

 -2.047013E-20 -4.358089E-11  1.408964E-07 -7.119836E-08

  4.357924E-11  8.444700E-16 -7.119836E-08  3.832441E-08

************************************************************************************************************************************
      6  Keyword, label(s) :  MATRIX      6                     


  Reference, before change of frame (part #     1)  : 
   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   3.00000000E+03   1.01403485E-01

           Frame for MATRIX calculation moved by :
            XC =    0.000 cm , YC =    0.000 cm ,   A =  0.00000 deg  ( = 0.000000 rad )


  Reference, after change of frame (part #     1)  : 
   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   3.00000000E+03   1.01403485E-01

  Reference particle (#     1), path length :   3000.0000     cm  relative momentum :    1.00000    


                  TRANSFER  MATRIX  ORDRE  1  (MKSA units)

          0.769172        0.441907        0.425152        0.244596         0.00000         0.00000    
         -0.394252        0.768688       -0.218234        0.426021         0.00000         0.00000    
         -0.425152       -0.244596        0.769172        0.441907         0.00000         0.00000    
          0.218234       -0.426021       -0.394252        0.768688         0.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         1.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         0.00000         1.00000    

          DetY-1 =      -0.2345240302,    DetZ-1 =      -0.2345240349

          R12=0 at  -0.5749     m,        R34=0 at  -0.5749     m

      First order symplectic conditions (expected values = 0) :
        -2.1123E-05   -2.1123E-05    8.6743E-04    4.8808E-04    4.8426E-04   -8.6742E-04

************************************************************************************************************************************
      7  Keyword, label(s) :  RESET       7                     


************************************************************************************************************************************
      8  Keyword, label(s) :  OBJET       8                     

                          MAGNETIC  RIGIDITY =        833.910 kG*cm

                                         CALCUL  DES  TRAJECTOIRES

                              OBJET  (5)  FORME  DE     11 POINTS 



                                Y (cm)         T (mrd)       Z (cm)        P (mrd)       S (cm)        dp/p 
               Sampling :          0.10          0.10          0.10          0.10           0.0          0.1000E-02
  Refrnce trajectry #      1 :      0.0           0.0           0.0           0.0           0.0           1.000    

************************************************************************************************************************************
      9  Keyword, label(s) :  SOLENOID    9                     


      -----  SOLENOID    : 
                Length  of  element  :    3000.      cm
                Inner radius  RO =   1.000      cm
                B-CNTRL  =   16.00         kG
                     XE =   20.00     cm,   XS =   20.00     cm
                MODL =  1 -> Solenoid model is axial field model   

                    Integration step :   1.000     cm

  A    1  1.0000     0.000     0.000     0.000     0.000         3040.000     0.000     0.000     0.000     0.000            1
  A    1  1.0000     0.100     0.000     0.000     0.000         3040.000     0.079    -0.000    -0.039     0.000            2
  A    1  1.0000    -0.100     0.000     0.000     0.000         3040.000    -0.079     0.000     0.039    -0.000            3
  A    1  1.0000     0.000     0.100     0.000     0.000         3040.000     0.004     0.000    -0.002    -0.000            4
  A    1  1.0000     0.000    -0.100     0.000     0.000         3040.000    -0.004    -0.000     0.002     0.000            5
  A    1  1.0000     0.000     0.000     0.100     0.000         3040.000     0.039    -0.000     0.079    -0.000            6
  A    1  1.0000     0.000     0.000    -0.100     0.000         3040.000    -0.039     0.000    -0.079     0.000            7
  A    1  1.0000     0.000     0.000     0.000     0.100         3040.000     0.002     0.000     0.004     0.000            8
  A    1  1.0000     0.000     0.000     0.000    -0.100         3040.000    -0.002    -0.000    -0.004    -0.000            9
  A    1  1.0010     0.000     0.000     0.000     0.000         3040.000     0.000     0.000     0.000     0.000           10
  A    1  0.9990     0.000     0.000     0.000     0.000         3040.000     0.000     0.000     0.000     0.000           11

 Cumulative length of optical axis =    30.0000000     m ;  Time  (for ref. rigidity & particle) =   3.886716E-07 s 

************************************************************************************************************************************
     10  Keyword, label(s) :  FAISCEAU    10                    

0                                             TRACE DU FAISCEAU

                                                 11 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

O  1   1.0000     0.000     0.000     0.000     0.000       0.0000    1.0000    0.000    0.000    0.000    0.000   3.000000E+03     1
A  1   1.0000     0.100     0.000     0.000     0.000       0.0000    1.0000    0.079   -0.394   -0.039    0.219   3.000001E+03     2
B  1   1.0000    -0.100     0.000     0.000     0.000       0.0000    1.0000   -0.079    0.394    0.039   -0.219   3.000001E+03     3
C  1   1.0000     0.000     0.100     0.000     0.000       0.0000    1.0000    0.004    0.075   -0.002   -0.046   3.000000E+03     4
D  1   1.0000     0.000    -0.100     0.000     0.000       0.0000    1.0000   -0.004   -0.075    0.002    0.046   3.000000E+03     5
E  1   1.0000     0.000     0.000     0.100     0.000       0.0000    1.0000    0.039   -0.219    0.079   -0.394   3.000001E+03     6
F  1   1.0000     0.000     0.000    -0.100     0.000       0.0000    1.0000   -0.039    0.219   -0.079    0.394   3.000001E+03     7
G  1   1.0000     0.000     0.000     0.000     0.100       0.0000    1.0000    0.002    0.046    0.004    0.075   3.000000E+03     8
H  1   1.0000     0.000     0.000     0.000    -0.100       0.0000    1.0000   -0.002   -0.046   -0.004   -0.075   3.000000E+03     9
I  1   1.0010     0.000     0.000     0.000     0.000       0.0000    1.0010    0.000    0.000    0.000    0.000   3.000000E+03    10
J  1   0.9990     0.000     0.000     0.000     0.000       0.0000    0.9990    0.000    0.000    0.000    0.000   3.000000E+03    11


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

   2.3170E-07    3.858        7.636        0.000        0.000          11      11    1.000    B-Dim     1      1
   2.3170E-07    3.858        7.636        0.000        0.000          11      11    1.000    B-Dim     2      1
   2.5227E-08  -1.8692E-14   1.7666E-07   0.1001        250.0          11      11    1.000    B-Dim     3      1


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  1.843846E-08  1.843846E-08
           AlpX, BetX = -3.858012E+00  7.635909E+00
           AlpZ, BetZ = -3.858012E+00  7.635909E+00

  1.407944E-07 -7.113580E-08  3.772041E-20  3.070005E-09

 -7.113580E-08  3.835578E-08 -3.070007E-09  8.426511E-16

  3.772041E-20 -3.070007E-09  1.407944E-07 -7.113580E-08

  3.070005E-09  8.426511E-16 -7.113580E-08  3.835578E-08

************************************************************************************************************************************
     11  Keyword, label(s) :  MATRIX      11                    


  Reference, before change of frame (part #     1)  : 
   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   3.00000000E+03   1.00069229E-01

           Frame for MATRIX calculation moved by :
            XC =    0.000 cm , YC =    0.000 cm ,   A =  0.00000 deg  ( = 0.000000 rad )


  Reference, after change of frame (part #     1)  : 
   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   0.00000000E+00   3.00000000E+03   1.00069229E-01

  Reference particle (#     1), path length :   3000.0000     cm  relative momentum :    1.00000    


                  TRANSFER  MATRIX  ORDRE  1  (MKSA units)

          0.785502        0.441849        0.393452        0.244692         0.00000         0.00000    
         -0.394003        0.751835       -0.219069        0.455477         0.00000         0.00000    
         -0.393452       -0.244692        0.785502        0.441849         0.00000         0.00000    
          0.219069       -0.455477       -0.394003        0.751835         0.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         1.00000         0.00000    
           0.00000         0.00000         0.00000         0.00000         0.00000         1.00000    

          DetY-1 =      -0.2353424906,    DetZ-1 =      -0.2353424954

          R12=0 at  -0.5877     m,        R34=0 at  -0.5877     m

      First order symplectic conditions (expected values = 0) :
        -2.5293E-03   -2.5293E-03    6.1581E-02    3.4568E-02    3.4116E-02   -6.1581E-02

************************************************************************************************************************************
     12  Keyword, label(s) :  END         12                    


                            11 particles have been launched
                     Made  it  to  the  end :     11

************************************************************************************************************************************

           MAIN PROGRAM : Execution ended upon key  END       

************************************************************************************************************************************
   
            ZGOUBI RUN COMPLETED. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  06-02-0015,  at  16:57:47 
  JOB  ENDED  ON    06-02-0015,  AT  16:57:47 

   CPU time, total :    4.099300000000000E-002
