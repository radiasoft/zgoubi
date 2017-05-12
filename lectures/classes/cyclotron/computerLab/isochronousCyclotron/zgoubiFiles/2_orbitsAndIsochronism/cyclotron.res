Test Half-Cell                                                                                                
 'OBJET'                                                                                                      1
1249.382414                                                                                                   
2                                                                                                             
6 1                                                                                                           
268.20121069999999       -15.580272709999999        0.  0.  0.         1.4000000059604645      'o'            
282.46123799999998       -26.403956470000001        0.  0.  0.         1.5000000074505806      'o'            
296.05739569999997       -38.103638570000001        0.  0.  0.         1.6000000089406967      'o'            
309.00757420000002       -50.260491320000000        0.  0.  0.         1.7000000104308128      'o'            
321.32946040000002       -62.556048799999999        0.  0.  0.         1.8000000119209290      'o'            
333.05392440000003       -74.952056110000001        0.  0.  0.         1.9000000134110451      'o'            
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1                                                                           
 'PARTICUL'                                                                                                   2
938.2723   1.602176487D-19  1. 0. 0.                                                                          
                                                                                                              
 'CYCLOTRON'                                                                                                  3
2   plot_spiral.H                                                                                             
1   45.0  2.760000E+02    1.0                                                                                 
0. 0. 0.99212277 51.4590015 .5 800. -0.476376328 2.27602517e-03 -4.81955890e-06 3.94715806e-09                
18.3000E+00  1.   28 -2.0                 g0,  k , G0, G1    Entrance face                                    
8 1.10243581,  3.12915071, -3.14287154,  3.0858059 , -1.43544992, 0.24047436 0. 0. 0.                         
11.0   3.5  35.E-3   0.E-4   3.E-8     1.  1.  1.       omega+, xi0,xi1,xi2,xi3,a,b,c                         
18.3000E+00  1.   28. -2.0                g0,  k , G0, G1    Exit face                                        
8  0.70490173, 4.16013051, -4.33095751,  3.54041582, -1.34727027, 0.18261076  0. 0. 0.                        
-8.5  2.  12.E-3   75.E-6   0.E-6     1.  1.  1.       omega-, xi0s,xi1s,xi2s,xi3s,aexit,bexit,cexit          
0. -1                                       g0,  k             Lateral face                                   
0  0.   0.   0.   0.   0. 0.  0.                        NC, C0...C5, shift                                    
0.  0.   0.    0.    0. 0.                              omega+, xi, 4                                         
2    10.   numerical field/derivatives,  flying mesh size is 0.415/10. (KIRD,RESOL)                           
1.0                                                                                                           
2                                                                                                             
0.  0. 0.  0.                                                                                                 
                                                                                                              
 'FAISCEAU'                                                                                                   4
                                                                                                              
 'END'                                                                                                        5

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       1249.382 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       6 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                          


     Particle  properties :

                     Mass          =    938.272        MeV/c2
                     Charge        =   1.602176E-19    C     
                     G  factor     =    1.00000              

              Reference  data :
                    mag. rigidity (kG.cm)   :   1249.3824      =p/q, such that dev.=B*L/rigidity
                    mass (MeV/c2)           :   938.27230    
                    momentum (MeV/c)        :   374.55542    
                    energy, total (MeV)     :   1010.2706    
                    energy, kinetic (MeV)   :   71.998295    
                    beta = v/c              :  0.3707476261    
                    gamma                   :   1.076734968    
                    beta*gamma              :  0.3991969334    
                    G*gamma                 :   1.076734968    
                    electric rigidity (MeV) :   138.8655346      =T[eV]*(gamma+1)/gamma, such that dev.=E*L/rigidity
  
 I, AMQ(1,I), AMQ(2,I)/QE, P/Pref, v/c, time :
  
     1   9.38272300E+02  1.00000000E+00  1.40000001E+00  4.87856067E-01  0.00000000E+00
     2   9.38272300E+02  1.00000000E+00  1.50000001E+00  5.13735640E-01  0.00000000E+00
     3   9.38272300E+02  1.00000000E+00  1.60000001E+00  5.38285263E-01  0.00000000E+00
     4   9.38272300E+02  1.00000000E+00  1.70000001E+00  5.61537317E-01  0.00000000E+00
     5   9.38272300E+02  1.00000000E+00  1.80000001E+00  5.83531351E-01  0.00000000E+00
     6   9.38272300E+02  1.00000000E+00  1.90000001E+00  6.04312387E-01  0.00000000E+00

************************************************************************************************************************************
      3  Keyword, label(s) :  CYCLOTRON                         


                OPEN FILE zgoubi.plt                                                                      
                FOR PRINTING TRAJECTORIES

                    FFAG  N-tuple,  number  of  dipoles  N :  1

            Total angular extent of the magnet :  45.00 Degres
            Reference geometrical radius R0  :     276.00 cm

     Dipole # 1
            Positionning  angle ACENT :    0.000     degrees
            Positionning  wrt.  R0  :    0.0     cm
            B0 =   51.46     kGauss,       K =  0.50000    

     Entrance  EFB
          Fringe  field  :  gap at R0 is  18.30 cm,    type is :  g_0(1-r**2) 1.00
           COEFFICIENTS :  8   1.10244   3.12915  -3.14287   3.08581  -1.43545   0.24047
           Shift  of  EFB  is            0 cm

          OMEGA =  11.00 deg.     Spiral  angle  =   3.50 deg.

         Exit  EFB
          Fringe  field  :  gap at R0 is  18.30 cm,    type is :  g_0(1-r**2) 1.00
           COEFFICIENTS :  8   0.70490   4.16013  -4.33096   3.54042  -1.34727   0.18261
           Shift  of  EFB  is            0 cm

          OMEGA =  -8.50 deg.     Spiral  angle  =   2.00 deg.
 
         Lateral face :  unused
 

      Field & deriv. calculation : interpolation
                    3*3-point  interpolation, size of flying mesh :  integration step /   10.00    

                    Integration step :   1.000     cm   (i.e.,   3.6232E-03 rad  at mean radius RM =    276.0    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.4000   268.201   -15.580     0.000     0.000            0.393   268.268    -0.016     0.000     0.000            1
  A    1  1.5000   282.461   -26.404     0.000     0.000            0.393   282.434    -0.030     0.000     0.000            2
  A    1  1.6000   296.057   -38.104     0.000     0.000            0.393   295.929    -0.043     0.000     0.000            3
  A    1  1.7000   309.008   -50.260     0.000     0.000            0.393   308.830    -0.056     0.000     0.000            4
  A    1  1.8000   321.329   -62.556     0.000     0.000            0.393   321.169    -0.068     0.000     0.000            5
  A    1  1.9000   333.054   -74.952     0.000     0.000            0.393   332.939    -0.079     0.000     0.000            6


                CONDITIONS  DE  MAXWELL  (     1463.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    53.8200000     m ;  Time  (for ref. rigidity & particle) =   4.842221E-07 s 

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  6 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   1.4000   268.201   -15.580     0.000     0.000       0.0000    1.4000  268.268  -15.802    0.000    0.000   2.166099E+02     1
               Time of flight (mus) :  1.48103721E-02 mass (MeV/c2) :   938.272    
o  1   1.5000   282.461   -26.404     0.000     0.000       0.0000    1.5000  282.434  -29.930    0.000    0.000   2.281330E+02     2
               Time of flight (mus) :  1.48124804E-02 mass (MeV/c2) :   938.272    
o  1   1.6000   296.057   -38.104     0.000     0.000       0.0000    1.6000  295.929  -43.489    0.000    0.000   2.390438E+02     3
               Time of flight (mus) :  1.48130438E-02 mass (MeV/c2) :   938.272    
o  1   1.7000   309.008   -50.260     0.000     0.000       0.0000    1.7000  308.830  -56.195    0.000    0.000   2.493744E+02     4
               Time of flight (mus) :  1.48133252E-02 mass (MeV/c2) :   938.272    
o  1   1.8000   321.329   -62.556     0.000     0.000       0.0000    1.8000  321.169  -67.995    0.000    0.000   2.591432E+02     5
               Time of flight (mus) :  1.48134077E-02 mass (MeV/c2) :   938.272    
o  1   1.9000   333.054   -74.952     0.000     0.000       0.0000    1.9000  332.939  -79.024    0.000    0.000   2.683616E+02     6
               Time of flight (mus) :  1.48128373E-02 mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

   7.9637E-04    75.42        770.0        3.016      -4.8739E-02       6       6    1.000    B-Dim     1      1
    0.000        75.42        770.0        0.000        0.000           6       6    1.000    B-Dim     2      1
   3.2319E-04  -0.9964       4.1519E-08   1.4813E-02    186.5           6       6    1.000    B-Dim     3      1


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  6.337310E-05  0.000000E+00
           AlpX, BetX = -7.542375E+01  7.700024E+02
           AlpZ, BetZ =           NaN  0.000000E+00

  4.879744E-02 -4.779837E-03  0.000000E+00  0.000000E+00

 -4.779837E-03  4.682798E-04  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      5  Keyword, label(s) :  END                               


                             6 particles have been launched
                     Made  it  to  the  end :      6

************************************************************************************************************************************

           MAIN PROGRAM : Execution ended upon key  END       

************************************************************************************************************************************
   
            Zgoubi run completed. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  17/09/14 ,  at  10:13:18 
  Job  ended  on    17/09/14 ,  at  10:13:18 

   CPU time, total :    0.204012000000000     
