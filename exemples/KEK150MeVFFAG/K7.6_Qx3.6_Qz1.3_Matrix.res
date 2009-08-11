FFAG triplet. 150MeV machine. CPU time, analyt. : 11:01:35->                                                  
  'OBJET'                                                                            1                             1
 1839.090113 150MeV                                                                                           
 5                                                                                                            
0.001 0.001 0.001 0.001 0.001 0.0001                                                                          
479.70   0.   0.00  0. 0.  .52   'o'     43MeV                                                                
 'PARTICUL'                                                                          2                             2
938.2723 1.60217733D-19 0. 0. 0.                                                                              
'FAISTORE'                                                                           3                             3
  b_zgoubi.fai #START                                                                                         
1                                                                                                             
 'FFAG'                                                                              4                             4
 0                                                                                                            
 3   30.     540.                        NMAG, AT=tetaF+2tetaD+2Atan(XFF/R0), R0                              
 6.465  0.  -1.21744691E+01  7.6             mag 1 : ACNT, dum, B0, K                                         
 6.3  03.                                EFB 1 : lambda, gap const/var=0/.ne.0                                
 4  .1455   2.2670  -.6395  1.1558  0. 0.  0.                                                                 
 1.715 0.    1.E6  -1.E6  1.E6  1.E6                                                                          
 6.3  03.                                EFB 2                                                                
 4  .1455   2.2670  -.6395  1.1558  0. 0.  0.                                                                 
-1.715  0.    1.E6  -1.E6  1.E6  1.E6                                                                         
 0. -1                                   EFB 3 : inhibited by iop=0                                           
 0  0.      0.      0.      0.      0. 0.  0.                                                                 
 0.  0.   0.    0.    0. 0.                                                                                   
 15.   0.   1.69055873E+01  7.6                  mag 2 : ACNT, dum, B0, K,dummie                              
 6.3  03.                                EFB 1                                                                
 4  .1455   2.2670  -.6395  1.1558  0. 0.  0.                                                                 
 5.12   0.    1.E6  -1.E6  1.E6  1.E6                                                                         
 6.3  03.                                EFB 2                                                                
 4  .1455   2.2670  -.6395  1.1558  0. 0.  0.                                                                 
-5.12   0.    1.E6  -1.E6  1.E6  1.E6                                                                         
 0. -1                                   EFB 3                                                                
 0  0.      0.      0.      0.      0. 0.  0.                                                                 
 0.  0.   0.    0.    0. 0.                                                                                   
23.535  0.  -1.21744691E+01  7.6             mag 3 : ACNT, dum, B0, K                                         
 6.3  03.                                EFB 1                                                                
 4  .1455   2.2670  -.6395  1.1558  0. 0.  0.                                                                 
 1.715 0.    1.E6  -1.E6  1.E6  1.E6                                                                          
 6.3  03.                                EFB 2                                                                
 4  .1455   2.2670  -.6395  1.1558  0. 0.  0.                                                                 
-1.715  0.    1.E6  -1.E6  1.E6  1.E6                                                                         
 0. -1                                   EFB 3                                                                
 0  0.      0.      0.      0.      0. 0.  0.                                                                 
 0.  0.   0.    0.    0. 0.                                                                                   
  0 2  125.                   KIRD anal/num (=0/2,25,4), resol(mesh=step/resol)                               
 .25                           integration step size (cm)                                                     
 2   0.  0. 0. 0.        (comp anal/num yield best num prec for resol=1e3*step)                               
'DRIFT'                                                                              5                             5
0.                                                                                                            
'MATRIX'                                                                             6                             6
1 11                                                                                                          
'END'                                                                                                              7

********************************************************************************************************************************
      1  OBJET       1                 
                          MAGNETIC  RIGIDITY =       1839.090 kG*cm

                                         CALCUL  DES  TRAJECTOIRES

                              OBJET  (5)  FORME  DE     11 POINTS 



                                Y (cm)         T (mrd)       Z (cm)        P (mrd)       S (cm)        dp/p 
               Sampling :          0.10E-02      0.10E-02      0.10E-02      0.10E-02      0.10E-02      0.1000E-03
  Reference trajectory #1 :        0.48E+03       0.0           0.0           0.0           0.0          0.5200    

********************************************************************************************************************************
      2  PARTICUL    2                 

     Particle  properties :

                     Mass          =    938.272        MeV/c2
                     Charge        =   1.602177E-19    C     

              Reference  data :
                    rigidity (kG.cm)      :   1839.09    
                    mass (MeV/c2)         :   938.272    
                    momentum (MeV/c)      :   551.346    
                    energy, total (MeV)   :   1088.27    
                    energy, kinetic (MeV) :   150.000    
                    beta = v/c            :  0.5066246390    
                    gamma                 :   1.159868459    
                    beta*gamma            :  0.5876179394    

********************************************************************************************************************************
      3  FAISTORE    3                 

                OPEN FILE b_zgoubi.fai                                                                    
                FOR PRINTING COORDINATES 

               Print will occur at element[s] labeled : 
                    i         
                    #START    

********************************************************************************************************************************
      4  FFAG        4                 
                    FFAG  N-tuple,  number  of  dipoles,  N :  3

            Total angular extent of the magnet :  30.0000 Degres
            Reference geometrical radius R0  :     540.0000 cm

     Dipole # 1
            Positionning  angle ACENT :    6.46500     degrees
            Positionning  wrt.  R0  :    0.00000     cm
            B0 =  -12.174469     kGauss,       K =   7.6000    

     Entrance  EFB
          Fringe  field  :  gap at R0 is   6.30 cm,    type is :  g_0(r0/r)^K
           COEFFICIENTS :  4   0.14550   2.26700  -0.63950   1.15580   0.00000   0.00000
           Shift  of  EFB  is    0.000     cm

          OMEGA =   1.72 deg.     Wedge  angle  =   0.00 deg.
           Radius 1 =  1.00E+06 CM
           Straight  segment 1 = -1.00E+06 CM
           Straight  segment 2 =  1.00E+06 CM
           Radius 2 =  1.00E+06 CM

         Exit  EFB
          Fringe  field  :  gap at R0 is   6.30 cm,    type is :  g_0(r0/r)^K
           COEFFICIENTS :  4   0.14550   2.26700  -0.63950   1.15580   0.00000   0.00000
           Shift  of  EFB  is    0.000     cm

          OMEGA =  -1.72 deg.     Wedge  angle  =   0.00 deg.
           Radius 1 =  1.00E+06 CM
           Straight  segment 1 = -1.00E+06 CM
           Straight  segment 2 =  1.00E+06 CM
           Radius 2 =  1.00E+06 CM

         Lateral face :  unused

     Dipole # 2
            Positionning  angle ACENT :    15.0000     degrees
            Positionning  wrt.  R0  :    0.00000     cm
            B0 =   16.905587     kGauss,       K =   7.6000    

     Entrance  EFB
          Fringe  field  :  gap at R0 is   6.30 cm,    type is :  g_0(r0/r)^K
           COEFFICIENTS :  4   0.14550   2.26700  -0.63950   1.15580   0.00000   0.00000
           Shift  of  EFB  is    0.000     cm

          OMEGA =   5.12 deg.     Wedge  angle  =   0.00 deg.
           Radius 1 =  1.00E+06 CM
           Straight  segment 1 = -1.00E+06 CM
           Straight  segment 2 =  1.00E+06 CM
           Radius 2 =  1.00E+06 CM

         Exit  EFB
          Fringe  field  :  gap at R0 is   6.30 cm,    type is :  g_0(r0/r)^K
           COEFFICIENTS :  4   0.14550   2.26700  -0.63950   1.15580   0.00000   0.00000
           Shift  of  EFB  is    0.000     cm

          OMEGA =  -5.12 deg.     Wedge  angle  =   0.00 deg.
           Radius 1 =  1.00E+06 CM
           Straight  segment 1 = -1.00E+06 CM
           Straight  segment 2 =  1.00E+06 CM
           Radius 2 =  1.00E+06 CM

         Lateral face :  unused

     Dipole # 3
            Positionning  angle ACENT :    23.5350     degrees
            Positionning  wrt.  R0  :    0.00000     cm
            B0 =  -12.174469     kGauss,       K =   7.6000    

     Entrance  EFB
          Fringe  field  :  gap at R0 is   6.30 cm,    type is :  g_0(r0/r)^K
           COEFFICIENTS :  4   0.14550   2.26700  -0.63950   1.15580   0.00000   0.00000
           Shift  of  EFB  is    0.000     cm

          OMEGA =   1.72 deg.     Wedge  angle  =   0.00 deg.
           Radius 1 =  1.00E+06 CM
           Straight  segment 1 = -1.00E+06 CM
           Straight  segment 2 =  1.00E+06 CM
           Radius 2 =  1.00E+06 CM

         Exit  EFB
          Fringe  field  :  gap at R0 is   6.30 cm,    type is :  g_0(r0/r)^K
           COEFFICIENTS :  4   0.14550   2.26700  -0.63950   1.15580   0.00000   0.00000
           Shift  of  EFB  is    0.000     cm

          OMEGA =  -1.72 deg.     Wedge  angle  =   0.00 deg.
           Radius 1 =  1.00E+06 CM
           Straight  segment 1 = -1.00E+06 CM
           Straight  segment 2 =  1.00E+06 CM
           Radius 2 =  1.00E+06 CM

         Lateral face :  unused


      Field & deriv. calculation : analytic     
                    Derivatives computed to order 2

                    Integration step :  0.2500     cm   (i.e.,   4.6296E-04 rad  at mean radius RM =    540.0    )

                               Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


 Cumulative length of optical axis =    2.93401     m ;   corresponding Time  (for ref. rigidity & particle) =   1.9318E-08 s 

********************************************************************************************************************************
      5  DRIFT       5                 

                              Drift,  length =     0.00000  cm

 TRAJ 1 IEX,D,Y,T,Z,P,S,time :  1  0.5200       479.7      6.8901E-03   0.000       0.000           259.90          2.96672E-02

 Cumulative length of optical axis =    2.93401     m ;   corresponding Time  (for ref. rigidity & particle) =   1.9318E-08 s 

********************************************************************************************************************************
      6  MATRIX      6                 

           Frame for MATRIX calculation moved by :
            XC =    0.000 cm , YC =  479.701 cm ,   A =  0.00039 deg  ( = 0.000007 rad )


Path length of particle  1 :   259.90391     cm


                  TRANSFER  MATRIX  ORDRE  1  (MKSA units)

         -0.393264        0.727864         0.00000         0.00000         0.00000        0.774519    
          -1.16140       -0.393267         0.00000         0.00000         0.00000        0.645625    
           0.00000         0.00000        0.771788         2.73618         0.00000         0.00000    
           0.00000         0.00000       -0.147777        0.771788         0.00000         0.00000    
          0.645626        0.774520         0.00000         0.00000         1.00000       -5.691068E-02
           0.00000         0.00000         0.00000         0.00000         0.00000         1.00000    

          DetY-1 =       0.0000000054,    DetZ-1 =       0.0000000447

          R12=0 at    1.851     m,        R34=0 at   -3.545     m

      First order sympletic conditions (expected values = 0) :
         5.4377E-09    4.4661E-08     0.000         0.000         0.000         0.000    


                TWISS  parameters,  periodicity  of   1  is  assumed :

       Beam  matrix  (beta/-alpha/-alpha/gamma) and  periodic  dispersion  (MKSA units)

           0.791652    -0.000002     0.000000     0.000000     0.000000     0.555903
          -0.000002     1.263182     0.000000     0.000000     0.000000     0.000000
           0.000000     0.000000     4.302980     0.000000     0.000000     0.000000
           0.000000     0.000000     0.000000     0.232397     0.000000     0.000000
           0.000000     0.000000     0.000000     0.000000     0.000000     0.000000
           0.000000     0.000000     0.000000     0.000000     0.000000     0.000000

                                   Betatron  tunes

                    NU_Y = 0.31432737         NU_Z = 0.10968143    

********************************************************************************************************************************
      7  END         7                 

                            11 particles have been launched
                     Made  it  to  the  end :     11

********************************************************************************************************************************

           MAIN PROGRAM : Execution ended upon key  END       

********************************************************************************************************************************
  Job  ended  on  11-Aug-09,  at  13:15:15 

   CPU time, total :    0.128008000000000     
