Warm snake, 3-D map with 'TOSCA'
 'OBJET'                                                                                                      1
7.20517829566983d3
2
1 1
-2.02391047    -0.00065841    -0.00222492     0.14641297     0.00000000     1. 'o'
1 1 1 1 1
 
 'FAISCEAU'                                                                                                   2
 
 'PARTICUL'                                                                                                   3
9.3827203E+02 1.602176487E-19 1.7928474 0 0
 'SPNTRK'                                                                                                     4
4.1
0. 0. 1.
 
 'DRIFT' DRIF      WSNK                                                                                       5
-50.
 'TOSCA'                                                                                                      6
0  2
1.e1   100. 100. 100.
HEADER_0 wsnake
801 29 29 1
Wsnk3D/table55_p.070
Wsnk3D/table55_p.065
Wsnk3D/table55_p.060
Wsnk3D/table55_p.055
Wsnk3D/table55_p.050
Wsnk3D/table55_p.045
Wsnk3D/table55_p.040
Wsnk3D/table55_p.035
Wsnk3D/table55_p.030
Wsnk3D/table55_p.025
Wsnk3D/table55_p.020
Wsnk3D/table55_p.015
Wsnk3D/table55_p.010
Wsnk3D/table55_p.005
Wsnk3D/table55_p.000
Wsnk3D/table55_m.005
Wsnk3D/table55_m.010
Wsnk3D/table55_m.015
Wsnk3D/table55_m.020
Wsnk3D/table55_m.025
Wsnk3D/table55_m.030
Wsnk3D/table55_m.035
Wsnk3D/table55_m.040
Wsnk3D/table55_m.045
Wsnk3D/table55_m.050
Wsnk3D/table55_m.055
Wsnk3D/table55_m.060
Wsnk3D/table55_m.065
Wsnk3D/table55_m.070
0 0 0 0
2
.1
2  0.  .00  0.  0.
 'DRIFT' DRIF      WSNK                                                                                       7
-50.
 
 'FAISCEAU'                                                                                                   8
 
 'FIT'                                                                                                        9
4
1 30 0  [-3,3]
1 31 0  [-3,3]
1 32 0  [-3,3]
1 33 0  [-3,3]
6
3.1 1 2 6 0. 1. 0
3.1 1 3 6 0. 1. 0
3.1 1 4 6 0. 1. 0
3.1 1 5 6 0. 1. 0
7.3 1 2 6 0. 100. 1 2
7.3 1 4 6 0. 100. 1 2
 
 'SPNPRT'   PRINT                                                                                            10
 
 'REBELOTE'         ! Will re-do the FIT for 4 different rigidities                                          11
4  0.1  0 1
1                   ! List of 4 successive values of D in OBJET follows
OBJET  35  1.3872739973  2.1368296674  4.8261190694  11.015241321
 
! Gg=6     6             9             20            45.5
 
 'END'                                                                                                       12

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      1)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   1.0000    -2.024    -0.001    -0.002     0.146       0.0000    1.0000   -2.024   -0.001   -0.002    0.146   0.000000E+00     1


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -2.023910E-02  -6.584100E-07        1        1    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00  -2.224920E-05   1.464130E-04        1        1    1.000      (Z,P)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   2.160058E+03        1        1    1.000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      3  Keyword, label(s) :  PARTICUL                          


     Particle  properties :

                     Mass          =    938.272        MeV/c2
                     Charge        =   1.602176E-19    C     
                     G  factor     =    1.79285              

              Reference  data :
                    mag. rigidity (kG.cm)   :   7205.1783      =p/q, such that dev.=B*L/rigidity
                    mass (MeV/c2)           :   938.27203    
                    momentum (MeV/c)        :   2160.0581    
                    energy, total (MeV)     :   2355.0383    
                    energy, kinetic (MeV)   :   1416.7663    
                    beta = v/c              :  0.9172072069    
                    gamma                   :   2.509973905    
                    beta*gamma              :   2.302166155    
                    G*gamma                 :   4.500000190    
                    electric rigidity (MeV) :   1981.220867      =T[eV]*(gamma+1)/gamma, such that dev.=E*L/rigidity
  
 I, AMQ(1,I), AMQ(2,I)/QE, P/Pref, v/c, time, s :
  
     1   9.38272030E+02  1.00000000E+00  1.00000000E+00  9.17207207E-01  0.00000000E+00  0.00000000E+00

************************************************************************************************************************************
      4  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              Same spin for all particles
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      5  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00 -2.023878E+00 -6.584100E-04 -9.545569E-03  1.464130E-01  -5.0000001E+01  -1.81837E-03
TRAJ #1 SX, SY, SZ, |S| :  1    0.000000E+00   0.000000E+00   1.000000E+00   1.000000E+00

 Cumulative length of optical axis =  -0.500000000     m   ;  Time  (for reference rigidity & particle) =  -1.818368E-09 s 

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


                OPEN FILE zgoubi.plt                                                                      
                FOR PRINTING TRAJECTORIES


     NDIM =   3 ;  Number of data file sets used is  29 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 1.0


           New field map(s) now used, cartesian mesh (MOD.le.19) ; 
           name(s) of map data file(s) : 

          Wsnk3D/table55_p.070
          Wsnk3D/table55_p.065
          Wsnk3D/table55_p.060
          Wsnk3D/table55_p.055
          Wsnk3D/table55_p.050
          Wsnk3D/table55_p.045
          Wsnk3D/table55_p.040
          Wsnk3D/table55_p.035
          Wsnk3D/table55_p.030
          Wsnk3D/table55_p.025
          Wsnk3D/table55_p.020
          Wsnk3D/table55_p.015
          Wsnk3D/table55_p.010
          Wsnk3D/table55_p.005
          Wsnk3D/table55_p.000
          Wsnk3D/table55_m.005
          Wsnk3D/table55_m.010
          Wsnk3D/table55_m.015
          Wsnk3D/table55_m.020
          Wsnk3D/table55_m.025
          Wsnk3D/table55_m.030
          Wsnk3D/table55_m.035
          Wsnk3D/table55_m.040
          Wsnk3D/table55_m.045
          Wsnk3D/table55_m.050
          Wsnk3D/table55_m.055
          Wsnk3D/table55_m.060
          Wsnk3D/table55_m.065
          Wsnk3D/table55_m.070

   ----
   Map file number    1 ( of 29) :

     Wsnk3D/table55_p.070 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number    2 ( of 29) :

     Wsnk3D/table55_p.065 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number    3 ( of 29) :

     Wsnk3D/table55_p.060 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number    4 ( of 29) :

     Wsnk3D/table55_p.055 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number    5 ( of 29) :

     Wsnk3D/table55_p.050 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number    6 ( of 29) :

     Wsnk3D/table55_p.045 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number    7 ( of 29) :

     Wsnk3D/table55_p.040 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number    8 ( of 29) :

     Wsnk3D/table55_p.035 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number    9 ( of 29) :

     Wsnk3D/table55_p.030 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   10 ( of 29) :

     Wsnk3D/table55_p.025 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   11 ( of 29) :

     Wsnk3D/table55_p.020 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   12 ( of 29) :

     Wsnk3D/table55_p.015 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   13 ( of 29) :

     Wsnk3D/table55_p.010 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   14 ( of 29) :

     Wsnk3D/table55_p.005 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   15 ( of 29) :

     Wsnk3D/table55_p.000 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   16 ( of 29) :

     Wsnk3D/table55_m.005 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   17 ( of 29) :

     Wsnk3D/table55_m.010 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   18 ( of 29) :

     Wsnk3D/table55_m.015 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   19 ( of 29) :

     Wsnk3D/table55_m.020 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   20 ( of 29) :

     Wsnk3D/table55_m.025 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   21 ( of 29) :

     Wsnk3D/table55_m.030 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   22 ( of 29) :

     Wsnk3D/table55_m.035 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   23 ( of 29) :

     Wsnk3D/table55_m.040 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   24 ( of 29) :

     Wsnk3D/table55_m.045 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   25 ( of 29) :

     Wsnk3D/table55_m.050 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   26 ( of 29) :

     Wsnk3D/table55_m.055 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   27 ( of 29) :

     Wsnk3D/table55_m.060 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   28 ( of 29) :

     Wsnk3D/table55_m.065 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.

   ----
   Map file number   29 ( of 29) :

     Wsnk3D/table55_m.070 map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.


     Min/max fields seen in map   :  -2.674110E-01                 /   3.119660E-01
       @  X,  Y, Z :  -1.04       6.00      -7.00                  /   1.05      -5.50      -7.00    

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.674110E+00                      3.119660E+00
       @  X (cm),  Y (cm), Z (cm) :  -104.       600.      -700.   /   105.      -550.      -700.    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =  -2.000000E+02 cm 
                                               to    XF =   2.000000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   29 ; Step in Z = -0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1  1.0000    -2.024    -0.001    -0.002     0.146          200.000    -2.023     0.000    -0.002     0.000            1


                CONDITIONS  DE  MAXWELL  (     4004.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -3.2026E-12      -1.0458E-13       -4.9555E-15
                                      1.6232E-12       -1.2802E-12
                                     -2.7349E-13       -8.2139E-12
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    3.50000000     m ;  Time  (for ref. rigidity & particle) =   1.272858E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00 -2.023284E+00  3.134711E-03 -9.675949E-03  1.463896E-01   3.0034191E+02   1.09226E-02
TRAJ #1 SX, SY, SZ, |S| :  1   -1.007678E-03   2.083835E-01   9.780467E-01   1.000000E+00

 Cumulative length of optical axis =    3.00000000     m   ;  Time  (for reference rigidity & particle) =   1.091021E-08 s 

************************************************************************************************************************************
      8  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      7)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   1.0000    -2.024    -0.001    -0.002     0.146       0.0000    1.0000   -2.023    0.003   -0.010    0.146   3.003419E+02     1
               Time of flight (mus) :  1.09226438E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -2.023284E-02   3.134711E-06        1        1    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00  -9.675949E-05   1.463896E-04        1        1    1.000      (Z,P)         1
   0.0000E+00   0.0000E+00   1.0000E+00   1.092264E-02   1.416766E+03        1        0    0.000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      9  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            1 ,   ok.
           variable #            2       IP =           31 ,   ok.
           variable #            3       IR =            1 ,   ok.
           variable #            3       IP =           32 ,   ok.
           variable #            4       IR =            1 ,   ok.
           variable #            4       IP =           33 ,   ok.
           constraint #            1       IR =            6 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            6 ,   ok.
           constraint #            2       I  =            1 ,   ok.
           constraint #            3       IR =            6 ,   ok.
           constraint #            3       I  =            1 ,   ok.
           constraint #            4       IR =            6 ,   ok.
           constraint #            4       I  =            1 ,   ok.
           constraint #            5       IR =            6 ,   ok.
           constraint #            6       IR =            6 ,   ok.

                    FIT  variables  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save [FileName]' command

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
  LMNT  VAR  PARAM   MINIMUM     INITIAL         FINAL         MAXIMUM     STEP     NAME       LBL1     LBL2
     1    1     30   -3.00       -2.02       -2.015381132       3.00      5.011E-42 OBJET      *          *         
     1    2     31   -3.00      -6.557E-04  -6.5566274669E-04   3.00      5.011E-42 OBJET      *          *         
     1    3     32   -3.00      -2.315E-03  -2.3147874791E-03   3.00      5.011E-42 OBJET      *          *         
     1    4     33   -3.00       0.146       0.1460889453       3.00      5.011E-42 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     6      0.0000000E+00     1.0000E+00    1.1876265E-05   8.0865E-06   TOSCA      *          *          *   0 : 
   3   1   3     6      0.0000000E+00     1.0000E+00    1.7894464E-04   1.8359E-03   TOSCA      *          *          *   0 : 
   3   1   4     6      0.0000000E+00     1.0000E+00    1.3662339E-05   1.0702E-05   TOSCA      *          *          *   0 : 
   3   1   5     6      0.0000000E+00     1.0000E+00    5.5387513E-06   1.7588E-06   TOSCA      *          *          *   0 : 
   7   1   2     6      0.0000000E+00     1.0000E+02   -4.1707491E-01   9.9731E-01   TOSCA      *          *          *   1 :   2.0E+00/
   7   1   4     6      0.0000000E+00     1.0000E+02   -1.2043324E-02   8.3156E-04   TOSCA      *          *          *   1 :   2.0E+00/
 Fit reached penalty value   1.7442E-05

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      1)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   1.0000    -2.015    -0.001    -0.002     0.146       0.0000    1.0000   -2.015   -0.001   -0.002    0.146   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -2.015381E-02  -6.556627E-07        1        1    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00  -2.314787E-05   1.460889E-04        1        1    1.000      (Z,P)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   1.416766E+03        1        0    0.000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      3  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      4  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              Same spin for all particles
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      5  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00 -2.015348E+00 -6.556627E-04 -9.619235E-03  1.460889E-01  -5.0000001E+01  -1.81837E-03
TRAJ #1 SX, SY, SZ, |S| :  1    0.000000E+00   0.000000E+00   1.000000E+00   1.000000E+00

 Cumulative length of optical axis =    3.00000000     m   ;  Time  (for reference rigidity & particle) =   1.091021E-08 s 

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is  29 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 1.0

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.070
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.065
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.060
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.055
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.050
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.045
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.040
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.035
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.030
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.025
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.020
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.015
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.010
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.005
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.000
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.005
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.010
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.015
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.020
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.025
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.030
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.035
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.040
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.045
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.050
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.055
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.060
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.065
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.070
 Pgm toscac,  restored mesh coordinates for field map #   1,  name : Wsnk3D/table55_m.070


     Min/max fields seen in map   :  -2.674110E-01                 /   3.119660E-01
       @  X,  Y, Z :  -1.04       6.00      -7.00                  /   1.05      -5.50      -7.00    

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.674110E+00                      3.119660E+00
       @  X (cm),  Y (cm), Z (cm) :  -104.       600.      -700.   /   105.      -550.      -700.    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =  -2.000000E+02 cm 
                                               to    XF =   2.000000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   29 ; Step in Z = -0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1  1.0000    -2.015    -0.001    -0.002     0.146          200.000    -2.015     0.000    -0.002     0.000            1


                CONDITIONS  DE  MAXWELL  (     4004.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -3.1901E-12      -1.0458E-13       -4.5069E-15
                                      1.6232E-12       -1.2802E-12
                                     -2.7349E-13       -8.2142E-12
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    7.00000000     m ;  Time  (for ref. rigidity & particle) =   2.545716E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  0.000000E+00 -2.015351E+00 -8.346074E-04 -9.632620E-03  1.460834E-01   3.0034191E+02   1.09226E-02
TRAJ #1 SX, SY, SZ, |S| :  1   -1.000116E-03   2.083811E-01   9.780472E-01   1.000000E+00

 Cumulative length of optical axis =    6.50000000     m   ;  Time  (for reference rigidity & particle) =   2.363879E-08 s 

************************************************************************************************************************************
      8  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      7)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   1.0000    -2.015    -0.001    -0.002     0.146       0.0000    1.0000   -2.015   -0.001   -0.010    0.146   3.003419E+02     1
               Time of flight (mus) :  1.09226441E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -2.015351E-02  -8.346074E-07        1        1    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00  -9.632620E-05   1.460834E-04        1        1    1.000      (Z,P)         1
   0.0000E+00   0.0000E+00   1.0000E+00   1.092264E-02   1.416766E+03        1        0    0.000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************

      9   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
     10  Keyword, label(s) :  SPNPRT      PRINT                 



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000        -0.001000   0.208381   0.978047   1.000000   4.500000  12.027635   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000          -0.001000  0.208381  0.978047  1.000000      2.5100   12.0276   12.0275    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -1.0001E-03 -1.0001E-03  2.0838E-01  2.0838E-01  9.7805E-01  9.7805E-01  1.0000E+00  1.0000E+00  1.00000E+00  2.50997E+00      1   1

************************************************************************************************************************************
     11  Keyword, label(s) :  REBELOTE                          


                                -----  REBELOTE  -----

     End of pass #        1 through the optical structure 

                     Total of          1 particles have been launched

     Multiple pass, 
          from element #     1 : OBJET     /label1=          /label2=           to REBELOTE /label1=          /label2=          
          ending at pass #       5 at element #    11 : REBELOTE  /label1=          /label2=          


     Parameter #   35 in element #    1 will be modified at each pass. 
     list of requested parameter values :
                   1   1.38727400E+00
                   2   2.13682967E+00
                   3   4.82611907E+00
                   4   1.10152413E+01

 Pgm rebel. At pass #    1/   5.  In element #    1,  parameter # 35  changed to    1.38727400E+00   (was    1.00000000E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
  LMNT  VAR  PARAM   MINIMUM     INITIAL         FINAL         MAXIMUM     STEP     NAME       LBL1     LBL2
     1    1     30   -3.00       -1.52       -1.522956209       3.00      1.218E-39 OBJET      *          *         
     1    2     31   -3.00      -3.381E-04  -3.3807340020E-04   3.00      1.218E-39 OBJET      *          *         
     1    3     32   -3.00      -4.954E-03  -4.9538386684E-03   3.00      1.218E-39 OBJET      *          *         
     1    4     33   -3.00       0.118       0.1176557231       3.00      1.218E-39 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     6      0.0000000E+00     1.0000E+00    1.2746671E-05   8.0848E-06   TOSCA      *          *          *   0 : 
   3   1   3     6      0.0000000E+00     1.0000E+00    3.6831641E-04   6.7502E-03   TOSCA      *          *          *   0 : 
   3   1   4     6      0.0000000E+00     1.0000E+00    1.4674377E-05   1.0715E-05   TOSCA      *          *          *   0 : 
   3   1   5     6      0.0000000E+00     1.0000E+00    1.4347029E-05   1.0242E-05   TOSCA      *          *          *   0 : 
   7   1   2     6      0.0000000E+00     1.0000E+02   -4.4649009E-01   9.9197E-01   TOSCA      *          *          *   1 :   2.0E+00/
   7   1   4     6      0.0000000E+00     1.0000E+02   -1.5873546E-02   1.2538E-03   TOSCA      *          *          *   1 :   2.0E+00/
 Fit reached penalty value   2.0097E-05

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      1)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   1.3873    -1.523     0.000    -0.005     0.118       0.0000    1.3873   -1.523    0.000   -0.005    0.118   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -1.522956E-02  -3.380734E-07        1        1    1.000      (Y,T)         2
   0.0000E+00   0.0000E+00   1.0000E+00  -4.953839E-05   1.176557E-04        1        1    1.000      (Z,P)         2
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   2.201779E+03        1        0    0.000      (t,K)         2

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      3  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      4  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              Same spin for all particles
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      5  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  3.872740E-01 -1.522939E+00 -3.380734E-04 -1.083662E-02  1.176557E-01  -5.0000000E+01  -1.74767E-03
TRAJ #1 SX, SY, SZ, |S| :  1    0.000000E+00   0.000000E+00   1.000000E+00   1.000000E+00

 Cumulative length of optical axis =    3.00000000     m   ;  Time  (for reference rigidity & particle) =   1.091021E-08 s 

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is  29 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 1.0

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.070
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.065
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.060
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.055
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.050
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.045
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.040
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.035
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.030
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.025
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.020
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.015
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.010
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.005
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.000
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.005
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.010
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.015
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.020
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.025
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.030
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.035
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.040
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.045
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.050
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.055
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.060
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.065
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.070
 Pgm toscac,  restored mesh coordinates for field map #   1,  name : Wsnk3D/table55_m.070


     Min/max fields seen in map   :  -2.674110E-01                 /   3.119660E-01
       @  X,  Y, Z :  -1.04       6.00      -7.00                  /   1.05      -5.50      -7.00    

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.674110E+00                      3.119660E+00
       @  X (cm),  Y (cm), Z (cm) :  -104.       600.      -700.   /   105.      -550.      -700.    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =  -2.000000E+02 cm 
                                               to    XF =   2.000000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   29 ; Step in Z = -0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1  1.3873    -1.523     0.000    -0.005     0.118          200.000    -1.523     0.000    -0.005     0.000            1


                CONDITIONS  DE  MAXWELL  (     4002.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -1.7382E-12      -5.6259E-14        2.2554E-15
                                      1.0599E-12       -9.1587E-13
                                     -1.5389E-13       -5.8796E-12
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    7.00000000     m ;  Time  (for ref. rigidity & particle) =   2.545716E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  3.872740E-01 -1.522934E+00 -7.063898E-04 -1.085058E-02  1.176414E-01   3.0017679E+02   1.04922E-02
TRAJ #1 SX, SY, SZ, |S| :  1   -7.375994E-04   1.961983E-01   9.805640E-01   1.000000E+00

 Cumulative length of optical axis =    6.50000000     m   ;  Time  (for reference rigidity & particle) =   2.363879E-08 s 

************************************************************************************************************************************
      8  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      7)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   1.3873    -1.523     0.000    -0.005     0.118       0.0000    1.3873   -1.523   -0.001   -0.011    0.118   3.001768E+02     1
               Time of flight (mus) :  1.04921727E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -1.522934E-02  -7.063898E-07        1        1    1.000      (Y,T)         2
   0.0000E+00   0.0000E+00   1.0000E+00  -1.085058E-04   1.176414E-04        1        1    1.000      (Z,P)         2
   0.0000E+00   0.0000E+00   1.0000E+00   1.049217E-02   2.201779E+03        1        0    0.000      (t,K)         2

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************

      9   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
     10  Keyword, label(s) :  SPNPRT      PRINT                 



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000        -0.000738   0.196198   0.980564   1.000000   6.000000  11.314812   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000          -0.000738  0.196198  0.980564  1.000000      3.3466   11.3148   11.3147    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -1.0001E-03 -7.3760E-04  1.9620E-01  2.0838E-01  9.7805E-01  9.8056E-01  1.0000E+00  1.0000E+00  1.38727E+00  3.34663E+00      1   1

************************************************************************************************************************************
     11  Keyword, label(s) :  REBELOTE                          


                                -----  REBELOTE  -----

     End of pass #        2 through the optical structure 

                     Total of          2 particles have been launched

 Pgm rebel. At pass #    2/   5.  In element #    1,  parameter # 35  changed to    2.13682967E+00   (was    1.38727400E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
  LMNT  VAR  PARAM   MINIMUM     INITIAL         FINAL         MAXIMUM     STEP     NAME       LBL1     LBL2
     1    1     30   -3.00       -1.13       -1.131290593       3.00      7.191E-35 OBJET      *          *         
     1    2     31   -3.00       2.910E-04   2.9097288651E-04   3.00      7.191E-35 OBJET      *          *         
     1    3     32   -3.00      -8.465E-03  -8.4654984757E-03   3.00      7.191E-35 OBJET      *          *         
     1    4     33   -3.00       8.288E-02   8.2877752955E-02   3.00      7.191E-35 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     6      0.0000000E+00     1.0000E+00    1.6406151E-05   7.7759E-06   TOSCA      *          *          *   0 : 
   3   1   3     6      0.0000000E+00     1.0000E+00    1.1206229E-03   3.6279E-02   TOSCA      *          *          *   0 : 
   3   1   4     6      0.0000000E+00     1.0000E+00    2.3110249E-05   1.5429E-05   TOSCA      *          *          *   0 : 
   3   1   5     6      0.0000000E+00     1.0000E+00    1.4192295E-05   5.8190E-06   TOSCA      *          *          *   0 : 
   7   1   2     6      0.0000000E+00     1.0000E+02   -5.7717611E-01   9.6240E-01   TOSCA      *          *          *   1 :   2.0E+00/
   7   1   4     6      0.0000000E+00     1.0000E+02   -2.1136380E-02   1.2906E-03   TOSCA      *          *          *   1 :   2.0E+00/
 Fit reached penalty value   3.4615E-05

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      1)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   2.1368    -1.131     0.000    -0.008     0.083       0.0000    2.1368   -1.131    0.000   -0.008    0.083   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -1.131291E-02   2.909729E-07        1        1    1.000      (Y,T)         3
   0.0000E+00   0.0000E+00   1.0000E+00  -8.465498E-05   8.287775E-05        1        1    1.000      (Z,P)         3
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   3.771805E+03        1        0    0.000      (t,K)         3

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      3  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      4  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              Same spin for all particles
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      5  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  1.136830E+00 -1.131305E+00  2.909729E-04 -1.260939E-02  8.287775E-02  -5.0000000E+01  -1.70193E-03
TRAJ #1 SX, SY, SZ, |S| :  1    0.000000E+00   0.000000E+00   1.000000E+00   1.000000E+00

 Cumulative length of optical axis =    3.00000000     m   ;  Time  (for reference rigidity & particle) =   1.091021E-08 s 

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is  29 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 1.0

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.070
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.065
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.060
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.055
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.050
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.045
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.040
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.035
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.030
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.025
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.020
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.015
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.010
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.005
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.000
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.005
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.010
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.015
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.020
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.025
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.030
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.035
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.040
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.045
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.050
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.055
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.060
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.065
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.070
 Pgm toscac,  restored mesh coordinates for field map #   1,  name : Wsnk3D/table55_m.070


     Min/max fields seen in map   :  -2.674110E-01                 /   3.119660E-01
       @  X,  Y, Z :  -1.04       6.00      -7.00                  /   1.05      -5.50      -7.00    

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.674110E+00                      3.119660E+00
       @  X (cm),  Y (cm), Z (cm) :  -104.       600.      -700.   /   105.      -550.      -700.    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =  -2.000000E+02 cm 
                                               to    XF =   2.000000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   29 ; Step in Z = -0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1  2.1368    -1.131     0.000    -0.008     0.083          200.000    -1.131     0.000    -0.008     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -8.4709E-13      -2.3751E-14        1.0627E-15
                                      7.4537E-13       -5.9519E-13
                                     -1.1184E-13       -3.8239E-12
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    7.00000000     m ;  Time  (for ref. rigidity & particle) =   2.545716E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  1.136830E+00 -1.131266E+00 -8.296500E-04 -1.263179E-02  8.286356E-02   3.0007428E+02   1.02141E-02
TRAJ #1 SX, SY, SZ, |S| :  1   -5.922127E-04   1.886199E-01   9.820500E-01   1.000000E+00

 Cumulative length of optical axis =    6.50000000     m   ;  Time  (for reference rigidity & particle) =   2.363879E-08 s 

************************************************************************************************************************************
      8  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      7)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   2.1368    -1.131     0.000    -0.008     0.083       0.0000    2.1368   -1.131   -0.001   -0.013    0.083   3.000743E+02     1
               Time of flight (mus) :  1.02141141E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -1.131266E-02  -8.296500E-07        1        1    1.000      (Y,T)         3
   0.0000E+00   0.0000E+00   1.0000E+00  -1.263179E-04   8.286356E-05        1        1    1.000      (Z,P)         3
   0.0000E+00   0.0000E+00   1.0000E+00   1.021411E-02   3.771805E+03        1        0    0.000      (t,K)         3

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************

      9   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
     10  Keyword, label(s) :  SPNPRT      PRINT                 



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000        -0.000592   0.188620   0.982050   1.000000   9.000000  10.872308   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000          -0.000592  0.188620  0.982050  1.000000      5.0199   10.8723   10.8723    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -1.0001E-03 -5.9221E-04  1.8862E-01  2.0838E-01  9.7805E-01  9.8205E-01  1.0000E+00  1.0000E+00  2.13683E+00  5.01995E+00      1   1

************************************************************************************************************************************
     11  Keyword, label(s) :  REBELOTE                          


                                -----  REBELOTE  -----

     End of pass #        3 through the optical structure 

                     Total of          3 particles have been launched

 Pgm rebel. At pass #    3/   5.  In element #    1,  parameter # 35  changed to    4.82611907E+00   (was    2.13682967E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
  LMNT  VAR  PARAM   MINIMUM     INITIAL         FINAL         MAXIMUM     STEP     NAME       LBL1     LBL2
     1    1     30   -3.00      -0.678      -0.6783291557       3.00      6.472E-34 OBJET      *          *         
     1    2     31   -3.00       3.313E-03   3.3129954003E-03   3.00      6.472E-34 OBJET      *          *         
     1    3     32   -3.00      -2.904E-02  -2.9041629732E-02   3.00      6.472E-34 OBJET      *          *         
     1    4     33   -3.00       3.730E-02   3.7296173821E-02   3.00      6.472E-34 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     6      0.0000000E+00     1.0000E+00    1.7546932E-05   4.0718E-06   TOSCA      *          *          *   0 : 
   3   1   3     6      0.0000000E+00     1.0000E+00    6.1636806E-03   5.0242E-01   TOSCA      *          *          *   0 : 
   3   1   4     6      0.0000000E+00     1.0000E+00    1.6450827E-05   3.5790E-06   TOSCA      *          *          *   0 : 
   3   1   5     6      0.0000000E+00     1.0000E+00    4.5084676E-06   2.6881E-07   TOSCA      *          *          *   0 : 
   7   1   2     6      0.0000000E+00     1.0000E+02   -6.1045022E-01   4.9282E-01   TOSCA      *          *          *   1 :   2.0E+00/
   7   1   4     6      0.0000000E+00     1.0000E+02   -5.9979939E-02   4.7577E-03   TOSCA      *          *          *   1 :   2.0E+00/
 Fit reached penalty value   7.5616E-05

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      1)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   4.8261    -0.678     0.003    -0.029     0.037       0.0000    4.8261   -0.678    0.003   -0.029    0.037   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -6.783292E-03   3.312995E-06        1        1    1.000      (Y,T)         4
   0.0000E+00   0.0000E+00   1.0000E+00  -2.904163E-04   3.729617E-05        1        1    1.000      (Z,P)         4
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   9.528565E+03        1        0    0.000      (t,K)         4

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      3  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      4  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              Same spin for all particles
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      5  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  3.826119E+00 -6.784948E-01  3.312995E-03 -3.090644E-02  3.729617E-02  -5.0000000E+01  -1.67456E-03
TRAJ #1 SX, SY, SZ, |S| :  1    0.000000E+00   0.000000E+00   1.000000E+00   1.000000E+00

 Cumulative length of optical axis =    3.00000000     m   ;  Time  (for reference rigidity & particle) =   1.091021E-08 s 

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is  29 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 1.0

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.070
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.065
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.060
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.055
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.050
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.045
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.040
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.035
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.030
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.025
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.020
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.015
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.010
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.005
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.000
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.005
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.010
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.015
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.020
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.025
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.030
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.035
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.040
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.045
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.050
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.055
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.060
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.065
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.070
 Pgm toscac,  restored mesh coordinates for field map #   1,  name : Wsnk3D/table55_m.070


     Min/max fields seen in map   :  -2.674110E-01                 /   3.119660E-01
       @  X,  Y, Z :  -1.04       6.00      -7.00                  /   1.05      -5.50      -7.00    

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.674110E+00                      3.119660E+00
       @  X (cm),  Y (cm), Z (cm) :  -104.       600.      -700.   /   105.      -550.      -700.    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =  -2.000000E+02 cm 
                                               to    XF =   2.000000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   29 ; Step in Z = -0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1  4.8261    -0.678     0.003    -0.029     0.037          200.000    -0.678     0.000    -0.029     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -2.2179E-13      -5.6051E-15       -1.7239E-15
                                      3.3947E-13       -2.6755E-13
                                     -5.3577E-14       -1.7217E-12
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    7.00000000     m ;  Time  (for ref. rigidity & particle) =   2.545716E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  3.826119E+00 -6.782042E-01 -2.850685E-03 -3.092266E-02  3.729167E-02   3.0001454E+02   1.00479E-02
TRAJ #1 SX, SY, SZ, |S| :  1   -4.837893E-04   1.841838E-01   9.828917E-01   1.000000E+00

 Cumulative length of optical axis =    6.50000000     m   ;  Time  (for reference rigidity & particle) =   2.363879E-08 s 

************************************************************************************************************************************
      8  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      7)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   4.8261    -0.678     0.003    -0.029     0.037       0.0000    4.8261   -0.678   -0.003   -0.031    0.037   3.000145E+02     1
               Time of flight (mus) :  1.00478603E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -6.782042E-03  -2.850685E-06        1        1    1.000      (Y,T)         4
   0.0000E+00   0.0000E+00   1.0000E+00  -3.092266E-04   3.729167E-05        1        1    1.000      (Z,P)         4
   0.0000E+00   0.0000E+00   1.0000E+00   1.004786E-02   9.528565E+03        1        0    0.000      (t,K)         4

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************

      9   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
     10  Keyword, label(s) :  SPNPRT      PRINT                 



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000        -0.000484   0.184184   0.982892   1.000000  20.000001  10.613586   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000          -0.000484  0.184184  0.982892  1.000000     11.1554   10.6136   10.6136    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -1.0001E-03 -4.8379E-04  1.8418E-01  2.0838E-01  9.7805E-01  9.8289E-01  1.0000E+00  1.0000E+00  4.82612E+00  1.11554E+01      1   1

************************************************************************************************************************************
     11  Keyword, label(s) :  REBELOTE                          


                                -----  REBELOTE  -----

     End of pass #        4 through the optical structure 

                     Total of          4 particles have been launched


      Next  pass  is  #     5 and  last  pass  through  the  optical  structure


 Pgm rebel. At pass #    4/   5.  In element #    1,  parameter # 35  changed to    1.10152413E+01   (was    4.82611907E+00)

     WRITE statements to zgoubi.res are re-established from now on.

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      1)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1  11.0152    -0.678     0.003    -0.029     0.037       0.0000   11.0152   -0.678    0.003   -0.029    0.037   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -6.783292E-03   3.312995E-06        1        1    1.000      (Y,T)         5
   0.0000E+00   0.0000E+00   1.0000E+00  -2.904163E-04   3.729617E-05        1        1    1.000      (Z,P)         5
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   2.287378E+04        1        0    0.000      (t,K)         5

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      3  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      4  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              Same spin for all particles
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      5  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  1.001524E+01 -6.784948E-01  3.312995E-03 -3.090644E-02  3.729617E-02  -5.0000000E+01  -1.66912E-03
TRAJ #1 SX, SY, SZ, |S| :  1    0.000000E+00   0.000000E+00   1.000000E+00   1.000000E+00

 Cumulative length of optical axis =  -0.500000000     m   ;  Time  (for reference rigidity & particle) =   2.182042E-08 s 

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is  29 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 1.0

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.070
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.065
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.060
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.055
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.050
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.045
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.040
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.035
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.030
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.025
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.020
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.015
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.010
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.005
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.000
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.005
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.010
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.015
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.020
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.025
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.030
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.035
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.040
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.045
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.050
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.055
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.060
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.065
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.070
 Pgm toscac,  restored mesh coordinates for field map #   1,  name : Wsnk3D/table55_m.070


     Min/max fields seen in map   :  -2.674110E-01                 /   3.119660E-01
       @  X,  Y, Z :  -1.04       6.00      -7.00                  /   1.05      -5.50      -7.00    

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.674110E+00                      3.119660E+00
       @  X (cm),  Y (cm), Z (cm) :  -104.       600.      -700.   /   105.      -550.      -700.    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =  -2.000000E+02 cm 
                                               to    XF =   2.000000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   29 ; Step in Z = -0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1 11.0152    -0.678     0.003    -0.029     0.037          200.000    -0.678     0.000    -0.022     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -9.7299E-14      -2.4558E-15       -7.3780E-16
                                      1.4873E-13       -1.1757E-13
                                     -2.3474E-14       -7.5667E-13
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    3.50000000     m ;  Time  (for ref. rigidity & particle) =   3.636737E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  1.001524E+01 -6.778349E-01 -1.439014E-04 -2.365568E-02  3.717610E-02   3.0000279E+02   1.00148E-02
TRAJ #1 SX, SY, SZ, |S| :  1   -4.385460E-04   1.833379E-01   9.830499E-01   1.000000E+00

 Cumulative length of optical axis =    3.00000000     m   ;  Time  (for reference rigidity & particle) =   3.454900E-08 s 

************************************************************************************************************************************
      8  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      7)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1  11.0152    -0.678     0.003    -0.029     0.037       0.0000   11.0152   -0.678    0.000   -0.024    0.037   3.000028E+02     1
               Time of flight (mus) :  1.00147935E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -6.778349E-03  -1.439014E-07        1        1    1.000      (Y,T)         5
   0.0000E+00   0.0000E+00   1.0000E+00  -2.365568E-04   3.717610E-05        1        1    1.000      (Z,P)         5
   0.0000E+00   0.0000E+00   1.0000E+00   1.001479E-02   2.287378E+04        1        0    0.000      (t,K)         5

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      9  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1


                    Saving new version of zgoubi.dat to zgoubi.FIT.out.dat, with variables updated.

                    Updated version of zgoubi.dat saved in  zgoubi.FIT.out.dat

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      1)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1  11.0152    -0.678     0.003    -0.029     0.037       0.0000   11.0152   -0.678    0.003   -0.029    0.037   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -6.783292E-03   3.312995E-06        1        1    1.000      (Y,T)         5
   0.0000E+00   0.0000E+00   1.0000E+00  -2.904163E-04   3.729617E-05        1        1    1.000      (Z,P)         5
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   2.287378E+04        1        0    0.000      (t,K)         5

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      3  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      4  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              Same spin for all particles
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      5  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  1.001524E+01 -6.784948E-01  3.312995E-03 -3.090644E-02  3.729617E-02  -5.0000000E+01  -1.66912E-03
TRAJ #1 SX, SY, SZ, |S| :  1    0.000000E+00   0.000000E+00   1.000000E+00   1.000000E+00

 Cumulative length of optical axis =  -0.500000000     m   ;  Time  (for reference rigidity & particle) =  -1.818368E-09 s 

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is  29 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 1.0

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.070
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.065
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.060
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.055
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.050
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.045
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.040
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.035
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.030
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.025
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.020
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.015
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.010
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.005
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.000
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.005
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.010
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.015
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.020
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.025
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.030
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.035
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.040
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.045
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.050
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.055
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.060
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.065
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.070
 Pgm toscac,  restored mesh coordinates for field map #   1,  name : Wsnk3D/table55_m.070


     Min/max fields seen in map   :  -2.674110E-01                 /   3.119660E-01
       @  X,  Y, Z :  -1.04       6.00      -7.00                  /   1.05      -5.50      -7.00    

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.674110E+00                      3.119660E+00
       @  X (cm),  Y (cm), Z (cm) :  -104.       600.      -700.   /   105.      -550.      -700.    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =  -2.000000E+02 cm 
                                               to    XF =   2.000000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   29 ; Step in Z = -0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1 11.0152    -0.678     0.003    -0.029     0.037          200.000    -0.678     0.000    -0.022     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -9.7299E-14      -2.4558E-15       -7.3780E-16
                                      1.4873E-13       -1.1757E-13
                                     -2.3474E-14       -7.5667E-13
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    3.50000000     m ;  Time  (for ref. rigidity & particle) =   1.272858E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  1.001524E+01 -6.778349E-01 -1.439014E-04 -2.365568E-02  3.717610E-02   3.0000279E+02   1.00148E-02
TRAJ #1 SX, SY, SZ, |S| :  1   -4.385460E-04   1.833379E-01   9.830499E-01   1.000000E+00

 Cumulative length of optical axis =    3.00000000     m   ;  Time  (for reference rigidity & particle) =   1.091021E-08 s 

************************************************************************************************************************************
      8  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      7)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1  11.0152    -0.678     0.003    -0.029     0.037       0.0000   11.0152   -0.678    0.000   -0.024    0.037   3.000028E+02     1
               Time of flight (mus) :  1.00147935E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -6.778349E-03  -1.439014E-07        1        1    1.000      (Y,T)         5
   0.0000E+00   0.0000E+00   1.0000E+00  -2.365568E-04   3.717610E-05        1        1    1.000      (Z,P)         5
   0.0000E+00   0.0000E+00   1.0000E+00   1.001479E-02   2.287378E+04        1        0    0.000      (t,K)         5

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      9  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            1 ,   ok.
           variable #            2       IP =           31 ,   ok.
           variable #            3       IR =            1 ,   ok.
           variable #            3       IP =           32 ,   ok.
           variable #            4       IR =            1 ,   ok.
           variable #            4       IP =           33 ,   ok.
           constraint #            1       IR =            6 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            6 ,   ok.
           constraint #            2       I  =            1 ,   ok.
           constraint #            3       IR =            6 ,   ok.
           constraint #            3       I  =            1 ,   ok.
           constraint #            4       IR =            6 ,   ok.
           constraint #            4       I  =            1 ,   ok.
           constraint #            5       IR =            6 ,   ok.
           constraint #            6       IR =            6 ,   ok.

                    FIT  variables  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save [FileName]' command

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
  LMNT  VAR  PARAM   MINIMUM     INITIAL         FINAL         MAXIMUM     STEP     NAME       LBL1     LBL2
     1    1     30   -3.00      -0.204      -0.2037655071       3.00      1.218E-39 OBJET      *          *         
     1    2     31   -3.00       2.758E-03   2.7576373690E-03   3.00      1.218E-39 OBJET      *          *         
     1    3     32   -3.00      -1.474E-03  -1.4740498503E-03   3.00      1.218E-39 OBJET      *          *         
     1    4     33   -3.00       1.615E-02   1.6153506996E-02   3.00      1.218E-39 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-10)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     6      0.0000000E+00     1.0000E+00    2.3175413E-06   2.1264E-07   TOSCA      *          *          *   0 : 
   3   1   3     6      0.0000000E+00     1.0000E+00    4.9567868E-03   9.7271E-01   TOSCA      *          *          *   0 : 
   3   1   4     6      0.0000000E+00     1.0000E+00    2.8905477E-06   3.3078E-07   TOSCA      *          *          *   0 : 
   3   1   5     6      0.0000000E+00     1.0000E+00    2.0342543E-04   1.6383E-03   TOSCA      *          *          *   0 : 
   7   1   2     6      0.0000000E+00     1.0000E+02   -8.0399559E-02   2.5591E-02   TOSCA      *          *          *   1 :   2.0E+00/
   7   1   4     6      0.0000000E+00     1.0000E+02   -3.7409358E-03   5.5405E-05   TOSCA      *          *          *   1 :   2.0E+00/
 Fit reached penalty value   2.5259E-05

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      1)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1  11.0152    -0.204     0.003    -0.001     0.016       0.0000   11.0152   -0.204    0.003   -0.001    0.016   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -2.037655E-03   2.757637E-06        1        1    1.000      (Y,T)         5
   0.0000E+00   0.0000E+00   1.0000E+00  -1.474050E-05   1.615351E-05        1        1    1.000      (Z,P)         5
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   2.287378E+04        1        0    0.000      (t,K)         5

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      3  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      4  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              Same spin for all particles
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      5  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  1.001524E+01 -2.039034E-01  2.757637E-03 -2.281725E-03  1.615351E-02  -5.0000000E+01  -1.66912E-03
TRAJ #1 SX, SY, SZ, |S| :  1    0.000000E+00   0.000000E+00   1.000000E+00   1.000000E+00

 Cumulative length of optical axis =    3.00000000     m   ;  Time  (for reference rigidity & particle) =   1.091021E-08 s 

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is  29 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 1.0

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.070
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.065
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.060
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.055
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.050
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.045
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.040
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.035
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.030
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.025
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.020
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.015
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.010
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.005
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_p.000
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.005
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.010
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.015
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.020
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.025
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.030
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.035
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.040
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.045
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.050
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.055
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.060
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.065
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           Wsnk3D/table55_m.070
 Pgm toscac,  restored mesh coordinates for field map #   1,  name : Wsnk3D/table55_m.070


     Min/max fields seen in map   :  -2.674110E-01                 /   3.119660E-01
       @  X,  Y, Z :  -1.04       6.00      -7.00                  /   1.05      -5.50      -7.00    

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.674110E+00                      3.119660E+00
       @  X (cm),  Y (cm), Z (cm) :  -104.       600.      -700.   /   105.      -550.      -700.    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =  -2.000000E+02 cm 
                                               to    XF =   2.000000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   29 ; Step in Z = -0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1 11.0152    -0.204     0.003    -0.001     0.016          200.000    -0.204     0.000    -0.001     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -2.6974E-14       5.2668E-17        2.1640E-16
                                      1.2915E-13       -1.1764E-13
                                     -1.9586E-14       -7.5696E-13
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    7.00000000     m ;  Time  (for ref. rigidity & particle) =   2.545716E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  DRIFT       DRIF        WSNK      


                              Drift,  length =   -50.00000  cm

TRAJ #1 IEX,D,Y,T,Z,P,S,time :  1  1.001524E+01 -2.036579E-01 -2.199149E-03 -2.274444E-03  1.595008E-02   3.0000279E+02   1.00148E-02
TRAJ #1 SX, SY, SZ, |S| :  1   -3.980044E-04   1.832817E-01   9.830604E-01   1.000000E+00

 Cumulative length of optical axis =    6.50000000     m   ;  Time  (for reference rigidity & particle) =   2.363879E-08 s 

************************************************************************************************************************************
      8  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU
                                           (follows element #      7)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1  11.0152    -0.204     0.003    -0.001     0.016       0.0000   11.0152   -0.204   -0.002   -0.002    0.016   3.000028E+02     1
               Time of flight (mus) :  1.00147935E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00  -2.036579E-03  -2.199149E-06        1        1    1.000      (Y,T)         5
   0.0000E+00   0.0000E+00   1.0000E+00  -2.274444E-05   1.595008E-05        1        1    1.000      (Z,P)         5
   0.0000E+00   0.0000E+00   1.0000E+00   1.001479E-02   2.287378E+04        1        0    0.000      (t,K)         5

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************

      9   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
     10  Keyword, label(s) :  SPNPRT      PRINT                 



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000        -0.000398   0.183282   0.983060   1.000000  45.500002  10.560991   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000          -0.000398  0.183282  0.983060  1.000000     25.3786   10.5610   10.5610    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -1.0001E-03 -3.9800E-04  1.8328E-01  2.0838E-01  9.7805E-01  9.8306E-01  1.0000E+00  1.0000E+00  1.10152E+01  2.53786E+01      1   1

************************************************************************************************************************************
     11  Keyword, label(s) :  REBELOTE                          


                         ****  End  of  'REBELOTE'  procedure  ****

      There  has  been          5  passes  through  the  optical  structure 

                     Total of          5 particles have been launched

************************************************************************************************************************************
     12  Keyword, label(s) :  END                               


************************************************************************************************************************************

           MAIN PROGRAM : Execution ended upon key  END       

************************************************************************************************************************************

                    Saving new version of zgoubi.dat to zgoubi.FIT.out.dat, with variables updated.

                    Updated version of zgoubi.dat saved in  zgoubi.FIT.out.dat
   
            ZGOUBI RUN COMPLETED. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  12-05-0016,  at  03:22:33 
  JOB  ENDED  ON    12-05-0016,  AT  03:24:03 

   CPU time, total :     89.3015800000000     
