Cyclotron, classical. Using DIPOLE, first, for reference, and then BEND.
 'OBJET'                                                                                                      1
64.62444403717985
2
5 1
12.9248888074 0. 0. 0. 0. 1.                          'm'   ! 200keV
28.9070891209 0. 0. 0. 0. 2.23654451125 'm'   ! 1 MeV
50.           0. 0. 0. 0. 3.86850523397 'o'   ! at RM
64.7070336799 0. 0. 0. 0. 5.0063899693  'M'   ! 5 MeV
91.6310723304 0. 0. 0. 0. 7.08950565809  'M'  ! 10 MeV
1 1 1 1 1
 
 'DIPOLE'                                                                                                     2
0
60. 50.
30.  5. 0. 0. 0.
0.  0.                                        ! EFB 1  hard-edge
4  .1455   2.2670  -.6395  1.1558  0. 0.  0.
30. 0.  1.E6  -1.E6  1.E6  1.E6
0.   0.                                       ! EFB 2
4  .1455   2.2670  -.6395  1.1558  0. 0.  0.
-30. 0.  1.E6  -1.E6  1.E6  1.E6
0. 0.                                         ! EFB 3
0  0.      0.      0.      0.      0. 0.  0.
0. 0.  1.E6  -1.E6  1.E6  1.E6 0.
4   10.
.1                 ! step size
2  0. 0. 0. 0.       ! could also be, e.g., 2 50. 0. 50. 0. with Y0 amended accordingly in OBJET
 'FAISCEAU'                                                                                                   3
 'FAISTORE'                                                                                                   4
zgoubi.fai
1
 
 
! Now using BEND
 
 'OBJET'                                                                                                      5
250.   ! Particle at 50cm (=RM in DIPOLE) taken as reference -> BORO =64.624444*3.86850523397 =250.
2      ! Thus, 50cm is reference radius in building BEND data
5 1
-42.8106508546 523.598775598 0. 0. 0. 0.258497776148 'm'   ! 200keV, Y =Brho/B-50.
-24.3559955481 523.598775598 0. 0. 0. 0.578141782415 'm'   ! 1 MeV,  Y =Brho/B-50.
0.             523.598775598 0. 0. 0. 1.             'o'   ! D=1 for reference particle, its Brho=250.
16.9822197081  523.598775598 0. 0. 0. 1.2941406736   'M'   ! 5 MeV,  Y =Brho/B-50.
41.6310723304  523.598775598 0. 0. 0. 1.83262144661  'M'   ! 10 MeV
1 1 1 1 1
 
! Starting coordinates for BEND, compared to DIPOLE case:
! (12.9248888074-50)/cos(theta/2) 523.598775598 0. 0. 0. 1./3.86850523397            'm'   ! 200keV
! (28.9070891209-50)/cos(theta/2) 523.598775598 0. 0. 0. 2.23654451125/3.86850523397 'm'   ! 1 MeV
! (50.-50)/cos(theta/2)           523.598775598 0. 0. 0. 3.86850523397/3.86850523397 'o'   ! at RM
! (64.7070336799-50)/cos(theta/2) 523.598775598 0. 0. 0. 5.0063899693/3.86850523397  'M'   ! 5 MeV
! (91.6310723304-50)/cos(theta/2) 523.598775598 0. 0. 0. 7.08950565809  'M'  ! 10 MeV
 
 'BEND'            ! Take field = B0 = 5 kG, as DIPOLE,                                                       6
0                  ! take ref particle #3, its rho=64.624444*3.86850523397/B0 =250/B0 =50
50.  0.  5.        ! thus, straight_L =2*50*sin(30deg)=50 (arc_L =rho*theta=50*pi/3=52.3598775598)
0.  0.   0.
4  .1455   2.2670  -.6395  1.1558  0. 0.  0.
0.  0.  0.
4  .1455   2.2670  -.6395  1.1558  0. 0.  0.
.2
2  0. 0. 0.
 'FAISCEAU'                                                                                                   7
 
! Particles do not end up here with the proper time of flight, this requires CHANGREF
! upstream and downstream of BEND, or equivalently KPOS=3 in BEND
 'FAISTORE'                                                                                                   8
zgoubi.fai
1
 
 
! Now, introduce KPOS=3 in BEND
 
 'OBJET'                                                                                                      9
250.   ! Particle at 50cm (=RM in DIPOLE) taken as reference -> BORO =64.624444*3.86850523397 =250.
2      ! Thus, 50cm is reference radius in building BEND data
5 1
-37.0751111926 0. 0. 0. 0. 0.258497776148 'm'   ! 200keV, Y =Brho/B-50.
-21.0929108791 0. 0. 0. 0. 0.578141782415 'm'   ! 1 MeV,  Y =Brho/B-50.
0.             0. 0. 0. 0. 1.             'o'   ! D=1 for reference particle, its Brho=250.
14.7070336799  0. 0. 0. 0. 1.2941406736   'M'   ! 5 MeV,  Y =Brho/B-50.
41.6310723304  523.598775598 0. 0. 0. 1.83262144661  'M'   ! 10 MeV
1 1 1 1 1
 
 'BEND'            ! Take field = B0 = 5 kG, as DIPOLE,                                                      10
2                  ! take ref particle #3, its rho=64.624444*3.86850523397/B0 =250/B0 =50
50.  0.  5.        ! thus, straight_L =2*50*sin(30deg)=50 (arc_L =rho*theta=50*pi/3=52.3598775598)
20.  0.   0.
4  .1455   2.2670  -.6395  1.1558  0. 0.  0.
20.  0.  0.
4  .1455   2.2670  -.6395  1.1558  0. 0.  0.
.2
3  0. 0. 0.
 'FAISCEAU'                                                                                                  11
! Arc lengths now should be : 12.9248888074*pi/3=13.5349119086, 28.9070891209*pi/3=30.2714329396,
! 50.*pi/3=52.3598775598 and 64.7070336799*pi/3=67.7610472148
 'FAISTORE'                                                                                                  12
zgoubi.fai
1
 
! Now, introduce short fringe fall-off (thus reduce integration step size)
 
 'OBJET'                                                                                                     13
250.   ! Particle at 50cm (=RM in DIPOLE) taken as reference -> BORO =64.624444*3.86850523397 =250.
2      ! Thus, 50cm is reference radius in building BEND data
5  1
-37.0751111926 0. 0. 0. 0. 0.258497776148 'm'   ! 200keV, Y =Brho/B-50.
-21.0929108791 0. 0. 0. 0. 0.578141782415 'm'   ! 1 MeV,  Y =Brho/B-50.
0.             0. 0. 0. 0. 1.             'o'   ! D=1 for reference particle, its Brho=250.
14.7070336799  0. 0. 0. 0. 1.2941406736   'M'   ! 5 MeV,  Y =Brho/B-50.
41.6310723304  523.598775598 0. 0. 0. 1.83262144661  'M'   ! 10 MeV
1 1 1 1 1
 
 'BEND'            ! Take field = B0 = 5 kG, as DIPOLE,                                                      14
2                  ! take ref particle #3, its rho=64.624444*3.86850523397/B0 =250/B0 =50
50.  0.  5.        ! thus, straight_L =2*50*sin(30deg)=50 (arc_L =rho*theta=50*pi/3=52.3598775598)
4.  1.   0.
4  .1455   2.2670  -.6395  1.1558  0. 0.  0.
4.  1.  0.
4  .1455   2.2670  -.6395  1.1558  0. 0.  0.
.01
3  0. 0. 0.
 'FAISCEAU'                                                                                                  15
 'FAISTORE'                                                                                                  16
zgoubi.fai
1
 
 'END'                                                                                                       17

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       5 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  DIPOLE                                                

                    Dipole  magnet

           ANGLES : A.TOTAL =   60.00     degrees     A.CENTRAL =   30.00     degrees
           RM =   50.00     cm
           HNORM =   5.000     kGauss     COEF.N =   0.000         COEF.B =   0.000         COEF.G=   0.000    

     Entrance  EFB 
          Fringe  field  : LAMBDA =   0.00 CM     QSI=   0.00
           COEFFICIENTS :  4   0.14550   2.26700  -0.63950   1.15580   0.00000   0.00000
           Shift  of  EFB  =    0.000     CM

          OMEGA =  30.00 deg.     Wedge  angle  =   0.00 deg.
           Radius 1 =  1.00E+06 CM
           Straight  segment 1 = -1.00E+06 CM
           Straight  segment 2 =  1.00E+06 CM
           Radius 2 =  1.00E+06 CM

     Exit  EFB     
          Fringe  field  : LAMBDA =   0.00 CM     QSI=   0.00
           COEFFICIENTS :  4   0.14550   2.26700  -0.63950   1.15580   0.00000   0.00000
           Shift  of  EFB  =    0.000     CM

          OMEGA = -30.00 deg.     Wedge  angle  =   0.00 deg.
           Radius 1 =  1.00E+06 CM
           Straight  segment 1 = -1.00E+06 CM
           Straight  segment 2 =  1.00E+06 CM
           Radius 2 =  1.00E+06 CM

     Lateral  EFB  
          Fringe  field  : LAMBDA =   0.00 CM     QSI=   0.00
           COEFFICIENTS :  0   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
           Shift  of  EFB  =    0.000     CM

                     Face centred on direction ACENT+OMEGA, A    0.000     CM

          OMEGA =   0.00 deg.     Wedge  angle  =   0.00 deg.
           Radius 1 =  1.00E+06 CM
           Straight  segment 1 = -1.00E+06 CM
           Straight  segment 2 =  1.00E+06 CM
           Radius 2 =  1.00E+06 CM

                     Interpolation  option : 4
                    25-point  interpolation, size of flying mesh :   STEP /   10.0

                    Integration step :  0.1000     cm   (i.e.,   2.0000E-03 rad  at mean radius RM =    50.00    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            1.047    12.925    -0.000     0.000     0.000            1
  A    1  2.2365    28.907     0.000     0.000     0.000            1.047    28.907    -0.000     0.000     0.000            2
  A    1  3.8685    50.000     0.000     0.000     0.000            1.047    50.000     0.000     0.000     0.000            3
  A    1  5.0064    64.707     0.000     0.000     0.000            1.047    64.707    -0.000     0.000     0.000            4
  A    1  7.0895    91.631     0.000     0.000     0.000            1.047    91.631     0.000     0.000     0.000            5

 Cumulative length of optical axis =   0.523598776     m ;  Time  (for ref. rigidity & particle) =   8.460221E-08 s 

************************************************************************************************************************************
      3  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      2)
                                                  5 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    0.0000   12.925   -0.000    0.000    0.000   1.353491E+01     1
m  1   2.2365    28.907     0.000     0.000     0.000       0.0000    1.2365   28.907   -0.000    0.000    0.000   3.027143E+01     2
o  1   3.8685    50.000     0.000     0.000     0.000       0.0000    2.8685   50.000    0.000    0.000    0.000   5.235988E+01     3
M  1   5.0064    64.707     0.000     0.000     0.000       0.0000    4.0064   64.707   -0.000    0.000    0.000   6.776105E+01     4
M  1   7.0895    91.631     0.000     0.000     0.000       0.0000    6.0895   91.631    0.000    0.000    0.000   9.595583E+01     5


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   3.0712E-11  -1.0012E+00   7.7091E+09   4.963402E-01  -2.714865E-11        5        4   0.8000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        5        5    1.000      (Z,P)         1
   1.4629E-09  -8.4742E+07   1.9748E+03   1.733753E-03   7.439952E+01        5        4   0.8000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   2.745268E-01
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   5.039191E-11

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   9.589427E-04
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   4.115053E+01


  Beam  sigma  matrix : 

   7.536496E-02   9.788070E-12   0.000000E+00   0.000000E+00
   9.788070E-12   2.539345E-21   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     9.776057E-12    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISTORE                                              


                OPEN FILE zgoubi.fai                                                                      
                FOR PRINTING COORDINATES 

               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      5  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =        250.000 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       5 POINTS 



************************************************************************************************************************************
      6  Keyword, label(s) :  BEND                                                  


      +++++        BEND  : 

                Length    =   5.000000E+01 cm
                Arc length    =   5.235988E+01 cm
                Deviation    =   6.000000E+01 deg.,    1.047198E+00 rad
                GAP   =   0.000000E+00 cm
                Gradient   =   0.000000E+00 kG/cm
                Grad-prime   =   0.000000E+00 kG/cm^2

                Field  =  5.0000000E+00  kG   (i.e.,   5.0000000E+00 * SCAL)
                Reference curvature radius (Brho/B) =   5.0000000E+01 cm
                Skew  angle  =   0.000000E+00  rad

               Entrance  face  
                DX =      0.000    LAMBDA =      0.000
                Wedge  angle  =  0.000000 RD

               Exit      face  
                DX =      2.887    LAMBDA =      0.000
                Wedge  angle  =  0.000000 RD

  ***  Warning : entrance sharp edge entails vertical wedge focusing approximated with first order kick, FINT values entr/exit :    0.000    

  ***  Warning : exit sharp edge entails vertical wedge focusing approximated with first order kick, FINT values entr/exit :    0.000    

                    Integration step :  0.2000     cm

  A    1  0.2585   -42.811   523.599     0.000     0.000           52.887   -42.811    -0.524     0.000     0.000            1
  A    1  0.5781   -24.356   523.599     0.000     0.000           52.887   -24.356    -0.524     0.000     0.000            2
  A    1  1.0000     0.000   523.599     0.000     0.000           52.887    -0.000    -0.524     0.000     0.000            3
  A    1  1.2941    16.982   523.599     0.000     0.000           52.887    16.467    -0.446     0.000     0.000            4
  A    1  1.8326    41.631   523.599     0.000     0.000           52.887    40.817    -0.277     0.000     0.000            5

     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    1.04719755     m ;  Time  (for ref. rigidity & particle) =   1.065367E-07 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      6)
                                                  5 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   0.2585   -42.811   523.599     0.000     0.000       0.0000   -0.7415  -42.811 -523.599    0.000    0.000   5.634556E+01     1
m  1   0.5781   -24.356   523.599     0.000     0.000       0.0000   -0.4219  -24.356 -523.599    0.000    0.000   5.462743E+01     2
o  1   1.0000     0.000   523.599     0.000     0.000       0.0000    0.0000   -0.000 -523.599    0.000    0.000   5.235988E+01     3
M  1   1.2941    16.982   523.599     0.000     0.000       0.0000    0.2941   16.467 -445.567    0.000    0.000   5.102163E+01     4
M  1   1.8326    41.631   523.599     0.000     0.000       0.0000    0.8326   40.817 -277.449    0.000    0.000   4.958375E+01     5


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   4.6589E-02  -1.6147E+00   5.8580E+00  -1.976457E-02  -4.587623E-01        5        4   0.8000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        5        5    1.000      (Z,P)         1
   2.0071E-03  -2.6771E+01   2.7076E-04   1.732276E-03   7.439952E+01        5        4   0.8000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   2.947425E-01
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   9.556159E-02

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   4.159079E-04
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   4.115053E+01


  Beam  sigma  matrix : 

   8.687317E-02   2.394581E-02   0.000000E+00   0.000000E+00
   2.394581E-02   9.132018E-03   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     1.482989E-02    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      8  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      9  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =        250.000 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       5 POINTS 



************************************************************************************************************************************
     10  Keyword, label(s) :  BEND                                                  


                OPEN FILE zgoubi.plt                                                                      
                FOR PRINTING TRAJECTORIES


      +++++        BEND  : 

                Length    =   5.000000E+01 cm
                Arc length    =   5.235988E+01 cm
                Deviation    =   6.000000E+01 deg.,    1.047198E+00 rad
                GAP   =   0.000000E+00 cm
                Gradient   =   0.000000E+00 kG/cm
                Grad-prime   =   0.000000E+00 kG/cm^2

                Field  =  5.0000000E+00  kG   (i.e.,   5.0000000E+00 * SCAL)
                Reference curvature radius (Brho/B) =   5.0000000E+01 cm
                Skew  angle  =   0.000000E+00  rad

               Entrance  face  
                DX =     20.000    LAMBDA =      0.000
                Wedge  angle  =  0.000000 RD

               Exit      face  
                DX =      2.887    LAMBDA =      0.000
                Wedge  angle  =  0.000000 RD

  ***  Warning : entrance sharp edge entails vertical wedge focusing approximated with first order kick, FINT values entr/exit :    20.00    

  ***  Warning : exit sharp edge entails vertical wedge focusing approximated with first order kick, FINT values entr/exit :    20.00    

     KPOS =  3 :  automatic positioning of element, 
        XCE, YCE, ALE =   0.000000000       0.000000000     -0.5235987756     cm/cm/rad

                    Integration step :  0.2000     cm

  A    1  0.2585   -37.075     0.000     0.000     0.000           72.887   -42.811    -0.524     0.000     0.000            1
  A    1  0.5781   -21.093     0.000     0.000     0.000           72.887   -24.356    -0.524     0.000     0.000            2
  A    1  1.0000     0.000     0.000     0.000     0.000           72.887    -0.000    -0.524     0.000     0.000            3
  A    1  1.2941    14.707     0.000     0.000     0.000           72.887    16.467    -0.446     0.000     0.000            4
  A    1  1.8326    41.631   523.599     0.000     0.000           72.887    81.516     0.062     0.000     0.000            5


                CONDITIONS  DE  MAXWELL  (     1248.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  3.  Automatic  positionning  of  element.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle = -0.523599     RAD


 Cumulative length of optical axis =    1.57079633     m ;  Time  (for ref. rigidity & particle) =   1.284712E-07 s 

************************************************************************************************************************************
     11  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #     10)
                                                  5 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   0.2585   -37.075     0.000     0.000     0.000       0.0000   -0.7415  -37.075   -0.000    0.000    0.000   1.353491E+01     1
m  1   0.5781   -21.093     0.000     0.000     0.000       0.0000   -0.4219  -21.093   -0.000    0.000    0.000   3.027143E+01     2
o  1   1.0000     0.000     0.000     0.000     0.000       0.0000    0.0000   -0.000   -0.000    0.000    0.000   5.235988E+01     3
M  1   1.2941    14.707     0.000     0.000     0.000       0.0000    0.2941   14.905   78.032    0.000    0.000   6.777132E+01     4
M  1   1.8326    41.631   523.599     0.000     0.000       0.0000    0.8326   97.610  585.326    0.000    0.000   1.363055E+02     5


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   9.9773E-02  -3.2171E+00   6.9137E+00   1.086927E-01   1.326716E-01        5        4   0.8000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        5        5    1.000      (Z,P)         1
   4.5586E-02  -5.4067E+00   2.5906E-04   1.230223E-03   7.439952E+01        5        4   0.8000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   4.685850E-01
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   2.283361E-01

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   1.938840E-03
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   4.115053E+01


  Beam  sigma  matrix : 

   2.195719E-01   1.021728E-01   0.000000E+00   0.000000E+00
   1.021728E-01   5.213737E-02   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     3.175879E-02    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
     12  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
     13  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =        250.000 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       5 POINTS 



************************************************************************************************************************************
     14  Keyword, label(s) :  BEND                                                  


     zgoubi.plt                                                                      
      already open...

      +++++        BEND  : 

                Length    =   5.000000E+01 cm
                Arc length    =   5.235988E+01 cm
                Deviation    =   6.000000E+01 deg.,    1.047198E+00 rad
                GAP   =   0.000000E+00 cm
                Gradient   =   0.000000E+00 kG/cm
                Grad-prime   =   0.000000E+00 kG/cm^2

                Field  =  5.0000000E+00  kG   (i.e.,   5.0000000E+00 * SCAL)
                Reference curvature radius (Brho/B) =   5.0000000E+01 cm
                Skew  angle  =   0.000000E+00  rad

               Entrance  face  
                DX =      4.000    LAMBDA =      1.000
                Wedge  angle  =  0.000000 RD
                Fringe  field  coefficients :
                  0.14550  2.26700 -0.63950  1.15580  0.00000  0.00000

               Exit      face  
                DX =      4.000    LAMBDA =      1.000
                Wedge  angle  =  0.000000 RD
                Fringe  field  coefficients :
                  0.14550  2.26700 -0.63950  1.15580  0.00000  0.00000

     KPOS =  3 :  automatic positioning of element, 
        XCE, YCE, ALE =   0.000000000       0.000000000     -0.5235987756     cm/cm/rad

                    Integration step :  1.0000E-02 cm

  A    1  0.2585   -37.075     0.000     0.000     0.000           58.000   -42.786    -0.523     0.000     0.000            1
  A    1  0.5781   -21.093     0.000     0.000     0.000           58.000   -24.352    -0.523     0.000     0.000            2
  A    1  1.0000     0.000     0.000     0.000     0.000           58.000     0.001    -0.524     0.000     0.000            3
  A    1  1.2941    14.707     0.000     0.000     0.000           58.000    19.769    -0.408     0.000     0.000            4
  A    1  1.8326    41.631   523.599     0.000     0.000           58.000   107.513     0.235     0.000     0.000            5


                CONDITIONS  DE  MAXWELL  (    32569.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  3.  Automatic  positionning  of  element.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle = -0.523599     RAD


 Cumulative length of optical axis =    2.09439510     m ;  Time  (for ref. rigidity & particle) =   1.504056E-07 s 

************************************************************************************************************************************
     15  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #     14)
                                                  5 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   0.2585   -37.075     0.000     0.000     0.000       0.0000   -0.7415  -37.070    0.736    0.000    0.000   1.352354E+01     1
m  1   0.5781   -21.093     0.000     0.000     0.000       0.0000   -0.4219  -21.091    0.115    0.000    0.000   3.026595E+01     2
o  1   1.0000     0.000     0.000     0.000     0.000       0.0000    0.0000    0.001    0.024    0.000    0.000   5.235645E+01     3
M  1   1.2941    14.707     0.000     0.000     0.000       0.0000    0.2941   18.273  116.063    0.000    0.000   6.971799E+01     4
M  1   1.8326    41.631   523.599     0.000     0.000       0.0000    0.8326  144.080  758.814    0.000    0.000   1.780006E+02     5


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   1.1914E-01  -4.9137E+00   1.0938E+01   2.083871E-01   1.751504E-01        5        4   0.8000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        5        5    1.000      (Z,P)         1
   8.8770E-02  -2.6266E+00   1.3181E-04   2.294018E-03   7.439952E+01        5        4   0.8000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   6.440619E-01
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   2.952564E-01

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   1.929856E-03
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   4.115053E+01


  Beam  sigma  matrix : 

   4.148157E-01   1.863436E-01   0.000000E+00   0.000000E+00
   1.863436E-01   8.717634E-02   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     3.792319E-02    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
     16  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
     17  Keyword, label(s) :  END                                                   


                             5 particles have been launched
                     Made  it  to  the  end :      5

************************************************************************************************************************************
 Pgm zgoubi : Execution ended normally, upon keyword END or FIN
   
            ZGOUBI RUN COMPLETED. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  06-04-2018,  at  10:45:42 
  JOB  ENDED  ON    06-04-2018,  AT  10:45:43 

   CPU time, total :    0.89738099999999998     
