Centering 5 helical orbits in the AGS warm helical snake 3-D OPERA map.
 'OBJET'            ! This data list may be copy-pasted and run, as is.                                       1
7.2051782956D3      ! Reference rigidity of the problem
2
1 1                 ! A single particle. Initial coordinates :
-2.2  0. 0. 0. 0.  1. 'o'           ! Yo, To, Zo, Po, So, p/po
1
 'PARTICUL'         ! proton data are necessary for spin tracking                                             2
938.27203 1.602176487E-19 1.7928474 0 0     ! M, Q, G factor
 
 'SPNTRK'                                                                                                     3
4.1                 ! Initial spin is positionned vertical
0. 0. 1.
 
 'FAISCEAU'                                                                                                   4
 'SPNPRT'                                                                                                     5
 
 'TOSCA'                                                                                                      6
0  20
1.e1   100. 100. 100.
HEADER_0 wsnake
801 29 29 12.1      ! The map is a 801x29x29 node 3-D mesh
warmSnake.map       ! AGS warm snake 3-D OPERA map
0 0 0 0
2
.1
2  0.  .0  0.  0.
 
 'FAISCEAU'                                                                                                   7
 
 'FIT'                                                                                                        8
2   save  nofinal   ! Two variables. Save FIT variables (in zgoubi.FITVALS.out).
1 30 0  [-3,3]          ! Vary initial coordinate Y_0 (horiz. position)
1 31 0  [-3,3]          ! Vary initial coordinate T_0 (horiz. angle)
3    1E-2           ! Three constraints (penalty 2E-3 requested) :
3.1 1 2 6 0. 1. 0       ! Y_0=Y after at exit of the magnet
3.1 1 3 6 0. 1. 0       ! T_0=T after at exit of the magnet
7.3 1 2 6 0. 1. 0       ! Y_min+Y_max=0 inside the OPERA field map
 
 'SPNPRT'  PRINT    ! Stack spin data (in zgoubi.SPNPRT.Out).                                                 9
 
 'SYSTEM'           ! Save zgoubi.FITVALS.out data following from successive REBELOTE.                       10
1
cat zgoubi.FITVALS.out >> zgoubi.FITVALS.out_cat
 
 'REBELOTE'         ! Will loop on re-doing the FIT for 4 additional particle rigidities                     11
4  0.1  0 1
1                   ! List of 4 successive values of D in OBJET follows
OBJET  35  1.3872739973  2.1368296674  4.8261190694  11.015241321
 
 'SYSTEM'           ! Save a copy of zgoubi.FITVALS.out_cat and of zgoubi.SPNPRT.Out_cat                     12
2
\cp zgoubi.FITVALS.out_cat zgoubi.FITVALS.out_cat_copy
\cp zgoubi.SPNPRT.Out      zgoubi.SPNPRT.Out_copy
 
 'END'                                                                                                       13

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


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
                    electric rigidity (MeV) :   1981.220867    =T[eV]*(gamma+1)/gamma, such that dev.=E*L/rigidity
  
 I, AMQ(1,I), AMQ(2,I)/QE, P/Pref, v/c, time, s :
  
     1   9.38272030E+02  1.00000000E+00  1.00000000E+00  9.17207207E-01  0.00000000E+00  0.00000000E+00

************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                                                


                SPIN  TRACKING  REQUESTED


                          Particle  mass          =    938.2720     MeV/c2
                          Gyromagnetic  factor  G =    1.792847    

                          Initial spin conditions type  4 :
                              Same spin for all particles
                              Particles # 1 to 1 may be subjected to spin matching using FIT procedure

                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.000000
                               <SY> =     0.000000
                               <SZ> =     1.000000
                               <S>  =     1.000000

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      3)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)      D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

o  1   1.0000    -2.200     0.000     0.000     0.000      0.0000    0.0000   -2.200    0.000    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


---------------  Concentration ellipses : 
   surface      alpha        beta         <X>            <XP>           numb. of prtcls   ratio      space      pass# 
                                                                        in ellips,  out 
   0.0000E+00   0.0000E+00   1.0000E+00  -2.200000E-02   0.000000E+00        1        1    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   1.416766E+03        1        1    1.000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      5  Keyword, label(s) :  SPNPRT                                                



                         Momentum  group  #1 ; average  over 1 particles at this pass : 

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000         0.000000   0.000000   1.000000   1.000000   4.500000   0.000000   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|         SX        SY        SZ        |S|       GAMMA   (Si,Sf)   (Z,Sf_X)  (Z,Sf)
                                                                                                      (deg.)    (deg.)    (deg.)
                                                                                           (Sf_X : projection of Sf on YZ plane)

 o  1  0.000000  0.000000  1.000000  1.000000     0.000000  0.000000  1.000000  1.000000      2.5100    0.000    0.000    0.000    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.00000E+00  2.50997E+00      1   1

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                                                 


                OPEN FILE zgoubi.plt                                                                      
                FOR PRINTING TRAJECTORIES


     Status of MOD.MOD2 is 12.1 ; NDIM = 3 ; number of field data files used is   1.

           New field map(s) opened, cartesian mesh (MOD.le.19) ; 
           name(s) of map data file(s) : 

          'warmSnake.map'

   ----
   Map file number    1 ( of 1) :

     warmSnake.map field map,
     FORMAT type : regular.

 HEADER  (0 lines) : 
  Sbr fmapw/fmapr3 : completed job of reading map.
     Field map scaled by factor  a(1) =   0.000000E+00

     Min/max fields found in map (series) read   :  -2.934230E-01  /   3.119730E-01
       @  X,  Y, Z : -1.045E-02 -6.000E-02  6.500E-02              /  1.050E-02  5.500E-02  7.000E-02
     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
     this yields min/max normalised fields (kG) : -2.934230E+00                      3.119730E+00
       @  X (cm),  Y (cm), Z (cm) :  -1.04      -6.00       6.50   /   1.05       5.50       7.00    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   29 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1  1.0000    -2.200     0.000     0.000     0.000          599.500    -2.182     0.000    -0.055     0.000            1


                CONDITIONS  DE  MAXWELL  (     4004.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -3.3328E-12      -1.0454E-13       -1.5080E-14
                                      1.5768E-12       -9.5122E-13
                                     -2.7052E-13       -6.0426E-12
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      6)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)      D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

o  1   1.0000    -2.200     0.000     0.000     0.000      0.0000    0.0000   -2.182    0.090   -0.055    0.009   4.003418E+02     1
               Time of flight (mus) :  1.45593751E-02 mass (MeV/c2) :   938.272    


---------------  Concentration ellipses : 
   surface      alpha        beta         <X>            <XP>           numb. of prtcls   ratio      space      pass# 
                                                                        in ellips,  out 
   0.0000E+00   0.0000E+00   1.0000E+00  -2.181677E-02   9.039558E-05        1        1    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00  -5.474115E-04   8.622527E-06        1        1    1.000      (Z,P)         1
   0.0000E+00   0.0000E+00   1.0000E+00   1.455938E-02   1.416766E+03        1        1    1.000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      8  Keyword, label(s) :  FIT                                                   

     FIT procedure launched.

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            1 ,   ok.
           variable #            2       IP =           31 ,   ok.
           constraint #            1       IR =            6 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            6 ,   ok.
           constraint #            2       I  =            1 ,   ok.
           constraint #            3       IR =            6 ,   ok.

                    FIT  variables  and  constraints  in  good  order,  FIT  will proceed. 

                    Final FIT status will be saved in zgoubi.FITVALS.out                                                              


 STATUS OF VARIABLES  (Iteration #     0 /    999 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30   -3.00       -1.84      -1.8400000       3.00      2.000E-02  OBJET      -                    -                   
   1   2    31   -3.00       6.000E-02  6.00000000E-02   3.00      2.000E-02  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-02)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     6    0.000000E+00    1.000E+00    8.141823E-03    7.26E-03 TOSCA      -                    -                    0
  3   1   3     6    0.000000E+00    1.000E+00    8.162918E-02    7.29E-01 TOSCA      -                    -                    0
  7   1   2     6    0.000000E+00    1.000E+00   -4.904629E-02    2.63E-01 TOSCA      -                    -                    0
 Fit reached penalty value   9.1352E-03

   Last run following FIT[2] is skipped, as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      9  Keyword, label(s) :  SPNPRT      PRINT                                     



                         Momentum  group  #1 ; average  over 1 particles at this pass : 

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000        -0.000966   0.208322   0.978060   1.000000   4.500000  12.024168   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|         SX        SY        SZ        |S|       GAMMA   (Si,Sf)   (Z,Sf_X)  (Z,Sf)
                                                                                                      (deg.)    (deg.)    (deg.)
                                                                                           (Sf_X : projection of Sf on YZ plane)

 o  1  0.000000  0.000000  1.000000  1.000000    -0.000966  0.208322  0.978060  1.000000      2.5100   12.024   12.024   12.024    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  1.00000E+00  2.50997E+00      1   1

************************************************************************************************************************************
     10  Keyword, label(s) :  SYSTEM                                                

     Number  of  commands :   1,  as  follows : 

 cat zgoubi.FITVALS.out >> zgoubi.FITVALS.out_cat

************************************************************************************************************************************
     11  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        1 through the optical structure 

                     Total of          1 particles have been launched

     Multiple pass, 
          from element #     1 : OBJET     /label1=                    /label2=                    
                             to  REBELOTE  /label1=                    /label2=                    
     ending at pass #       5 at element #    11 : REBELOTE  /label1=                    /label2=                    


 Pgm rebel. At pass #    1/   5.  In element #    1,  parameter # 35  changed to    1.38727400E+00   (was    1.00000000E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

     An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.


************************************************************************************************************************************

 STATUS OF VARIABLES  (Iteration #     0 /    999 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30   -3.00       -1.30      -1.3000000       3.00      2.000E-02  OBJET      -                    -                   
   1   2    31   -3.00       0.120      0.12000000       3.00      2.000E-02  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-02)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     6    0.000000E+00    1.000E+00    3.656197E-02    2.20E-01 TOSCA      -                    -                    0
  3   1   3     6    0.000000E+00    1.000E+00    5.916558E-02    5.77E-01 TOSCA      -                    -                    0
  7   1   2     6    0.000000E+00    1.000E+00    3.510739E-02    2.03E-01 TOSCA      -                    -                    0
 Fit reached penalty value   6.0699E-03

   Last run following FIT[2] is skipped, as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      9  Keyword, label(s) :  SPNPRT      PRINT                                     



                         Momentum  group  #1 ; average  over 1 particles at this pass : 

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000        -0.000671   0.196147   0.980574   1.000000   6.000000  11.311791   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|         SX        SY        SZ        |S|       GAMMA   (Si,Sf)   (Z,Sf_X)  (Z,Sf)
                                                                                                      (deg.)    (deg.)    (deg.)
                                                                                           (Sf_X : projection of Sf on YZ plane)

 o  1  0.000000  0.000000  1.000000  1.000000    -0.000671  0.196147  0.980574  1.000000      3.3466   11.312   11.312   11.312    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  1.38727E+00  3.34663E+00      1   1

************************************************************************************************************************************
     10  Keyword, label(s) :  SYSTEM                                                

     Number  of  commands :   1,  as  follows : 

 cat zgoubi.FITVALS.out >> zgoubi.FITVALS.out_cat

************************************************************************************************************************************
     11  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        2 through the optical structure 

                     Total of          2 particles have been launched

 Pgm rebel. At pass #    2/   5.  In element #    1,  parameter # 35  changed to    2.13682967E+00   (was    1.38727400E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

     An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.


************************************************************************************************************************************

 STATUS OF VARIABLES  (Iteration #     0 /    999 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30   -3.00      -0.880     -0.88000000       3.00      2.000E-02  OBJET      -                    -                   
   1   2    31   -3.00       0.120      0.12000000       3.00      2.000E-02  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-02)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     6    0.000000E+00    1.000E+00    4.256638E-02    4.41E-01 TOSCA      -                    -                    0
  3   1   3     6    0.000000E+00    1.000E+00    2.834126E-02    1.96E-01 TOSCA      -                    -                    0
  7   1   2     6    0.000000E+00    1.000E+00   -3.860018E-02    3.63E-01 TOSCA      -                    -                    0
 Fit reached penalty value   4.1051E-03

   Last run following FIT[2] is skipped, as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      9  Keyword, label(s) :  SPNPRT      PRINT                                     



                         Momentum  group  #1 ; average  over 1 particles at this pass : 

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000        -0.000530   0.188579   0.982058   1.000000   9.000000  10.869907   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|         SX        SY        SZ        |S|       GAMMA   (Si,Sf)   (Z,Sf_X)  (Z,Sf)
                                                                                                      (deg.)    (deg.)    (deg.)
                                                                                           (Sf_X : projection of Sf on YZ plane)

 o  1  0.000000  0.000000  1.000000  1.000000    -0.000530  0.188579  0.982058  1.000000      5.0199   10.870   10.870   10.870    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  2.13683E+00  5.01995E+00      1   1

************************************************************************************************************************************
     10  Keyword, label(s) :  SYSTEM                                                

     Number  of  commands :   1,  as  follows : 

 cat zgoubi.FITVALS.out >> zgoubi.FITVALS.out_cat

************************************************************************************************************************************
     11  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        3 through the optical structure 

                     Total of          3 particles have been launched

 Pgm rebel. At pass #    3/   5.  In element #    1,  parameter # 35  changed to    4.82611907E+00   (was    2.13682967E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

     An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.


************************************************************************************************************************************

 STATUS OF VARIABLES  (Iteration #     0 /    999 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30   -3.00      -0.400     -0.40000000       3.00      2.000E-02  OBJET      -                    -                   
   1   2    31   -3.00       0.120      0.12000000       3.00      2.000E-02  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-02)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     6    0.000000E+00    1.000E+00    4.565989E-02    7.55E-01 TOSCA      -                    -                    0
  3   1   3     6    0.000000E+00    1.000E+00    1.223178E-02    5.42E-02 TOSCA      -                    -                    0
  7   1   2     6    0.000000E+00    1.000E+00   -2.292763E-02    1.90E-01 TOSCA      -                    -                    0
 Fit reached penalty value   2.7601E-03

   Last run following FIT[2] is skipped, as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      9  Keyword, label(s) :  SPNPRT      PRINT                                     



                         Momentum  group  #1 ; average  over 1 particles at this pass : 

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000        -0.000434   0.184148   0.982898   1.000000  20.000001  10.611513   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|         SX        SY        SZ        |S|       GAMMA   (Si,Sf)   (Z,Sf_X)  (Z,Sf)
                                                                                                      (deg.)    (deg.)    (deg.)
                                                                                           (Sf_X : projection of Sf on YZ plane)

 o  1  0.000000  0.000000  1.000000  1.000000    -0.000434  0.184148  0.982898  1.000000     11.1554   10.612   10.611   10.612    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  4.82612E+00  1.11554E+01      1   1

************************************************************************************************************************************
     10  Keyword, label(s) :  SYSTEM                                                

     Number  of  commands :   1,  as  follows : 

 cat zgoubi.FITVALS.out >> zgoubi.FITVALS.out_cat

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
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                                                


                SPIN  TRACKING  REQUESTED


                          Particle  mass          =    938.2720     MeV/c2
                          Gyromagnetic  factor  G =    1.792847    

                          Initial spin conditions type  4 :
                              Same spin for all particles
                              Particles # 1 to 1 may be subjected to spin matching using FIT procedure

                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.000000
                               <SY> =     0.000000
                               <SZ> =     1.000000
                               <S>  =     1.000000

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      3)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)      D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

o  1  11.0152    -0.400     0.120     0.000     0.000      0.0000   10.0152   -0.400    0.120    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


---------------  Concentration ellipses : 
   surface      alpha        beta         <X>            <XP>           numb. of prtcls   ratio      space      pass# 
                                                                        in ellips,  out 
   0.0000E+00   0.0000E+00   1.0000E+00  -4.000000E-03   1.200000E-04        1        1    1.000      (Y,T)         5
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         5
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   2.287378E+04        1        1    1.000      (t,K)         5

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      5  Keyword, label(s) :  SPNPRT                                                



                         Momentum  group  #1 ; average  over 1 particles at this pass : 

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000         0.000000   0.000000   1.000000   1.000000  45.500002   0.000000   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|         SX        SY        SZ        |S|       GAMMA   (Si,Sf)   (Z,Sf_X)  (Z,Sf)
                                                                                                      (deg.)    (deg.)    (deg.)
                                                                                           (Sf_X : projection of Sf on YZ plane)

 o  1  0.000000  0.000000  1.000000  1.000000     0.000000  0.000000  1.000000  1.000000     25.3786    0.000    0.000    0.000    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  1.10152E+01  2.53786E+01      1   1

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     Status of MOD.MOD2 is 12.1 ; NDIM = 3 ; number of field data files used is   1.

     No  map  file  opened :  was  already  stored.
     Skip  reading  field  map  file :           warmSnake.map
     Restored mesh coordinates for field map #IMAP=   1,  name : warmSnake.map
     Field map scaled by factor  a(1) =   0.000000E+00

     Min/max fields found in map (series) read   :  -2.934230E-01  /   3.119730E-01
       @  X,  Y, Z : -1.045E-02 -6.000E-02  6.500E-02              /  1.050E-02  5.500E-02  7.000E-02
     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
     this yields min/max normalised fields (kG) : -2.934230E+00                      3.119730E+00
       @  X (cm),  Y (cm), Z (cm) :  -1.04      -6.00       6.50   /   1.05       5.50       7.00    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   29 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1 11.0152    -0.400     0.120     0.000     0.000          599.500    -0.353     0.000    -0.006    -0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -4.7510E-14      -2.4439E-15       -2.3047E-16
                                      1.4452E-13       -8.7736E-14
                                     -2.3242E-14       -5.5918E-13
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   2.909389E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      6)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)      D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

o  1  11.0152    -0.400     0.120     0.000     0.000      0.0000   10.0152   -0.353    0.115   -0.006   -0.000   4.000028E+02     1
               Time of flight (mus) :  1.33530270E-02 mass (MeV/c2) :   938.272    


---------------  Concentration ellipses : 
   surface      alpha        beta         <X>            <XP>           numb. of prtcls   ratio      space      pass# 
                                                                        in ellips,  out 
   0.0000E+00   0.0000E+00   1.0000E+00  -3.528743E-03   1.154151E-04        1        1    1.000      (Y,T)         5
   0.0000E+00   0.0000E+00   1.0000E+00  -5.722836E-05  -2.200187E-07        1        1    1.000      (Z,P)         5
   0.0000E+00   0.0000E+00   1.0000E+00   1.335303E-02   2.287378E+04        1        1    1.000      (t,K)         5

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      8  Keyword, label(s) :  FIT                                                   


     An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.


************************************************************************************************************************************

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                                                


                SPIN  TRACKING  REQUESTED


                          Particle  mass          =    938.2720     MeV/c2
                          Gyromagnetic  factor  G =    1.792847    

                          Initial spin conditions type  4 :
                              Same spin for all particles
                              Particles # 1 to 1 may be subjected to spin matching using FIT procedure

                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.000000
                               <SY> =     0.000000
                               <SZ> =     1.000000
                               <S>  =     1.000000

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      3)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)      D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

o  1  11.0152    -0.400     0.120     0.000     0.000      0.0000   10.0152   -0.400    0.120    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


---------------  Concentration ellipses : 
   surface      alpha        beta         <X>            <XP>           numb. of prtcls   ratio      space      pass# 
                                                                        in ellips,  out 
   0.0000E+00   0.0000E+00   1.0000E+00  -4.000000E-03   1.200000E-04        1        1    1.000      (Y,T)         5
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         5
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   2.287378E+04        1        1    1.000      (t,K)         5

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      5  Keyword, label(s) :  SPNPRT                                                



                         Momentum  group  #1 ; average  over 1 particles at this pass : 

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000         0.000000   0.000000   1.000000   1.000000  45.500002   0.000000   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|         SX        SY        SZ        |S|       GAMMA   (Si,Sf)   (Z,Sf_X)  (Z,Sf)
                                                                                                      (deg.)    (deg.)    (deg.)
                                                                                           (Sf_X : projection of Sf on YZ plane)

 o  1  0.000000  0.000000  1.000000  1.000000     0.000000  0.000000  1.000000  1.000000     25.3786    0.000    0.000    0.000    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  1.10152E+01  2.53786E+01      1   1

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     Status of MOD.MOD2 is 12.1 ; NDIM = 3 ; number of field data files used is   1.

     No  map  file  opened :  was  already  stored.
     Skip  reading  field  map  file :           warmSnake.map
     Restored mesh coordinates for field map #IMAP=   1,  name : warmSnake.map
     Field map scaled by factor  a(1) =   0.000000E+00

     Min/max fields found in map (series) read   :  -2.934230E-01  /   3.119730E-01
       @  X,  Y, Z : -1.045E-02 -6.000E-02  6.500E-02              /  1.050E-02  5.500E-02  7.000E-02
     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
     this yields min/max normalised fields (kG) : -2.934230E+00                      3.119730E+00
       @  X (cm),  Y (cm), Z (cm) :  -1.04      -6.00       6.50   /   1.05       5.50       7.00    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   29 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1 11.0152    -0.400     0.120     0.000     0.000          599.500    -0.353     0.000    -0.006    -0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -4.7510E-14      -2.4439E-15       -2.3047E-16
                                      1.4452E-13       -8.7736E-14
                                     -2.3242E-14       -5.5918E-13
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      6)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)      D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

o  1  11.0152    -0.400     0.120     0.000     0.000      0.0000   10.0152   -0.353    0.115   -0.006   -0.000   4.000028E+02     1
               Time of flight (mus) :  1.33530270E-02 mass (MeV/c2) :   938.272    


---------------  Concentration ellipses : 
   surface      alpha        beta         <X>            <XP>           numb. of prtcls   ratio      space      pass# 
                                                                        in ellips,  out 
   0.0000E+00   0.0000E+00   1.0000E+00  -3.528743E-03   1.154151E-04        1        1    1.000      (Y,T)         5
   0.0000E+00   0.0000E+00   1.0000E+00  -5.722836E-05  -2.200187E-07        1        1    1.000      (Z,P)         5
   0.0000E+00   0.0000E+00   1.0000E+00   1.335303E-02   2.287378E+04        1        1    1.000      (t,K)         5

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALF^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      8  Keyword, label(s) :  FIT                                                   

     FIT procedure launched.

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            1 ,   ok.
           variable #            2       IP =           31 ,   ok.
           constraint #            1       IR =            6 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            6 ,   ok.
           constraint #            2       I  =            1 ,   ok.
           constraint #            3       IR =            6 ,   ok.

                    FIT  variables  and  constraints  in  good  order,  FIT  will proceed. 

                    Final FIT status will be saved in zgoubi.FITVALS.out                                                              


 STATUS OF VARIABLES  (Iteration #     0 /    999 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30   -3.00      -0.220     -0.22000000       3.00      2.000E-02  OBJET      -                    -                   
   1   2    31   -3.00       0.120      0.12000000       3.00      2.000E-02  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-02)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     6    0.000000E+00    1.000E+00    4.704789E-02    2.30E-01 TOSCA      -                    -                    0
  3   1   3     6    0.000000E+00    1.000E+00    4.978871E-03    2.57E-03 TOSCA      -                    -                    0
  7   1   2     6    0.000000E+00    1.000E+00   -8.599176E-02    7.68E-01 TOSCA      -                    -                    0
 Fit reached penalty value   9.6329E-03

   Last run following FIT[2] is skipped, as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      9  Keyword, label(s) :  SPNPRT      PRINT                                     



                         Momentum  group  #1 ; average  over 1 particles at this pass : 

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000        -0.000419   0.183281   0.983060   1.000000  45.500002  10.560969   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|         SX        SY        SZ        |S|       GAMMA   (Si,Sf)   (Z,Sf_X)  (Z,Sf)
                                                                                                      (deg.)    (deg.)    (deg.)
                                                                                           (Sf_X : projection of Sf on YZ plane)

 o  1  0.000000  0.000000  1.000000  1.000000    -0.000419  0.183281  0.983060  1.000000     25.3786   10.561   10.561   10.561    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  1.10152E+01  2.53786E+01      1   1

************************************************************************************************************************************
     10  Keyword, label(s) :  SYSTEM                                                

     Number  of  commands :   1,  as  follows : 

 cat zgoubi.FITVALS.out >> zgoubi.FITVALS.out_cat

************************************************************************************************************************************
     11  Keyword, label(s) :  REBELOTE                                              


                         ****  End  of  'REBELOTE'  procedure  ****

      There  has  been          5  passes  through  the  optical  structure 

                     Total of          5 particles have been launched

************************************************************************************************************************************
     12  Keyword, label(s) :  SYSTEM                                                

     Number  of  commands :   2,  as  follows : 

 \cp zgoubi.FITVALS.out_cat zgoubi.FITVALS.out_cat_copy
 \cp zgoubi.SPNPRT.Out      zgoubi.SPNPRT.Out_copy

************************************************************************************************************************************
     13  Keyword, label(s) :  END                                                   


************************************************************************************************************************************
 Pgm zgoubi : Execution ended normally, upon keyword END or FIN
   
            ZGOUBI RUN COMPLETED. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  12-07-2018,  at  16:50:16 
  JOB  ENDED  ON    12-07-2018,  AT  16:50:21 

   CPU time, total :     4.2244560000000000     

     An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.


************************************************************************************************************************************
