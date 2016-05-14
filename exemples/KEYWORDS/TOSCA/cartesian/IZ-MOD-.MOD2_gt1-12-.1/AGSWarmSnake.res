Centering 5 helical orbits in the AGS warm helical snake 3-D OPERA map.
 'OBJET'                                                                                                      1
7.2051782956d3      ! Reference rigidity of the problem
2
1 1                 ! A single particle. Initial coordinates :
-2.2  0. 0. 0. 0.  1. 'o'           ! Yo, To, Zo, Po, So, p/po
1
 'PARTICUL'         ! proton / data used for spin tracking :                                                  2
938.27203 1.602176487E-19 1.7928474 0 0     ! M, Q, G factor
 
 'SPNTRK'                                                                                                     3
4.1                 ! Initial spin is positionned vertical
0. 0. 1.
 
 'TOSCA'                                                                                                      4
0  20
1.e1   100. 100. 100.
HEADER_0 wsnake
801 29 29 12.1
table55.tab       ! AGS warm snake 3-D OPERA map
0 0 0 0
2
.1
2  0.  .00  0.  0.
 
 'FIT'                                                                                                        5
2                   ! Two variables :
1 30 0  [-3,3]          ! Vary initial coordinate Y_0 (horiz. position)
1 31 0  [-3,3]          ! Vary initial coordinate T_0 (horiz. angle)
3    2E-3           ! Three constraints (penalty 2E-3 requested) :
3.1 1 2 4 0. 1. 0       ! Y_0=Y after at exit of the magnet
3.1 1 3 4 0. 1. 0       ! T_0=T after at exit of the magnet
7.3 1 2 4 0. 1. 0       ! Y_min+Y_max=0 inside the OPERA field map
 
 'SPNPRT'  PRINT    ! Print spin coordinates in zgoubi.SPNPRT.Out                                             6
 
 'REBELOTE'         ! Will re-do the FIT for 4 different rigidities                                           7
4  0.1  0 1
1                   ! List of 4 successive values of D in OBJET follows
OBJET  35  1.3872739973  2.1368296674  4.8261190694  11.015241321
 
 'END'                                                                                                        8

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
                    electric rigidity (MeV) :   1981.220867      =T[eV]*(gamma+1)/gamma, such that dev.=E*L/rigidity
  
 I, AMQ(1,I), AMQ(2,I)/QE, P/Pref, v/c, time, s :
  
     1   9.38272030E+02  1.00000000E+00  1.00000000E+00  9.17207207E-01  0.00000000E+00  0.00000000E+00

************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                            


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
      4  Keyword, label(s) :  TOSCA                             


                OPEN FILE zgoubi.plt                                                                      
                FOR PRINTING TRAJECTORIES


     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 12.1


           New field map(s) now used, cartesian mesh (MOD.le.19) ; 
           name(s) of map data file(s) : 

          table55.tab

   ----
   Map file number    1 ( of 1) :

     table55.tab map,  FORMAT type : regular
 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.


     Min/max fields seen in map   :  -2.934230E-01                 /   3.119730E-01
       @  X,  Y, Z : -1.045E-02 -6.000E-02  6.500E-02              /  1.050E-02  5.500E-02  7.000E-02

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.934230E+00                      3.119730E+00
       @  X (cm),  Y (cm), Z (cm) :  -1.04      -6.00       6.50   /   1.05       5.50       7.00    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

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
      5  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            1 ,   ok.
           variable #            2       IP =           31 ,   ok.
           constraint #            1       IR =            4 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            4 ,   ok.
           constraint #            2       I  =            1 ,   ok.
           constraint #            3       IR =            4 ,   ok.

                    FIT  variables  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save [FileName]' command

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
  LMNT  VAR  PARAM   MINIMUM     INITIAL         FINAL         MAXIMUM     STEP     NAME       LBL1     LBL2
     1    1     30   -3.00       -1.82       -1.822930101       3.00      4.002E-92 OBJET      *          *         
     1    2     31   -3.00       4.255E-02   4.2551440329E-02   3.00      4.002E-92 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-03)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     4      0.0000000E+00     1.0000E+00    1.1962572E-04   1.7560E-06   TOSCA      *          *          *   0 : 
   3   1   3     4      0.0000000E+00     1.0000E+00    8.7960089E-02   9.4940E-01   TOSCA      *          *          *   0 : 
   7   1   2     4      0.0000000E+00     1.0000E+00   -2.0305509E-02   5.0595E-02   TOSCA      *          *          *   0 : 
 Fit reached penalty value   8.1493E-03

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

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
      4  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 12.1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           table55.tab
 Pgm toscac,  restored mesh coordinates for field map #   1,  name : table55.tab


     Min/max fields seen in map   :  -2.934230E-01                 /   3.119730E-01
       @  X,  Y, Z : -1.045E-02 -6.000E-02  6.500E-02              /  1.050E-02  5.500E-02  7.000E-02

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.934230E+00                      3.119730E+00
       @  X (cm),  Y (cm), Z (cm) :  -1.04      -6.00       6.50   /   1.05       5.50       7.00    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1  1.0000    -1.823     0.043     0.000     0.000          599.500    -1.823     0.000    -0.045     0.000            1


                CONDITIONS  DE  MAXWELL  (     4004.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -2.7700E-12      -1.0454E-13        7.1092E-15
                                      1.5768E-12       -9.5295E-13
                                     -2.7052E-13       -6.0520E-12
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    8.00000000     m ;  Time  (for ref. rigidity & particle) =   2.909389E-08 s 

************************************************************************************************************************************

      5   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
      6  Keyword, label(s) :  SPNPRT      PRINT                 



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000        -0.000949   0.208319   0.978060   1.000000   4.500000  12.023986   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000          -0.000949  0.208319  0.978060  1.000000      2.5100   12.0240   12.0239    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.4869E-04 -9.4869E-04  2.0832E-01  2.0832E-01  9.7806E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.00000E+00  2.50997E+00      1   1

************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                          


                                -----  REBELOTE  -----

     End of pass #        1 through the optical structure 

                     Total of          1 particles have been launched

     Multiple pass, 
          from element #     1 : OBJET     /label1=          /label2=           to REBELOTE /label1=          /label2=          
          ending at pass #       5 at element #     7 : REBELOTE  /label1=          /label2=          


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
     1    1     30   -3.00       -1.31       -1.306463750       3.00      1.395E-82 OBJET      *          *         
     1    2     31   -3.00       2.538E-02   2.5382980954E-02   3.00      1.395E-82 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-03)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     4      0.0000000E+00     1.0000E+00    9.5422494E-05   3.2061E-06   TOSCA      *          *          *   0 : 
   3   1   3     4      0.0000000E+00     1.0000E+00    5.2904384E-02   9.8549E-01   TOSCA      *          *          *   0 : 
   7   1   2     4      0.0000000E+00     1.0000E+00   -6.4182407E-03   1.4504E-02   TOSCA      *          *          *   0 : 
 Fit reached penalty value   2.8401E-03

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

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
      4  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 12.1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           table55.tab
 Pgm toscac,  restored mesh coordinates for field map #   1,  name : table55.tab


     Min/max fields seen in map   :  -2.934230E-01                 /   3.119730E-01
       @  X,  Y, Z : -1.045E-02 -6.000E-02  6.500E-02              /  1.050E-02  5.500E-02  7.000E-02

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.934230E+00                      3.119730E+00
       @  X (cm),  Y (cm), Z (cm) :  -1.04      -6.00       6.50   /   1.05       5.50       7.00    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1  1.3873    -1.306     0.025     0.000     0.000          599.500    -1.307     0.000    -0.038     0.000            1


                CONDITIONS  DE  MAXWELL  (     4002.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -1.4558E-12      -5.6219E-14       -2.7384E-15
                                      1.0172E-12       -6.7886E-13
                                     -1.4808E-13       -4.3142E-12
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    8.00000000     m ;  Time  (for ref. rigidity & particle) =   2.909389E-08 s 

************************************************************************************************************************************

      5   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
      6  Keyword, label(s) :  SPNPRT      PRINT                 



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000        -0.000669   0.196151   0.980574   1.000000   6.000000  11.312016   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000          -0.000669  0.196151  0.980574  1.000000      3.3466   11.3120   11.3120    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.4869E-04 -6.6853E-04  1.9615E-01  2.0832E-01  9.7806E-01  9.8057E-01  1.0000E+00  1.0000E+00  1.38727E+00  3.34663E+00      1   1

************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                          


                                -----  REBELOTE  -----

     End of pass #        2 through the optical structure 

                     Total of          2 particles have been launched

 Pgm rebel. At pass #    2/   5.  In element #    1,  parameter # 35  changed to    2.13682967E+00   (was    1.38727400E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
  LMNT  VAR  PARAM   MINIMUM     INITIAL         FINAL         MAXIMUM     STEP     NAME       LBL1     LBL2
     1    1     30   -3.00      -0.826      -0.8264637503       3.00      2.000E-02 OBJET      *          *         
     1    2     31   -3.00      -3.462E-02  -3.4617019046E-02   3.00      2.000E-02 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-03)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     4      0.0000000E+00     1.0000E+00    1.9760875E-02   2.6646E-01   TOSCA      *          *          *   0 : 
   3   1   3     4      0.0000000E+00     1.0000E+00    3.0581482E-02   6.3816E-01   TOSCA      *          *          *   0 : 
   7   1   2     4      0.0000000E+00     1.0000E+00    1.1822931E-02   9.5382E-02   TOSCA      *          *          *   0 : 
 Fit reached penalty value   1.4655E-03

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

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
      4  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 12.1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           table55.tab
 Pgm toscac,  restored mesh coordinates for field map #   1,  name : table55.tab


     Min/max fields seen in map   :  -2.934230E-01                 /   3.119730E-01
       @  X,  Y, Z : -1.045E-02 -6.000E-02  6.500E-02              /  1.050E-02  5.500E-02  7.000E-02

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.934230E+00                      3.119730E+00
       @  X (cm),  Y (cm), Z (cm) :  -1.04      -6.00       6.50   /   1.05       5.50       7.00    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1  2.1368    -0.826    -0.035     0.000     0.000          599.500    -0.846     0.000    -0.027     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -6.1228E-13      -2.3726E-14        1.3210E-16
                                      7.2217E-13       -4.4136E-13
                                     -1.0949E-13       -2.8078E-12
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    8.00000000     m ;  Time  (for ref. rigidity & particle) =   2.909389E-08 s 

************************************************************************************************************************************

      5   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
      6  Keyword, label(s) :  SPNPRT      PRINT                 



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000        -0.000494   0.188577   0.982058   1.000000   9.000000  10.869780   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000          -0.000494  0.188577  0.982058  1.000000      5.0199   10.8698   10.8697    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.4869E-04 -4.9358E-04  1.8858E-01  2.0832E-01  9.7806E-01  9.8206E-01  1.0000E+00  1.0000E+00  2.13683E+00  5.01995E+00      1   1

************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                          


                                -----  REBELOTE  -----

     End of pass #        3 through the optical structure 

                     Total of          3 particles have been launched

 Pgm rebel. At pass #    3/   5.  In element #    1,  parameter # 35  changed to    4.82611907E+00   (was    2.13682967E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
  LMNT  VAR  PARAM   MINIMUM     INITIAL         FINAL         MAXIMUM     STEP     NAME       LBL1     LBL2
     1    1     30   -3.00      -0.386      -0.3864637503       3.00      6.667E-03 OBJET      *          *         
     1    2     31   -3.00       6.538E-02   6.5382980954E-02   3.00      6.667E-03 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-03)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     4      0.0000000E+00     1.0000E+00    2.3799058E-02   7.3020E-01   TOSCA      *          *          *   0 : 
   3   1   3     4      0.0000000E+00     1.0000E+00    1.2288927E-02   1.9469E-01   TOSCA      *          *          *   0 : 
   7   1   2     4      0.0000000E+00     1.0000E+00   -7.6326512E-03   7.5106E-02   TOSCA      *          *          *   0 : 
 Fit reached penalty value   7.7567E-04

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

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
      4  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 12.1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           table55.tab
 Pgm toscac,  restored mesh coordinates for field map #   1,  name : table55.tab


     Min/max fields seen in map   :  -2.934230E-01                 /   3.119730E-01
       @  X,  Y, Z : -1.045E-02 -6.000E-02  6.500E-02              /  1.050E-02  5.500E-02  7.000E-02

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.934230E+00                      3.119730E+00
       @  X (cm),  Y (cm), Z (cm) :  -1.04      -6.00       6.50   /   1.05       5.50       7.00    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1  4.8261    -0.386     0.065     0.000     0.000          599.500    -0.363     0.000    -0.013     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -1.1123E-13      -5.5780E-15       -6.0172E-16
                                      3.2984E-13       -1.9945E-13
                                     -5.3047E-14       -1.2710E-12
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    8.00000000     m ;  Time  (for ref. rigidity & particle) =   2.909389E-08 s 

************************************************************************************************************************************

      5   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
      6  Keyword, label(s) :  SPNPRT      PRINT                 



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000        -0.000423   0.184148   0.982898   1.000000  20.000001  10.611502   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000          -0.000423  0.184148  0.982898  1.000000     11.1554   10.6115   10.6115    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.4869E-04 -4.2316E-04  1.8415E-01  2.0832E-01  9.7806E-01  9.8290E-01  1.0000E+00  1.0000E+00  4.82612E+00  1.11554E+01      1   1

************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                          


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
      4  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 12.1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           table55.tab
 Pgm toscac,  restored mesh coordinates for field map #   1,  name : table55.tab


     Min/max fields seen in map   :  -2.934230E-01                 /   3.119730E-01
       @  X,  Y, Z : -1.045E-02 -6.000E-02  6.500E-02              /  1.050E-02  5.500E-02  7.000E-02

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.934230E+00                      3.119730E+00
       @  X (cm),  Y (cm), Z (cm) :  -1.04      -6.00       6.50   /   1.05       5.50       7.00    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1 11.0152    -0.386     0.065     0.000     0.000          599.500    -0.361     0.000    -0.006     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -4.8711E-14      -2.4439E-15       -2.4364E-16
                                      1.4452E-13       -8.7736E-14
                                     -2.3242E-14       -5.5918E-13
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   4.364084E-08 s 

************************************************************************************************************************************
      5  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1


                    Saving new version of zgoubi.dat to zgoubi.FIT.out.dat, with variables updated.

                    Updated version of zgoubi.dat saved in  zgoubi.FIT.out.dat

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
      4  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 12.1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           table55.tab
 Pgm toscac,  restored mesh coordinates for field map #   1,  name : table55.tab


     Min/max fields seen in map   :  -2.934230E-01                 /   3.119730E-01
       @  X,  Y, Z : -1.045E-02 -6.000E-02  6.500E-02              /  1.050E-02  5.500E-02  7.000E-02

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.934230E+00                      3.119730E+00
       @  X (cm),  Y (cm), Z (cm) :  -1.04      -6.00       6.50   /   1.05       5.50       7.00    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1 11.0152    -0.386     0.065     0.000     0.000          599.500    -0.361     0.000    -0.006     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -4.8711E-14      -2.4439E-15       -2.4364E-16
                                      1.4452E-13       -8.7736E-14
                                     -2.3242E-14       -5.5918E-13
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

************************************************************************************************************************************
      5  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            1 ,   ok.
           variable #            2       IP =           31 ,   ok.
           constraint #            1       IR =            4 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            4 ,   ok.
           constraint #            2       I  =            1 ,   ok.
           constraint #            3       IR =            4 ,   ok.

                    FIT  variables  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save [FileName]' command

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
  LMNT  VAR  PARAM   MINIMUM     INITIAL         FINAL         MAXIMUM     STEP     NAME       LBL1     LBL2
     1    1     30   -3.00      -0.146      -0.1464637503       3.00      2.000E-02 OBJET      *          *         
     1    2     31   -3.00       5.383E-03   5.3829809538E-03   3.00      2.000E-02 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-03)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     4      0.0000000E+00     1.0000E+00    1.1677191E-03   1.0710E-03   TOSCA      *          *          *   0 : 
   3   1   3     4      0.0000000E+00     1.0000E+00    5.1455260E-03   2.0795E-02   TOSCA      *          *          *   0 : 
   7   1   2     4      0.0000000E+00     1.0000E+00    3.5290075E-02   9.7813E-01   TOSCA      *          *          *   0 : 
 Fit reached penalty value   1.2732E-03

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

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
      4  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 12.1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           table55.tab
 Pgm toscac,  restored mesh coordinates for field map #   1,  name : table55.tab


     Min/max fields seen in map   :  -2.934230E-01                 /   3.119730E-01
       @  X,  Y, Z : -1.045E-02 -6.000E-02  6.500E-02              /  1.050E-02  5.500E-02  7.000E-02

     Given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
       min/max normalised fields (kG) : -2.934230E+00                      3.119730E+00
       @  X (cm),  Y (cm), Z (cm) :  -1.04      -6.00       6.50   /   1.05       5.50       7.00    

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of mesh nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm

  A    1 11.0152    -0.146     0.005     0.000     0.000          599.500    -0.145     0.000    -0.006     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -1.8560E-14       4.3083E-17        1.4809E-16
                                      1.2338E-13       -8.7761E-14
                                     -1.9032E-14       -5.5921E-13
                       LAPLACIEN SCALAIRE =   0.000    



     KPOS =  2.  Element  is  mis-aligned  wrt.  the  optical  axis.
          X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD


 Cumulative length of optical axis =    8.00000000     m ;  Time  (for ref. rigidity & particle) =   2.909389E-08 s 

************************************************************************************************************************************

      5   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
      6  Keyword, label(s) :  SPNPRT      PRINT                 



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

      <SX>       <SY>       <SZ>       <|S|>            <SX>       <SY>       <SZ>      <|S|>    <G.gma>    <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
   0.000000   0.000000   1.000000   1.000000        -0.000395   0.183280   0.983061   1.000000  45.500002  10.560903   0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000          -0.000395  0.183280  0.983061  1.000000     25.3786   10.5609   10.5609    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.4869E-04 -3.9519E-04  1.8328E-01  2.0832E-01  9.7806E-01  9.8306E-01  1.0000E+00  1.0000E+00  1.10152E+01  2.53786E+01      1   1

************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                          


                         ****  End  of  'REBELOTE'  procedure  ****

      There  has  been          5  passes  through  the  optical  structure 

                     Total of          5 particles have been launched

************************************************************************************************************************************
      8  Keyword, label(s) :  END                               


************************************************************************************************************************************

           MAIN PROGRAM : Execution ended upon key  END       

************************************************************************************************************************************

                    Saving new version of zgoubi.dat to zgoubi.FIT.out.dat, with variables updated.

                    Updated version of zgoubi.dat saved in  zgoubi.FIT.out.dat
   
            ZGOUBI RUN COMPLETED. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  12-05-0016,  at  03:56:17 
  JOB  ENDED  ON    12-05-0016,  AT  03:56:54 

   CPU time, total :     36.8623030000000     
