Uniform field sector
 'OBJET'                                                                                                      1
64.62444403717985                                   ! 200keV
2
1 1
12.9248888074 0.  0.   0.   0.  1. 'm'              ! 200keV. R=Brho/B=*/.5
1 1 1 1
 
 'PARTICUL'   ! This is required only because we want to get the time-of-flight,                              2
PROTON       ! otherwise zgoubi only requires rigidity.
 
 'TOSCA'                                                                                                      3
0  0
1.  1. 1. 1.
HEADER_8
106 121 1 22.1 1.       IZ=1 -> 2D ; MOD=22 -> polar map ; .MOD2=.1 -> one map file
geneSectorMap.out
0 0 0 0
2
.01
2
0. 0. 0. 0.
 
 'FAISCEAU'                                                                                                   4
 
 'FIT'                                                                                                        5
1   nofinal
1 30 0 [12.,65.]     ! Variable : Y_0
2  2e-8  199         ! Penalty; max numb of calls to function
3.1 1 2 #End 0. 1. 0    ! Constraint :  Y_final=Y_0
3.1 1 3 #End 0. 1. 0    ! Constraint :  T_final=T_0
 
 'FAISTORE'                                                                                                   6
zgoubi.fai
1
 
 'REBELOTE'                                                                                                   7
20 0.2  0 1        ! 20 different rigidities ;  ; coordinates as found in OBJET ; change parameter(s)
1
OBJET 35     1:5.0063899693    ! 0.2 MeV to 5 MeV
 
 'END'                                                                                                        8

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


     Particle  properties :
     PROTON
                     Mass          =    938.272        MeV/c2
                     Charge        =   1.602176E-19    C     
                     G  factor     =    1.79285              
                     COM life-time =   1.000000E+99    s     

              Reference  data :
                    mag. rigidity (kG.cm)   :   64.624444      =p/q, such that dev.=B*L/rigidity
                    mass (MeV/c2)           :   938.27203    
                    momentum (MeV/c)        :   19.373921    
                    energy, total (MeV)     :   938.47203    
                    energy, kinetic (MeV)   :  0.20000000    
                    beta = v/c              :  2.0644111178E-02
                    gamma                   :   1.000213158    
                    beta*gamma              :  2.0648511631E-02
                    G*gamma                 :   1.793229509    
                    electric rigidity (MeV) :  0.3999573775      =T[eV]*(gamma+1)/gamma, such that dev.=E*L/rigidity
  
 I, AMQ(1,I), AMQ(2,I)/QE, P/Pref, v/c, time, s :
  
     1   9.38272030E+02  1.00000000E+00  1.00000000E+00  2.06441112E-02  0.00000000E+00  0.00000000E+00

************************************************************************************************************************************
      3  Keyword, label(s) :  TOSCA                                                 


     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00

          New field map(s) now used, polar mesh (MOD .ge. 20) ;  name(s) of map data file(s) are : 
          geneSectorMap.out

   ----
   Map file number    1 ( of 1) :

     geneSectorMap.out map,  FORMAT type : regular.  Field multiplication factor :  1.00000000E+00

 HEADER  (8 lines) : 
        10.000000000000000       0.50000000000000000       0.57142857142857140      -4.97067035767729219E-042      !  Rmi/cm,
      # Field map generated using geneSectorMap.f                                                                            
     # AT/rd,  AT/deg, Rmi/cm, Rma/cm, RM/cm, NR, dR/cm, NX, dX/cm, dA/rd :                                                  
     #   1.04719755E+00   6.00000000E+01   1.00000000E+01   7.00000000E+01   5.00000000E+01 121   5.00000000E-01 106   4.9866
      #                                                                                                                      
      # R, Z, X, BY, BZ, BX                                                                                                  
      # cm cm rd kG  kG  kG                                                                                                  
      #                                                                                                                      
 
R_min (cm), DR (cm), DTTA (deg), DZ (cm) :   1.000000E+01  5.000000E-01  5.714286E-01 -4.970670E-42



     Field map limits, angle :  min, max, max-min (rad) :  -0.52858543E+00   0.51861212E+00   0.10471976E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+02   0.70000000E+02   0.60000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      10.0     -0.00      /   0.00      10.0     -0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  106/ 121/   1
     Node  distance  in   X/Y/Z :   9.973310E-03/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :  1.0000E-02 cm   (i.e.,   2.5000E-04 rad  at mean radius RM =    40.00    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            0.519    12.925     0.000     0.000     0.000            1

 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      3)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    0.0000   12.925    0.000    0.000    0.000   1.353491E+01     1
               Time of flight (mus) :  2.18694843E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01   2.388097E-12        1        1    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         1
   0.0000E+00   0.0000E+00   1.0000E+00   2.186948E-02   2.000000E-01        1        1    1.000      (t,K)         1

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
      5  Keyword, label(s) :  FIT                                                   

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           constraint #            1       IR =            4 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            4 ,   ok.
           constraint #            2       I  =            1 ,   ok.

                    FIT  variables  and  constraints  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save [FileName]' command

 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        12.9       12.924889       65.0      0.177      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.818456E-11    5.80E-05 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    2.388097E-09    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   5.7033E-18

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


                OPEN FILE zgoubi.fai                                                                      
                FOR PRINTING COORDINATES 

               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        1 through the optical structure 

                     Total of          1 particles have been launched

     Multiple pass, 
          from element #     1 : OBJET     /label1=                    /label2=                    
                             to  REBELOTE  /label1=                    /label2=                    
     ending at pass #      21 at element #     7 : REBELOTE  /label1=                    /label2=                    


     Parameter #   35 in element #    1 will be modified at each pass. 
     list of requested parameter values :
                   1   1.00000000E+00
                   2   1.21086263E+00
                   3   1.42172526E+00
                   4   1.63258789E+00
                   5   1.84345052E+00
                   6   2.05431315E+00
                   7   2.26517578E+00
                   8   2.47603841E+00
                   9   2.68690104E+00
                  10   2.89776367E+00
                  11   3.10862630E+00
                  12   3.31948893E+00
                  13   3.53035156E+00
                  14   3.74121419E+00
                  15   3.95207682E+00
                  16   4.16293945E+00
                  17   4.37380208E+00
                  18   4.58466471E+00
                  19   4.79552734E+00
                  20   5.00638997E+00

 Pgm rebel. At pass #    1/  21.  In element #    1,  parameter # 35  changed to    1.00000000E+00   (was    1.00000000E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        12.9       12.924889       65.0      0.177      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.818456E-11    5.80E-05 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    2.388097E-09    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   5.7033E-18

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        2 through the optical structure 

                     Total of          2 particles have been launched

 Pgm rebel. At pass #    2/  21.  In element #    1,  parameter # 35  changed to    1.21086263E+00   (was    1.00000000E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        15.7       15.650266       65.0      2.992E-06  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    5.068885E-07    8.16E-05 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    5.609851E-05    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   3.1473E-09

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        3 through the optical structure 

                     Total of          3 particles have been launched

 Pgm rebel. At pass #    3/  21.  In element #    1,  parameter # 35  changed to    1.42172526E+00   (was    1.21086263E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        18.4       18.375643       65.0      2.992E-06  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.013794E-06    1.13E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    9.555826E-05    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   9.1324E-09

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        4 through the optical structure 

                     Total of          4 particles have been launched

 Pgm rebel. At pass #    4/  21.  In element #    1,  parameter # 35  changed to    1.63258789E+00   (was    1.42172526E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        21.1       21.101020       65.0      2.992E-06  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.520700E-06    1.48E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    1.248249E-04    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.5584E-08

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        5 through the optical structure 

                     Total of          5 particles have been launched

 Pgm rebel. At pass #    5/  21.  In element #    1,  parameter # 35  changed to    1.84345052E+00   (was    1.63258789E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        23.8       23.826394       65.0      9.973E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    5.316727E-07    1.89E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    3.864994E-05    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.4941E-09

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        6 through the optical structure 

                     Total of          6 particles have been launched

 Pgm rebel. At pass #    6/  21.  In element #    1,  parameter # 35  changed to    2.05431315E+00   (was    1.84345052E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        26.6       26.551771       65.0      2.992E-06  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.038579E-06    2.35E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    6.774977E-05    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   4.5911E-09

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        7 through the optical structure 

                     Total of          7 particles have been launched

 Pgm rebel. At pass #    7/  21.  In element #    1,  parameter # 35  changed to    2.26517578E+00   (was    2.05431315E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        29.3       29.277148       65.0      2.992E-06  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.545488E-06    2.86E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    9.143182E-05    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   8.3622E-09

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        8 through the optical structure 

                     Total of          8 particles have been launched

 Pgm rebel. At pass #    8/  21.  In element #    1,  parameter # 35  changed to    2.47603841E+00   (was    2.26517578E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        32.0       32.002516       65.0      8.078E-05  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    2.435406E-06    3.41E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    1.318097E-04    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.7380E-08

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        9 through the optical structure 

                     Total of          9 particles have been launched

 Pgm rebel. At pass #    9/  21.  In element #    1,  parameter # 35  changed to    2.68690104E+00   (was    2.47603841E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        34.7       34.727893       65.0      2.992E-06  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.928499E-06    4.02E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    9.618365E-05    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   9.2550E-09

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       10 through the optical structure 

                     Total of         10 particles have been launched

 Pgm rebel. At pass #   10/  21.  In element #    1,  parameter # 35  changed to    2.89776367E+00   (was    2.68690104E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        37.5       37.453270       65.0      2.992E-06  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.421595E-06    4.67E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    6.574239E-05    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   4.3241E-09

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       11 through the optical structure 

                     Total of         11 particles have been launched

 Pgm rebel. At pass #   11/  21.  In element #    1,  parameter # 35  changed to    3.10862630E+00   (was    2.89776367E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        40.2       40.178647       65.0      2.992E-06  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    9.146849E-07    5.38E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    3.943095E-05    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.5556E-09

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       12 through the optical structure 

                     Total of         12 particles have been launched

 Pgm rebel. At pass #   12/  21.  In element #    1,  parameter # 35  changed to    3.31948893E+00   (was    3.10862630E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        42.9       42.904024       65.0      2.992E-06  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    4.077765E-07    6.13E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    1.646222E-05    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.7117E-10

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       13 through the optical structure 

                     Total of         13 particles have been launched

 Pgm rebel. At pass #   13/  21.  In element #    1,  parameter # 35  changed to    3.53035156E+00   (was    3.31948893E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        45.6       45.629402       65.0      2.992E-06  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    9.912580E-08    6.94E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    3.762796E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.4168E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       14 through the optical structure 

                     Total of         14 particles have been launched

 Pgm rebel. At pass #   14/  21.  In element #    1,  parameter # 35  changed to    3.74121419E+00   (was    3.53035156E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        48.4       48.354770       65.0      8.078E-05  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    3.881756E-06    7.79E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    1.390435E-04    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.9348E-08

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       15 through the optical structure 

                     Total of         15 particles have been launched

 Pgm rebel. At pass #   15/  21.  In element #    1,  parameter # 35  changed to    3.95207682E+00   (was    3.74121419E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        51.1       51.080156       65.0      2.992E-06  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.112935E-06    8.69E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    3.773819E-05    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.4254E-09

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       16 through the optical structure 

                     Total of         16 particles have been launched

 Pgm rebel. At pass #   16/  21.  In element #    1,  parameter # 35  changed to    4.16293945E+00   (was    3.95207682E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        53.8       53.805524       65.0      8.078E-05  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    2.867951E-06    9.64E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    9.232207E-05    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   8.5316E-09

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       17 through the optical structure 

                     Total of         17 particles have been launched

 Pgm rebel. At pass #   17/  21.  In element #    1,  parameter # 35  changed to    4.37380208E+00   (was    4.16293945E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        56.5       56.530910       65.0      2.992E-06  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    2.126758E-06    1.06E-03 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    6.516153E-05    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   4.2505E-09

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       18 through the optical structure 

                     Total of         18 particles have been launched

 Pgm rebel. At pass #   18/  21.  In element #    1,  parameter # 35  changed to    4.58466471E+00   (was    4.37380208E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        59.3       59.256278       65.0      8.078E-05  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.854130E-06    1.17E-03 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    5.419623E-05    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.9407E-09

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       19 through the optical structure 

                     Total of         19 particles have been launched

 Pgm rebel. At pass #   19/  21.  In element #    1,  parameter # 35  changed to    4.79552734E+00   (was    4.58466471E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        62.0       61.981655       65.0      2.992E-06  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.347240E-06    1.28E-03 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    3.764777E-05    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.4192E-09

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       20 through the optical structure 

                     Total of         20 particles have been launched


      Next  pass  is  #    21 and  last  pass  through  the  optical  structure


 Pgm rebel. At pass #   20/  21.  In element #    1,  parameter # 35  changed to    5.00638997E+00   (was    4.79552734E+00)

     WRITE statements to zgoubi.res are re-established from now on.

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  TOSCA                                                 


     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.52858543E+00   0.51861212E+00   0.10471976E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+02   0.70000000E+02   0.60000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      10.0     -0.00      /   0.00      10.0     -0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  106/ 121/   1
     Node  distance  in   X/Y/Z :   9.973310E-03/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :  1.0000E-02 cm   (i.e.,   2.5000E-04 rad  at mean radius RM =    40.00    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  5.0064    61.982     0.000     0.000     0.000            0.519    63.301     0.036     0.000     0.000            1

 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      3)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   5.0064    61.982     0.000     0.000     0.000       0.0000    4.0064   63.301   36.484    0.000    0.000   6.540028E+01     1
               Time of flight (mus) :  2.12155179E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   6.330128E-01   3.648399E-02        1        1    1.000      (Y,T)        21
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        21
   0.0000E+00   0.0000E+00   1.0000E+00   2.121552E-02   5.000000E+00        1        1    1.000      (t,K)        21

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
      5  Keyword, label(s) :  FIT                                                   

     FIT procedure launched. Method is 1


                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  TOSCA                                                 


     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.52858543E+00   0.51861212E+00   0.10471976E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+02   0.70000000E+02   0.60000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      10.0     -0.00      /   0.00      10.0     -0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  106/ 121/   1
     Node  distance  in   X/Y/Z :   9.973310E-03/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :  1.0000E-02 cm   (i.e.,   2.5000E-04 rad  at mean radius RM =    40.00    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  5.0064    61.982     0.000     0.000     0.000            0.519    63.301     0.036     0.000     0.000            1

 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      3)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   5.0064    61.982     0.000     0.000     0.000       0.0000    4.0064   63.301   36.484    0.000    0.000   6.540028E+01     1
               Time of flight (mus) :  2.12155179E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   6.330128E-01   3.648399E-02        1        1    1.000      (Y,T)        21
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        21
   0.0000E+00   0.0000E+00   1.0000E+00   2.121552E-02   5.000000E+00        1        1    1.000      (t,K)        21

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
      5  Keyword, label(s) :  FIT                                                   

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           constraint #            1       IR =            4 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            4 ,   ok.
           constraint #            2       I  =            1 ,   ok.

                    FIT  variables  and  constraints  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save [FileName]' command

 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        64.7       64.707032       65.0      2.992E-06  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-08)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    8.403265E-07    1.39E-03 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    2.249346E-05    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   5.0666E-10

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                         ****  End  of  'REBELOTE'  procedure  ****

      There  has  been         21  passes  through  the  optical  structure 

                     Total of         21 particles have been launched

************************************************************************************************************************************
      8  Keyword, label(s) :  END                                                   


************************************************************************************************************************************
   
            ZGOUBI RUN COMPLETED. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  05-07-2017,  at  14:49:37 
  JOB  ENDED  ON    05-07-2017,  AT  14:49:38 

   CPU time, total :     1.2600000000000000     

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.


