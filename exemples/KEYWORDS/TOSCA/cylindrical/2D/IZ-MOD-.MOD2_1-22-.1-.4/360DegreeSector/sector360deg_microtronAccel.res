Uniform field sector
 'OBJET'                                                                                                      1
64.62444403717985                       ! 200keV proton
2
1 1
12.9248888074 0.  0.   0.   0.  1. 'm'
1
 
 'PARTICUL'              ! This is required only because we want to get the time-of-flight                    2
938.27203D0 1.602176487D-19  1.79284735D0 0. 0.  ! otherwise zgoubi only requires rigidity.
! 'PARTICUL'
! PROTON
 
 'FAISTORE'                                                                                                   3
zgoubi.fai   #End
1
 'TOSCA'                                                                                                      4
0  2
1.  1. 1. 1.
HEADER_8
629 151 1 22.1 1.       IZ=1 -> 2D ; MOD=22 -> polar map ; .MOD2=.1 -> one map file
geneSectorMap.out
0 0 0 0
2
1.
2
0. 0. 0. 0.
 
 'CAVITE'                                                                                                     5
3
0. 0.
100e3  1.57079632679
 
 'FAISCEAU'                                                                                                   6
 
 'REBELOTE'      ! K = 99 : coordinates at end of previous pass are used as initial                           7
60  1.1 99       ! coordinates for the next pass ; idem for spin components.
!                 Note that, (i) Y0 remains constant (due to the "microtron configuration"),
!                 (ii) is updated by passage through CAVITE
 'FAISCEAU'                                                                                                   8
 
 'FAISCEAU'  #End                                                                                             9
 'END'                                                                                                       10

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


     Particle  properties :
      
                     Mass          =    938.272        MeV/c2
                     Charge        =   1.602176E-19    C     
                     G  factor     =    1.79285              

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
                    electric rigidity (MeV) :  0.3999573775    =T[eV]*(gamma+1)/gamma, such that dev.=E*L/rigidity
  
 I, AMQ(1,I), AMQ(2,I)/QE, P/Pref, v/c, time, s :
  
     1   9.38272030E+02  1.00000000E+00  1.00000000E+00  2.06441112E-02  0.00000000E+00  0.00000000E+00

************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


                OPEN FILE zgoubi.fai                                                                      
                FOR PRINTING COORDINATES 

               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


                OPEN FILE zgoubi.plt                                                                      
                FOR PRINTING TRAJECTORIES


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
        1.0000000000000000       0.50000000000000000       0.57324840764331209       2.12199579145933796E-314      !  Rmi/cm,
      # Field map generated using geneSectorMap.f                                                                            
     # AT/rd,  AT/deg, Rmi/cm, Rma/cm, RM/cm, NR, dR/cm, NX, dX/cm, dA/rd :                                                  
     #   6.28318531E+00   3.60000000E+02   1.00000000E+00   7.60000000E+01   5.00000000E+01 151   5.00000000E-01 629   5.0025
      # For TOSCA:          629         151  1 22.1 1.  !IZ=1 -> 2D ; MOD=22 -> polar map ; .MOD2=.1 -> one map file         
      # R*cosA (A:0->360), Z==0, R*sinA, BY, BZ, BX                                                                          
      # cm                 cm    cm      kG  kG  kG                                                                          
      #                                                                                                                      
 
R_min (cm), DR (cm), DTTA (deg), DZ (cm) :   1.000000E+00  5.000000E-01  5.732484E-01  2.121996-314



     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925     0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (       82.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    0.2248   12.925    0.000    0.000    0.000   8.120937E+01     1
               Time of flight (mus) :  0.13121675     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01   1.969356E-08        1        1    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         1
   0.0000E+00   0.0000E+00   1.0000E+00   1.312167E-01   3.000000E-01        1        1    1.000      (t,K)         1

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        1 through the optical structure 

                     Total of          1 particles have been launched

     Multiple pass, 
          from element #     1 : OBJET     /label1=                    /label2=                    
                             to  REBELOTE  /label1=                    /label2=                    
     ending at pass #      61 at element #     7 : REBELOTE  /label1=                    /label2=                    


************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      101.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    0.4143   12.925   -0.000    0.000    0.000   1.806729E+02     1
               Time of flight (mus) :  0.26244757     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -7.831703E-10        1        1    1.000      (Y,T)         2
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         2
   0.0000E+00   0.0000E+00   1.0000E+00   2.624476E-01   4.000000E-01        1        1    1.000      (t,K)         2

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        2 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      116.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    0.5813   12.925   -0.000    0.000    0.000   2.955265E+02     1
               Time of flight (mus) :  0.39369240     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -3.370547E-09        1        1    1.000      (Y,T)         3
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         3
   0.0000E+00   0.0000E+00   1.0000E+00   3.936924E-01   5.000000E-01        1        1    1.000      (t,K)         3

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        3 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      130.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    0.7322   12.925   -0.000    0.000    0.000   4.239402E+02     1
               Time of flight (mus) :  0.52495122     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -8.587041E-09        1        1    1.000      (Y,T)         4
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         4
   0.0000E+00   0.0000E+00   1.0000E+00   5.249512E-01   6.000000E-01        1        1    1.000      (t,K)         4

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        4 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      142.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    0.8711   12.925   -0.000    0.000    0.000   5.646141E+02     1
               Time of flight (mus) :  0.65622404     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -1.015956E-08        1        1    1.000      (Y,T)         5
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         5
   0.0000E+00   0.0000E+00   1.0000E+00   6.562240E-01   7.000000E-01        1        1    1.000      (t,K)         5

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        5 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      153.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    1.0003   12.925   -0.000    0.000    0.000   7.165633E+02     1
               Time of flight (mus) :  0.78751084     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -9.742285E-09        1        1    1.000      (Y,T)         6
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         6
   0.0000E+00   0.0000E+00   1.0000E+00   7.875108E-01   8.000000E-01        1        1    1.000      (t,K)         6

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        6 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      164.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    1.1217   12.925   -0.000    0.000    0.000   8.790082E+02     1
               Time of flight (mus) :  0.91881163     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -1.089252E-08        1        1    1.000      (Y,T)         7
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         7
   0.0000E+00   0.0000E+00   1.0000E+00   9.188116E-01   9.000000E-01        1        1    1.000      (t,K)         7

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        7 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      173.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    1.2365   12.925   -0.000    0.000    0.000   1.051312E+03     1
               Time of flight (mus) :   1.0501264     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -9.566627E-09        1        1    1.000      (Y,T)         8
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         8
   0.0000E+00   0.0000E+00   1.0000E+00   1.050126E+00   1.000000E+00        1        1    1.000      (t,K)         8

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        8 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      183.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    1.3458   12.925   -0.000    0.000    0.000   1.232940E+03     1
               Time of flight (mus) :   1.1814552     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -9.838052E-09        1        1    1.000      (Y,T)         9
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         9
   0.0000E+00   0.0000E+00   1.0000E+00   1.181455E+00   1.100000E+00        1        1    1.000      (t,K)         9

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        9 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      192.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    1.4501   12.925   -0.000    0.000    0.000   1.423439E+03     1
               Time of flight (mus) :   1.3127979     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -1.014996E-08        1        1    1.000      (Y,T)        10
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        10
   0.0000E+00   0.0000E+00   1.0000E+00   1.312798E+00   1.200000E+00        1        1    1.000      (t,K)        10

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       10 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      200.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    1.5503   12.925   -0.000    0.000    0.000   1.622414E+03     1
               Time of flight (mus) :   1.4441546     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -9.783237E-09        1        1    1.000      (Y,T)        11
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        11
   0.0000E+00   0.0000E+00   1.0000E+00   1.444155E+00   1.300000E+00        1        1    1.000      (t,K)        11

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       11 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      208.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    1.6466   12.925   -0.000    0.000    0.000   1.829519E+03     1
               Time of flight (mus) :   1.5755253     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -9.318792E-09        1        1    1.000      (Y,T)        12
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        12
   0.0000E+00   0.0000E+00   1.0000E+00   1.575525E+00   1.400000E+00        1        1    1.000      (t,K)        12

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       12 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      216.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    1.7396   12.925   -0.000    0.000    0.000   2.044447E+03     1
               Time of flight (mus) :   1.7069100     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -9.063765E-09        1        1    1.000      (Y,T)        13
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        13
   0.0000E+00   0.0000E+00   1.0000E+00   1.706910E+00   1.500000E+00        1        1    1.000      (t,K)        13

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       13 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      224.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    1.8295   12.925   -0.000    0.000    0.000   2.266926E+03     1
               Time of flight (mus) :   1.8383087     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -9.125750E-09        1        1    1.000      (Y,T)        14
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        14
   0.0000E+00   0.0000E+00   1.0000E+00   1.838309E+00   1.600000E+00        1        1    1.000      (t,K)        14

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       14 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      231.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    1.9166   12.925   -0.000    0.000    0.000   2.496706E+03     1
               Time of flight (mus) :   1.9697213     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -8.988208E-09        1        1    1.000      (Y,T)        15
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        15
   0.0000E+00   0.0000E+00   1.0000E+00   1.969721E+00   1.700000E+00        1        1    1.000      (t,K)        15

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       15 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      238.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    2.0013   12.925   -0.000    0.000    0.000   2.733565E+03     1
               Time of flight (mus) :   2.1011480     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -8.809769E-09        1        1    1.000      (Y,T)        16
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        16
   0.0000E+00   0.0000E+00   1.0000E+00   2.101148E+00   1.800000E+00        1        1    1.000      (t,K)        16

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       16 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      245.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    2.0836   12.925   -0.000    0.000    0.000   2.977298E+03     1
               Time of flight (mus) :   2.2325886     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -8.697546E-09        1        1    1.000      (Y,T)        17
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        17
   0.0000E+00   0.0000E+00   1.0000E+00   2.232589E+00   1.900000E+00        1        1    1.000      (t,K)        17

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       17 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      252.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    2.1638   12.925   -0.000    0.000    0.000   3.227715E+03     1
               Time of flight (mus) :   2.3640432     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -8.685164E-09        1        1    1.000      (Y,T)        18
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        18
   0.0000E+00   0.0000E+00   1.0000E+00   2.364043E+00   2.000000E+00        1        1    1.000      (t,K)        18

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       18 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      258.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    2.2420   12.925   -0.000    0.000    0.000   3.484645E+03     1
               Time of flight (mus) :   2.4955117     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -8.500478E-09        1        1    1.000      (Y,T)        19
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        19
   0.0000E+00   0.0000E+00   1.0000E+00   2.495512E+00   2.100000E+00        1        1    1.000      (t,K)        19

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       19 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      264.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    2.3184   12.925   -0.000    0.000    0.000   3.747927E+03     1
               Time of flight (mus) :   2.6269943     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -8.216542E-09        1        1    1.000      (Y,T)        20
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        20
   0.0000E+00   0.0000E+00   1.0000E+00   2.626994E+00   2.200000E+00        1        1    1.000      (t,K)        20

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       20 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      271.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    2.3931   12.925   -0.000    0.000    0.000   4.017412E+03     1
               Time of flight (mus) :   2.7584909     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -8.170208E-09        1        1    1.000      (Y,T)        21
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        21
   0.0000E+00   0.0000E+00   1.0000E+00   2.758491E+00   2.300000E+00        1        1    1.000      (t,K)        21

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       21 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      277.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    2.4661   12.925   -0.000    0.000    0.000   4.292961E+03     1
               Time of flight (mus) :   2.8900014     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -8.106736E-09        1        1    1.000      (Y,T)        22
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        22
   0.0000E+00   0.0000E+00   1.0000E+00   2.890001E+00   2.400000E+00        1        1    1.000      (t,K)        22

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       22 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.32 X =   -0.11 T = 0.081 Trajectory #       1 went out of field region.

    Y =   76.19 X =    0.13 T =-0.095 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      283.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    2.5377   12.925   -0.000    0.000    0.000   4.574444E+03     1
               Time of flight (mus) :   3.0215259     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -8.051415E-09        1        1    1.000      (Y,T)        23
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        23
   0.0000E+00   0.0000E+00   1.0000E+00   3.021526E+00   2.500000E+00        1        1    1.000      (t,K)        23

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       23 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.39 X =   -0.28 T = 0.197 Trajectory #       1 went out of field region.

    Y =   76.05 X =    0.30 T =-0.212 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      288.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    2.6079   12.925   -0.000    0.000    0.000   4.861738E+03     1
               Time of flight (mus) :   3.1530644     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -7.835065E-09        1        1    1.000      (Y,T)        24
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        24
   0.0000E+00   0.0000E+00   1.0000E+00   3.153064E+00   2.600000E+00        1        1    1.000      (t,K)        24

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       24 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.46 X =   -0.37 T = 0.264 Trajectory #       1 went out of field region.

    Y =   76.19 X =    0.38 T =-0.273 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      294.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    2.6767   12.925   -0.000    0.000    0.000   5.154730E+03     1
               Time of flight (mus) :   3.2846169     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -7.689869E-09        1        1    1.000      (Y,T)        25
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        25
   0.0000E+00   0.0000E+00   1.0000E+00   3.284617E+00   2.700000E+00        1        1    1.000      (t,K)        25

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       25 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.35 X =   -0.44 T = 0.319 Trajectory #       1 went out of field region.

    Y =   76.22 X =    0.45 T =-0.322 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      300.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    2.7442   12.925   -0.000    0.000    0.000   5.453312E+03     1
               Time of flight (mus) :   3.4161833     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -7.619619E-09        1        1    1.000      (Y,T)        26
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        26
   0.0000E+00   0.0000E+00   1.0000E+00   3.416183E+00   2.800000E+00        1        1    1.000      (t,K)        26

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       26 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.27 X =   -0.51 T = 0.363 Trajectory #       1 went out of field region.

    Y =   75.94 X =    0.52 T =-0.371 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      305.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    2.8106   12.925   -0.000    0.000    0.000   5.757380E+03     1
               Time of flight (mus) :   3.5477637     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -7.476369E-09        1        1    1.000      (Y,T)        27
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        27
   0.0000E+00   0.0000E+00   1.0000E+00   3.547764E+00   2.900000E+00        1        1    1.000      (t,K)        27

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       27 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.36 X =   -0.55 T = 0.397 Trajectory #       1 went out of field region.

    Y =   76.15 X =    0.56 T =-0.401 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      311.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    2.8759   12.925   -0.000    0.000    0.000   6.066839E+03     1
               Time of flight (mus) :   3.6793582     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -7.420487E-09        1        1    1.000      (Y,T)        28
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        28
   0.0000E+00   0.0000E+00   1.0000E+00   3.679358E+00   3.000000E+00        1        1    1.000      (t,K)        28

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       28 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.30 X =   -0.60 T = 0.430 Trajectory #       1 went out of field region.

    Y =   76.20 X =    0.60 T =-0.432 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      316.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    2.9400   12.925   -0.000    0.000    0.000   6.381596E+03     1
               Time of flight (mus) :   3.8109666     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -7.331025E-09        1        1    1.000      (Y,T)        29
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        29
   0.0000E+00   0.0000E+00   1.0000E+00   3.810967E+00   3.100000E+00        1        1    1.000      (t,K)        29

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       29 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.57 X =   -0.63 T = 0.453 Trajectory #       1 went out of field region.

    Y =   76.11 X =    0.64 T =-0.461 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      321.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.0032   12.925   -0.000    0.000    0.000   6.701565E+03     1
               Time of flight (mus) :   3.9425889     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -7.219096E-09        1        1    1.000      (Y,T)        30
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        30
   0.0000E+00   0.0000E+00   1.0000E+00   3.942589E+00   3.200000E+00        1        1    1.000      (t,K)        30

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       30 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.30 X =   -0.67 T = 0.483 Trajectory #       1 went out of field region.

    Y =   75.88 X =    0.68 T =-0.490 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      326.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.0654   12.925   -0.000    0.000    0.000   7.026663E+03     1
               Time of flight (mus) :   4.0742253     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -7.097913E-09        1        1    1.000      (Y,T)        31
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        31
   0.0000E+00   0.0000E+00   1.0000E+00   4.074225E+00   3.300000E+00        1        1    1.000      (t,K)        31

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       31 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.41 X =   -0.70 T = 0.505 Trajectory #       1 went out of field region.

    Y =   75.99 X =    0.71 T =-0.511 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      331.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.1266   12.925   -0.000    0.000    0.000   7.356809E+03     1
               Time of flight (mus) :   4.2058757     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -6.977535E-09        1        1    1.000      (Y,T)        32
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        32
   0.0000E+00   0.0000E+00   1.0000E+00   4.205876E+00   3.400000E+00        1        1    1.000      (t,K)        32

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       32 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.46 X =   -0.72 T = 0.525 Trajectory #       1 went out of field region.

    Y =   76.01 X =    0.73 T =-0.532 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      336.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.1870   12.925   -0.000    0.000    0.000   7.691930E+03     1
               Time of flight (mus) :   4.3375400     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -6.865485E-09        1        1    1.000      (Y,T)        33
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        33
   0.0000E+00   0.0000E+00   1.0000E+00   4.337540E+00   3.500000E+00        1        1    1.000      (t,K)        33

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       33 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.45 X =   -0.75 T = 0.545 Trajectory #       1 went out of field region.

    Y =   75.94 X =    0.76 T =-0.552 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      341.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.2465   12.925   -0.000    0.000    0.000   8.031952E+03     1
               Time of flight (mus) :   4.4692183     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -6.767320E-09        1        1    1.000      (Y,T)        34
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        34
   0.0000E+00   0.0000E+00   1.0000E+00   4.469218E+00   3.600000E+00        1        1    1.000      (t,K)        34

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       34 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.40 X =   -0.77 T = 0.564 Trajectory #       1 went out of field region.

    Y =   75.78 X =    0.79 T =-0.572 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      346.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.3052   12.925   -0.000    0.000    0.000   8.376807E+03     1
               Time of flight (mus) :   4.6009106     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -6.687063E-09        1        1    1.000      (Y,T)        35
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        35
   0.0000E+00   0.0000E+00   1.0000E+00   4.600911E+00   3.700000E+00        1        1    1.000      (t,K)        35

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       35 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.29 X =   -0.80 T = 0.582 Trajectory #       1 went out of field region.

    Y =   76.09 X =    0.80 T =-0.585 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      351.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.3631   12.925   -0.000    0.000    0.000   8.726427E+03     1
               Time of flight (mus) :   4.7326169     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -6.626943E-09        1        1    1.000      (Y,T)        36
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        36
   0.0000E+00   0.0000E+00   1.0000E+00   4.732617E+00   3.800000E+00        1        1    1.000      (t,K)        36

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       36 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.71 X =   -0.81 T = 0.593 Trajectory #       1 went out of field region.

    Y =   75.77 X =    0.83 T =-0.604 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      355.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.4202   12.925   -0.000    0.000    0.000   9.080751E+03     1
               Time of flight (mus) :   4.8643371     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -6.518202E-09        1        1    1.000      (Y,T)        37
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        37
   0.0000E+00   0.0000E+00   1.0000E+00   4.864337E+00   3.900000E+00        1        1    1.000      (t,K)        37

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       37 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.54 X =   -0.83 T = 0.610 Trajectory #       1 went out of field region.

    Y =   75.94 X =    0.84 T =-0.617 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      360.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.4767   12.925   -0.000    0.000    0.000   9.439715E+03     1
               Time of flight (mus) :   4.9960713     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -6.438428E-09        1        1    1.000      (Y,T)        38
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        38
   0.0000E+00   0.0000E+00   1.0000E+00   4.996071E+00   4.000000E+00        1        1    1.000      (t,K)        38

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       38 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.33 X =   -0.86 T = 0.627 Trajectory #       1 went out of field region.

    Y =   76.06 X =    0.86 T =-0.630 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      365.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.5324   12.925   -0.000    0.000    0.000   9.803263E+03     1
               Time of flight (mus) :   5.1278196     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -6.386784E-09        1        1    1.000      (Y,T)        39
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        39
   0.0000E+00   0.0000E+00   1.0000E+00   5.127820E+00   4.100000E+00        1        1    1.000      (t,K)        39

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       39 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.68 X =   -0.87 T = 0.636 Trajectory #       1 went out of field region.

    Y =   76.12 X =    0.88 T =-0.642 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      369.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.5875   12.925   -0.000    0.000    0.000   1.017134E+04     1
               Time of flight (mus) :   5.2595818     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -6.305785E-09        1        1    1.000      (Y,T)        40
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        40
   0.0000E+00   0.0000E+00   1.0000E+00   5.259582E+00   4.200000E+00        1        1    1.000      (t,K)        40

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       40 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.41 X =   -0.89 T = 0.652 Trajectory #       1 went out of field region.

    Y =   76.13 X =    0.89 T =-0.655 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      374.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.6419   12.925   -0.000    0.000    0.000   1.054388E+04     1
               Time of flight (mus) :   5.3913580     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -6.256348E-09        1        1    1.000      (Y,T)        41
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        41
   0.0000E+00   0.0000E+00   1.0000E+00   5.391358E+00   4.300000E+00        1        1    1.000      (t,K)        41

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       41 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.72 X =   -0.90 T = 0.660 Trajectory #       1 went out of field region.

    Y =   76.09 X =    0.91 T =-0.667 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      378.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.6957   12.925   -0.000    0.000    0.000   1.092084E+04     1
               Time of flight (mus) :   5.5231481     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -6.186589E-09        1        1    1.000      (Y,T)        42
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        42
   0.0000E+00   0.0000E+00   1.0000E+00   5.523148E+00   4.400000E+00        1        1    1.000      (t,K)        42

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       42 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.40 X =   -0.92 T = 0.675 Trajectory #       1 went out of field region.

    Y =   75.99 X =    0.92 T =-0.679 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      382.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.7488   12.925   -0.000    0.000    0.000   1.130218E+04     1
               Time of flight (mus) :   5.6549523     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -6.100883E-09        1        1    1.000      (Y,T)        43
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        43
   0.0000E+00   0.0000E+00   1.0000E+00   5.654952E+00   4.500000E+00        1        1    1.000      (t,K)        43

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       43 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.69 X =   -0.93 T = 0.683 Trajectory #       1 went out of field region.

    Y =   75.83 X =    0.94 T =-0.692 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      387.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.8015   12.925   -0.000    0.000    0.000   1.168783E+04     1
               Time of flight (mus) :   5.7867704     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -6.050292E-09        1        1    1.000      (Y,T)        44
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        44
   0.0000E+00   0.0000E+00   1.0000E+00   5.786770E+00   4.600000E+00        1        1    1.000      (t,K)        44

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       44 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.32 X =   -0.95 T = 0.697 Trajectory #       1 went out of field region.

    Y =   75.63 X =    0.96 T =-0.704 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      391.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.8535   12.925   -0.000    0.000    0.000   1.207775E+04     1
               Time of flight (mus) :   5.9186025     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.988823E-09        1        1    1.000      (Y,T)        45
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        45
   0.0000E+00   0.0000E+00   1.0000E+00   5.918603E+00   4.700000E+00        1        1    1.000      (t,K)        45

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       45 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.58 X =   -0.95 T = 0.705 Trajectory #       1 went out of field region.

    Y =   76.03 X =    0.96 T =-0.710 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      395.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.9050   12.925   -0.000    0.000    0.000   1.247190E+04     1
               Time of flight (mus) :   6.0504486     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.919626E-09        1        1    1.000      (Y,T)        46
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        46
   0.0000E+00   0.0000E+00   1.0000E+00   6.050449E+00   4.800000E+00        1        1    1.000      (t,K)        46

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       46 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.83 X =   -0.96 T = 0.712 Trajectory #       1 went out of field region.

    Y =   75.74 X =    0.98 T =-0.722 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      399.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    3.9559   12.925   -0.000    0.000    0.000   1.287023E+04     1
               Time of flight (mus) :   6.1823087     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.845592E-09        1        1    1.000      (Y,T)        47
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        47
   0.0000E+00   0.0000E+00   1.0000E+00   6.182309E+00   4.900000E+00        1        1    1.000      (t,K)        47

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       47 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.41 X =   -0.98 T = 0.725 Trajectory #       1 went out of field region.

    Y =   76.06 X =    0.99 T =-0.728 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      404.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.0064   12.925   -0.000    0.000    0.000   1.327270E+04     1
               Time of flight (mus) :   6.3141827     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.805392E-09        1        1    1.000      (Y,T)        48
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        48
   0.0000E+00   0.0000E+00   1.0000E+00   6.314183E+00   5.000000E+00        1        1    1.000      (t,K)        48

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       48 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.64 X =   -0.99 T = 0.732 Trajectory #       1 went out of field region.

    Y =   75.68 X =    1.00 T =-0.740 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      408.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.0563   12.925   -0.000    0.000    0.000   1.367927E+04     1
               Time of flight (mus) :   6.4460708     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.763260E-09        1        1    1.000      (Y,T)        49
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        49
   0.0000E+00   0.0000E+00   1.0000E+00   6.446071E+00   5.100000E+00        1        1    1.000      (t,K)        49

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       49 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.86 X =   -1.00 T = 0.738 Trajectory #       1 went out of field region.

    Y =   75.93 X =    1.01 T =-0.746 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      412.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.1058   12.925   -0.000    0.000    0.000   1.408989E+04     1
               Time of flight (mus) :   6.5779728     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.720220E-09        1        1    1.000      (Y,T)        50
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        50
   0.0000E+00   0.0000E+00   1.0000E+00   6.577973E+00   5.200000E+00        1        1    1.000      (t,K)        50

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       50 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.39 X =   -1.01 T = 0.750 Trajectory #       1 went out of field region.

    Y =   76.15 X =    1.02 T =-0.752 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      416.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.1548   12.925   -0.000    0.000    0.000   1.450453E+04     1
               Time of flight (mus) :   6.7098888     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.677522E-09        1        1    1.000      (Y,T)        51
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        51
   0.0000E+00   0.0000E+00   1.0000E+00   6.709889E+00   5.300000E+00        1        1    1.000      (t,K)        51

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       51 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.59 X =   -1.02 T = 0.757 Trajectory #       1 went out of field region.

    Y =   75.64 X =    1.03 T =-0.764 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      420.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.2033   12.925   -0.000    0.000    0.000   1.492315E+04     1
               Time of flight (mus) :   6.8418188     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.636288E-09        1        1    1.000      (Y,T)        52
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        52
   0.0000E+00   0.0000E+00   1.0000E+00   6.841819E+00   5.400000E+00        1        1    1.000      (t,K)        52

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       52 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.79 X =   -1.03 T = 0.763 Trajectory #       1 went out of field region.

    Y =   75.79 X =    1.04 T =-0.770 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      424.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.2514   12.925   -0.000    0.000    0.000   1.534571E+04     1
               Time of flight (mus) :   6.9737628     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.597355E-09        1        1    1.000      (Y,T)        53
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        53
   0.0000E+00   0.0000E+00   1.0000E+00   6.973763E+00   5.500000E+00        1        1    1.000      (t,K)        53

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       53 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.28 X =   -1.04 T = 0.774 Trajectory #       1 went out of field region.

    Y =   75.91 X =    1.05 T =-0.777 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      428.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.2991   12.925   -0.000    0.000    0.000   1.577218E+04     1
               Time of flight (mus) :   7.1057207     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.561067E-09        1        1    1.000      (Y,T)        54
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        54
   0.0000E+00   0.0000E+00   1.0000E+00   7.105721E+00   5.600000E+00        1        1    1.000      (t,K)        54

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       54 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.46 X =   -1.05 T = 0.780 Trajectory #       1 went out of field region.

    Y =   75.99 X =    1.05 T =-0.783 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      432.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.3464   12.925   -0.000    0.000    0.000   1.620251E+04     1
               Time of flight (mus) :   7.2376926     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.526579E-09        1        1    1.000      (Y,T)        55
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        55
   0.0000E+00   0.0000E+00   1.0000E+00   7.237693E+00   5.700000E+00        1        1    1.000      (t,K)        55

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       55 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.63 X =   -1.05 T = 0.785 Trajectory #       1 went out of field region.

    Y =   76.05 X =    1.06 T =-0.789 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      435.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.3932   12.925   -0.000    0.000    0.000   1.663669E+04     1
               Time of flight (mus) :   7.3696786     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.474005E-09        1        1    1.000      (Y,T)        56
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        56
   0.0000E+00   0.0000E+00   1.0000E+00   7.369679E+00   5.800000E+00        1        1    1.000      (t,K)        56

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       56 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.80 X =   -1.06 T = 0.790 Trajectory #       1 went out of field region.

    Y =   76.07 X =    1.07 T =-0.796 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      439.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.4396   12.925   -0.000    0.000    0.000   1.707467E+04     1
               Time of flight (mus) :   7.5016785     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.427819E-09        1        1    1.000      (Y,T)        57
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        57
   0.0000E+00   0.0000E+00   1.0000E+00   7.501678E+00   5.900000E+00        1        1    1.000      (t,K)        57

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       57 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.96 X =   -1.07 T = 0.796 Trajectory #       1 went out of field region.

    Y =   76.07 X =    1.08 T =-0.802 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      443.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.4857   12.925   -0.000    0.000    0.000   1.751642E+04     1
               Time of flight (mus) :   7.6336923     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.388253E-09        1        1    1.000      (Y,T)        58
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        58
   0.0000E+00   0.0000E+00   1.0000E+00   7.633692E+00   6.000000E+00        1        1    1.000      (t,K)        58

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       58 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.40 X =   -1.08 T = 0.806 Trajectory #       1 went out of field region.

    Y =   76.03 X =    1.09 T =-0.808 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      447.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.5314   12.925   -0.000    0.000    0.000   1.796191E+04     1
               Time of flight (mus) :   7.7657202     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.354749E-09        1        1    1.000      (Y,T)        59
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        59
   0.0000E+00   0.0000E+00   1.0000E+00   7.765720E+00   6.100000E+00        1        1    1.000      (t,K)        59

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       59 through the optical structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.55 X =   -1.09 T = 0.811 Trajectory #       1 went out of field region.

    Y =   75.97 X =    1.09 T =-0.815 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      450.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.5767   12.925   -0.000    0.000    0.000   1.841111E+04     1
               Time of flight (mus) :   7.8977620     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.306868E-09        1        1    1.000      (Y,T)        60
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        60
   0.0000E+00   0.0000E+00   1.0000E+00   7.897762E+00   6.200000E+00        1        1    1.000      (t,K)        60

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
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       60 through the optical structure 

                     Total of          1 particles have been launched


      Next  pass  is  #    61 and  last  pass  through  the  optical  structure


************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 



               Final  coordinates  of  previous  run  taken  as  initial  coordinates ;         1 particles

************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 
                    #End                

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                                                 


     zgoubi.plt                                                                      
      already open...

     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+01   0.76000000E+02   0.75000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 151/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad


    Y =   76.70 X =   -1.09 T = 0.816 Trajectory #       1 went out of field region.

    Y =   75.88 X =    1.10 T =-0.821 Trajectory #       1 came back in. field region.
  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925    -0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (      454.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  CAVITE                                                

                Accelerating cavity. Type is :   OPTION 3  


                    Synchronous  phase                 =     0.1571E+01 rad
                    Synchronous energy  gain           =     0.1000E+00 MeV
                    Cumulated  distance  from  origin  =     0.0000E+00 m
                    Synchronous  time                  =     0.0000E+00 s
                    Particle mass           =    0.93827E+03 MeV/c2
                             charge         =   0.160218E-18 C

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.6216   12.925   -0.000    0.000    0.000   1.886398E+04     1
               Time of flight (mus) :   8.0298179     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.267047E-09        1        1    1.000      (Y,T)        61
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        61
   0.0000E+00   0.0000E+00   1.0000E+00   8.029818E+00   6.300000E+00        1        1    1.000      (t,K)        61

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
      7  Keyword, label(s) :  REBELOTE                                              


                         ****  End  of  'REBELOTE'  procedure  ****

      There  has  been         61  passes  through  the  optical  structure 

                     Total of          1 particles have been launched

************************************************************************************************************************************
      8  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      7)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.6216   12.925   -0.000    0.000    0.000   1.886398E+04     1
               Time of flight (mus) :   8.0298179     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.267047E-09        1        1    1.000      (Y,T)        62
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        62
   0.0000E+00   0.0000E+00   1.0000E+00   8.029818E+00   6.300000E+00        1        1    1.000      (t,K)        62

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
      9  Keyword, label(s) :  FAISCEAU    #End                                      

0                                             TRACE DU FAISCEAU
                                           (follows element #      8)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    4.6216   12.925   -0.000    0.000    0.000   1.886398E+04     1
               Time of flight (mus) :   8.0298179     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01  -5.267047E-09        1        1    1.000      (Y,T)        62
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        62
   0.0000E+00   0.0000E+00   1.0000E+00   8.029818E+00   6.300000E+00        1        1    1.000      (t,K)        62

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
     10  Keyword, label(s) :  END                                                   


************************************************************************************************************************************
 Pgm zgoubi : Execution ended normally, upon keyword END or FIN
   
            ZGOUBI RUN COMPLETED. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  01-03-2018,  at  13:58:25 
  JOB  ENDED  ON    01-03-2018,  AT  13:58:26 

   CPU time, total :    0.97199699999999989     
