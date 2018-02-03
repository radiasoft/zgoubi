Uniform field sector
 'OBJET'                                                                                                      1
64.62444403717985                       ! 200keV proton
2
1 1
12.9248888074 0.  0.   0.   0.  1. 'm'
1
 
 'PARTICUL'              ! This is required only because we want to get the time-of-flight                    2
938.27203D0 1.602176487D-19  1.79284735D0 0. 0.  ! otherwise zgoubi only requires rigidity.
 
 'FAISTORE'                                                                                                   3
zgoubi.fai   #End
1
 'TOSCA'                                                                                                      4
0  2
1.  1. 1. 1.
HEADER_8
629 121 1 22.1 1.       IZ=1 -> 2D ; MOD=22 -> polar map ; .MOD2=.1 -> one map file
geneSectorMap.out
0 0 0 0
2
1.
2
0. 0. 0. 0.
 
 'FAISCEAU'  #End                                                                                             5
 'END'                                                                                                        6

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
        10.000000000000000       0.50000000000000000       0.57324840764331209      -5.49755858038609913E-042      !  Rmi/cm,
      # Field map generated using geneSectorMap.f                                                                            
     # AT/rd,  AT/deg, Rmi/cm, Rma/cm, RM/cm, NR, dR/cm, NX, dX/cm, dA/rd :                                                  
     #   6.28318531E+00   3.60000000E+02   1.00000000E+01   7.00000000E+01   5.00000000E+01 121   5.00000000E-01 629   5.0025
      # For TOSCA:  629 121 1 22.1 1.  !IZ=1 -> 2D ; MOD=22 -> polar map ; .MOD2=.1 -> one map file                          
      # R*cosA (A:0->360), Z==0, R*sinA, BY, BZ, BX                                                                          
      # cm                 cm    cm      kG  kG  kG                                                                          
      #                                                                                                                      
 
R_min (cm), DR (cm), DTTA (deg), DZ (cm) :   1.000000E+01  5.000000E-01  5.732484E-01 -5.497559E-42



     Field map limits, angle :  min, max, max-min (rad) :  -0.31415927E+01   0.31415927E+01   0.62831853E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+02   0.70000000E+02   0.60000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      10.0     -0.00      /   0.00      10.0     -0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  629/ 121/   1
     Node  distance  in   X/Y/Z :   1.000507E-02/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5000E-02 rad  at mean radius RM =    40.00    )

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
      5  Keyword, label(s) :  FAISCEAU    #End                                      

0                                             TRACE DU FAISCEAU
                                           (follows element #      4)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    0.0000   12.925    0.000    0.000    0.000   8.120937E+01     1
               Time of flight (mus) :  0.13121675     mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01   2.412022E-08        1        0    0.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         1
   0.0000E+00   0.0000E+00   1.0000E+00   1.312167E-01   2.000000E-01        1        1    1.000      (t,K)         1

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

   4.057292E-37   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      6  Keyword, label(s) :  END                                                   


                             1 particles have been launched
                     Made  it  to  the  end :      1

************************************************************************************************************************************
 Pgm zgoubi : Execution ended normally, upon keyword END or FIN
   
            ZGOUBI RUN COMPLETED. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  03-02-2018,  at  12:08:30 
  JOB  ENDED  ON    03-02-2018,  AT  12:08:30 

   CPU time, total :    0.86005299999999996     
