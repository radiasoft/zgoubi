#!/usr/local/bin python2.7
from string import *
import commands
import sys
from numpy import matrix
from numpy import linalg
from numpy import matlib
from numpy import math
from copy import deepcopy

#def factorial(n):
#    fac=1
#    for i in range(n):
#        fac=fac*(i+1)
#    return fac

#-------------------------------------------------------------
#   this program to read the magnet information into and 
#   create an elegant input files 
#   Written by Fanglei Lin, Jefferson Laboratory, 2019
#-------------------------------------------------------------

print "-----------------------------------------"
print "This code convert madx output to zgoubi"
print "It accepts: "
print "SBEND, RBEND, QUAD, SEXT, MULT, SOLENOID, DRIFT, MARKER, MATRIX, "
print "CAVITY, MONITOR, BPM, HKICKER, VKICKER, HMONITOR, VMONITOR."
print "-----------------------------------------"
print "To run the program: "
print "0) input file from MADX: madzg.in"
print "   output file         : zgoubi.dat" 
print "1) create parameters_input with madx: "
print "use, period=...;  "
print "select,flag=twiss,clear;"
print "select,flag=twiss,column=NAME,KEYWORD,S,L,ANGLE,tilt,E1,E2,h1,h2,hkic,vkic,"
print "K1L,K2L,K3L,K4L,KSI,betx,bety,alfx,alfy,dx,dpx,dy,dpy,mux,muy,x,y,"
print "RE11,RE12,RE13,RE14,RE15,RE16,RE21,RE22,RE23,RE24,RE25,RE26,"
print "RE31,RE32,RE33,RE34,RE35,RE36,RE41,RE42,RE43,RE44,RE45,RE46,"
print "RE51,RE52,RE53,RE54,RE55,RE56,RE61,RE62,RE63,RE64,RE65,RE66;"
print "twiss,rmatrix,betx=1,alfx=0,bety=1,alfy=0,file=madzg.in;"
print "stop; "
print "------------------------------------------"

#----read into the whole file
f1name="./madzg.in"
f1=open(f1name,"r")
holder=f1.readlines()
length=len(holder)
f1.close()

#----write into a ZGOUBI input file
f2name="./zgoubi.dat"
f2=open(f2name,"w")
#f2.write('Q: CHARGE,TOTAL=16e-09;        \n')
#f2.write('!! ZGOUBI input for MEIC \n')
#f2.write('   \n')

temp=holder[7]
temp=temp.split()
pgevoc=float(temp[3])
brho=pgevoc/0.299792458
temp=holder[8]
temp=temp.split()
gamma=float(temp[3])
ge=1.15965218076e-3
temp=holder[17]
temp=temp.split()
circ=float(temp[3])

f2.write("Generated by MADX  -> Zgoubi translator \n") 
f2.write("'OBJET'    \n")
f2.write("  "+str(brho*1000.0)+"    "+"! reference rigidity (kG.cm),"+"   "+" momentum p (GeV/c) = "+str(pgevoc)+",   "+"G.gamma = "+str(ge*gamma)+".\n")

kobj=2
imax=1
idmax=1
coy=0.0
cot=0.0
coz=0.0
cop=0.0
cos=0.0
cod=1.0
iex=1
f2.write("  "+str(kobj)+"\n")  
f2.write("  "+str(imax)+"   "+str(idmax)+"\n")  # total number of particles and number of distinct momenta
f2.write("  "+str(coy)+"  "+str(cot)+"  "+str(coz)+"  "+str(cop)+"  "+str(cos)+"  "+str(cod)+"  "+"'o'"+"\n")  # coordinates and tagging of the imax particles
f2.write("  "+str(iex)+"\n") 

f2.write("'PARTICUL'    \n") 
f2.write("  "+"0.51099892   1.60217653e-19   1.15965218076e-3   0.0   0.0"+"\n")  # electron mass, charge (positive), gyromagnetic factor, COM life-time, unused

f2.write("'SRLOSS'    \n") 
ksr=0     # switch off 0 or on 1 radiation loss. If 1.1, output into zgoubi.SRLOSS.out 
f2.write("  "+str(ksr)+"   !  .srloss"+"\n")  # electron mass, charge, gyromagnetic factor, COM life-time, unused
f2.write("  "+"BEND"+"\n")
f2.write("  "+"1"+"   "+"123456"+"\n")

f2.write("'SPNTRK'"+"\n") 
f2.write("  "+"3"+"\n")

f2.write("'FAISCEAU'"+"\n")   # print particle coordinates at the location where the keyword is introduced
f2.write("'FAISTORE'"+"\n")
ip=1
f2.write("  "+"zgoubi.fai"+"\n")
f2.write("  "+str(ip)+"\n")  # store every ip other pass
f2.write("'FAISCEAU'"+"\n") 

f2.write("'SCALING'"+"\n")
f2.write("  "+"1"+"   "+"3"+"\n")
f2.write("  "+"BEND"+"\n")
f2.write("  "+"2"+"\n")
f2.write("  "+str(brho)+"  "+str(brho)+"\n")
f2.write("  "+"1"+"  "+"999999"+"\n")
f2.write("  "+"MULTIPOL"+"\n")
f2.write("  "+"2"+"\n")
f2.write("  "+str(brho)+"  "+str(brho)+"\n")
f2.write("  "+"1"+"  "+"999999"+"\n")
f2.write("  "+"SOLENOID"+"\n")
f2.write("  "+"-1"+"\n")
f2.write("  "+"1"+"\n")
f2.write("  "+"1"+"\n")
#f2.write("  "+"SOLENOID"+"\n")
#f2.write("  "+"2"+"\n")
#f2.write("  "+str(brho)+"  "+str(brho)+"\n")
#f2.write("  "+"1"+"  "+"999999"+"\n")

f2.write("'OPTIONS'"+"\n")
f2.write("  "+"1"+"   "+"1"+"\n")
f2.write("  "+"WRITE ON"+"\n")

f2.write("'MARKER'"+"    "+"MARK"+"    "+"COMPLETE_RING$START"+"\n")

x2=2.0          #---radius at solenoids
x10=10.0        #---radius at pole tip: cm
total_l=0.
total_ang=0
nameseq=[]      #---the element names in the beam line 
namelist=[]     #---individual different elements
rmf=matlib.eye(6)

for i in range(length-49):        # All counting starts from 0
    temp=holder[i+48]
    temp=temp.split()
    name=temp[0][1:(len(temp[0])-1)]
    type=temp[1][1:(len(temp[1])-1)]   
    s=temp[2]
    l=temp[3]
    total_l=total_l+float(l)
    angle=temp[4]
    total_ang=total_ang+float(angle)
    tilt=temp[5]
    e1=temp[6]
    e2=temp[7]
    h1=temp[8]
    h2=temp[9]
    hkic=temp[10]
    vkic=temp[11]
    k1l=temp[12]
    k2l=temp[13]
    k3l=temp[14]
    k4l=temp[15]
    ksi=temp[16]
    betx=temp[17]
    bety=temp[18]
    alfx=temp[19]
    alfy=temp[20]
    dx=temp[21]
    dpx=temp[22]
    dy=temp[23]
    dpy=temp[24]
    mux=temp[25]
    muy=temp[26]
    x=temp[27]
    y=temp[28]
    rmi=deepcopy(rmf)
    for i in range(6):
        for j in range(6):
            rmf[i,j]=temp[29+i*6+j]
    rm=rmf*rmi.I
    nameseq.append(name)
    if name not in namelist:
         namelist.append(name)
    if type=="DRIFT":
       temp_write="'DRIFT'" +"    "+"DRIF"+"      "+name+"\n" 
       f2.write(temp_write)
       temp_write="  "+str(float(l)*100.0)+"\n"
       f2.write(temp_write)
    elif type=="MONITOR":
       temp_write="'DRIFT'" +"    "+"MONI"+"      "+name+"\n"
       f2.write(temp_write)
       temp_write="  "+str(float(l)*100.0)+"\n"
       f2.write(temp_write)
    elif type=="KICKER":
       temp_write="'DRIFT'" +"    "+"KICK"+"      "+name+"\n"
       f2.write(temp_write)
       temp_write="  "+str(float(l)*100.0)+"\n"
       f2.write(temp_write)
#       temp_write="'MULTIPOL'" +"    "+"KICK"+"      "+name+"\n"
#       f2.write(temp_write)
#       temp_write="0"+"  "+".Kicker"+"\n"
#       f2.write(temp_write)
##       temp_write="  "+str(float(l)/float(angle)*2.0*math.sin(float(angle)/2.0)*100.0)+"  "+str(x10)+"  "+str(float(angle)/float(l)*10.0)+"  "+str(float(k1l)/float(l)*0.1*10.0)+"  "+str(float(k2l)/float(l)*(x10/100.0)*(x10/100.0)/2.0*10.0)+"  0.0  0.0  0.0  0.0  0.0  0.0  0.0   "+"\n"
#       temp_write="  "+str(float(l)*100.0)+"  "+str(x10)+"  "+str(float(angle)/float(l)*10.0)+"  "+str(float(k1l)/float(l)*0.1*10.0)+"  "+str(float(k2l)/float(l)*(x10/100.0)*(x10/100.0)/2.0*10.0)+"  0.0  0.0  0.0  0.0  0.0  0.0  0.0   "+"\n"
#       f2.write(temp_write)
#       temp_write="  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0"+"\n"
#       f2.write(temp_write)
#       temp_write="  4  .1455   2.2670  -.6395  1.1558  0. 0.  0.  \n"
#       f2.write(temp_write)
#       temp_write="  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0"+"\n"
#       f2.write(temp_write)
#       temp_write="  4  .1455   2.2670  -.6395  1.1558  0. 0.  0.  \n"
#       f2.write(temp_write)
#       temp_write="  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  \n"
#       f2.write(temp_write)
##       temp_write="#20|4|20"+"  "+"Kick"+"\n"
#       temp_write="  1.  cm"+"\n"
#       f2.write(temp_write)
#       temp_write="  1  0.0  0.0  0.0  "+"\n"
#       f2.write(temp_write)
    elif type=="SBEND":
       temp_write="'BEND'"+"     "+"SBEN"+"      "+name+"\n"
       f2.write(temp_write)
       temp_write="0"+"  "+".Bend"+"\n"
       f2.write(temp_write)
       temp_write="  "+str(2.0*float(l)/float(angle)*math.sin(float(angle)/2.0)*100.0)+"  "+tilt+"  "+str(float(angle)/float(l)*10.0)+"\n"
       f2.write(temp_write)
       temp_write="  0.0  0.0  0.0"+"\n"
       f2.write(temp_write)
       temp_write="  4 .2401  1.8639  -.5572  .3904 0. 0. 0.  \n"
       f2.write(temp_write)
       temp_write="  0.0  0.0  0.0"+"\n"
       f2.write(temp_write)
       temp_write="  4 .2401  1.8639  -.5572  .3904 0. 0. 0.  \n"
       f2.write(temp_write)
#       temp_write="  #200|"+str(int(float(l)*100.0))+"|200"+"  "+"Bend"+"  "+name+"\n"
       temp_write="  1.  cm"+"\n"
       f2.write(temp_write)
       temp_write="  3  0.0  0.0  "+str(-0.5*float(angle))+"\n"
#       temp_write="  3  0.0  0.0  0.0"+"\n"
       f2.write(temp_write)
    elif type=="RBEND":
       temp_write="'MULTIPOL'"+"    "+"RBEN"+"      "+name+"\n"
       f2.write(temp_write)
       temp_write="0"+"  "+".Dip"+"\n"
       f2.write(temp_write)
       temp_write="  "+str(float(l)/float(angle)*2.0*math.sin(float(angle)/2.0)*100.0)+"  "+str(x10)+"  "+str(float(angle)/float(l)*10.0)+"  "+str(float(k1l)/float(l)*0.1*10.0)+"  "+str(float(k2l)/float(l)*(x10/100.0)*(x10/100.0)/2.0*10.0)+"  0.0  0.0  0.0  0.0  0.0  0.0  0.0   "+"\n"
       f2.write(temp_write)
       temp_write="  0.0  0.0  10.0  4.0  0.8  0.0  0.0  0.0  0.0  0.0  0.0  \n"
       f2.write(temp_write)
       temp_write="  4  .1455   2.2670  -.6395  1.1558  0. 0.   \n"
       f2.write(temp_write)
       temp_write="  0.0  0.0  10.0  4.0  0.8  0.0  0.0  0.0  0.0  0.0  0.0  \n"
       f2.write(temp_write)
       temp_write="  4  .1455   2.2670  -.6395  1.1558  0. 0.   \n"
       f2.write(temp_write)
       temp_write="  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  \n"
       f2.write(temp_write)
#       temp_write="  #30|"+str(int(float(l)*100.0))+"|30"+"  "+"Dip"+"  "+name+"\n"
       temp_write="  1.  cm"+"\n"
       f2.write(temp_write)
       temp_write="  3  0.0  0.0  "+str(-0.5*float(angle))+"\n"
#       temp_write="  3  0.0  0.0  0.0"+"\n"
       f2.write(temp_write)
    elif type=="QUADRUPOLE" and float(l)!=0.:
       temp_write="'MULTIPOL'"+"    "+"QUAD"+"      "+name+"\n"
       f2.write(temp_write)
       temp_write="0"+"  "+".Quad"+"\n"
       f2.write(temp_write)
       temp_write="  "+str(float(l)*100.0)+"  "+str(x10)+"  "+str(float(angle)/float(l)*10.0)+"  "+str(float(k1l)/float(l)*0.1*10.0)+"  "+str(float(k2l)/float(l)*(x10/100.0)*(x10/100.0)/2.0*10.0)+"  0.0  0.0  0.0  0.0  0.0  0.0  0.0   "+"\n"
       f2.write(temp_write)
       temp_write="  0.0  0.0  6.00  3.00 1.00 0.00 0.00 0.00 0.00 0. 0.  \n"
       f2.write(temp_write)
       temp_write="  6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723   \n"
       f2.write(temp_write)
       temp_write="  0.0  0.0  6.00  3.00 1.00 0.00 0.00 0.00 0.00 0. 0.  \n"
       f2.write(temp_write)
       temp_write="  6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723   \n"
       f2.write(temp_write)
       temp_write="  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  \n"
       f2.write(temp_write)
#       temp_write="  #30|"+str(int(float(l)*100.0)*4)+"|30"+"  "+"Quad"+"  "+name+"\n"
       temp_write="  1.  cm"+"\n"
       f2.write(temp_write)
       temp_write="  1  0.0  0.0  0.0  "+"\n"
       f2.write(temp_write)
    elif type=="SEXTUPOLE" and float(l)!=0:
       temp_write="'MULTIPOL'"+"    "+"SEXT"+"      "+name+"\n"
       f2.write(temp_write)
       temp_write="0"+"  "+".Sext"+"\n"
       f2.write(temp_write)
       temp_write="  "+str(float(l)*100.0)+"  "+str(x10)+"  "+str(float(angle)/float(l)*10.0)+"  "+str(float(k1l)/float(l)*0.1*10.0)+"  "+str(float(k2l)/float(l)*(x10/100.0)*(x10/100.0)/2.0*10.0)+"  0.0  0.0  0.0  0.0  0.0  0.0  0.0   "+"\n"
       f2.write(temp_write)
       temp_write="  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  \n"
       f2.write(temp_write)
       temp_write="  6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723   \n"
       f2.write(temp_write)
       temp_write="  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  \n"
       f2.write(temp_write)
       temp_write="  6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723   \n"
       f2.write(temp_write)
       temp_write="  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  \n"
       f2.write(temp_write)
#       temp_write="  #30|"+str(int(float(l)*100.0))+"|30"+"  "+"Sext"+"  "+name+"\n"
       temp_write="  1.  cm"+"\n"
       f2.write(temp_write)
       temp_write="  1  0.0  0.0  0.0  "+"\n"
       f2.write(temp_write)
    elif type=="SOLENOID" and float(l)!=0.:
       temp_write="'SOLENOID'"+"    "+"SOLE"+"      "+name+"\n"
       f2.write(temp_write)
       temp_write="0"+"  "+".sole"+"\n"
       f2.write(temp_write)
       temp_write="  "+str(float(l)*100.0)+"  "+str(x2)+"  "+str(float(ksi)/float(l)*brho*10.0)+"\n"
       f2.write(temp_write)
       temp_write="  "+"2.5"+"  "+"2.5"+"  "+"\n"
       f2.write(temp_write)
 #      temp_write="  #30|"+str(int(float(l)*100.0))+"|30"+"  "+"Sole"+"  "+name+"\n"
       temp_write="  1.  cm"+"\n"
       f2.write(temp_write)
       temp_write="  1  0.0  0.0  0.0  "+"\n"
       f2.write(temp_write)
    elif type=="MARKER" and float(l)==0.:
       temp_write="'MARKER'"+"    "+"MARK"+"      "+name+"\n"
       f2.write(temp_write)
    elif type=="MATRIX":
       temp_write="'TRANSMAT'"+"    "+"MATX"+"      "+name+"\n"
       f2.write(temp_write)
       temp_write="  "+"1"+"\n"
       f2.write(temp_write)
       temp_write="  "+str(float(l)*1.0)+"\n"
       f2.write(temp_write)
       temp_write=("  "+str(rm[0,0])+"  "+str(rm[0,1])+"  "+str(rm[0,2])+
                   "  "+str(rm[0,3])+"  "+str(rm[0,4])+"  "+str(rm[0,5])+"\n")
       f2.write(temp_write)
       temp_write=("  "+str(rm[1,0])+"  "+str(rm[1,1])+"  "+str(rm[1,2])+
                   "  "+str(rm[1,3])+"  "+str(rm[1,4])+"  "+str(rm[1,5])+"\n")
       f2.write(temp_write)
       temp_write=("  "+str(rm[2,0])+"  "+str(rm[2,1])+"  "+str(rm[2,2])+
                   "  "+str(rm[2,3])+"  "+str(rm[2,4])+"  "+str(rm[2,5])+"\n")
       f2.write(temp_write)
       temp_write=("  "+str(rm[3,0])+"  "+str(rm[3,1])+"  "+str(rm[3,2])+
                   "  "+str(rm[3,3])+"  "+str(rm[3,4])+"  "+str(rm[3,5])+"\n")
       f2.write(temp_write)
       temp_write=("  "+str(rm[4,0])+"  "+str(rm[4,1])+"  "+str(rm[4,2])+
                   "  "+str(rm[4,3])+"  "+str(rm[4,4])+"  "+str(rm[4,5])+"\n")
       f2.write(temp_write)
       temp_write=("  "+str(rm[5,0])+"  "+str(rm[5,1])+"  "+str(rm[5,2])+
                   "  "+str(rm[5,3])+"  "+str(rm[5,4])+"  "+str(rm[5,5])+"\n")
       f2.write(temp_write)
#    elif type=="HMONITOR" and float(l)==0.:
#       temp_write= name+":"+"MARKER; \n"
#       f2.write(temp_write)
#    elif type=="VMONITOR" and float(l)==0.:
#       temp_write= name+":"+"MARKER; \n"
#       f2.write(temp_write)
#    elif type=="MONITOR" and float(l)==0.:
#       temp_write= name+":"+"MARKER; \n"
#       f2.write(temp_write)
#    elif type=="HKICKER" or type=="VKICKER":
#       temp_write= name+":"+ "KICKER" +",L="+l+ "; \n"
#       f2.write(temp_write)
#    elif type=="RFCAVITY" and rf_first==1:
#       rf_first=0
#       temp_write= name+":"+ "RFCA,  FREQ = 78196.2811281, VOLT = 5e6;\n"
#       f2.write(temp_write)
    else:
       temp_write="'DRIFT'" +"    "+"DRIF"+"      "+name+"\n" 
       f2.write(temp_write)
       temp_write="  "+str(float(l)*100.0)+"\n"
       f2.write(temp_write) 

harm=3416
vpeak=3432000.000             # peak voltage, V
phis=2.88153859504264         # synchronous phase, rad (pi=3.1415926535897932384626433832795)
f2.write("'CAVITE'"+"\n")
f2.write("  "+"2"+"\n")
f2.write("  "+str(circ)+"  "+str(harm)+"\n")
f2.write("  "+str(vpeak)+"  "+str(phis)+"\n")

f2.write("'MARKER'"+"    "+"MARK"+"    "+"COMPLETE_RING$END"+"\n")
f2.write("'OPTIONS'"+"\n")
f2.write("  "+"1"+"   "+"1"+"\n")
f2.write("  "+"WRITE ON"+"\n")
f2.write("'FAISCEAU'"+"\n") 
f2.write("'SPNPRT'"+"\n") 
f2.write("'TWISS'"+"\n")
f2.write("  "+"2"+"   "+"1."+"   "+"1."+"   "+"PRINT"+"\n")
f2.write("'REBELOTE'"+"\n")
f2.write("  "+"299"+"   "+"0.2"+"   "+"99"+"\n")
f2.write("'END'"+"\n") 

print "total length is: ",total_l
print "total number of elements =  ", len(nameseq)
print "number of different elements =  ",len(namelist)
print "total angle = ",total_ang

elem_per_line=10
mline=len(nameseq)/elem_per_line

#for i in range(mline):
#    for j in range(elem_per_line):
#        f2.write(nameseq[elem_per_line*i+j]+",")
#    f2.write("&\n")

#for i in range(len(nameseq)-elem_per_line*mline-1):
#    f2.write(nameseq[elem_per_line*mline+i]+",")
#f2.write(nameseq[len(nameseq)-1]+")\n")

f2.close()
sys.exit()
