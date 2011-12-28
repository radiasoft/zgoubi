        subroutine bblmap(linearmap,tunex,tuney,tunez,
     >              alphax,betax,alphay,betay,alfmom,ga,bet,clight,circ)
        implicit double precision (a-h,o-z)
        double precision linearmap(6,6) 
        double precision  tunex,tuney,tunez,alphax,betax,
     >                                   alphay,betay
        double precision  twopi,cx,sx,gx,cy,sy,gy,cz,sz
        double precision  alfmom,ga,clight,circ,szspz,bet

        twopi = 4*dasin(1.0d0)
        gx = (1+alphax*alphax)/betax
	cx=dcos(twopi*tunex)
	sx=dsin(twopi*tunex)
        linearmap(1,1) = cx + alphax*sx
        linearmap(2,1) = -gx*sx
        linearmap(1,2) = betax*sx
        linearmap(2,2) = cx - alphax*sx
        gy = (1+alphay*alphay)/betay
	cy=dcos(twopi*tuney)
	sy=dsin(twopi*tuney)
        linearmap(3,3) = cy + alphay*sy
        linearmap(4,3) = -gy*sy
        linearmap(3,4) = betay*sy
        linearmap(4,4) = cy - alphay*sy
	cz=dcos(twopi*tunez)
	sz=dsin(twopi*tunez)
        szspz = (alfmom-1/(ga*ga))*bet*clight/
     >            (tunez*bet*clight*twopi/circ)
        linearmap(5,5) = cz
        linearmap(6,5) = sz/szspz
        linearmap(5,6) = -sz*szspz
        linearmap(6,6) = cz

!         Set coupling terms to 0

        linearmap(1,3)=0
        linearmap(1,4)=0
        linearmap(1,5)=0
        linearmap(1,6)=0
        linearmap(2,3)=0
        linearmap(2,4)=0
        linearmap(2,5)=0
        linearmap(2,6)=0

        linearmap(3,1)=0
        linearmap(3,2)=0
        linearmap(3,5)=0
        linearmap(3,6)=0
        linearmap(4,1)=0
        linearmap(4,2)=0
        linearmap(4,5)=0
        linearmap(4,6)=0

        linearmap(5,1)=0
        linearmap(5,2)=0
        linearmap(5,3)=0
        linearmap(5,4)=0
        linearmap(6,1)=0
        linearmap(6,2)=0
        linearmap(6,3)=0
        linearmap(6,4)=0

        return
        end
