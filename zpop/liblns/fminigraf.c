/*
  ZGOUBI, a program for computing the trajectories of charged particles
  in electric and magnetic fields
  Copyright (C) 1988-2007  François Méot

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor,
  Boston, MA  02110-1301  USA

  François Meot <fmeot@bnl.gov>
  Brookhaven National Laboratory      
  C-AD, Bldg 911
  Upton, NY, 11973
*/
/*                                                                 */
/* AUTEUR : J.L.HAMEL                                              */
/*                                                                 */
/*                                                                 */
/* *****LBMITV:BIBLIOTHEQUE DE TRACE DE COURBES EN FORTRAN   ***** */
/* -----CE MODULE CONTIENT:                                        */
/*                         TRAXE                                   */
/*                         TRDEF                                   */
/*                         TRTEXT                                  */
/*                         GRADAX                                  */
/*                         TRGRAD                                  */
/*                         TRGRAL                                  */
/*                         TRDIV                                   */
/*                         LABAXE                                  */
/*                         LABLOG                                  */
/*                         INIVCF                                  */
/*                         FINVCF                                  */
/*                         VECTPH                                  */
/*                         TRSYMB                                  */
/*                         INIVPH                                  */
/*                         VECTF                                   */
/*                         TXTFBG                                  */
/*                         FBGTT                                   */
/*                         INISCA                                  */
/*                         SAVECR                                  */
/*                         DEFMKR                                  */
/*                         TRGETD                                  */
/*                                                                 */

/* GENERATION CONDITIONNELLE SELON SYSTEME */

#if defined(SUN) || defined(LINUX) || defined(GFORTRAN4)
#define TRAXE  traxe_  
#define TRDEF  trdef_
#define TRTEXT trtext_
#define INIVCF inivcf_
#define FINVCF finvcf_
#define VECTPH vectph_
#define TRSYMB trsymb_
#define INIVPH inivph_
#define VECTF  vectf_
#define TXTFBG txtfbg_
#define FBGTXT fbgtxt_
#define INISCA inisca_
#define SAVECR savecr_
#define DEFMKR defmkr_
#define DEFCAR defcar_
#define TRGETD trgetd_
#endif

#ifdef HP
#define TRAXE  traxe  
#define TRDEF  trdef
#define TRTEXT trtext
#define INIVCF inivcf
#define FINVCF finvcf
#define VECTPH vectph
#define TRSYMB trsymb
#define INIVPH inivph
#define VECTF  vectf
#define TXTFBG txtfbg
#define FBGTXT fbgtxt
#define INISCA inisca
#define SAVECR savecr
#define DEFMKR defmkr
#define DEFCAR defcar
#define TRGETD trgetd
#endif

/* DECLARATIONS GLOBALES */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define FALSE	0
#define TRUE	1
#define MAXXTK	1023
#define MAXYTK	779
#define MAXX	maxx
#define MAXY	maxy
#define HCARTK	19
#define LCARTK	12
#define HCAR	hcar
#define LCAR	lcar
#define M	1

#define SCALX(x)	x = (x * MAXXTK)/maxx
#define SCALY(y)        y = (y * MAXYTK)/maxy
/* Sam Tygier, Oct 2014 - #define XTERM		( !strcmp(getenv("TERM"),"xterm") ||	\ */
#define XTERM		( !strncmp(getenv("TERM"),"xterm",5) ||	\
			((getenv("DISPLAY") &&\
				 !strcmp(getenv("TERM"),"vt100"))) )

static int v_init,v_ix,v_iy;
static int v_ix1,v_ix2,v_iy1,v_iy2;
static char v_typlin,g_typlin; 
static int g_motr,g_ntic,g_logx,g_logy,g_point;
static int g_ix0,g_ix1,g_ix2,g_iy0,g_iy1,g_iy2;
static double g_x1,g_x2,g_y1,g_y2,g_dx,g_dy,g_cx,g_cy,g_div;
static double g_sincr[10];

static FILE *g_pps, *g_ppsc;
static int g_ps = FALSE;

static int g_tcar = 2, g_pscar = 2;
/* static int g_hcar[4] = {12,14,19,21}; */
/* static int g_hcar[4] = {21,26,31,36}; */
/* static int g_hcar[4] = {21,26,36,46}; */
static int g_hcar[4] = {27,26,30,42};
static int g_lcar[4] = {8,9,12,13};
static char g_ccar[4] = {';',':','9','8'};

static int g_imark = 0, g_xmark, g_ymark;
static char g_tmk[5][30] = {
	{2,0,-M,0,M, 2,-M,0,M,0, 0},
	{2,0,-M,0,M, 2,-M,0,M,0, 2,-M,-M,M,M, 2,-M,M,M,-M, 0},
	{5,-M,M,M,M,M,-M,-M,-M,-M,M, 0},
	{2,-M,-M,M,M, 2,-M,M,M,-M, 0},
	{5,-M,0,0,M,M,0,0,-M,-M,0}
};

static int maxx = MAXXTK, maxy = MAXYTK;
static int hcar = HCARTK, lcar = LCARTK;
static int g_mode = FALSE;
static int g_nmod[7] = {0,1,2,4,6,9,17};
static int g_ldiv[7] = {0,0x1FF,0x1FF,0x1FF,0x109,0x108,0x100};
static int g_lpoint[7] = {0,0x1FF,0x109,0x108,0x100,0x100,0x100};
static int g_lbit[10] = {0,0x1,0x2,0x4,0x8,0x10,0x20,0x40,0x80,0x100};
static double g_divlog[10] = {0.0,.30103,.47712,.60206,.69897,
				 .77815,.84510,.90309,.95425,1.};
void FBGTXT(),TXTFBG();
void _traxe(),_trdef(),_trtext(),_vectph(),_trsymb();
void _inivph(),_vectf();
static void gradax(),trgrad(),trgral(),trdiv(),labaxe(),lablog(),gclip();
static void gpoint(),xytek(),postext();
static int gadjust();
static double flmod(),sgn();
static void inipost(),finpost(),vecpost(),dashpost(),textpost();

/*
   *****TRGETD: LECTURE PARAMETRES DES AXES
*/

void TRGETD(c1,c2,div,ntic)
	float *c1,*c2,*div;
	int *ntic;

{
	(void)gradax(*c2-*c1);
	*div = g_div;
	*ntic = g_ntic;
}

/*
   *****TRAXE:TRACE D'AXES AUTOMATIQUE AVEC OPTIONS LOG ET 1/2 LOG
*/

void TRAXE(xmin,xmax,ymin,ymax,mode)
	float *xmin,*xmax,*ymin,*ymax;
	int *mode;
{
	_traxe(*xmin,*xmax,*ymin,*ymax,*mode);
}

void _traxe(xmin,xmax,ymin,ymax,mode)
	double xmin,xmax,ymin,ymax;
	int mode;

{
	int ier,modlog,erlog,lx1,lx2,ly1,ly2,idx,idy,x0,y0,savtcar;

        savtcar = g_tcar;
        g_tcar = 2;

	/* TEST VALIDITE DES PARAMETRES PHYSIQUES */
	_inivph(1);
	g_motr = (mode % 10)-1;
	modlog = mode / 10;
	ier = 0;

	if( (g_motr < 0) || (g_motr > 2) || (modlog  < 0) || (modlog > 3) ) {
		ier%=1;
	}

	if( ((xmax <= xmin) || (ymax <= ymin)) && !ier) {
		ier = 2;
	}

	if(!ier) {
		g_logx = ( (modlog == 1) || (modlog == 3) );
		g_logy = ( (modlog == 2) || (modlog == 3) );
		if(g_logx && ( (xmin <= 0) || (xmax <= 0) ) ) {
			ier = 3;
		}
		if(g_logy && ( (ymin <= 0) || (ymax <=0) ) ) {
			ier = 3;
		}
	}

	if(!ier) {
		g_x1 = xmin;
		g_x2 = xmax;
		if(g_logx) {
			g_x1 = log10((double)g_x1);
			g_x2 = log10((double)g_x2);
			lx1 = g_x1 + 0.0001 * sgn(g_x1);
			lx2 = g_x2 + 0.0001 * sgn(g_x2);
			if(fabs((double)(g_x1 - lx1)) >= 0.0001) {
				lx1 = lx1 - 0.5 + 0.5 * sgn(g_x1);
			}
			if(fabs((double)(g_x2 - lx2)) >= 0.0001) {
				lx2 = lx2 + 0.5 + 0.5 * sgn(g_x2);
			}
			g_x1 = lx1;
			g_x2 = lx2;
		}
		g_y1 = ymin;
		g_y2 = ymax;
		if(g_logy) {
			g_y1 = log10((double)(g_y1));
			g_y2 = log10((double)(g_y2));
			ly1 = g_y1 + 0.0001 * sgn(g_y1);
			ly2 = g_y2 + 0.0001 * sgn(g_y2);
			if(fabs((double)(g_y1 - ly1)) >= 0.0001) {
				ly1 = ly1 - 0.5 + 0.5 * sgn(g_y1);
			}
			if(fabs((double)(g_y2 - ly2)) >= 0.0001) {
				ly2 = ly2 + 0.5 + 0.5 * sgn(g_y2);
			}
			g_y1 = ly1;
			g_y2 = ly2;
		}
		g_dx = g_x2 - g_x1;
		g_dy = g_y2 - g_y1;

		/* PARAMETRES DE L'ECRAN */
		g_ix1 = 0.1  * MAXX + 0.5;
		g_ix2 = 0.95 * MAXX + 0.5;
		_vectf(g_ix1,g_ix2,6);
		g_iy1 = 0.2  * MAXY + 0.5;
		g_iy2 = 0.95 * MAXY + 0.5;
		_vectf(g_iy1,g_iy2,7);
		idx = g_ix2  - g_ix1;
		idy = g_iy2 - g_iy1;
		g_cx = idx / g_dx;
		g_cy = idy / g_dy;

		/* POSITION DES DROITES X=0 ET Y=0 SUR L'ECRAN */
		x0 = -g_x1 * g_cx + g_ix1;
		y0 = -g_y1 * g_cy + g_iy1;
		g_ix0 = -1;
		g_iy0 = -1;
		if(!g_logx && (x0 > g_ix1) && (x0 < g_ix2) ) {
			g_ix0 = x0 + 0.5;
		}
		if(!g_logy && (y0 > g_iy1) && (y0 < g_iy2) ) {
			g_iy0 = y0 + 0.5;
		}

		/* TRACE DES AXES S'IL Y A LIEU */
		if(g_motr != 2) {
				/* GRADUATIONS SUR LES AXES HORIZONTAUX */
			if(!g_logx) {
				(void)gradax(g_dx);
				trgrad(0,g_x1,g_dx,g_cx,g_ix1,g_iy1,
					 g_iy0,g_iy2);
				}
			else {
				trgral(0,lx1,lx2,g_ix1,g_ix2,g_iy1,
					 g_iy2,g_motr,&erlog);
				if(erlog) {
					ier = 4;
				}
			}
				/* GRADUATIONS SUR LES AXES VERTICAUX */
			if(!g_logy) {
				(void)gradax(g_dy);
					trgrad(1,g_y1,g_dy,g_cy,g_iy1,g_ix1,
					 g_ix0,g_ix2);
				}
			else {
				trgral(1,ly1,ly2,g_iy1,g_iy2,g_ix1,
				 g_ix2,g_motr,&erlog);
				if(erlog) {
					ier = 4;
				}
			}
			/* TRACE DU CADRE ET DES DROITES X=0,Y=0 */
			/* S'IL Y A LIEU                         */
			_inivph(1);
			_vectf(g_ix1,g_iy1,4);
			_vectf(g_ix2,g_iy1,2);
			_vectf(g_ix2,g_iy2,2);
			_vectf(g_ix1,g_iy2,2);
			_vectf(g_ix1,g_iy1,2);
			if(g_ix0 != -1) {
				_vectf(g_ix0,g_iy1,4);
				_vectf(g_ix0,g_iy2,2);
			}
			if(g_iy0 != -1) {
				_vectf(g_ix1,g_iy0,4);
				_vectf(g_ix2,g_iy0,2);
			}
		}
	}

	/* SORTIE EN ERREUR */
	if(ier) {
		printf("*ERREUR TRAXE:%d*\n",ier);
	}
	g_tcar = savtcar;
}

/*
   *****TRDEF:DEFINITION DE L'ESPACE OBJET
*/

void TRDEF(xmin,xmax,ymin,ymax,imin,imax,jmin,jmax)
	float *xmin,*xmax,*ymin,*ymax;
	int   *imin,*imax,*jmin,*jmax;

{
	_trdef(*xmin,*xmax,*ymin,*ymax,*imin,*imax,*jmin,*jmax);
}

void _trdef(xmin,xmax,ymin,ymax,imin,imax,jmin,jmax)
	double xmin,xmax,ymin,ymax;
	int   imin,imax,jmin,jmax;

{
	int ier,idx,idy;

	/* RAZ G_LOGX,G_LOGY */
	g_logx = FALSE;
	g_logy = FALSE;

	/* TEST VALIDITE DES PARAMETRES PHYSIQUES */
	ier = 0;
	if( (xmax <= xmin) || (ymax <= ymin) ) {
		ier = 1;
	}

	/* TEST VALIDITE DES PARAMETRES ECRAN */
	if(!ier) {
		if( (imax <= imin) || (jmax <= jmin) || (imin < 0) ||
		    (imax > MAXX) || (jmin < 0) || (jmax > MAXY) ) {
			ier = 2;
		    }
	}

	if(!ier) {

		/* PARAMETRES ESPACE SUJET */
		g_x1 = xmin;
		g_x2 = xmax;
		g_y1 = ymin;
		g_y2 = ymax;
		g_dx = g_x2 - g_x1;
		g_dy = g_y2 - g_y1;

		/* PARAMETRES ESPACE OBJET */
		g_ix1 = imin;
		g_ix2 = imax;
		_vectf(g_ix1,g_ix2,6);
		g_iy1 = jmin;
		g_iy2 = jmax;
		_vectf(g_iy1,g_iy2,7);
		idx = g_ix2 - g_ix1;
		idy = g_iy2 - g_iy1;
		g_cx = idx / g_dx;
		g_cy = idy / g_dy;
	}
	else {
		/* SORTIE EN ERREUR */
		printf("*ERREUR TRDEF:%d*",ier);
	}
}

/*
   *****TRTEXT:ECRITURE DE TEXTE
*/

void TRTEXT(x,y,texte,ityp,n)
	float *x,*y;
	char *texte;
	int *ityp,n;

{
	char txt[200];
	int i;

	if(n > 199) return;
	for(i = 0; i < n; i++) txt[i] = texte[i];
	txt[n] = 0;
	_trtext(*x,*y,txt,*ityp);
}

void _trtext(x,y,texte,ityp)
	double x,y;
	char *texte;
	int ityp;

{
	double xx,yy;
	int ier,jx,jy;

	xx = x;
	yy = y;
	ier = FALSE;
	if(ityp != 0) {
		if(g_logx) {
			xx = log10((double)xx);
		}
		if(g_logy) {
			yy = log10((double)yy);
		}
	xx = (xx - g_x1) * g_cx + 0.5 + g_ix1;
	yy = (yy - g_y1) * g_cy + 0.5 + g_iy1;
	ier = ( (xx < 0) || (xx > MAXX) || (yy < 0) || (yy > MAXY) );
	}
	if(!ier) {
		jx = xx;
		jy = yy;
		postext(jx,jy,texte);
	}
}

/*
   *****GRADAX:CALCUL PARAMETRES GRADUATIONS D'AXES
*/

static void gradax(delta)
	double delta;

{
	double s,xn;
	int n;

	s = 0;
	if( (fabs((double)(delta-1)) >= 1.0E-06) ) {
		s = log10(delta);
	}
	n = s;
	s = s - n;
	if(s < 0) {
	s = 1 + s;
	n = n - 1;
	}
	xn = pow(10.0,(double)n);
	if(s > 0.69897) {
		g_div = xn;
		g_ntic = 5;
	}
	else {
		if(s <= 0.30103) {
			g_div = 0.2 * xn;
			g_ntic = 4;
		}
		else {
			g_div = 0.5 * xn;
			g_ntic = 5;
		}
	}
}

/*
   *****TRGRAD:TRACE DE GRADUATION D'AXES
*/

static void trgrad(nat,c,dc,cc,ic,j1,j0,j2)
	int nat,ic,j1,j0,j2;
	double c,dc,cc;

{
	int unsur2,occup,premdv,typ1,typ2,typ;
	double tic,tic1,epsdc,dce,epsdiv,epstic;
	double tica,ticb,diva,divb,div1;
	int nbtic,nbtica,nbtic1,nbdiv,nbdiva,i,m,n;
	int nparz,ioftic;
	double x,z;
	int ix,llong;

	/* PREPARATION DU TRACE DES GRADUATIONS */

	/* LARGEUR D'UN TIC */
	tic = g_div / g_ntic;

	/* APPROXIMATIONS */
	epsdc = 0.00001 * dc;
	dce = dc + epsdc;
	epstic = 0.0001 * tic;
	epsdiv = 0.0001 * g_div;

	/* COORDONNEE DU 1ER TIC */
	tica = flmod(c,tic);
	tica = fabs((double)tica);
	ticb = tic - tica;
	if( (tica < epstic) || (ticb < epstic) ) {
		tica = 0;
	}
	if( (c >= 0) && (tica != 0) ) {
		tica = ticb;
	}
	tic1 = tica + c;

	/* COORDONNEE DE LA 1ERE DIVISION */
	diva = flmod(c,g_div);
	diva = fabs((double)diva);
	divb = g_div - diva;
	if( (diva < epsdiv) || (divb < epsdiv) ) {
		diva = 0;
	}
	if( (c >= 0) && (diva != 0) ) {
		diva = divb;
	}
	div1 = diva + c;

	/* NOMBRE DE TICS JUSQU'A LA 1ERE DIVISION */
	nbtic1 = (div1 - tic1 + epstic) / tic + 1.0;

	/* NOMBRE TOTAL DE TICS (Y COMPRIS LES DIVISIONS) */
	nbtic = dce / tic + 1.0;
	nbtica = nbtic - 1;
	if( (tica + nbtica * tic) > dce) {
	nbtic = nbtica;
	}

	/* NOMBRE TOTAL DE DIVISIONS */
	nbdiv = dce / g_div + 1.0;
	nbdiva = nbdiv - 1;
	if( (diva + nbdiva * g_div) > dce) {
	nbdiv = nbdiva;
	}

	/* PARITE DU NB DE DIVISIONS JUSQU'A ZERO */
	nparz = 0;
	x = (fabs(c) + epsdiv) / g_div;
	if(x <= 32767.0) {
	nparz = x;
	}
	nparz = nparz % 2;

	/* OFFSET SUR LE COMPTE INITIAL DE TICS */
	ioftic = g_ntic - nbtic1;

	/* BOUCLE SUR LES TICS POUR TRACER LES GRADUATIONS */

	unsur2 = (nbdiv >= 8);
	typ = (g_div >= 1.0) && (g_div <= 9999.0);
	typ1 = (g_x1 >= -999.0) && (g_x2 <= 9999.0) && typ;
	typ2 = (g_y1 >= -9999.0) && (g_y2 <= 99999.0) && typ;
	m = 0;
	for(i = 1; i <= nbtic; i++) {
		n = i - 1;

		/* COORDONNEE DU TIC COURANT */
		x = tic1 + n * tic;
		if(fabs((double)x) < epsdiv) {
			x = 0;
		}
		ix = (int)((x - c) * cc + 0.5) + ic;

		/* TRAITEMENT D'UNE DIVISION */
		/* TIC ORDINAIRE*/
		if( ( (i + ioftic) % g_ntic) != 0) {
			llong = MAXX/200;
		}
		else {
			llong = MAXX/100;
			premdv = (fabs((double)(x - div1)) < tic);
			m = m + 1;
			if(g_motr == 1) {
				_inivph(2);
				trdiv(2,ix,j1,0,j2 - j1,nat);
				_inivph(1);
			}
			if(nat == 0) {
				if(!(unsur2 && (((m + nparz) % 2) == 0) && 
				   (!typ1)) ) {
					labaxe(ix,j1 - HCAR/2,x,typ1);
					if(premdv) {
						occup = (ix < (g_ix1 + 10) );
					}
				   }
			}
			else {
				z = fabs((double)x);
				if(!(premdv && occup && (ix < (g_iy1 + 20)) &&
				 ((z < 0.1) || (z >= 10000.0)) && (!typ2) )) {
					labaxe(MAXX/20,ix + HCAR/2,x,typ2);
			     	 }
			}
		}

		/* TRAITEMENT COMMUN:TIC ET DIVISION */
		trdiv(2,ix,j1,0,llong,nat);
		if(j0 != -1) {
			trdiv(2,ix,j0,-llong,llong,nat);
		}
		trdiv(2,ix,j2,-llong,0,nat);
	}
}

/*
   *****TRGRAL:TRACE D'AXES LOGARITHMIQUES
*/

static void trgral(nat,l1,l2,k1,k2,j1,j2,motr,perlog)
	int nat,l1,l2,k1,k2,j1,j2,motr,*perlog;

{
	int unsur2,unsur5,labint,maindv,i,j,l;
	double sx,smodul;
	int nmodul,erlog,im,im1,ix,lx,llong,lg,lg1,lg2;

	/* PREPARATION DES PARAMETRES */

	/* CALCUL DU NOMBRE DE MODULES */
	erlog = FALSE;
	nmodul = l2 - l1;
	erlog = (nmodul > 40);
	if(erlog) goto err_trgral;

	/* DETERMINER LA CLASSE DU TRACE */
	for(im = 2; im <= 6; im++) {
		lg = (g_nmod[im] > nmodul);
		if(lg) break;
	}
	if(!lg) {
		im = 7;
	}
	im = im - 1;

	/* PARAMETRES DE LA CLASSE */
	unsur2 = (im == 5);
	unsur5 = (im == 6);
	labint = (im <= 2);

	/* LONGUEUR D'UN MODULE SUR L'ECRAN */
	smodul = (double)(k2 - k1) / (double)nmodul;

	/* INCREMENTS */
	for(i = 1; i <= 9; i++) {
		g_sincr[i] = smodul * g_divlog[i];
	}

	/* BOUCLE DE TRACE */

	/* TRACE DU 1ER LABEL */
	if(!(unsur2 && ((l1 % 2) != 0)) || (unsur5 && ((l1 % 5) != 0)) ) {
		if(nat == 0) {
			lablog(k1,j1-HCAR/2,l1,TRUE);
		}
		else {
			lablog(MAXX/20,k1 + HCAR/2,l1,TRUE);
		}
	}

	/* BOUCLE SUR LES MODULES */
	for(i = 1; i <= nmodul; i++) {
		im1 = i - 1;
		sx = im1 * smodul;

		/* BOUCLE SUR LES INCREMENTS */
		for(j = 1; j <= 9; j++) {
			maindv = (j == 9);

			/* COORDONNEE DU TIC COURANT */
			ix = sx + g_sincr[j] + k1;
			l = l1 + i;

			/* LONGUEUR D'UN TIC */
			if (maindv && (im < 6) ) {
				llong = MAXX/100;
			}
			else {
				llong = MAXX/200;
			}

			/* TEST TRACE DES POINTILLES ET DES LABELS */
			lg1 = ((g_ldiv[im]   & g_lbit[j]) == 0);
			lg2 = ((g_lpoint[im] & g_lbit[j]) == 0);
			lg  = ( (unsur2 && ((l % 2) != 0)) || lg2);
			lg  = ( (unsur5 && ((l % 5) != 0)) || lg );

			/* LABELS */
			if ( (maindv || labint) && (!((j == 8) && 
			     (nat == 0))) && (!lg) ) {
				lx = MAXX/20;
				if (!maindv) {
					l = j + 1;
				}
				if(nat == 0) {
					lablog(ix,j1-HCAR/2,l,maindv);
				}
				else {
					lablog(lx,ix + HCAR/2,l,maindv);
				}
				if(im >= 6) {
					llong = 4;
				}
			}
			/* TRACE DES TICS ET POINTILLES */
			if(!(lg1 || ((i == nmodul) && maindv)) ) {
				if(g_motr == 1) {

					/* POINTILLES */
					_inivph(2);
					trdiv(2,ix,j1,0,j2-j1,nat);
					_inivph(1);
				}

				/* TICS */
				trdiv(2,ix,j1,0,llong,nat);
				trdiv(2,ix,j2,-llong,0,nat);
			}

			/* FIN DES BOUCLES */
		}
	}

	err_trgral:

	*perlog = erlog;
}

/*
   *****TRDIV:TRACE D'UN TIC OU D'UNE DIVISION SUR UN AXE
*/

static void trdiv(mode,ix,iy,ideb,ifin,nat)
	int mode,ix,iy,ideb,ifin,nat;

{
	if(nat != 0) {

		/* AXE VERTICAL:PERMUTER LES COORDONNEES */
		_vectf(iy + ideb,ix,4);
		_vectf(iy + ifin,ix,mode);
	}
	else {

		/* AXE HORIZONTAL */
		_vectf(ix,iy + ideb,4);
		_vectf(ix,iy + ifin,mode);
	}
}

/*
   *****LABAXE:TRACE DE LABEL DE GRADUATION DES AXES
*/

static void labaxe(ix,iy,val,typi)
	int ix,iy,typi;
	double val;

{
	int jx,jy,n,e;
	double m;
	char ns[30],ms[30],es[30],aux[30];

	jx = ix;
	jy = iy - HCAR;

	/* TEST TYPE DE FORMAT */
	if(!typi) {

		/* FORMAT FLOTTANT */
		e = 0;
		if(val == 0) {
		strcpy(ms,"0.0");
		}
		else {
			sprintf(ms,"%.4g",val);
			if(!strncmp(ms,"-0",2)) {
				strcpy(aux,"-");
				strcat(aux,&ms[2]);
				strcpy(ms,aux);
			}
			if(strchr(ms,'e')) {
				e = (int)(log10(fabs(val * 1.001)));
			        m = (val / pow(10.0,(double)e));
				sprintf(ms,"%g",m);
				sprintf(es,"%d",e);
				if(e >= 0.0) {
					strcpy(aux,"E+");
					strcat(aux,es);
					strcpy(es,aux);
				}
				else {
					strcpy(aux,"E");
					strcat(aux,es);
					strcpy(es,aux);
				}
			}
			if(!strchr(ms,'.')) {
				strcat(ms,".");
			}
		}
		jx = ix - strlen(ms) * LCAR / 2;
		postext(jx,jy,ms);
		if(e != 0) {
			jx = ix - strlen(es) * LCAR / 2;
			jy = jy - HCAR;
			postext(jx,jy,es);
		}
	}
	else {

		/* FORMAT FIXE */
		sprintf(ns,"%d",(int)(val + 0.5 * sgn(val)));
		jx = ix - strlen(ns) * LCAR / 2;
		postext(jx,jy,ns);
	}
}

/*
   *****LABLOG:LABEL POUR AXE LOGARITHMIQUE
*/

static void lablog(ix,iy,ival,typi)
	int ix,iy,ival,typi;

{
	int jx,jy;

	char ns[30],aux[30];
	jy = iy - HCAR;
	sprintf(ns,"%d",ival);

	/* TEST FORMAT */
	if(typi) {

		/* LABEL PRINCIPAL */
		if(ival >= 0) {
			strcpy(aux,"1E+");
			strcat(aux,ns);
			strcpy(ns,aux);
		}
		else {
			strcpy(aux,"1E");
			strcat(aux,ns);
			strcpy(ns,aux);
		}
	}
	jx = ix - strlen(ns) * LCAR / 2;
	postext(jx,jy,ns);
}

/*
   *****INISCA:INIT TAILLE ECRAN
*/

void INISCA(ex,ey)
	int *ex,*ey;

{
	maxx = *ex;
	maxy = *ey;
	hcar = (HCARTK*maxy)/MAXYTK;
	lcar = (LCARTK*maxx)/MAXXTK;
}


/*
   *****INIVCF:INIT PARAMETRES DE VECTF
*/

void INIVCF()

{
	inipost();
	if(!v_init) {
	v_init = TRUE;
	}
	g_logx = FALSE;
	g_logy = FALSE;
	_inivph(1);
/* Seems to be _vectf that causes graphic window opening... :	*/
	_vectf(0,MAXX,6);
	_vectf(0,MAXY,7);
	_vectf(0,0,0);
}

/*
   *****FINVCF: ARRET DU GRAPHIQUE
*/

void FINVCF()
{
	finpost();
	v_init = FALSE;
	if(g_mode) FBGTXT();
}

/*
   *****VECTPH:TRACE DE VECTEUR DANS L'ESPACE PHYSIQUE
*/

void VECTPH(x,y,mode)
	float *x,*y;
	int *mode;

{
	_vectph(*x,*y,*mode);
}

void _vectph(x,y,mode)
	double x,y;
	int mode;

{
	double xx,yy;
	int   ix,iy;

	xx = x;
	yy = y;
	if(g_logx) {
		xx = log10((double)xx);
	}
	if(g_logy) {
		yy = log10((double)yy);
	}
	xx = (xx - g_x1) * g_cx + 0.5 + g_ix1;
	yy = (yy - g_y1) * g_cy + 0.5 + g_iy1;
	if(fabs((double)xx) >= 32767.0) {
		xx = sgn(xx) * 32767.0;
	}
	if(fabs((double)yy) >= 32767.0) {
		yy = sgn(yy) * 32767.0;
	}
	ix = xx;
	iy = yy;
	_vectf(ix,iy,mode);
}

/*
   *****TRSYMB:TRACE DE SYMBOLE SUR MITV
*/

void TRSYMB(x,y,mode,isymb,itypx,itypy,coefx,coefy)
	float *x,*y,*coefx,*coefy;
	int *mode,*itypx,*itypy;
	short isymb[];

{
	_trsymb(*x,*y,*mode,isymb,*itypx,*itypy,*coefx,*coefy);
}

void _trsymb(x,y,mode,isymb,itypx,itypy,coefx,coefy)
	double x,y,coefx,coefy;
	int mode,itypx,itypy;
	short isymb[];

{
	int n,i,j,mod,ix,iy,savx,savy;
	double xc,yc,xx,yy;

	savx = v_ix;
	savy = v_iy;
	i = 0;

	/* INIT D'UNE SUITE DE POINTS ET TEST FIN */
	n = 1;
	while(n != 0) {
		n = isymb[i++];

		/* BOUCLE SUR LES POINTS D'UNE SUITE*/
		for(j = 1; j <= n; j++) {
			xc = isymb[i++];
			if(xc > 32767) xc=xc - 65536;
			yc = isymb[i++];
			if(yc > 32767) yc=yc - 65536;
			xc = xc * coefx;
			yc = yc * coefy;
			if(itypx != 0) {

				/* ESPACE SUJET X */
				xc = xc + x;
				if(g_logx) xc = log10(xc);
				xx = (xc - g_x1) * g_cx + 0.5 + g_ix1;
			}
			else {

				/* ESPACE OBJET X */
				xx = x;
				if(g_logx) xx = log10(xx);
				xx = (xx - g_x1) * g_cx + 0.5 + xc + g_ix1;
			}

			if(itypy != 0) {

				/* ESPACE SUJET Y */
				yc = yc + y;
				if (g_logy) yc = log10(yc);
				yy = (yc - g_y1) * g_cy + 0.5 + g_iy1;
			}
			else {

				/* ESPACE OBJET Y */
				yy = y;
				if(g_logy) yy = log10(yy);
				yy = (yy - g_y1) * g_cy + 0.5 + yc + g_iy1;
			}

			/* TEST LIMITES */
			if(fabs((double)(xx)) >= 32767.) xx = sgn(xx) * 32767.;
			if(fabs((double)(yy)) >= 32767.) yy = sgn(yy) * 32767.;
			ix = xx;
			iy = yy;

			/* TRACE */
			mod = mode;
			if(j == 1) mod = 4;
			_vectf(ix,iy,mod);
		}
	}
	v_ix = savx;
	v_iy = savy;
}


/*
   *****INIVPH:INIT TYPE DE LIGNE POUR TRACE DE VECTEUR
*/

void INIVPH(n)
	int *n;

{
	_inivph(*n);
}

void _inivph(n)
	int n;

{
	int typlin;

	if( (n >= 1) && (n < 9) ) {
		typlin = 0x60 + ((n - 1) % 5);
		dashpost((char)typlin);
		_vectf(typlin,0,5);
	}
	g_point = (n == 9);
}

/*
   *****DEFMKR: DEFINITION DES MARQUEURS
*/

void DEFMKR(haut,type)
	int *haut,*type;

{
	if((*type >= 1) && (*type <= 6)) g_imark = *type - 1;
	else g_imark = 0;
	g_xmark = *haut;
	g_ymark = *haut;
	SCALX(g_xmark);
	SCALY(g_ymark);
}

/*
   *****DEFCAR: TAILLE DES CARACTERES
*/

void DEFCAR(haut,form,effa)
	int *haut,*form,*effa;

{
	g_tcar = *haut - 1;
	if((g_tcar < 0) || (g_tcar > 3)) g_tcar = 2;
}

/*
   *****VECTF:TRACE DE POINTS ET DE VECTEURS
*/

void VECTF(ix,iy,mode)
	int *ix,*iy,*mode;


{
	_vectf(*ix,*iy,*mode);
}

void _vectf(ix,iy,mode)
	int ix,iy,mode;


{
	int savix,saviy;

	if(!g_mode) TXTFBG();

	if(mode == 0) {
		/* Mise a l'origine, effacement de l'ecran */
		if(g_ps) {
			fclose(g_pps);
			inipost();
		}
		printf("\035\033\014");
		v_ix = 0;
		v_iy = 0;
	}
	else if(mode > 0) {

		switch(mode) {
			/* Trace avec effacement: impossible en mode TEKTRO */
			case 1 : break;

			/* Trace normal (pas d'inverse video) */
			case 2:
			case 3: {
				SCALX(ix); SCALY(iy);
				g_typlin = 
				   (mode ==2) ? v_typlin : v_typlin | 0x70;
				if(!g_point) {
					/* LIGNE */
					gclip(ix,iy);
				}
				else {
					/* POINT */
					gpoint(ix,iy,g_imark);
				}
				v_ix = ix;
				v_iy = iy;
				break;
			}

			/* Deplacement sans trace */
			case 4: {
				SCALX(ix); SCALY(iy);
				v_ix = ix;
				v_iy = iy;
				break;
			}

			/* Changer type de ligne */
			case 5: {
				v_typlin = (char)ix;
				break;
			}

			/* Clipping en X */
			case 6: {
				SCALX(ix); SCALX(iy);
				v_ix1 = ix;
				v_ix2 = iy;
				break;
			}

			/* Clipping en Y */
			case 7: {
				SCALY(ix); SCALY(iy);
				v_iy1 = ix;
				v_iy2 = iy;
				break;
			}
		}
	}
}


/* Routine elementaire de trace de texte */

static void postext(x,y,s)
	int x,y;
	char *s;

{
	char sxy[6];

	SCALX(x); SCALY(y);

	textpost(x,y,s);
	if(!g_mode) TXTFBG();
	xytek(x,y,sxy);
	printf("\033%c\035%s\037%s",g_ccar[g_tcar],sxy,s);
	fflush(stdout);
}

/* Routine xttek de conversion de coordonnees pour TEKTRO*/

static void xytek(x,y,s)
	int x,y;
	char s[6];

{
	s[0] = (y >> 5) + 0x20;
	s[1] = (y & 0x1F) + 0x60;
	s[2] = (x >> 5) + 0x20;
	s[3] = (x & 0x1F) + 0x40;
	s[4] = 0;
}

/* Routine de trace de point */

static void gpoint(ix,iy,g_imark)
	int ix,iy,g_imark;
{
	char s1[6],s2[6],*mk;
	int  savix,saviy,n,i,j,xc,yc;

	if(!g_imark) {
		if( (ix < v_ix1) || (ix > v_ix2) ||
		    (iy < v_iy1) || (iy > v_iy2) ) 
			return;
		vecpost(ix,iy,ix,iy+1);
		xytek(ix,iy,s1); xytek(ix+1,iy+1,s2);
		printf("\035\033\001\035%s%s\037",s1,s2);
		fflush(stdout);
	}
	else {
		savix = v_ix; saviy = v_iy;

		i = 0;
		mk = &g_tmk[g_imark - 1][0];

		n = 1;
		while(n != 0) {
			n = mk[i++];
	                for(j = 1; j <= n; j++) {
 	                        xc = mk[i++]*g_xmark + savix;
	                        yc = mk[i++]*g_ymark + saviy;
				if(j != 1) gclip(xc,yc);
				v_ix = xc; v_iy = yc;
			}
		}

		v_ix =savix;  v_iy  =saviy;
	}
}

/* Routine de clipping avant trace */

static void gclip(ix,iy)
	int ix,iy;

{
	int ix1,iy1,ix2,iy2;
	int lx1,lx2,ly1,ly2;
	char s1[6],s2[6];

	ix1 = v_ix; ix2 = ix; iy1 = v_iy; iy2 = iy;

	/* variables logiques: -1 = inferieur; 0 = dedans; 1 = superieur */  
	lx1 = ( (ix1 > v_ix2) - (ix1 < v_ix1) );
	ly1 = ( (iy1 > v_iy2) - (iy1 < v_iy1) );
	lx2 = ( (ix2 > v_ix2) - (ix2 < v_ix1) );
	ly2 = ( (iy2 > v_iy2) - (iy2 < v_iy1) );

	/* Les deux points tous deux inf. ou sup.: on sort */
	if( ((lx1 == lx2) && lx1) || ((ly1 == ly2) && ly1 ) ) return;

	/* On ajuste les points hors bornes et on sort s'il y a lieu */
	if (lx1 || ly1) if(gadjust(&ix1,&iy1,ix2,iy2)) return;
	if (lx2 || ly2) if(gadjust(&ix2,&iy2,ix1,iy1)) return;

	/* Trace des points */
	vecpost(ix1,iy1,ix2,iy2);
	xytek(ix1,iy1,s1);
	xytek(ix2,iy2,s2);
	printf("\035\033%c\035%s%s\037",g_typlin,s1,s2);
	fflush(stdout);
}

/* Routine d'ajustage des points hors bornes */

static int gadjust(ax,ay,x,y)
	int *ax,*ay,x,y;

{
	int xi,yi,xf,yf,i,notout;

	xi = *ax; yi = *ay; xf = x; yf = y;
	do {
		if ( (*ax < v_ix1) || (*ax > v_ix2) ||
		     (*ay < v_iy1) || (*ay > v_iy2) ) {
			xi = *ax; yi = *ay;
 		}
		else {
			xf =*ax; yf = *ay;
		}
		*ax = (xi + xf) >> 1; *ay = (yi + yf) >> 1;
		notout = ( (abs(*ax - x) > 1) || (abs(*ay - y) > 1) );
	}
	while( (notout) && ((*ax != v_ix1) && (*ax != v_ix2) &&
			  (*ay != v_iy1) && (*ay != v_iy2)) );
	return(notout ? 0 : 1);
}

/*
 ***** MODULO FLOTTANT
*/

static double flmod(x,y)
	double x,y;

{
	return( (double)(x - ((int)(x / y) * y)) );
}

/*
 ***** SIGNE
*/

static double sgn(x)
	double x;

{
	return( (x >= 0.0) ? 1.0 : -1.0);
}

/*
    ***** PASSAGE EN MODE GRAPHIQUE ET VICE VERSA
*/

void FBGTXT()

{
	char *s;

	if(g_mode) {
		if(XTERM)
			printf("\033\003");
		else
			printf("\033[?38l");
		g_mode = FALSE;
	}
}

void TXTFBG()

{
        char *s;

	if(!g_mode) {
		if(XTERM)
			printf("\033[?38h");
		else
			printf("\033[?38h");
		g_mode = TRUE;
	}
}

/*
    ***** FONCTIONS POSTSCRIPT
*/

/* Init fichier PostScript */

static void inipost()

{
	g_ps = FALSE;
	if ( (g_pps = tmpfile()) == 0) {
		perror("Erreur ouverture fichier temporaire PostScript:");
		exit(-1);
	}
}

/* Hardcopy PostScript */

void SAVECR(savfic,nsavfic)
	char *savfic;
	int nsavfic;

{
	char buf[100], g_pname[60], namfic[50], *ncpr, *ncbl;
	int i;

	if(!v_init) return;

	if(nsavfic >= sizeof(namfic)) nsavfic = sizeof(namfic) - 1;

	strncpy(namfic,savfic,nsavfic);
	namfic[nsavfic] = 0;

	if(ncpr = (char *)strchr(namfic,':'))
	  {
	    *ncpr = 0;
	    snprintf(g_pname,sizeof(g_pname),"%s%04X",namfic,
		     ((unsigned)time(0))&0xFFFF);
	  }
        if(ncbl = (char *)strchr(namfic,' ')) *ncbl = 0;

       	if ( (g_ppsc = fopen(ncpr ? g_pname : namfic,"w")) == 0) {
                perror("Erreur ouverture fichier PostScript:");
		exit(-1);
       	}

	fputs("%!PS-Adobe-2.0 EPSF-2.0\n",g_ppsc);
	fprintf(g_ppsc,"%%%%Title: %s\n",namfic);
	fprintf(g_ppsc,"%%%%Titre: %s\n",namfic);
	fputs("%%Creator: Fminigraf Library\n",g_ppsc);
	fputs("%%BoundingBox: 20 40 567 457\n",g_ppsc);
	fputs("%%EndComments\n",g_ppsc);
	fputs("/offx 20 def\n",g_ppsc);
	fputs("/offy 40 def\n",g_ppsc);
	fputs("/lx 595 def\n",g_ppsc);
	fputs("/ly 842 def\n",g_ppsc);
	fputs("/maxx 1023 def\n",g_ppsc);
	fputs("/maxy 779 def\n",g_ppsc);
	fputs("/reduc 0.92 def\n",g_ppsc);
	fputs("/coeff lx maxx div reduc mul def\n",g_ppsc);
	fputs("offx offy translate\n",g_ppsc);
	fputs("coeff coeff scale\n",g_ppsc);
	fputs("1 setlinewidth newpath 0 0 moveto\n",g_ppsc);
/*	fputs("/Courier-Bold findfont 19 scalefont setfont\n",g_ppsc); */
	fputs("/Courier-Bold findfont 32 scalefont setfont\n",g_ppsc);
	fputs("/M {moveto} def\n",g_ppsc);
	fputs("/L {lineto} def\n",g_ppsc);
	fputs("/S {stroke} def\n",g_ppsc);
	fputs("/H {show} def\n",g_ppsc);
	fputs("/D {setdash} def\n",g_ppsc);
	fputs("0 0 0 setrgbcolor\n",g_ppsc);

	rewind(g_pps);
	while(!feof(g_pps)) {
		fgets(buf,82,g_pps);
		if(feof(g_pps)) break;
		fputs(buf,g_ppsc);
	}
	fputs("showpage\n",g_ppsc);
	fclose(g_ppsc);

	if(ncpr) {
		/* namfic[nsavfic-1] = 0; */
#if defined(SUN) || defined(LINUX) || defined(GFORTRAN4)
		sprintf(buf,"lpr -P%s %s;rm -f %s",namfic,g_pname,g_pname);
#endif
#ifdef HP
		sprintf(buf,"lp  -d%s %s;rm -f %s",namfic,g_pname,g_pname);
#endif

		system(buf);
		free(g_pname);
	}
}

/* Fin fichier PostScript*/

static void finpost()

{
	fclose(g_pps);
}

/* Trace de vecteurs */

static void vecpost(ix1,iy1,ix2,iy2)
	int ix1,iy1,ix2,iy2;

{
	g_ps = TRUE;
	fprintf(g_pps,"%d %d M %d %d L S\n",
		ix1,iy1,ix2,iy2);
}

/* Type de ligne */

static void dashpost(typlin)
	char typlin;

{
	char s[20];
	switch(typlin) {
		case '`' : { strcpy(s,"[] 0"); break; }
		case 'a' : { strcpy(s,"[2 3] 0"); break; }
		case 'b' : { strcpy(s,"[2 3 6 3] 0"); break; }
		case 'c' : { strcpy(s,"[6 3] 0"); break; }
		case 'd' : { strcpy(s,"[10 3] 0"); break; }
		default  : { strcpy(s,"[] 0"); break; }
	}
	fprintf(g_pps,"%s D\n",s);
}

/* Trace de caracteres */

static void textpost(x,y,s)
	int x,y;
	char *s;

{
	g_ps = TRUE;
	if(g_tcar != g_pscar
) {
		fprintf(g_pps,
		   "/Courier-Bold findfont %d scalefont setfont\n",g_hcar[g_tcar]); 
		g_pscar = g_tcar;
	}
	fprintf(g_pps,"%d %d M (%s) H\n",x,y,s);
}
