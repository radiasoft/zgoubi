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

  François Méot <meot@lpsc.in2p3.fr>
  Service Accélerateurs
  LPSC Grenoble
  53 Avenue des Martyrs
  38026 Grenoble Cedex
  France
*/
/*          ROUTINES DE GESTION DU TERMINAL EN FORTRAN           */

/* INITRM: init du terminal pour utiliser les autres routines    */

/* WRSTAT: ecriture en Fortran dans la derniere ligne de l'ecran */
/* Ce programme ne marche qu'avec les terminaux de type VTXXX    */
/* ou "xterm"                                                    */
/* Il suppose que l'on a appele au moins une fois le prog. tterm */
/* Ceci est fait normalement par le .login standard du L.N.S.    */
/* ou par appel a INITRM                                         */

/* GETPAG: obtention du nombre de lignes de la page              */

#include <stdlib.h>
#include <curses.h>
#include <stdio.h>

#if defined(LINUX) || defined(GFORTRAN4)
#define WRSTAT wrstat_
#define INITRM initrm_
#define GETPAG getpag_
#define TTERM  "clear"
#endif

#ifdef SUN
#define WRSTAT wrstat_
#define INITRM initrm_
#define GETPAG getpag_
#define TTERM  "clear;/home/bin/tterm"
#endif

#ifdef HP
#define WRSTAT wrstat
#define INITRM initrm
#define GETPAG getpag
#define TTERM  "clear;/home/bin/tterm"
#endif


/***** INITRM *****/
void INITRM()

{
	(void)system(TTERM);
}

/***** WRSTAT *****/
void WRSTAT(pbuf,lg)
char pbuf[];
int  lg;

{
	char term_bp[1024], *pterm;
	int nli, i, nc1, nc2;
	FILE *cons;

	/* Trouver le terminal */
	pterm = getenv("TERM");
	if (tgetent(term_bp, pterm) <= 0)  return;
	if (!strncmp(pterm,"vt",2) || !strcmp(pterm,"xterm")) {
		/* Traitement terminal VTxxx */

		/* Trouver le nb de lignes */
		if ((nli = tgetnum("lines")) <= 0) return;
		nli++;

		/* Sauver et positionner le curseur, puis inverse video */
		printf("\0337\033[%d;1H\033[7m",nli);

		/* Ecrire les caracteres et centrer s'il y a lieu */
		if(lg >= 80) {
			lg = 80;
			for(i = 0; i < lg; i++) putchar(pbuf[i]);
		}
		else {
			nc1 = (80 - lg) / 2;
			nc2 = (80 - lg - nc1);
			for(i = 0; i < nc1; i++) putchar(' ');
			for(i = 0; i <  lg; i++) putchar(pbuf[i]);
			for(i = 0; i < nc2; i++) putchar(' ');
		}

		/* Restaurer la video et le curseur */
		printf("\033[0m\0338");
	}
}

/***** GETPAG *****/
GETPAG(nline)
int *nline;

{
        char term_bp[1024];
        char *p_nom_term;
        int nli;

        /* acquisition terminal */
        p_nom_term = getenv("TERM");

        tgetent(term_bp, p_nom_term);

        *nline = ((nli = tgetnum("lines")) <= 0) ? 24 : nli;
}
