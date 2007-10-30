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
/********************************************************/
/*          Test si un terminal est "FALCO"             */
/********************************************************/

#include <fcntl.h>
#include <termios.h>
#include <stdio.h>
#include <stdlib.h>

static int inpstr();

static int console;

int tfalco()

{
	char	buf1[80];
	int	n;
	struct termios argp, args;

        /* Ouverture du terminal */
        if ( (console = open("/dev/tty", O_RDWR)) < 0) {
                perror("open");
                exit(-1);
        };

	/* Lecture et sauvegarde des caracteristiques */
	if ( tcgetattr(console, &argp) ) {
		perror("tcgetattr");
		exit(-1);
	};
	args = argp;

	/* Sucrer echo et mode ligne, mettre un time-out */
	argp.c_lflag = (argp.c_lflag & ~ICANON & ~ECHO);
	argp.c_cc[VMIN] = 0;
	argp.c_cc[VTIME] = 5;
	if ( tcsetattr(console, TCSAFLUSH,&argp) ) {
		perror("tcsetattr");
		exit(-1);
	};

	/* Demande a un terminal: id. FALCO ? avec sauvegarde curseur */
	write(console, "\0337\033[<I",6);
	n = inpstr(buf1, 80, 'I');
	buf1[n] = 0;
	if (!strncmp(buf1,"\033[<",3)) buf1[0] = 'E'; else n = 0;
	write(console,"\0338");

	/* Restituer les caracteristiques de la console */
	if ( tcsetattr(console, TCSAFLUSH, &args) ) {
		perror("tcsetattr");
		exit(-1);
	}

	return(n);
}


/* Entree de chaine avec time-out*/
static int inpstr(str,nmax,charfin)

char    str[];
int     nmax;
char    charfin;

{
        char    buf;
        int     n;
        int     compt = 0;

        str[0] = 0;
        do {
                buf = 0;
                n = read(console, &buf, 1);
                if ( (n > 0) && (buf != 0) ) 
                        str[compt++] = buf;
        }        while ( (n > 0) && (compt < nmax) && (buf != charfin) );
        str[compt] = 0;
        if (n <= 0) {
                sleep(1);
                tcflush(console,TCIOFLUSH);
        }
        return( compt );
}

