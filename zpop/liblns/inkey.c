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
/* Fonctions FORTRAN de traitement du clavier*/
#include <errno.h>
#include <fcntl.h>
#include <termios.h>
#include <stdio.h>
#include <signal.h>

#ifdef LINUX

#include <stdlib.h>

#define	FTSTTY		ftstty_
#define	OPEN_KEY	open_key_
#define	TEST_KEY	test_key_
#define	GET_C		get_c_
#define	PUT_C		put_c_
#define	CLOSE_KEY	close_key_
#define	CTRLC_KEY	ctrlc_key_
#define SEND_INT	send_int_
#define FSLEEP          fsleep_
#endif

/* Case LINUX in version originale */
#ifdef XUNIL

#include <stdlib.h>

#define	FTSTTY		ftstty_
#define	OPEN_KEY	open_key__
#define	TEST_KEY	test_key__
#define	GET_C		get_c__
#define	PUT_C		put_c__
#define	CLOSE_KEY	close_key__
#define	CTRLC_KEY	ctrlc_key__
#define SEND_INT	send_int__
#define FSLEEP          fsleep_
#endif

#ifdef GFORTRAN4

#include <stdlib.h>

#define	FTSTTY		ftstty_
#define	OPEN_KEY	open_key_
#define	TEST_KEY	test_key_
#define	GET_C		get_c_
#define	PUT_C		put_c_
#define	CLOSE_KEY	close_key_
#define	CTRLC_KEY	ctrlc_key_
#define SEND_INT	send_int_
#define FSLEEP          fsleep_
#endif

#ifdef SUN

#include <stdlib.h>

#define	FTSTTY		ftstty_
#define	OPEN_KEY	open_key_
#define	TEST_KEY	test_key_
#define	GET_C		get_c_
#define	PUT_C		put_c_
#define	CLOSE_KEY	close_key_
#define	CTRLC_KEY	ctrlc_key_
#define SEND_INT	send_int_
#define FSLEEP          fsleep_
#endif

#ifdef HP
#define	FTSTTY		ftstty
#define	OPEN_KEY	open_key
#define	TEST_KEY	test_key
#define	GET_C		get_c
#define	PUT_C		put_c
#define	CLOSE_KEY	close_key
#define	CTRLC_KEY	ctrlc_key
#define SEND_INT        send_int
#define FSLEEP          fsleep
#endif

static int userterm_;
static int uin = 5;
static struct termios save_userterm_,save_tempterm_;
static int *ptrevt;

/* Test si vrai terminal */
int FTSTTY(id)
int *id;

{
	struct termios argp;

	if ( tcgetattr(*id, &argp) ) {
		if( errno != ENOTTY ) {
			perror("ioctl GET");
			exit(1);
		}
		else {
			return(0); 
		};
	}
	else return(1);
}


/* Ouverture du clavier virtuel*/
void OPEN_KEY()

{
	char term_bp[1024];
	char *p_nom_term;
	struct termios argp;

	/* Recuperer le fd du terminal (stdin)*/
	userterm_ = 0;
        /*userterm_ = getfd_(&uin);*/

	/* Lecture et sauvegarde des caracteristiques */
	if ( tcgetattr(userterm_, &argp) ) {
			perror("ioctl GET");
			exit(1);
	};
	save_userterm_ = argp;

	/* Test du terminal */
	if ( (p_nom_term = getenv("TERM")) == NULL  ) {
		printf("*Terminal inconnu*\n");
		exit(1);
	};
	if (strncmp(p_nom_term,"vt",2) &&
	      strcmp(p_nom_term,"xterm")) {
		printf("* Le terminal doit etre \"xterm\" ou \"VTxx\"*\n");
		exit(1);
	};
	if (tgetent(term_bp, p_nom_term) <= 0) {
		printf("*Terminal inconnu*\n");
		exit(1);
	};
	if (tgetnum("co") < 80) {
		printf("*Le terminal doit avoir au moins 80 col.*\n");
		exit(1);
	};

	/* Sucrer echo et mode ligne, lire 1 car a la fois */
	argp.c_lflag = (argp.c_lflag & ~ICANON & ~ECHO);
	argp.c_iflag = (argp.c_oflag & ~ICRNL); 
	argp.c_cc[VMIN]  = 1;
	argp.c_cc[VTIME] = 0;
        argp.c_cc[VINTR] = 0;
        argp.c_cc[VSUSP] = 0;
        argp.c_cc[VQUIT] = 0;
	save_tempterm_ =argp;
	if ( tcsetattr(userterm_, TCSADRAIN, &argp) ) {
		perror("ioctl SET");
		exit(1);
	};
}

/* Test et lecture d'un caractere */
int TEST_KEY(ptcar)
char *ptcar;

{
	struct termios argp;
	int    n,st;

	/* Mettre le time-out*/
	argp = save_tempterm_;
	argp.c_cc[VMIN]  = 0;
	argp.c_cc[VTIME] = 1; 
	(void) tcsetattr(userterm_, TCSADRAIN, &argp);

	/* Lire le caractere s'il est la */
	n = getc(stdin) ;

	/* Restaurer le terminal */
	(void) tcsetattr(userterm_, TCSADRAIN, &save_tempterm_);

	if (n == -1) st = fseek(stdin,0L,0); else *ptcar=(char)n;
 
	return( (n == -1) ? 0 : 1);
}


/* Lecture caractere*/
void GET_C(ptcar)
char *ptcar;

{
	*ptcar = (char)getc(stdin);
}

/* Ecriture caractere*/
void PUT_C(ptcar)
char *ptcar;

{
	(void) putc(*ptcar,stdout);
}

/* Cloture du clavier */
void CLOSE_KEY()

{
	if ( tcsetattr(userterm_, TCSADRAIN, &save_userterm_)) {
		perror("ioctl SET");
		exit(1);
	};
}

/* Handler du ctrl/c */
void ctrlc_han()

{
	*ptrevt = 1;
	signal(SIGINT, ctrlc_han);
}

/* Armement du ctrl/c */
void CTRLC_KEY(adrevt)
int *adrevt;

{
	ptrevt = adrevt;
	signal(SIGINT, ctrlc_han);
}

/* Envoi du signal INT */
void SEND_INT()

{
	(void)kill(getpid(),SIGINT);
}

/* Interface Fortran pour sleep */
void FSLEEP(time)
	int *time;

{
	sleep(*time);
}
