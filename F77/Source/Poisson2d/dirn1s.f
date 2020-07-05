C
C   Copyright (c) 1997 Silvano Bonazzola
C
C    This file is part of LORENE.
C
C    LORENE is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    LORENE is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with LORENE; if not, write to the Free Software
C    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
	SUBROUTINE DIRN1S(N,NDIM,NEQ,IPAR,MULT,CS,CC)

	IMPLICIT double PRECISION(A-H,O-Z)

C	
C		ROUTINE POUR LA MULTIPLICATION OU DIVISION PAR r DE NEQ
C		FONCTIONS SIMULTANEMENT AVEC ECHANTLLONNAGE RAREFIE' A
C	 	L'ORIGINE.(FONCTION PAIRE OU IMPAIRE),(DI(vision) M(ultiple)
C		R2(echantillonage rarefie') S(imple precision).
C
C		ARGUMENTS DE LA ROUTINE:
C
C		N	=NOMBRE DE DEGRES DE LIBERTE-1
C		NDIM	=DIMENSIONS DES TABLEAUX COMME DECLAREES DANS LE PRO-
C		         GRAMME APPELANT.
C		NEQ	=NOMBRE D'EQUATIONS A TRAITER SIMULTANEMENT
C		IPAR	=PARAMETRE: SI LA FONCTION QUI DOIT ETRE DIVISEE
C			 OU MULTPLIEE EST PAIRE IPAR=0, IPAR=1 DANS LE
C			 CAS CONTRAIRE.
C		MULT	=PARAMETRE: SI MULT=1 LA MULTIPLICATION PAR x EST
C			 EFFECTUEE SI MULT=-1 LA DIVISION EST EFFECTUEE.
C		CS	=TABLEAU CONTENANT LE N+1 COEFFICIENTS DE TCHEBY-
C			 TCHEV DE LA FONCTION QUI DOIT ETRE DIVISEE OU
C			 MULTIPLIEE PAR x.
C		CC	=TABLEAU CONTENANT LES COEFF. DE LA FONCTION DIVISEE.
C
C		(Routine testee le 28/9/1986)
C
C
C $Id: dirn1s.f,v 1.2 2012/03/30 12:12:43 j_novak Exp $
C $Log: dirn1s.f,v $
C Revision 1.2  2012/03/30 12:12:43  j_novak
C Cleaning of fortran files
C
C Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
C LORENE
C
c Revision 1.1  1997/10/21  14:03:48  eric
c Initial revision
c
C
C $Header: /cvsroot/Lorene/F77/Source/Poisson2d/dirn1s.f,v 1.2 2012/03/30 12:12:43 j_novak Exp $
C
C
	character*120 header
	data header/'$Header: /cvsroot/Lorene/F77/Source/Poisson2d/dirn1s.f,v 1.2 2012/03/30 12:12:43 j_novak Exp $'/


	DIMENSION CS(NDIM,*),CC(NDIM,*)
	N1=N+1
	N0=N-1
	X1=1
C
C
C		CAS FONCTIONS FONCTIONS PAIRES
C
	IF(IPAR.EQ.0) THEN
C
C		MULTIPLICATION FONCTIONS PAIRES (OK)
C
	X5=-.5
	IF(MULT.EQ.1) THEN
	DO 1 LEQ=1,NEQ
C	CC(1,LEQ)=-CS(1,LEQ)
	DO 2 L=1,N0
	CC(L,LEQ)=(CS(L,LEQ)+CS(L+1,LEQ))*X5
   2	CONTINUE
	CC(N,NEQ)=CS(N,LEQ)*X5
   1	CONTINUE
	RETURN
	ENDIF
C
C		DIVISION FONCTIONS PAIRES (OK)
C
	IF(MULT.EQ.-1)THEN
C
	DO 24 LEQ=1,NEQ
	CC(1,LEQ)=-CS(1,LEQ)
	DO 5 L=2,N
	CC(L,LEQ)=-(2*CS(L,LEQ)+CC(L-1,LEQ))
   5	CONTINUE
	CC(N1,LEQ)=0
 24	CONTINUE
	RETURN
	ENDIF
	ENDIF
C
C		CAS FONCTIONS IPAIRES.
C
	IF(IPAR.EQ.1) THEN
C
C		MULTIPLICATION FONCTIONS IMPAIRES (OK)
C
	IF(MULT.EQ.1) THEN
C
	X5=-.5
	DO 23 LEQ=1,NEQ
	CC(1,LEQ)=-CS(1,LEQ)
	DO 4 L=2,N
	CC(L,LEQ)=(CS(L-1,LEQ)+CS(L,LEQ))*X5
   4	CONTINUE
	CC(N1,LEQ)=-CS(N,LEQ)
  23	CONTINUE
	RETURN
	ENDIF
C
C
C		DIVISION FONCTIONS IMPAIRES (OK)
C
	IF(MULT.EQ.-1) THEN
C
	DO 3 LEQ=1,NEQ
	CC(N,LEQ)=-CS(N,LEQ)*2
	DO 10 L=N0,1,-1
	CC(L,LEQ)=-(CC(L+1,LEQ)+2*CS(L,LEQ))
   10	CONTINUE
	CC(N1,LEQ)=0
   3	CONTINUE
	RETURN
	ENDIF
	ENDIF
C
	END

