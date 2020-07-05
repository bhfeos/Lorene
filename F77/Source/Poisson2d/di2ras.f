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

	SUBROUTINE DI2RAS(N,NDR,IPAR,MULT,NEQ,CS,CC)

	IMPLICIT NONE
C
C	
C		ROUTINE POUR LA MULTIPLICATION PAR x D' UNE FONCTION
C		AVEC ECHANTLLONAGE RAREFIE' A L'ORIGINE.(FONCTION PAIRE OU
C		IMPAIRE),(DI(vision) M(ultiple) RA(echantillonage rarefie') 
C		S(imple precision).
C
C		ARGUMENTS DE LA ROUTINE:
C
C		N	=NOMBREDE DEGRES DE LIBERTE-1
C		NDR	=DIMENSION DES TABLEAUX COMME DECLAREES DANS LE
C			 PROGRAMME APPELLANT LA ROUTINE
C		IPAR	=PARAMETRE: SI LA FONCTION QUI DOIT ETRE DIVISEE
C			 OU MULTPLIEE EST PAIIRE IPAR=0, IPAR=1 DANS LE
C			 CAS CONTRAIRE.
C		MULT	=PARAMETRE: SI MULT=1 LA MULTIPLICATION PAR x  EST
C			 EFFECTUEE SI MULT=-1 LA DIVISION EST EFFECTUEE.
C		NEQ	=NEQ EST LE NOMBRE DE FONCTIONS QUI DOIVENT ETRE
C			 DIVISEES (OU MULTIPLIEES) PAR x
C		CS(NDR,*)=TABLEAU CONTENANT LE N+1 COEFFICIENTS DE TCHEBY-
C			 TCHEV DE LA FONCTION QUI DOIT ETRE DIVISEE OU
C			 MULTIPLIEE PAR x.
C		CC(NDR,*)	=TABLEAU CONTENANT LES COEFF. DE LA FONCTION DIVISEE.
C			 LES DIMENSIONS DES TABLEAUX CC ET CS DOIBENT
C			 ETRE > (NEQ+1)*N1
C
C
C $Id: di2ras.f,v 1.2 2012/03/30 12:12:43 j_novak Exp $
C $Log: di2ras.f,v $
C Revision 1.2  2012/03/30 12:12:43  j_novak
C Cleaning of fortran files
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.2  1997/05/23  11:30:05  hyc
c *** empty log message ***
c
C Revision 1.1  1997/05/08 07:37:26  hyc
C Initial revision
C
C
C $Header: /cvsroot/Lorene/F77/Source/Poisson2d/di2ras.f,v 1.2 2012/03/30 12:12:43 j_novak Exp $
C
C
	character*120 header
	data header/'$Header: /cvsroot/Lorene/F77/Source/Poisson2d/di2ras.f,v 1.2 2012/03/30 12:12:43 j_novak Exp $'/

	INTEGER NDR,N1,N,IPAR,MULT,J,L,NEQ
C
	double PRECISION CS,CC
C
	DIMENSION CS(NDR,*),CC(NDR,*)
	N1=N+1
C
C		CAS FONCTIONS PAIRES.
C
	IF(IPAR.EQ.0) THEN
C
C		MULTIPLICATION FONCTIONS PAIRES
C
	IF(MULT.EQ.1) THEN
	DO 2 J=1,NEQ
	DO 1 L=1,N
	CC(L,J)=-(CS(L,J)+CS(L+1,J))*.5
   1	CONTINUE
   2	CONTINUE
C
C		LE DERNIER CORFFICIENT ESTPOSE'=0.
C
	DO J=1,NEQ
	CC(N1,J)=0
	ENDDO
C
	RETURN
	ENDIF
C
C		DIVISION FONCTIONS PAIRES
C
	IF(MULT.EQ.-1) THEN
	DO 4 J=1,NEQ
	CC(1,J)=-CS(1,J)
   4	CONTINUE	
	DO 6 J=1,NEQ
	DO 5 L=2,N
	CS(L,J)=-(CS(L,J)*2+CC(L-1,J))
	CC(L,J)=CS(L,J)
   5	CONTINUE
   6	CONTINUE
C
	DO 8 J=1,NEQ
	CC(N1,J)=0
   8	CONTINUE
	ENDIF
	RETURN
	ENDIF
C
C		CAS FONCTIONS FONCTIONS IMPAIRES
C
	IF(IPAR.EQ.1) THEN
C
C		MULTIPLICATION FONCTIONS IMPAIRES
C
	IF(MULT.EQ.1) THEN
	DO 9 J=1,NEQ
	CC(1,J)=-CS(1,J)
   9	CONTINUE
C
	DO 11 J=1,NEQ
	DO 10 L=2,N
	CC(L,J)=-.5*(CS(L-1,J)+CS(L,J))
  10	CONTINUE
  11	CONTINUE
C
	DO 12 J=1,NEQ
	CC(N1,J)=-CS(N,J)
  12	CONTINUE
	RETURN
	ENDIF
C
C		DIVISION FONCTIONS IMPAIRES
C
	IF(MULT.EQ.-1)THEN
	DO 13 J=1,NEQ
	CC(N1,J)=0
  13	CONTINUE
C
	DO 14 J=1,NEQ
	CC(N,J)=-2*CS(N,J)
  14	CONTINUE
C
	DO 15 J=1,NEQ
	DO 16 L=N-1,1,-1
	CS(L,J)=-(CC(L+1,J)+2*CS(L,J))
	CC(L,J)=CS(L,J)
  16	CONTINUE	
  15	CONTINUE
C
	ENDIF
	ENDIF
	RETURN
  100	FORMAT(1X,10E10.3)
  101	FORMAT(1X,' ')
	END
C	
