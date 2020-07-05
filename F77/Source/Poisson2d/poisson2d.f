C
C   Copyright (c) 1998 Eric Gourgoulhon
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
	SUBROUTINE POISS2D(NDL1,NDR,NDT,NDF,INDD,ERRE,SOUMAT,SOUQUAD,
     1			     ALAMB,POT)

C
C $Id: poisson2d.f,v 1.8 2013/09/04 14:12:10 j_novak Exp $
C $Log: poisson2d.f,v $
C Revision 1.8  2013/09/04 14:12:10  j_novak
C Removed some comments.
C
C Revision 1.7  2013/09/04 13:56:53  j_novak
C Dynamical memory allocation for working arrays.
C
C Revision 1.6  2012/06/08 12:08:35  j_novak
C Increase of NDR0
C
C Revision 1.5  2012/03/30 12:12:44  j_novak
C Cleaning of fortran files
C
C Revision 1.4  2010/12/20 10:05:28  m_bejger
C Increase of NDZ0 from 5 to 8
C
C Revision 1.3  2002/09/09 13:50:28  e_gourgoulhon
C
C Change of subroutine names:
C 	POISSON2D -> POISS2D
C 	POISSON2DI -> POISS2DI
C to avoid any clash with Map::poisson2d and Map::poisson2di.
C
C Revision 1.2  2002/03/25 09:16:59  m_bejger
C Increased the number of domains (NZOE) from 4 to 5
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.2  1998/07/20  12:43:52  eric
c Augmentation NDR0, NDT0
c
C Revision 1.1  1998/07/01 13:27:06  eric
C Initial revision
C
C
C $Header: /cvsroot/Lorene/F77/Source/Poisson2d/poisson2d.f,v 1.8 2013/09/04 14:12:10 j_novak Exp $
C
C

		IMPLICIT NONE

	INTEGER NDL1(*), NDR, NDT, NDF, INDD(*)
	REAL*8 ERRE(NDR,*), SOUMAT(NDR,NDT,NDF,*)
	REAL*8 SOUQUAD(NDR,NDT,NDF,*), ALAMB, POT(NDR,NDT,NDF,*)

	character*120 header
	data header/'$Header: /cvsroot/Lorene/F77/Source/Poisson2d/poisson2d.f,v 1.8 2013/09/04 14:12:10 j_novak Exp $'/

	INTEGER N64, ND64Q, ND2Z, NDEQ

	REAL*8,allocatable::  TRA0(:,:),TRA1(:,:),TRA2(:,:),TRA3(:,:)
	REAL*8,allocatable::  TRAB0(:,:,:), TRAB1(:,:,:)
	REAL*8,allocatable::  DEN1(:,:,:), DEN2(:,:,:)
	REAL*8,allocatable::  BB(:,:), ERRE0(:,:), SOLHH(:,:,:,:)
	REAL*8,allocatable::  CC(:), CS(:), C64(:)

	INTEGER,ALLOCATABLE:: NDL(:)

	INTEGER IND, LZON, NR1, LR, LY, NZOE, NZON, NY1, NF, LZ, LT, LF
	REAL*8 AMAQ, AMAT, X0, X1, C1

C******************************************************************************

C 
C... Allocation memoire
C
	ND64Q = (NDR+2)*(NDT+2)*NDF
	NZOE = NDL1(1)
	ND2Z = MAX(NZOE, NDF, 8)
	NDEQ = NZOE + 8
	allocate(CC(ND64Q), CS(ND64Q), C64(ND64Q) )
	allocate(TRA0(NDR,NDT))
	allocate(TRA1(NDR,NDT))
	allocate(TRA2(NDR,NDT))
	allocate(TRA3(NDR,NDT))
	allocate(TRAB0(NDR,NDT,ND2Z),TRAB1(NDR,NDT,ND2Z))
	allocate(DEN1(NDR,NDT,ND2Z), DEN2(NDR,NDT,ND2Z) )
	allocate(BB(NDR,12), ERRE0(NDR, ND2Z) )
	allocate(SOLHH(NDR, NDT, 8, NZOE) )
	N64 = 20 
	ALLOCATE(NDL(NDEQ))

	NZON = NZOE - 1
	NY1 = NDL1(NZOE+2)
	NF = NDL1(NZOE+3)

C
C... Tableau NDL decrivant les zones non compactifies:
C
	NDL(1) = NZON 
		DO LZON = 1, NZON
	NDL(1+LZON) = NDL1(1+LZON)
		ENDDO
	NDL(NZON+2)=NY1
	NDL(NZON+3)=NF	
	NDL(NZON+4)=NDL1(NZOE+4)

C
C... Tableau ERRE0 decrivant les points de collocation en r dimensionne comme
C		(NDR0,*), contrairement a ERRE qui est dimensionne comme (NDR,*)
C
		DO LZON=1,NZOE
	NR1=NDL1(LZON+1)
			DO LR=1,NR1
	ERRE0(LR,LZON)=ERRE(LR,LZON)
			ENDDO
	ENDDO

C
C Mise a zero des tableaux de travail
C
	DO LR = 1, ND64Q
		CC(LR) = 0
		CS(LR) = 0 
		C64(LR) = 0
	ENDDO

	DO LT = 1, NDT
		DO LR = 1, NDR
			TRA0(LR,LT) = 0
			TRA1(LR,LT) = 0
			TRA2(LR,LT) = 0
			TRA3(LR,LT) = 0
		ENDDO
	ENDDO
	
	DO LZ = 1, ND2Z
		DO LT = 1, NDT
			DO LR = 1, NDR
				TRAB0(LR,LT,LZ) = 0
				TRAB1(LR,LT,LZ) = 0
				DEN1(LR,LT,LZ) = 0
				DEN2(LR,LT,LZ) = 0
			ENDDO
		ENDDO
	ENDDO

	DO LZ = 1, 12
		DO LR = 1, NDR
			BB(LR,LZ) = 0 
		ENDDO
	ENDDO

	DO LZ = 1, NZOE
		DO LF = 1, 8
			DO LT = 1, NDT
				DO LR = 1, NDR
					SOLHH(LR,LT,LF,LZ) = 0
				ENDDO
			ENDDO
		ENDDO
	ENDDO


C
C Les termes quadratiques sont contenus dans SOUQUAD; on les passe dans l'espace
C  des coefficients en theta et en r:
C
C
	CALL FCEZ3S(NDL1,NDR,NDT,NDF,N64,2,0,C64,CC,CS,DEN2,DEN1,SOUQUAD)
	CALL FCEZ3S(NDL1,NDR,NDT,NDF,N64,1,0,C64,CC,CS,DEN2,DEN1,SOUQUAD)
C
C Terme "monopolaire" l=0 du developpement de SOUQUAD en cos(l*theta) ---> TRA0:
C 
	DO LZON=1,NZOE
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	TRA0(LR,LZON)=SOUQUAD(LR,1,1,LZON)
	ENDDO
	ENDDO

C  Calcul de la masse generee par le termes quadratiques
C
C Calcul de   AMAQ = int_0^{+infini} SOUQUAD(r)(l=0) r dr   [cf. Eq.(4.6)] 
C
	IND=2
	CALL GPAR2S(NDL1,NDR,IND,C64,ERRE0,TRA0,TRA1)

C	PRINT*,'TRA1=',TRA1(1,1),TRA1(2,1)
	AMAQ=TRA1(1,1)+TRA1(2,1)

C
C
C Les termes matiere sont contenus dans SOUMAT
C  Passage en Tchebyshev en theta et en r:
C
	CALL FCEZ3S(NDL,NDR,NDT,NDF,N64,2,0,C64,CC,CS,DEN2,DEN1,SOUMAT)
	CALL FCEZ3S(NDL,NDR,NDT,NDF,N64,1,0,C64,CC,CS,DEN2,DEN1,SOUMAT)
C
C Terme "monopolaire" l=0 du developpement de SOUMAT en cos(l*theta) ---> TRA2:
C 
	DO LZON=1,NZON
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	TRA2(LR,LZON)=SOUMAT(LR,1,1,LZON)
	ENDDO
	ENDDO
C
C Calcul de   AMAT = int_0^{+infini} SOUMAT(r)(l=0) r dr   [cf. Eq.(4.6)] 

	IND=1
	CALL GPAR2S(NDL,NDR,IND,C64,ERRE0,TRA2,TRA3)
	AMAT=TRA3(1,1)

C
C----------------------------------------
C Identite du viriel GRV2 [ Eq.(4.6), Eq.(6.36) ] 
C----------------------------------------
C   ALAMB = Abs(lambda), lambda etant defini par l'Eq.(6.38)
C    ALAMB doit etre egal a 1 pour une solution exacte
C
	X0=0
	X1=1 
	ALAMB=0

C##	IF(AMAQ.NE.X0) ALAMB=ABS(AMAT/AMAQ)
	IF(AMAQ.NE.X0) ALAMB = - AMAT / AMAQ 
C
C Preparation de la source totale pour dzeta suivant l'Eq.(6.37):
C
	DO LZON=1,NZOE
	NR1=NDL1(LZON+1)
	DO LY=1,NY1
	DO LR=1,NR1
		POT(LR,LY,1,LZON) = ALAMB * SOUQUAD(LR,LY,1,LZON)
	ENDDO
	ENDDO
	ENDDO
C
C ... On ajoute la matiere dans les zones non compactifiees : 
	DO LZON=1,NZON
	NR1=NDL1(LZON+1)
	DO LY=1,NY1
	DO LR=1,NR1
		POT(LR,LY,1,LZON) = POT(LR,LY,1,LZON) + SOUMAT(LR,LY,1,LZON)
	ENDDO
	ENDDO
	ENDDO

	IND=4
C
	DO LZON=1,NZOE
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	TRA0(LR,LZON)=POT(LR,1,1,LZON)
	ENDDO
	ENDDO
C
	CALL GR2P1S(NDL1,NDR,2,ERRE0,TRA2,TRA0)
C
	CALL GR2P3S(NDL1,NDR,NDT,NDF,INDD,0,CC,C64,BB,DEN1,
     1		    DEN2,TRAB0,TRAB1,ERRE,SOLHH,POT)
C
	DO LZON=1,NZOE
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	POT(LR,1,1,LZON)=TRA0(LR,LZON)
	ENDDO
	ENDDO
C
C
	CALL FCIQ3S(NDL1,NDR,NDT,NDF,INDD,1,0,0,2,C64,CC,CS,DEN1, 
     1	DEN2,ERRE,POT)
C
C
C Constante d'integration t.q. valeur nulle a l'infini
C
	NR1=NDL1(NZOE+1)
	C1=POT(NR1,1,1,NZOE) ! valeur du coef l=0 de dzeta a l'infini 
C
	DO LZON=1,NZOE
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	POT(LR,1,1,LZON)=POT(LR,1,1,LZON)-C1
	ENDDO
	ENDDO

CC
C Retour dans l'espace des configurations en theta (on y etait deja en r):
C
	CALL FCIQ3S(NDL1,NDR,NDT,NDF,INDD,2,0,0,2,C64,CC,CS,DEN1,
     1	DEN2,ERRE,POT)

	deallocate(TRA0, TRA1, TRA2, TRA3)
	DEALLOCATE(TRAB0, TRAB1, DEN1, DEN2)
	DEALLOCATE(BB, ERRE0, SOLHH)
	DEALLOCATE(CC, CS, C64)
	DEALLOCATE(NDL)

	RETURN 

	END

