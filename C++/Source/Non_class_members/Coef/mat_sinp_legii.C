/*
 *   Copyright (c) 2003 Jerome Novak
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */


 

/*
 * Fournit la matrice de passage pour la transformation des coefficients du
 * developpement en sin(2j*theta) 
 * dans les coefficients du developpement en fonctions associees de Legendre 
 * P_l^m(cos(theta)) avec l pair et m impair.
 * 
 * Cette routine n'effectue le calcul de la matrice que si celui-ci n'a pas 
 * deja ete fait, sinon elle renvoie le pointeur sur une valeur precedemment 
 * calculee.
 * 
 * Entree:
 * -------
 *  int np  :	Nombre de degres de liberte en phi 
 *  int nt  :	Nombre de degres de liberte en theta
 *
 * Sortie (valeur de retour) :
 * ---------------------------
 *  double* mat_sinp_legii : pointeur sur le tableau contenant l'ensemble
 *			        (pour les np/2 valeurs de m: m=1,3,...,np-1) des 
 *				matrices de passage.
 *				La dimension du tableau est (np/2+1)*nt^2
 *				Le stokage est le suivant: 
 *
 *		mat_sinp_legii[ nt*nt* m/2 + nt*l + j] = A_{mlj}
 *			    
 *				ou A_{mlj} est defini par
 *
 *   sin(2*j*theta) = som_{l=(m+1)/2}^{nt-2} A_{mlj} P_{2l}^m( cos(theta) )
 *						pour 1 <= j <= nt-2
 *
 *  ou P_n^m(x) represente la fonction de Legendre associee de degre n et 
 *     d'ordre m normalisee de facon a ce que
 *
 *	 int_0^pi [ P_n^m(cos(theta)) ]^2  sin(theta) dtheta = 1
 *
 * 
 */

/*
 * $Id: mat_sinp_legii.C,v 1.5 2016/12/05 16:18:02 j_novak Exp $
 * $Log: mat_sinp_legii.C,v $
 * Revision 1.5  2016/12/05 16:18:02  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:14  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:16:03  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2005/02/18 13:14:15  j_novak
 * Changing of malloc/free to new/delete + suppression of some unused variables
 * (trying to avoid compilation warnings).
 *
 * Revision 1.1  2003/09/16 08:58:01  j_novak
 * New functions for the T_LEG_II base
 *
 *
 * $Header: /cvsroot/Lorene/C++/Source/Non_class_members/Coef/mat_sinp_legii.C,v 1.5 2016/12/05 16:18:02 j_novak Exp $
 *
 */

// headers du C
#include <cstdlib>
#include <cmath>
#include <cassert>

// Prototypage
#include "headcpp.h"
#include "proto.h"

namespace Lorene {
//******************************************************************************

double* mat_sinp_legii(int np, int nt) {

#define NMAX	30		// Nombre maximun de couples(np,nt) differents 
static	double*	tab[NMAX] ;	// Tableau des pointeurs sur les tableaux 
static	int	nb_dejafait = 0 ;	// Nombre de tableaux deja initialises 
static	int	np_dejafait[NMAX] ;    // Valeurs de np pour lesquelles le
				       // calcul a deja ete fait
static	int	nt_dejafait[NMAX] ;    // Valeurs de np pour lesquelles le
				       // calcul a deja ete fait

int i, indice,  j,  j2,  m,  l ;
    
    {
	
    // Les matrices A_{mlj} pour ce couple (np,nt) ont-elles deja ete calculees ? 
    indice = -1 ;
    for ( i=0 ; i < nb_dejafait ; i++ ) {
	if ( (np_dejafait[i] == np) && (nt_dejafait[i] == nt) ) indice = i ;
    }


// Si le calcul n'a pas deja ete fait, il faut le faire : 
    if (indice == -1) {		   
	if ( nb_dejafait >= NMAX ) {
	    cout << "mat_sinp_legii: nb_dejafait >= NMAX : "
		<< nb_dejafait << " <-> " << NMAX << endl ;
	    abort () ;	
	    exit(-1) ;
	    }
	indice = nb_dejafait ; 
	nb_dejafait++ ; 
	np_dejafait[indice] = np ;
	nt_dejafait[indice] = nt ;

	tab[indice] = new double[(np/2+1)*nt*nt] ;//(double *) malloc( sizeof(double) * (np/2+1)*nt*nt ) ;

//-----------------------
// Preparation du calcul 
//-----------------------

// Sur-echantillonnage pour calculer les produits scalaires sans aliasing:
	int nt2 = 2*nt - 1 ; 
	int nt2m1 = nt2 - 1 ; 

	int deg[3] ;
	deg[0] = 1 ;
	deg[1] = 1 ;
	deg[2] = nt2 ;

// Tableaux de travail
	double* yy = new double[nt2] ;//(double*)( malloc( nt2*sizeof(double) ) ) ;
	double* sint = new double[nt*nt2];//(double*)( malloc( nt*nt2*sizeof(double) ) ) ;
   
// Calcul des sin( 2j theta)  aux points de collocation
//  de l'echantillonnage double : 

	double dt = M_PI / double(2*(nt2-1)) ;
	for (j=0; j<nt-1; j++) {
	    for (j2=0; j2<nt2; j2++) {
		double theta = j2*dt ;
		sint[nt2*j + j2] = sin( 2*j * theta ) ;
	    }
	}	


//-------------------
// Boucle sur m
//-------------------

	int m_max = (np == 1) ? 1 : np-1 ; 

	for (m=1; m <= m_max ; m+=2) {

	    // Recherche des fonctions de Legendre associees d'ordre m :

	    double* leg = legendre_norm(m, nt) ;	
	
	    for (l=(m+1)/2; l<nt-1; l++) {	// boucle sur les P_{2l}^m
	    
		int ll = 2*l ;	// degre des fonctions de Legendre
		    
		for (j=0; j<nt-1; j++) {  // boucle sur les sin((2j+1)theta)
			
		//... produit scalaire de sin(2j theta) par 
		//				    P_{2l}^m(cos(theta)) 

		    for (j2=0; j2<nt2; j2++) {
			yy[nt2m1-j2] = sint[nt2*j + j2] 
						    * leg[nt2* (ll-m) + j2] ;
		    }

//....... on passe en Tchebyshev vis-a-vis de x=cos(theta) pour	calculer
//	    l'integrale (routine int1d_chebp) : 		
		    cfrchebp(deg, deg, yy, deg, yy) ;
		    tab[indice][ nt*nt* ((m-1)/2) + nt*l + j] = 
						2.*int1d_chebp(nt2, yy) ;


		}	// fin de la boucle sur j  (indice de  sin((2j)theta) )
	    
	    }  // fin de la boucle sur l (indice de P_{2l}^m)
	    	    
	delete [] leg ;
	   	    
	}  // fin de la boucle sur m

// Liberation espace memoire
// -------------------------

	delete [] yy ;
	delete [] sint ;

    } // fin du cas ou le calcul etait necessaire

    } //Fin de zone critique
    
    return tab[indice] ;

}


}
