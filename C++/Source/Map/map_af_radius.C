/*
 *  Methods of the class Map_af relative to the function
 *	    r = R_l(xi, theta', phi')
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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
 * $Id: map_af_radius.C,v 1.10 2016/12/05 16:17:57 j_novak Exp $
 * $Log: map_af_radius.C,v $
 * Revision 1.10  2016/12/05 16:17:57  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.9  2014/10/13 08:53:03  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.8  2014/10/06 15:13:12  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.7  2013/06/05 15:10:42  j_novak
 * Suppression of FINJAC sampling in r. This Jacobi(0,2) base is now
 * available by setting colloc_r to BASE_JAC02 in the Mg3d constructor.
 *
 * Revision 1.6  2008/09/01 08:12:03  j_novak
 * Improved test on the [rmin, rmax] interval.
 *
 * Revision 1.5  2007/12/11 15:28:14  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.4  2006/09/13 13:59:21  j_novak
 * Higher tolerance thereshold for Map_af::val_lx
 *
 * Revision 1.3  2006/07/10 07:44:51  j_novak
 * Correction of the comparison between rmin and rr (now has to be greater than
 * some threshold).
 *
 * Revision 1.2  2006/07/05 12:36:51  n_vasset
 * Added a test on rmin to see whether the point lies in the computational domains.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.5  1999/12/16  14:19:08  eric
 * Introduction de l'argument const Param& par dans val_lx et val_lx_jk.
 * (en remplacement de l'argument Tbl& param).
 *
 * Revision 1.4  1999/12/07  14:51:37  eric
 * val_r_kj --> val_r_jk
 * val_lx_kj -->val_lx_jk
 * Changement ordre des arguments val_r, val_lx
 *
 * Revision 1.3  1999/12/06  16:47:21  eric
 * Surcharge de val_lx avec la version sans param.
 *
 * Revision 1.2  1999/12/06  15:34:06  eric
 * Ajout des fonctions val_r_kj et val_lx_kj.
 *
 * Revision 1.1  1999/12/06  13:12:16  eric
 * Initial revision
 *
 *
 * $Header: /cvsroot/Lorene/C++/Source/Map/map_af_radius.C,v 1.10 2016/12/05 16:17:57 j_novak Exp $
 *
 */

#include <cmath>

// Headers Lorene
#include "map.h"

			//------------------------------// 
			//	    val_r		//
			//------------------------------// 

 
namespace Lorene {
double Map_af::val_r(int l, double xi, double, double) const {

    assert( l>=0 ) ; 
    assert( l<mg->get_nzone() ) ; 
    
    double resu ; 

    switch( mg->get_type_r(l) ) {

	case FIN: case RARE: {
	    resu = alpha[l] * xi + beta[l] ;
	    break ;
	}
	
	case UNSURR: {
	    resu = double(1) / ( alpha[l] * xi + beta[l] ) ;
	    break ;
	}

	default: {
	    cout << "Map_af::val_r: unknown type_r ! " << endl ;
	    abort () ;
	}	   
    }
             
    return resu ;    
}
			
			//------------------------------// 
			//	    val_lx		//
			//------------------------------// 

void Map_af::val_lx(double rr, double, double, int& lz, double& xi) const {
			   
    // In which domain is located r ? 
    // ----------------------------
    int nz = mg->get_nzone() ;
    lz = - 1 ;
    
    for (int l=0; l<nz; l++) {
        
	double rmax = alpha[l] + beta[l] ;
	double rmin = beta[l] - alpha[l] ;
	if (mg->get_type_r(l) == RARE) rmin = 0. ;
	if (mg->get_type_r(l) == UNSURR) {
	    rmin = double(1)/rmin ; 
	    rmax = double(1)/rmax ; 
	}		
	if ((rr - rmin >= -1.e-14*fabs(rmin)) && ( rr <= rmax )) { 
	    lz = l ;
	    break ; 
	}	
    }		// fin de la boucle sur les zones
    
    if (lz == -1) {		    // On n'a pas trouve la zone 
	cout.precision(16);
	cout.setf(ios::showpoint);
	cout << "Map_af::val_lx: the domain containing r = " << rr <<
		 " has not been found ! " 
	      << endl ;
	for (int l=0; l<nz; l++) {
	    double rmin = -alpha[l] + beta[l] ;	
	    if (mg->get_type_r(l) == UNSURR) rmin = double(1)/rmin ; 
	    if (mg->get_type_r(l) == RARE) rmin = 0. ;
	    cout << "domain " << l << " :  r_min = " << rmin ; 
	    double rmax = alpha[l] + beta[l] ;	
	    if (mg->get_type_r(l) == UNSURR) rmax = double(1)/rmax ; 
	    cout << " :  r_max = " << rmax << endl ; 
	}
	abort () ;
    }

    // Computation of xi
    // ----------------- 

    switch( mg->get_type_r(lz) ) {

	case FIN: case RARE: {
	    xi = ( rr - beta[lz] ) / alpha[lz]  ;
	    break ;
	}
	
	case UNSURR: {
	    xi = ( double(1)/rr - beta[lz] ) / alpha[lz]  ;
	    break ;
	}

	default: {
	    cout << "Map_af::val_lx: unknown type_r ! " << endl ;
	    abort () ;
	}	   
    }

} 


void Map_af::val_lx(double rr, double, double, const Param&,  
			    int& lz, double& xi) const {

    val_lx(rr, 0., 0., lz, xi) ;

}


			//------------------------------// 
			//	    val_r_jk		//
			//------------------------------// 

 
double Map_af::val_r_jk(int l, double xi, int, int) const {

    return val_r(l, xi, 0., 0.) ; 
    
}
			
			//------------------------------// 
			//	    val_lx_jk		//
			//------------------------------// 

void Map_af::val_lx_jk(double rr, int, int, const Param& par, 
			    int& l, double& xi) const {
			   
    val_lx(rr, 0., 0., par, l, xi) ; 

} 


}
