/*
 *  Methods of class Sym_tensor
 *
 *   (see file tensor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 1999-2001 Philippe Grandclement (Cmp version)
 *   Copyright (c) 2000-2001 Eric Gourgoulhon (Cmp version)
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
 * $Id: sym_tensor.C,v 1.25 2016/12/05 16:18:17 j_novak Exp $
 * $Log: sym_tensor.C,v $
 * Revision 1.25  2016/12/05 16:18:17  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.24  2014/10/13 08:53:43  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.23  2014/10/06 15:13:19  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.22  2007/12/21 16:07:08  j_novak
 * Methods to filter Tensor, Vector and Sym_tensor objects.
 *
 * Revision 1.21  2007/12/03 13:00:00  n_vasset
 * Adjusting memory management for new member p_tilde_c
 *
 * Revision 1.20  2006/06/12 07:27:20  j_novak
 * New members concerning A and tilde{B}, dealing with the transverse part of the
 * Sym_tensor.
 *
 * Revision 1.19  2005/04/04 15:25:24  j_novak
 * Added new members www, xxx, ttt and the associated methods.
 *
 * Revision 1.18  2005/04/01 14:28:32  j_novak
 * Members p_eta and p_mu are now defined in class Sym_tensor.
 *
 * Revision 1.17  2004/03/30 14:01:19  j_novak
 * Copy constructors and operator= now copy the "derived" members.
 *
 * Revision 1.16  2004/02/26 22:48:50  e_gourgoulhon
 * -- Method divergence: call to Tensor::divergence and cast of the
 *    result.
 * -- Added method derive_lie.
 *
 * Revision 1.15  2004/01/04 20:54:00  e_gourgoulhon
 * Sym_tensor is now a derived class of Tensor_sym.
 * Methods indices and position have been suppressed (they are now
 * implemented at the Tensor_sym level).
 *
 * Revision 1.14  2003/12/30 23:09:47  e_gourgoulhon
 * Change in methods derive_cov() and divergence() to take into account
 *  the change of name: Metric::get_connect() --> Metric::connect().
 *
 * Revision 1.13  2003/11/26 21:58:15  e_gourgoulhon
 * Added new data member p_transverse and p_longit_pot.
 * Modified the memory management consequently.
 *
 * Revision 1.12  2003/10/28 12:34:08  e_gourgoulhon
 * Corrected bug in the copy constructor and constructor from Tensor:
 * the cmp have already been created by the (special) Tensor constructor called
 * by these constructors.
 *
 * Revision 1.11  2003/10/20 14:26:03  j_novak
 * New assignement operators.
 *
 * Revision 1.10  2003/10/16 14:21:36  j_novak
 * The calculation of the divergence of a Tensor is now possible.
 *
 * Revision 1.9  2003/10/13 13:52:39  j_novak
 * Better managment of derived quantities.
 *
 * Revision 1.8  2003/10/11 16:47:10  e_gourgoulhon
 * Suppressed the call to Ibtl::set_etat_qcq() after the construction
 * of the Itbl's, thanks to the new property of the Itbl class.
 *
 * Revision 1.7  2003/10/07 09:56:59  j_novak
 * method Sym_tensor::inverse() implemented (but not tested!)
 *
 * Revision 1.6  2003/10/06 13:58:48  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.5  2003/10/03 11:21:48  j_novak
 * More methods for the class Metric
 *
 * Revision 1.4  2003/10/02 15:45:51  j_novak
 * New class Metric
 *
 * Revision 1.3  2003/10/01 15:39:43  e_gourgoulhon
 * Added assert to insure that both indices have the same type.
 *
 * Revision 1.2  2003/09/26 08:05:31  j_novak
 * New class Vector.
 *
 * Revision 1.1  2003/09/25 13:37:40  j_novak
 * Symmetric tensors of valence 2 are now implemented (not tested yet).
 *
 *
 * $Header: /cvsroot/Lorene/C++/Source/Tensor/sym_tensor.C,v 1.25 2016/12/05 16:18:17 j_novak Exp $
 *
 */

// Headers C
#include <cstdlib>
#include <cassert>
#include <cmath>

// Headers Lorene
#include "metric.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor 
// --------------------
namespace Lorene {
Sym_tensor::Sym_tensor(const Map& map, const Itbl& tipe, 
                       const Base_vect& triad_i) 
            : Tensor_sym(map, 2, tipe, triad_i, 0, 1) {
		
    set_der_0x0() ;

}

// Standard constructor when all the indices are of the same type
// --------------------------------------------------------------
Sym_tensor::Sym_tensor(const Map& map, int tipe, const Base_vect& triad_i)  
            : Tensor_sym(map, 2, tipe, triad_i, 0, 1) {

    set_der_0x0() ;
}

// Copy constructor
// ----------------
Sym_tensor::Sym_tensor(const Sym_tensor& source) 
            : Tensor_sym( source ) {

    set_der_0x0() ;

    for (int i_met = 0; i_met < N_MET_MAX; i_met++) {

      if ( source.p_transverse[i_met] != 0x0 ) {
	set_dependance( *source.met_depend[i_met] ) ;
	int jp = get_place_met( *source.met_depend[i_met] ) ;
	assert ((jp>=0) && (jp<N_MET_MAX)) ;
	p_transverse[jp] = 
	new Sym_tensor_trans ( *source.p_transverse[i_met] ) ;
      }

      if ( source.p_longit_pot[i_met] != 0x0 ) {
	set_dependance( *source.met_depend[i_met] ) ;
	int jp = get_place_met( *source.met_depend[i_met] ) ;
	assert ((jp>=0) && (jp<N_MET_MAX)) ;
	p_longit_pot[jp] = 
	new Vector ( *source.p_longit_pot[i_met] ) ;
      }
    }
    if (source.p_eta != 0x0) p_eta = new Scalar( *(source.p_eta) ) ; 
    if (source.p_mu != 0x0) p_mu = new Scalar( *(source.p_mu) ) ; 
    if (source.p_www != 0x0) p_www = new Scalar( *(source.p_www) ) ; 
    if (source.p_xxx != 0x0) p_xxx = new Scalar( *(source.p_xxx) ) ; 
  
}   


// Constructor from a Tensor
// --------------------------
Sym_tensor::Sym_tensor(const Tensor& source) 
            : Tensor_sym(*source.mp, 2, source.type_indice, *(source.triad), 
                         0, 1) {
	
    assert(source.valence == 2) ;
	
    for (int ic=0 ; ic<n_comp ; ic++) {
        int posi = source.position(indices(ic)) ;
        *(cmp[ic]) = *(source.cmp[posi]) ;
    }
	    
    set_der_0x0() ;
}   

	
// Constructor from a file
// -----------------------
Sym_tensor::Sym_tensor(const Map& map, const Base_vect& triad_i, FILE* fd)
			: Tensor_sym(map, triad_i, fd) {
	
	assert (valence == 2) ;
	assert (n_comp == 6) ;
	set_der_0x0() ;
}

			//--------------//
			//  Destructor  //
			//--------------//

Sym_tensor::~Sym_tensor() {

  Sym_tensor::del_deriv() ;

}



			//--------------//
			//  Assignment  //
			//--------------//

void Sym_tensor::operator=(const Sym_tensor& source) {
    
    Tensor_sym::operator=(source) ; 
    
    del_deriv() ;

    for (int i_met = 0; i_met < N_MET_MAX; i_met++) {

      if ( source.p_transverse[i_met] != 0x0 ) {
	set_dependance( *source.met_depend[i_met] ) ;
	int jp = get_place_met( *source.met_depend[i_met] ) ;
	assert ((jp>=0) && (jp<N_MET_MAX)) ;
	p_transverse[jp] = 
	new Sym_tensor_trans ( *source.p_transverse[i_met] ) ;
      }

      if ( source.p_longit_pot[i_met] != 0x0 ) {
	set_dependance( *source.met_depend[i_met] ) ;
	int jp = get_place_met( *source.met_depend[i_met] ) ;
	assert ((jp>=0) && (jp<N_MET_MAX)) ;
	p_longit_pot[jp] = 
	new Vector ( *source.p_longit_pot[i_met] ) ;
      }

    }
    if (source.p_eta != 0x0) p_eta = new Scalar( *(source.p_eta) ) ; 
    if (source.p_mu != 0x0) p_mu = new Scalar( *(source.p_mu) ) ; 
    if (source.p_www != 0x0) p_www = new Scalar( *(source.p_www) ) ; 
    if (source.p_xxx != 0x0) p_xxx = new Scalar( *(source.p_xxx) ) ; 
    
}


void Sym_tensor::operator=(const Tensor_sym& tt) {
    
    Tensor_sym::operator=(tt) ; 

    del_deriv() ;
}


void Sym_tensor::operator=(const Tensor& tt) {
    
    Tensor_sym::operator=(tt) ; 

    del_deriv() ;
}

			//-------------------//
			// Memory managment  //
			//-------------------//

void Sym_tensor::del_deriv() const {

	for (int i=0; i<N_MET_MAX; i++) 
    	del_derive_met(i) ;
	
	if (p_eta != 0x0) delete p_eta ; 
	if (p_mu != 0x0) delete p_mu ; 
	if (p_www != 0x0) delete p_www ; 
	if (p_xxx != 0x0) delete p_xxx ; 
	if (p_ttt != 0x0) delete p_ttt ; 
	if (p_aaa != 0x0) delete p_aaa ;
	if (p_tilde_b != 0x0) delete p_tilde_b ;
	if (p_tilde_c != 0x0) delete p_tilde_c ;
        
	set_der_0x0() ;
	Tensor::del_deriv() ;

}

void Sym_tensor::set_der_0x0() const {

  for (int i=0; i<N_MET_MAX; i++) 
    set_der_met_0x0(i) ;
  p_eta = 0x0 ; 
  p_mu = 0x0 ; 
  p_www = 0x0 ; 
  p_xxx = 0x0 ; 
  p_ttt = 0x0 ; 
  p_aaa = 0x0 ;
  p_tilde_b = 0x0 ;
  p_tilde_c = 0x0 ;
}


void Sym_tensor::del_derive_met(int j) const {

	assert( (j>=0) && (j<N_MET_MAX) ) ;

	if (met_depend[j] != 0x0) {
    	if ( p_transverse[j] != 0x0) delete p_transverse[j] ;
    	if ( p_longit_pot[j] != 0x0) delete p_longit_pot[j] ;
    
    	set_der_met_0x0(j) ;
    
		Tensor::del_derive_met(j) ;
	}
}


void Sym_tensor::set_der_met_0x0(int i) const {

  assert( (i>=0) && (i<N_MET_MAX) ) ;

  p_transverse[i] = 0x0 ;
  p_longit_pot[i] = 0x0 ;

}


                //----------------------------------//
                //  Computation of derived members  //
                //----------------------------------//
	
const Vector& Sym_tensor::divergence(const Metric& gam) const {
  
    const Vector* pvect = 
        dynamic_cast<const Vector*>( &(Tensor::divergence(gam)) ) ;

    assert(pvect != 0x0) ;

    return *pvect ;
}


Sym_tensor Sym_tensor::derive_lie(const Vector& vv) const {

    Sym_tensor resu(*mp, type_indice, *triad) ; 
    
    compute_derive_lie(vv, resu) ;
    
    return resu ; 
    
}



Sym_tensor* Sym_tensor::inverse() const {

  //Le resultat :
  Sym_tensor* res = 
    new Sym_tensor(*mp, -type_indice(0), *triad) ;
    
  // le determinant :
  Scalar determ1(*mp) ;
  determ1 = double(1)/
    (operator()(1, 1)*operator()(2, 2)*operator()(3, 3) 
     + operator()(1, 2)*operator()(2, 3)*operator()(1, 3)
     + operator()(1, 3)*operator()(1, 2)*operator()(2, 3) 
     - operator()(1, 3)*operator()(2, 2)*operator()(1, 3)
     - operator()(2, 3)*operator()(2, 3)*operator()(1, 1) 
     - operator()(3, 3)*operator()(1, 2)*operator()(1, 2) ) ;
    
  int sgn ;	// Le signe du co-facteur ...
  int l_up, l_down, c_left, c_right ;	    // Coordonnees du cofacteur :
    
  Scalar cofacteur(*mp) ;
    
  for (int i=1 ; i<=3 ; i++) {
    sgn = 1 ;
    for (int j=i ; j<=3 ; j++) {
	    
      switch (j) {
		
      case 1 : {
	c_left = 2 ;
	c_right = 3 ;
	break ;
      }
      case 2 : {
	c_left = 1 ;
	c_right = 3 ;
	break ;
      }
      default : {
	c_left = 1 ;
	c_right = 2 ;
	break ;
      }
      }
	    
      switch (i) {
		
      case 1 : {
	l_up = 2 ;
	l_down = 3 ;
	break ;
      }
      case 2 : {
	l_up = 1 ;
	l_down = 3 ;
	break ;
      }
      default : {
	l_up = 1 ;
	l_down = 2 ;
	break ;
      } 
      }
	    
      cofacteur = sgn*(operator()(l_up, c_left)*operator()(l_down, c_right)-
		       operator()(l_up, c_right)*operator()(l_down, c_left))*determ1 ;
	    
      res->set(i, j) = cofacteur ;
      sgn *= -1 ;
    }
  }
  return res ;

}

void Sym_tensor::exponential_filter_r(int lzmin, int lzmax, int p, 
			  double alpha) {
    if( triad->identify() == (mp->get_bvect_cart()).identify() ) 
	for (int i=0; i<n_comp; i++)
	    cmp[i]->exponential_filter_r(lzmin, lzmax, p, alpha) ;
    else {
	assert( triad->identify() == (mp->get_bvect_spher()).identify()) ;
	Scalar srr_tmp = operator()(1,1) ; 
	srr_tmp.exponential_filter_r(lzmin, lzmax, p, alpha) ;
	Scalar eta_tmp = eta() ;
	eta_tmp.div_r() ; //## to change one day...
	eta_tmp.exponential_filter_r(lzmin, lzmax, p, alpha) ;
	Scalar mu_tmp = mu() ;
	mu_tmp.div_r() ; //## to change one day...
	mu_tmp.exponential_filter_r(lzmin, lzmax, p, alpha) ;
	Scalar w_tmp = www() ;
	w_tmp.exponential_filter_r(lzmin, lzmax, p, alpha) ;
	Scalar x_tmp = xxx() ;
	x_tmp.exponential_filter_r(lzmin, lzmax, p, alpha) ;
	Scalar t_tmp = ttt() ;
	t_tmp.exponential_filter_r(lzmin, lzmax, p, alpha) ;
	set_auxiliary(srr_tmp, eta_tmp, mu_tmp, w_tmp, x_tmp, t_tmp) ;
    }
}

void Sym_tensor::exponential_filter_ylm(int lzmin, int lzmax, int p, 
			  double alpha) {
    if( triad->identify() == (mp->get_bvect_cart()).identify() ) 
	for (int i=0; i<n_comp; i++)
	    cmp[i]->exponential_filter_ylm(lzmin, lzmax, p, alpha) ;
    else {
	assert( triad->identify() == (mp->get_bvect_spher()).identify()) ;
	Scalar srr_tmp = operator()(1,1) ; 
	srr_tmp.exponential_filter_ylm(lzmin, lzmax, p, alpha) ;
	Scalar eta_tmp = eta() ;
	eta_tmp.div_r() ; //## to change one day...
	eta_tmp.exponential_filter_ylm(lzmin, lzmax, p, alpha) ;
	Scalar mu_tmp = mu() ;
	mu_tmp.div_r() ; //## to change one day...
	mu_tmp.exponential_filter_ylm(lzmin, lzmax, p, alpha) ;
	Scalar w_tmp = www() ;
	w_tmp.exponential_filter_ylm(lzmin, lzmax, p, alpha) ;
	Scalar x_tmp = xxx() ;
	x_tmp.exponential_filter_ylm(lzmin, lzmax, p, alpha) ;
	Scalar t_tmp = ttt() ;
	t_tmp.exponential_filter_ylm(lzmin, lzmax, p, alpha) ;
	set_auxiliary(srr_tmp, eta_tmp, mu_tmp, w_tmp, x_tmp, t_tmp) ;
    }
}
}
