/*
 *  Methods of class Connection_fcart.
 *
 *	(see file connection.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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
 * $Id: connection_fcart.C,v 1.15 2016/12/05 16:17:50 j_novak Exp $
 * $Log: connection_fcart.C,v $
 * Revision 1.15  2016/12/05 16:17:50  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.14  2014/10/13 08:52:50  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.13  2014/10/06 15:13:04  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.12  2004/01/28 13:25:40  j_novak
 * The ced_mult_r arguments have been suppressed from the Scalar::*dsd* methods.
 * In the div/mult _r_dzpuis, there is no more default value.
 *
 * Revision 1.11  2004/01/04 21:00:50  e_gourgoulhon
 * Better handling of tensor symmetries in methods p_derive_cov() and
 * p_divergence() (thanks to the new class Tensor_sym).
 *
 * Revision 1.10  2004/01/01 11:24:04  e_gourgoulhon
 * Full reorganization of method p_derive_cov: the main loop is now
 * on the indices of the *output* tensor (to take into account
 * symmetries in the input and output tensors).
 *
 * Revision 1.9  2003/12/27 14:59:52  e_gourgoulhon
 * -- Method derive_cov() suppressed.
 * -- Change of the position of the derivation index from the first one
 *    to the last one in methods p_derive_cov() and p_divergence().
 *
 * Revision 1.8  2003/10/17 13:46:15  j_novak
 * The argument is now between 1 and 3 (instead of 0->2)
 *
 * Revision 1.7  2003/10/16 21:37:08  e_gourgoulhon
 * Corrected deriv index in divergence.
 *
 * Revision 1.6  2003/10/16 15:26:03  e_gourgoulhon
 * Suppressed unsued variable
 *
 * Revision 1.5  2003/10/16 14:21:36  j_novak
 * The calculation of the divergence of a Tensor is now possible.
 *
 * Revision 1.4  2003/10/11 16:45:43  e_gourgoulhon
 * Suppressed the call to Itbl::set_etat_qcq() after
 * the construction of the Itbl's.
 *
 * Revision 1.3  2003/10/11 14:39:50  e_gourgoulhon
 * Suppressed declaration of unusued arguments in some methods.
 *
 * Revision 1.2  2003/10/06 13:58:46  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.1  2003/10/03 14:11:48  e_gourgoulhon
 * Methods of class Connection_fcart.
 *
 * 
 *
 * $Header: /cvsroot/Lorene/C++/Source/Connection/connection_fcart.C,v 1.15 2016/12/05 16:17:50 j_novak Exp $
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <cstdlib>

// Lorene headers
#include "connection.h"


        //------------------------------//
        //          Constructors        //
        //------------------------------//



// Contructor from a Cartesian flat-metric-orthonormal basis

namespace Lorene {
Connection_fcart::Connection_fcart(const Map& mpi, const Base_vect_cart& bi) 
  : Connection_flat(mpi, bi) {

}		

// Copy constructor
Connection_fcart::Connection_fcart(const Connection_fcart& ci) 
  : Connection_flat(ci) {

}		

	
        //----------------------------//
        //          Destructor        //
        //----------------------------//


Connection_fcart::~Connection_fcart(){
	
}


        //--------------------------------//
        //      Mutators / assignment     //
        //--------------------------------//


void Connection_fcart::operator=(const Connection_fcart& ) {
	
  cout << "Connection_fcart::operator= : not implemented yet !" << endl ; 
  abort() ; 

}	



        //-----------------------------//
        //    Computational methods    //
        //-----------------------------//

// Covariant derivative, returning a pointer.
//-------------------------------------------

Tensor* Connection_fcart::p_derive_cov(const Tensor& uu) const {

    // Notations: suffix 0 in name <=> input tensor
    //            suffix 1 in name <=> output tensor

    int valence0 = uu.get_valence() ; 
    int valence1 = valence0 + 1 ; 
    int valence1m1 = valence1 - 1 ; // same as valence0, but introduced for 
                                    // the sake of clarity
	
    // Protections
    // -----------
    if (valence0 >= 1) {
        assert(uu.get_triad() == triad) ; 
    }

    // Creation of the result (pointer)
    // --------------------------------
    Tensor* resu ;

    // If uu is a Scalar, the result is a vector
    if (valence0 == 0) 
        resu = new Vector(*mp, COV, triad) ;
    else {

        // Type of indices of the result :
        Itbl tipe(valence1) ; 
        const Itbl& tipeuu = uu.get_index_type() ;  
        for (int id = 0; id<valence0; id++) {
            tipe.set(id) = tipeuu(id) ;   // First indices = same as uu
        }
        tipe.set(valence1m1) = COV ;  // last index is the derivation index

        // if uu is a Tensor_sym, the result is also a Tensor_sym:
        const Tensor* puu = &uu ; 
        const Tensor_sym* puus = dynamic_cast<const Tensor_sym*>(puu) ; 
        if ( puus != 0x0 ) {    // the input tensor is symmetric
            resu = new Tensor_sym(*mp, valence1, tipe, *triad,
                                  puus->sym_index1(), puus->sym_index2()) ;
        }
        else {  
            resu = new Tensor(*mp, valence1, tipe, *triad) ;  // no symmetry  
        }

    }

    int ncomp1 = resu->get_n_comp() ; 
	
    Itbl ind1(valence1) ; // working Itbl to store the indices of resu
    Itbl ind0(valence0) ; // working Itbl to store the indices of uu
	
    // Loop on all the components of the output tensor
    // -----------------------------------------------
    for (int ic=0; ic<ncomp1; ic++) {
    
        // indices corresponding to the component no. ic in the output tensor
        ind1 = resu->indices(ic) ; 
    
        // Component no. ic:
        Scalar& cresu = resu->set(ind1) ; 
		
        // Indices of the input tensor
        for (int id = 0; id < valence0; id++) {
            ind0.set(id) = ind1(id) ; 
        }
 
        // Value of last index (derivation index)
        int k = ind1(valence1m1) ; 
        
        // Partial derivation with respect to x^k:

        cresu = (uu(ind0)).deriv(k) ; 	

    }

    // C'est fini !
    // -----------
    return resu ; 

}



// Divergence, returning a pointer.
//---------------------------------

Tensor* Connection_fcart::p_divergence(const Tensor& uu) const {

    // Notations: suffix 0 in name <=> input tensor
    //            suffix 1 in name <=> output tensor

  int valence0 = uu.get_valence() ; 
  int valence1 = valence0 - 1 ; 
	
  // Protections
  // -----------
  assert (valence0 >= 1) ;
  assert (uu.get_triad() == triad) ; 
  
  // Last index must be contravariant:
  assert (uu.get_index_type(valence0-1) == CON) ;


  // Creation of the pointer on the result tensor
  // --------------------------------------------
    Tensor* resu ;

    if (valence0 == 1)      // if u is a Vector, the result is a Scalar
        resu = new Scalar(*mp) ;
    else {
    
        // Type of indices of the result :
        Itbl tipe(valence1) ; 
        const Itbl& tipeuu = uu.get_index_type() ;  
        for (int id = 0; id<valence1; id++) {
            tipe.set(id) = tipeuu(id) ;     // type of remaining indices = 
        }                                   //  same as uu indices

        if (valence0 == 2) {  // if u is a rank 2 tensor, the result is a Vector
            resu = new Vector(*mp, tipe(0), *triad) ;
        }
        else {
            // if uu is a Tensor_sym, the result might be also a Tensor_sym:
            const Tensor* puu = &uu ; 
            const Tensor_sym* puus = dynamic_cast<const Tensor_sym*>(puu) ; 
            if ( puus != 0x0 ) {    // the input tensor is symmetric

                if (puus->sym_index2() != valence0 - 1) {
                 
                    // the symmetry is preserved: 

                    if (valence1 == 2) {
                        resu = new Sym_tensor(*mp, tipe, *triad) ;
                    }
                    else {
                        resu = new Tensor_sym(*mp, valence1, tipe, *triad,
                                  puus->sym_index1(), puus->sym_index2()) ;
                    }
                }
                else { // the symmetry is lost: 
                
	            resu = new Tensor(*mp, valence1, tipe, *triad) ;
                }
            }
            else { // no symmetry in the input tensor:
	        resu = new Tensor(*mp, valence1, tipe, *triad) ;
            }
        }
    }


  int ncomp1 = resu->get_n_comp() ;
	
  Itbl ind0(valence0) ; // working Itbl to store the indices of uu
	
  Itbl ind1(valence1) ; // working Itbl to store the indices of resu
	
  // Loop on all the components of the output tensor
  for (int ic=0; ic<ncomp1; ic++) {
	
    ind1 = resu->indices(ic) ; 
    Scalar& cresu = resu->set(ind1) ;
    cresu.set_etat_zero() ;

    for (int k=1; k<=3; k++) {
      
      // indices (ind1,k) in the input tensor
      for (int id = 0; id < valence1; id++) {
	ind0.set(id) = ind1(id) ; 
      }
      ind0.set(valence0-1) = k ; 

      cresu += uu(ind0).deriv(k) ; //Addition of dT^i/dx^i
    }

  }

  // C'est fini !
  // -----------
  return resu ; 

}

}
