/*
 *  Methods of class Sym_tensor_tt related to eta and mu
 *
 *   (see file sym_tensor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
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
 * $Id: sym_tensor_tt_etamu.C,v 1.19 2016/12/05 16:18:17 j_novak Exp $
 * $Log: sym_tensor_tt_etamu.C,v $
 * Revision 1.19  2016/12/05 16:18:17  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.18  2014/10/13 08:53:44  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.17  2014/10/06 15:13:19  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.16  2006/10/24 13:03:19  j_novak
 * New methods for the solution of the tensor wave equation. Perhaps, first
 * operational version...
 *
 * Revision 1.15  2005/04/01 14:28:32  j_novak
 * Members p_eta and p_mu are now defined in class Sym_tensor.
 *
 * Revision 1.14  2004/06/04 09:25:58  e_gourgoulhon
 * Method eta(): eta is no longer computed from h^rr but from khi (in the
 *   case where khi is known).
 *
 * Revision 1.13  2004/05/25 15:08:44  f_limousin
 * Add parameters in argument of the functions update, eta and mu for
 * the case of a Map_et.
 *
 * Revision 1.12  2004/05/24 13:45:29  e_gourgoulhon
 * Added parameter dzp to method Sym_tensor_tt::update.
 *
 * Revision 1.11  2004/05/05 14:24:54  e_gourgoulhon
 * Corrected a bug in method set_khi_mu: the division of khi by r^2
 * was ommitted in the case dzp=2 !!!
 *
 * Revision 1.10  2004/04/08 16:38:43  e_gourgoulhon
 * Sym_tensor_tt::set_khi_mu: added argument dzp (dzpuis of resulting h^{ij}).
 *
 * Revision 1.9  2004/03/04 09:53:04  e_gourgoulhon
 * Methods eta(), mu() and upate(): use of Scalar::mult_r_dzpuis and
 * change of dzpuis behavior of eta and mu.
 *
 * Revision 1.8  2004/03/03 13:16:21  j_novak
 * New potential khi (p_khi) and the functions manipulating it.
 *
 * Revision 1.7  2004/02/05 13:44:50  e_gourgoulhon
 * Major modif. of methods eta(), mu() and update() to treat
 * any value of dzpuis, thanks to the new definitions of
 * Scalar::mult_r(), Scalar::dsdr(), etc...
 *
 * Revision 1.6  2004/01/28 13:25:41  j_novak
 * The ced_mult_r arguments have been suppressed from the Scalar::*dsd* methods.
 * In the div/mult _r_dzpuis, there is no more default value.
 *
 * Revision 1.5  2003/11/05 15:28:31  e_gourgoulhon
 * Corrected error in update.
 *
 * Revision 1.4  2003/11/04 23:03:34  e_gourgoulhon
 * First full version of method update().
 * Add method set_rr_mu.
 * Method set_eta_mu ---> set_rr_eta_mu.
 *
 * Revision 1.3  2003/11/04 09:35:27  e_gourgoulhon
 * First operational version of update_tp().
 *
 * Revision 1.2  2003/11/03 22:33:36  e_gourgoulhon
 * Added methods update_tp and set_eta_mu.
 *
 * Revision 1.1  2003/11/03 17:08:37  e_gourgoulhon
 * First version
 *
 *
 * $Header: /cvsroot/Lorene/C++/Source/Tensor/sym_tensor_tt_etamu.C,v 1.19 2016/12/05 16:18:17 j_novak Exp $
 *
 */

// Headers C
#include <cstdlib>
#include <cassert>

// Headers Lorene
#include "tensor.h"

			//--------------//
			//     khi      //
			//--------------//

namespace Lorene {
const Scalar& Sym_tensor_tt::khi() const {

  if (p_khi == 0x0) {   // a new computation is necessary
		
    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

    // khi is computed from $h^{rr}$ component

    p_khi = new Scalar(operator()(1,1)) ;
    p_khi->mult_r() ;
    p_khi->mult_r() ;
  }

  return *p_khi ; 

}


			//--------------//
			//     eta      //
			//--------------//
			
			
const Scalar& Sym_tensor_tt::eta(Param* par) const {


    if (p_eta == 0x0) {   // a new computation is necessary
	

	    // All this has a meaning only for spherical components:
	    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

        // eta is computed from the divergence-free condition:

        int dzp = operator()(1,1).get_dzpuis() ; 
        int dzp_resu = ((dzp == 0) ? 0 : dzp-1) ;

        Scalar source_eta(*mp) ; 
        
        if (p_khi == 0x0) { //  eta is computed from h^rr
                            //  -------------------------
        
	        source_eta = - operator()(1,1).dsdr() ; 	
        
        // dhrr contains - dh^{rr}/dr in all domains but the CED,                                           
        // in the CED:   - r^2 dh^{rr}/dr        if dzp = 0          (1)
        //               - r^(dzp+1) dh^{rr}/dr  if dzp > 0          (2)
                                                    
		        
	        // Multiplication by r of (-d h^{rr}/dr) (with the same dzpuis as h^{rr})
	        source_eta.mult_r_dzpuis( dzp ) ;                           

            // Substraction of the h^rr part and multiplication by r :
            source_eta -= 3. * operator()(1,1) ;                          

            source_eta.mult_r_dzpuis(dzp_resu) ;
        }
        else {      // eta is computed from khi
                    // ------------------------

            source_eta = - p_khi->dsdr() ;
            int diff_dzp = source_eta.get_dzpuis() - dzp_resu ; 
            assert( diff_dzp >= 0 ) ;
            source_eta.dec_dzpuis(diff_dzp) ;
            
            Scalar tmp(*p_khi) ; 
            tmp.div_r_dzpuis(dzp_resu) ; 
            
            source_eta -= tmp ;  
            
        }

        
	    // Resolution of the angular Poisson equation for eta
	    // --------------------------------------------------
	    if (dynamic_cast<const Map_af*>(mp) != 0x0) {
	        p_eta = new Scalar( source_eta.poisson_angu() ) ; 
	    }
	    else {
	        Scalar resu (*mp) ;
	        resu = 0. ;
	        mp->poisson_angu(source_eta, *par, resu) ;
	        p_eta = new Scalar( resu ) ;  	    
	    }
	
    }

    return *p_eta ; 

}

			

			//-------------------//
			//  set_rr_eta_mu    //
			//-------------------//
			

void Sym_tensor_tt::set_rr_eta_mu(const Scalar& hrr, const Scalar& eta_i, 
		const Scalar& mu_i) {

		// All this has a meaning only for spherical components:
		assert( dynamic_cast<const Base_vect_spher*>(triad) != 0x0 ) ; 
						
		set(1,1) = hrr ; 	// h^{rr}
							// calls del_deriv() and therefore delete previous
							// p_eta and p_mu
		
		p_eta = new Scalar( eta_i ) ; 	// eta

		p_mu = new Scalar( mu_i ) ; 	// mu 
		
		update( hrr.get_dzpuis() ) ; // all h^{ij}, except for h^{rr}
		
}
			
			//---------------//
			//  set_rr_mu    //
			//---------------//
			

void Sym_tensor_tt::set_rr_mu(const Scalar& hrr, const Scalar& mu_i) {

		// All this has a meaning only for spherical components:
		assert( dynamic_cast<const Base_vect_spher*>(triad) != 0x0 ) ; 
						
		set(1,1) = hrr ; 	// h^{rr}
							// calls del_deriv() and therefore delete previous
							// p_eta and p_mu
		
		p_mu = new Scalar( mu_i ) ; 	// mu 
		
		eta() ; // computes eta form the divergence-free condition
		
		update( hrr.get_dzpuis() ) ; // all h^{ij}, except for h^{rr}
		
}
			
			//-------------------//
			//  set_khi_eta_mu    //
			//-------------------//
			

void Sym_tensor_tt::set_khi_eta_mu(const Scalar& khi_i, const Scalar& eta_i, 
		const Scalar& mu_i) {

  // All this has a meaning only for spherical components:
  assert( dynamic_cast<const Base_vect_spher*>(triad) != 0x0 ) ; 
			
  set(1,1) = khi_i ;
  set(1,1).div_r() ;
  set(1,1).div_r() ;     // h^{rr}

  // calls del_deriv() and therefore delete previous
  // p_khi, p_eta and p_mu
		
  p_khi = new Scalar( khi_i ) ;        // khi

  p_eta = new Scalar( eta_i ) ; 	// eta
  
  p_mu = new Scalar( mu_i ) ; 	// mu 
  
  update( khi_i.get_dzpuis() ) ; // all h^{ij}, except for h^{rr}
		
}
			
			//---------------//
			//  set_khi_mu    //
			//---------------//
			

void Sym_tensor_tt::set_khi_mu(const Scalar& khi_i, const Scalar& mu_i, 
                              int dzp, Param* par1, Param* par2, Param* par3) {

  // All this has a meaning only for spherical components:
  assert( dynamic_cast<const Base_vect_spher*>(triad) != 0x0 ) ; 
						
  set(1,1) = khi_i ; 
                        // calls del_deriv() and therefore delete previous
                        // p_eta and p_mu

  assert( khi_i.check_dzpuis(0) ) ; 
  if (dzp == 0) {
    set(1,1).div_r() ;
    set(1,1).div_r() ;	// h^{rr}
  }
  else {
    assert(dzp == 2) ; //## temporary: the other cases are not treated yet
    set(1,1).div_r_dzpuis(1) ;
    set(1,1).div_r_dzpuis(2) ;
  }		
  
  p_khi = new Scalar ( khi_i ) ;  // khi

  p_mu = new Scalar( mu_i ) ; 	// mu 

  if (dynamic_cast<const Map_af*>(mp) != 0x0) {		
      eta() ; // computes eta form the divergence-free condition
      // dzp = 0 ==> eta.dzpuis = 0
      // dzp = 2 ==> eta.dzpuis = 1
      
      update(dzp) ; // all h^{ij}, except for h^{rr}
  }
  else {
      eta(par1) ; // computes eta form the divergence-free condition
      // dzp = 0 ==> eta.dzpuis = 0
      // dzp = 2 ==> eta.dzpuis = 1
      
      update(dzp, par2, par3) ; // all h^{ij}, except for h^{rr}
  }

}

			

			//-------------//
			//   update    //
			//-------------//
			


void Sym_tensor_tt::update(int dzp, Param* par1, Param* par2) {

    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

    assert( (p_eta != 0x0) && (p_mu != 0x0) ) ; 
	        
    Itbl idx(2) ;
    idx.set(0) = 1 ;	// r index
	
	// h^{r theta} : 
	// ------------
	idx.set(1) = 2 ;	// theta index
	*cmp[position(idx)] = p_eta->srdsdt() - p_mu->srstdsdp() ; 
    
    if (dzp == 0) {
        assert( cmp[position(idx)]->check_dzpuis(2) ) ; 
        cmp[position(idx)]->dec_dzpuis(2) ;
    }
    
    assert( cmp[position(idx)]->check_dzpuis(dzp) ) ; 
    
	// h^{r phi} :
	// ------------
	idx.set(1) = 3 ;	// phi index
	*cmp[position(idx)] = p_eta->srstdsdp() + p_mu->srdsdt() ; 

    if (dzp == 0) {
        assert( cmp[position(idx)]->check_dzpuis(2) ) ; 
        cmp[position(idx)]->dec_dzpuis(2) ;
    }
	
    assert( cmp[position(idx)]->check_dzpuis(dzp) ) ; 
	
	// h^{theta phi} and h^{phi phi}
	// -----------------------------
	
	//--------------  Computation of T^theta   --> taut : 
    
    Scalar tautst = operator()(1,2).dsdr() ; 

    // dhrr contains  dh^{rt}/dr in all domains but the CED,                                           
    // in the CED:    r^2 dh^{rt}/dr        if dzp = 0          (1)
    //                r^(dzp+1) dh^{rt}/dr  if dzp > 0          (2)
                                                    
	// Multiplication by r of dh^{rt}/dr (with the same dzpuis than h^{rt})
    tautst.mult_r_dzpuis( operator()(1,2).get_dzpuis() ) ; 	
    
    	        
    // Addition of the remaining parts :	
	tautst += 3 * operator()(1,2) - operator()(1,1).dsdt() ; 
	tautst.mult_sint() ; 
	
	Scalar tmp = operator()(1,1) ;
	tmp.mult_cost() ; 		// h^{rr} cos(th)
	
	tautst -= tmp ; 	// T^th / sin(th)
	
	Scalar taut = tautst ; 
	taut.mult_sint() ; 	// T^th
	

	//----------- Computation of T^phi   --> taup : 
    
	Scalar taupst = - operator()(1,3).dsdr() ; 

    // dhrr contains  - dh^{rp}/dr in all domains but the CED,                                           
    // in the CED:    - r^2 dh^{rp}/dr        if dzp = 0          (3)
    //                - r^(dzp+1) dh^{rp}/dr  if dzp > 0          (4)
                                                    	        
	// Multiplication by r of -dh^{rp}/dr  (with the same dzpuis than h^{rp})
    taupst.mult_r_dzpuis( operator()(1,3).get_dzpuis() ) ; 	

                          
    // Addition of the remaining part :	
	
	taupst -= 3 * operator()(1,3) ; 
	taupst.mult_sint() ; 	// T^ph / sin(th)
	
	Scalar taup = taupst ; 
	taup.mult_sint() ; 		// T^ph 
	
    
    //------------------- Computation of F and h^[ph ph}
    
	tmp = tautst ; 
	tmp.mult_cost() ; 	// T^th / tan(th)
	
	// dT^th/dth + T^th / tan(th) + 1/sin(th) dT^ph/dph :
	tmp = taut.dsdt() + tmp + taup.stdsdp() ;

	Scalar tmp2 (*mp) ;		
	if (dynamic_cast<const Map_af*>(mp) != 0x0) {		
	    tmp2 = tmp.poisson_angu() ;  // F
	}
	else {
	    tmp2 = 0. ;
	    mp->poisson_angu(tmp, *par1, tmp2) ; // F
	}
	    

	tmp2.div_sint() ; 
	tmp2.div_sint() ; // h^{ph ph}
	
	idx.set(0) = 3 ;	// phi index
	idx.set(1) = 3 ;	// phi index
	*cmp[position(idx)] = tmp2 ; 		// h^{ph ph} is updated
	
    
    //------------------- Computation of G and h^[th ph}
    
	tmp = taupst ; 
	tmp.mult_cost() ; // T^ph / tan(th)
	
	// - 1/sin(th) dT^th/dph + dT^ph/dth + T^ph / tan(th) :
	tmp = - taut.stdsdp() + taup.dsdt() + tmp ; 
	
	if (dynamic_cast<const Map_af*>(mp) != 0x0) {		
	    tmp2 = tmp.poisson_angu() ;  // G
	}
	else {
	    tmp2 = 0. ;
	    mp->poisson_angu(tmp, *par2, tmp2) ;  // G
	}
	
	tmp2.div_sint() ; 
	tmp2.div_sint() ; // h^{th ph}
	
	idx.set(0) = 2 ;	// theta index
	idx.set(1) = 3 ;	// phi index
	*cmp[position(idx)] = tmp2 ; 		// h^{th ph} is updated
	
	// h^{th th}  (from the trace-free condition)
	// ---------
	idx.set(1) = 2 ;	// theta index
	*cmp[position(idx)] = - operator()(1,1) - operator()(3,3) ; 
	

	Sym_tensor_trans::del_deriv() ; //## in order not to delete p_eta and p_mu
	


}			


			//-----------------//
			//  set_A_tilde_B  //
			//-----------------//

void Sym_tensor_tt::set_A_tildeB(const Scalar& a_in, const Scalar& tb_in, 
				 Param* par_bc, Param* par_mat) {

    Scalar zero(*mp) ;
    zero.set_etat_zero() ;
    set_AtB_trace(a_in, tb_in, zero, par_bc, par_mat) ;
    return ;
}






}
