/*
 *  Method of class Star_bhns compute the derivative of metric quantities
 *   from the companion black-hole
 *
 *    (see file star_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2006-2007 Keisuke Taniguchi
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
 * $Id: star_bhns_upmetr_der.C,v 1.4 2016/12/05 16:18:16 j_novak Exp $
 * $Log: star_bhns_upmetr_der.C,v $
 * Revision 1.4  2016/12/05 16:18:16  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:41  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2008/05/15 19:18:17  k_taniguchi
 * Change of some parameters.
 *
 * Revision 1.1  2007/06/22 01:32:55  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header: /cvsroot/Lorene/C++/Source/Star_bhns/star_bhns_upmetr_der.C,v 1.4 2016/12/05 16:18:16 j_novak Exp $
 *
 */

// C++ headers
//#include <>

// C headers
//#include <>

// Lorene headers
#include "star_bhns.h"
#include "hole_bhns.h"
#include "utilitaires.h"

namespace Lorene {
void Star_bhns::update_met_der_comp_bhns(const Hole_bhns& hole) {

    // Computation of d_lapconf_comp
    // -----------------------------

    if ( (hole.get_d_lapconf_auto())(1).get_etat() == ETATZERO ) {
        assert( (hole.get_d_lapconf_auto())(2).get_etat() == ETATZERO ) ;
	assert( (hole.get_d_lapconf_auto())(3).get_etat() == ETATZERO ) ;

	d_lapconf_comp.set_etat_zero() ;
    }
    else {
        d_lapconf_comp.set_etat_qcq() ;
	Vector comp_dlapconf( hole.get_d_lapconf_auto() ) ;
	comp_dlapconf.dec_dzpuis(2) ; // dzpuis : 2 -> 0 for import

	(d_lapconf_comp.set(1)).import( comp_dlapconf(1) ) ;
	(d_lapconf_comp.set(2)).import( comp_dlapconf(2) ) ;
	(d_lapconf_comp.set(3)).import( comp_dlapconf(3) ) ;

	d_lapconf_comp.std_spectral_base() ;
	d_lapconf_comp.inc_dzpuis(2) ;  // dzpuis : 0 -> 2
    }


    // Computation of d_shift_comp
    // ---------------------------

    if ( (hole.get_d_shift_auto())(1,2).get_etat() == ETATZERO ) {
        assert( (hole.get_d_shift_auto())(1,1).get_etat() == ETATZERO ) ;
	assert( (hole.get_d_shift_auto())(1,3).get_etat() == ETATZERO ) ;

        d_shift_comp.set_etat_zero() ;
    }
    else {

        d_shift_comp.set_etat_qcq() ;
	Tensor comp_dshift( hole.get_d_shift_auto() ) ;
	comp_dshift.dec_dzpuis(2) ;  // dzpuis : 2 -> 0 for import

	(d_shift_comp.set(1,1)).import( comp_dshift(1,1) ) ;
	(d_shift_comp.set(1,2)).import( comp_dshift(1,2) ) ;
	(d_shift_comp.set(1,3)).import( comp_dshift(1,3) ) ;
	(d_shift_comp.set(2,1)).import( comp_dshift(2,1) ) ;
	(d_shift_comp.set(2,2)).import( comp_dshift(2,2) ) ;
	(d_shift_comp.set(2,3)).import( comp_dshift(2,3) ) ;
	(d_shift_comp.set(3,1)).import( comp_dshift(3,1) ) ;
	(d_shift_comp.set(3,2)).import( comp_dshift(3,2) ) ;
	(d_shift_comp.set(3,3)).import( comp_dshift(3,3) ) ;

	d_shift_comp.std_spectral_base() ;
	d_shift_comp.inc_dzpuis(2) ;  // dzpuis : 0 -> 2
    }


    // Computation of d_confo_comp
    // ---------------------------

    if ( (hole.get_d_confo_auto())(1).get_etat() == ETATZERO ) {
        assert( (hole.get_d_confo_auto())(2).get_etat() == ETATZERO ) ;
	assert( (hole.get_d_confo_auto())(3).get_etat() == ETATZERO ) ;

        d_confo_comp.set_etat_zero() ;
    }
    else {
        d_confo_comp.set_etat_qcq() ;
	Vector comp_dconfo( hole.get_d_confo_auto() ) ;
	comp_dconfo.dec_dzpuis(2) ;  // dzpuis : 2 -> 0 for import

	(d_confo_comp.set(1)).import( comp_dconfo(1) ) ;
	(d_confo_comp.set(2)).import( comp_dconfo(2) ) ;
	(d_confo_comp.set(3)).import( comp_dconfo(3) ) ;

	d_confo_comp.std_spectral_base() ;
	d_confo_comp.inc_dzpuis(2) ;  // dzpuis : 0 -> 2
    }


    // The derived quantities are obsolete
    // -----------------------------------

    del_deriv() ;

}
}
