/*
 *  Methods of the class Cmp for partial differential equations
 *   with a multipole falloff condition at the outer boundary
 *
 *    (see file cmp.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Joshua A. Faber
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
 * $Id: cmp_pde_ylm.C,v 1.3 2016/12/05 16:17:49 j_novak Exp $
 * $Log: cmp_pde_ylm.C,v $
 * Revision 1.3  2016/12/05 16:17:49  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.2  2014/10/13 08:52:48  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1  2004/12/29 16:27:48  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header: /cvsroot/Lorene/C++/Source/Cmp/cmp_pde_ylm.C,v 1.3 2016/12/05 16:17:49 j_novak Exp $
 *
 */

// Lorene headers
#include "map.h"
#include "cmp.h"
#include "param.h"

		    //-----------------------------------//
		    //      Scalar Poisson equation	 //
		    //-----------------------------------//

// Version without parameters
// --------------------------

namespace Lorene {
Cmp Cmp::poisson_ylm(int nylm, double* intvec) const {
    
    Param bidon ;
    Cmp resu(*mp) ; 
    
    mp->poisson_ylm(*this, bidon, resu, nylm, intvec) ; 

    return resu ;          
}

// Version with parameters
// -----------------------

void Cmp::poisson_ylm(Param& par, Cmp& uu, int nylm, double* intvec) const {
    
    mp->poisson_ylm(*this, par, uu, nylm, intvec) ;     
    
}
}
