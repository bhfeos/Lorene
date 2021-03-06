/*
 *  Main code for reading a time slice Sigma_t stored in file. 
 *
 */

/*
 *   Copyright (c) 2004  Eric Gourgoulhon & Jerome Novak
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
 * $Id: read_tslice.C,v 1.7 2016/12/05 16:18:24 j_novak Exp $
 * $Log: read_tslice.C,v $
 * Revision 1.7  2016/12/05 16:18:24  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:53:56  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:09:44  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2004/06/24 15:42:16  j_novak
 * Display of log(n) and log(psi) instead of N/psi.
 *
 * Revision 1.3  2004/06/15 10:25:27  e_gourgoulhon
 * Added plot of hrt and hrp.
 *
 * Revision 1.2  2004/05/31 20:32:50  e_gourgoulhon
 * Added graphical outputs.
 *
 * Revision 1.1  2004/05/27 15:25:38  e_gourgoulhon
 * First version.
 *
 *
 * $Header: /cvsroot/Lorene/Codes/Einstein/read_tslice.C,v 1.7 2016/12/05 16:18:24 j_novak Exp $
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <cstdlib>
#include <cmath>
#include <cstring>

// Lorene headers
#include "time_slice.h"
#include "param.h"
#include "utilitaires.h"
#include "graphique.h"

using namespace Lorene ;

int main(int argc, char** argv){
    
    if (argc < 2) {
		cout << 
		"read_tslice : the name of a file containing a configuration"
		<< endl << " must be given in argument !" << endl ; 
		abort() ; 
    }
    
    char* nomresu = argv[1] ; 
    cout << "Name of the file to be read : " << nomresu << endl ;         
    
    FILE* fich = fopen(nomresu, "r") ; 
    if (fich == 0x0) {
    	cout << "Problem in opening the file " << nomresu << " ! " << endl ; 
		perror(" reason") ; 
		abort() ; 
    }
    
    Mg3d mgrid(fich) ;
    Map_af map(mgrid, fich) ;
    
    Base_vect* ptriad_file = Base_vect::bvect_from_file(fich) ; 
    
    cout << "Computational grid :\n" 
         << "------------------ \n" 
         << "  " << mgrid << endl ; 

    cout << "Mapping computational grid --> physical space :\n" 
         << "---------------------------------------------\n" 
         << "  " << map << endl ;  
    
    const Base_vect_spher* triad_s = 
              dynamic_cast<const Base_vect_spher*>(ptriad_file) ; 
    const Base_vect_cart* triad_c = 
              dynamic_cast<const Base_vect_cart*>(ptriad_file) ; 
    
    const Base_vect* ptriad ; 
    
    if ( triad_s != 0x0 ) ptriad = &(map.get_bvect_spher()) ;
    else {
        assert( triad_c != 0x0) ; 
        ptriad = &(map.get_bvect_cart()) ;
    }
     
    delete ptriad_file ; 

   
    // Flat metric f
    // -------------

    const Metric_flat& ff = map.flat_met_spher() ; 
    
    int depth ; 
	fread_be(&depth, sizeof(int), 1, fich) ;	
    
    // Time slice Sigma_t
    // -------------------
    Tslice_dirac_max sigmat(map, *ptriad, ff, fich, false, depth) ;  
    
    fclose(fich) ;     

    Scalar lnn = log(sigmat.nn()) ;
    lnn.std_spectral_base() ;
    
    cout << sigmat << endl ; 
    
    // For graphical outputs:
    char graph_device[40] ; 
    strcpy(graph_device, "/xwin") ;
    int ngraph0 = 20 ;  // index of the first graphic device to be used
    int nz = mgrid.get_nzone() ; 
    double ray_des = 1.25 * map.val_r(nz-2, 1., 0., 0.) ; // outermost radius
                                                          // for plots
    des_meridian(lnn, 0., ray_des, "log(N)", ngraph0,
                     graph_device) ; 
    des_meridian(sigmat.ln_psi(), 0., ray_des, "log(\\gQ)", ngraph0+1,
                     graph_device) ; 
    des_meridian(sigmat.beta()(1), 0., ray_des, "\\gb\\ur\\d", ngraph0+6,
                     graph_device) ; 
    des_meridian(sigmat.beta()(2), 0., ray_des, "\\gb\\u\\gh\\d", ngraph0+7,
                     graph_device) ; 
    des_meridian(sigmat.beta()(3), 0., ray_des, "\\gb\\u\\gf\\d", ngraph0+8,
                     graph_device) ; 
    des_meridian(sigmat.khi(), 0., ray_des, "\\gx", ngraph0+9,
                     graph_device) ; 
    des_meridian(sigmat.mu(), 0., ray_des, "\\gm", ngraph0+10,
                     graph_device) ;
    des_meridian(sigmat.trh(), 0., ray_des, "tr h", ngraph0+11,
                     graph_device) ; 
    des_meridian(sigmat.hh()(1,1), 0., ray_des, "h\\urr\\d", ngraph0+12,
                     graph_device) ; 
    des_meridian(sigmat.hh()(1,2), 0., ray_des, "h\\ur\\gh\\d", ngraph0+13,
                     graph_device) ; 
    des_meridian(sigmat.hh()(1,3), 0., ray_des, "h\\ur\\gf\\d", ngraph0+14,
                     graph_device) ; 
    des_meridian(sigmat.hh()(2,3), 0., ray_des, "h\\u\\gh\\gf\\d", ngraph0+15,
                     graph_device) ; 
    des_meridian(sigmat.hh()(3,3), 0., ray_des, "h\\u\\gf\\gf\\d", ngraph0+16,
                     graph_device) ; 
    des_meridian(sigmat.aa()(1,1), 0., ray_des, "A\\urr\\d", ngraph0+17, 
                         graph_device) ; 
    des_meridian(sigmat.aa()(2,3), 0., ray_des, "A\\u\\gh\\gf\\d", ngraph0+18,
                         graph_device) ; 
    des_meridian(sigmat.aa()(3,3), 0., ray_des, "A\\u\\gf\\gf\\d", ngraph0+19,
                          graph_device) ; 
    arrete() ; 

    return EXIT_SUCCESS ; 
}
