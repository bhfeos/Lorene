/*
 *  Prototypes of non class-member functions
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 1999-2002 Eric Gourgoulhon
 *   Copyright (c) 2002-2003 Jerome Novak
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


#ifndef	__PROTO_H_
#define	__PROTO_H_


/*
 * $Id: proto.h,v 1.51 2014/10/13 08:52:36 j_novak Exp $
 * $Log: proto.h,v $
 * Revision 1.51  2014/10/13 08:52:36  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.50  2013/06/13 14:18:18  j_novak
 * Inclusion of new bases R_LEG, R_LEGP and R_LEGI.
 *
 * Revision 1.49  2013/06/06 15:31:31  j_novak
 * Functions to compute Legendre coefficients (not fully tested yet).
 *
 * Revision 1.48  2013/06/05 15:06:10  j_novak
 * Legendre bases are treated as standard bases, when the multi-grid
 * (Mg3d) is built with BASE_LEG.
 *
 * Revision 1.47  2010/10/22 08:08:40  j_novak
 * Removal of the method Star_rot_dirac::lambda_grv2() and call to the C++ version of integrale2d.
 *
 * Revision 1.46  2010/01/20 14:53:50  n_vasset
 * Adding spectral cutoff functions for use in elliptic tensor equations.
 *
 * Revision 1.45  2009/10/23 12:55:46  j_novak
 * New base T_LEG_MI
 *
 * Revision 1.44  2009/10/13 13:50:39  j_novak
 * New base T_LEG_MP.
 *
 * Revision 1.43  2009/08/31 19:39:07  n_vasset
 * removal of obsolete function get_kerr()
 *
 * Revision 1.42  2008/11/27 12:12:38  j_novak
 * New function to initialize parameters for wave equation.
 *
 * Revision 1.41  2008/08/20 11:51:25  n_vasset
 * new functions to solve the Kerr problem, using degenerate elliptic operators
 *
 * Revision 1.40  2008/08/19 06:41:59  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.39  2008/07/18 12:28:41  j_novak
 * Corrected some mistakes.
 *
 * Revision 1.38  2008/07/18 09:17:35  j_novak
 * New function tilde_laplacian().
 *
 * Revision 1.37  2008/07/11 13:20:08  j_novak
 * Miscellaneous functions for the system wave equation.
 *
 * Revision 1.36  2008/07/10 10:34:35  p_grandclement
 * forgot this one
 *
 * Revision 1.35  2007/12/11 15:28:05  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.34  2007/05/06 10:48:08  p_grandclement
 * Modification of a few operators for the vorton project
 *
 * Revision 1.33  2007/01/23 17:08:43  j_novak
 * New function pois_vect_r0.C to solve the l=0 part of the vector Poisson
 * equation, which involves only the r-component.
 *
 * Revision 1.32  2006/04/27 09:12:29  p_grandclement
 * First try at irrotational black holes
 *
 * Revision 1.31  2006/04/10 15:19:18  j_novak
 * New definition of 1D operators dsdx and sx in the nucleus (bases R_CHEBP and
 * R_CHEBI).
 *
 * Revision 1.30  2005/11/30 11:09:03  p_grandclement
 * Changes for the Bin_ns_bh project
 *
 * Revision 1.29  2005/08/26 14:02:38  p_grandclement
 * Modification of the elliptic solver that matches with an oscillatory exterior solution
 * small correction in Poisson tau also...
 *
 * Revision 1.28  2005/06/09 07:56:25  f_limousin
 * Implement a new function sol_elliptic_boundary() and
 * Vector::poisson_boundary(...) which solve the vectorial poisson
 * equation (method 6) with an inner boundary condition.
 *
 * Revision 1.27  2005/05/13 08:50:29  j_novak
 * Added the function int1d_chebi.
 *
 * Revision 1.26  2005/02/16 15:04:07  m_forot
 * Add int1D_cheb function
 *
 * Revision 1.25  2005/02/08 10:08:57  f_limousin
 * Add neumann_binaire(...), dirichlet_binaire(...) and
 * poisson_vect_binaire(...) with Scalars and Vectors in argument.
 *
 * Revision 1.24  2004/12/17 13:35:00  m_forot
 * Add the case T_LEG
 *
 * Revision 1.23  2004/11/23 15:05:40  m_forot
 * Added prototypes for all new functions in the case there is no
 * symmetry in theta.
 *
 * Revision 1.22  2004/09/28 15:59:47  f_limousin
 * Add function poisson_vect_boundary which is the same as
 * poisson_vect_frontiere but for the new classes Tensor and Scalar.
 *
 * Revision 1.21  2004/08/24 09:14:40  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.20  2004/06/22 08:49:57  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * Revision 1.19  2004/03/17 15:58:47  p_grandclement
 * Slight modification of sol_elliptic_no_zec
 *
 * Revision 1.18  2004/02/17 09:21:38  j_novak
 * New functions for calculating values of the derivatives of a function
 * using its Chebyshev coefficients.
 *
 * Revision 1.17  2004/02/11 09:47:44  p_grandclement
 * Addition of a new elliptic solver, matching with the homogeneous solution
 * at the outer shell and not solving in the external domain (more details
 * coming soon ; check your local Lorene dealer...)
 *
 * Revision 1.16  2004/02/09 08:55:30  j_novak
 * Corrected error in the arguments of _solp_r_chebu_cinq
 *
 * Revision 1.15  2004/02/06 10:53:51  j_novak
 * New dzpuis = 5 -> dzpuis = 3 case (not ready yet).
 *
 * Revision 1.14  2004/01/28 16:46:22  p_grandclement
 * Addition of the sol_elliptic_fixe_der_zero stuff
 *
 * Revision 1.13  2004/01/15 09:15:36  p_grandclement
 * Modification and addition of the Helmholtz operators
 *
 * Revision 1.12  2003/12/11 14:48:47  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * Revision 1.11  2003/09/16 13:07:40  j_novak
 * New files for coefficient trnasformation to/from the T_LEG_II base.
 *
 * Revision 1.10  2003/09/16 08:53:05  j_novak
 * Addition of the T_LEG_II base (odd in theta, only for odd m) and the
 * transformation functions from and to T_SIN_P.
 *
 * Revision 1.9  2003/06/18 08:45:26  j_novak
 * In class Mg3d: added the member get_radial, returning only a radial grid
 * For dAlembert solver: the way the coefficients of the operator are defined has been changed.
 *
 * Revision 1.8  2003/02/13 16:40:24  p_grandclement
 * Addition of various things for the Bin_ns_bh project, non of them being
 * completely tested
 *
 * Revision 1.7  2002/11/12 17:45:19  j_novak
 * Added transformation function for T_COS basis.
 *
 * Revision 1.6  2002/09/09 13:00:39  e_gourgoulhon
 * Modification of declaration of Fortran 77 prototypes for
 * a better portability (in particular on IBM AIX systems):
 * All Fortran subroutine names are now written F77_* and are
 * defined in the new file C++/Include/proto_f77.h.
 *
 * Revision 1.5  2002/05/11 12:39:08  e_gourgoulhon
 * Added declaration of som_tet_cossin_si.
 *
 * Revision 1.4  2002/05/05 16:24:48  e_gourgoulhon
 * Added som_tet_cossin_sp
 *
 * Revision 1.3  2002/01/03 15:30:27  j_novak
 * Some comments modified.
 *
 * Revision 1.2  2002/01/02 14:07:56  j_novak
 * Dalembert equation is now solved in the shells. However, the number of
 * points in theta and phi must be the same in each domain. The solver is not
 * completely tested (beta version!).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.69  2001/05/07  09:11:26  phil
 * *** empty log message ***
 *
 * Revision 2.68  2001/04/03  12:41:23  phil
 * modification de itemax dans separation
 *
 * Revision 2.67  2001/03/22  10:40:13  phil
 * modification prototypage se separation
 *
 * Revision 2.66  2001/03/02  10:21:26  phil
 * *** empty log message ***
 *
 * Revision 2.65  2001/03/02  10:18:47  phil
 * modif parametrage separation
 *
 * Revision 2.64  2001/02/28  11:23:30  phil
 * ajout separation
 *
 * Revision 2.63  2001/01/29  14:30:10  phil
 * ajout type rotation
 *
 * Revision 2.62  2000/12/13  15:42:14  phil
 * ajout des trucs relatifs a Lindquist
 *
 * Revision 2.61  2000/12/04  13:33:28  novak
 * *** empty log message ***
 *
 * Revision 2.60  2000/12/04 13:29:08  novak
 * Added prototypes for the dalembertian
 *
 * Revision 2.59  2000/10/19 10:07:49  phil
 * ajout de regle
 *
 * Revision 2.58  2000/10/19  09:35:44  phil
 * *** empty log message ***
 *
 * Revision 2.57  2000/10/04  14:40:34  eric
 * *** empty log message ***
 *
 * Revision 2.56  2000/09/29  14:02:06  eric
 * *** empty log message ***
 *
 * Revision 2.55  2000/09/28  10:01:36  eric
 * *** empty log message ***
 *
 * Revision 2.54  2000/09/08  16:26:26  eric
 * *** empty log message ***
 *
 * Revision 2.53  2000/09/08  16:07:43  eric
 * *** empty log message ***
 *
 * Revision 2.52  2000/09/07  15:07:40  eric
 * *** empty log message ***
 *
 * Revision 2.51  2000/09/06  13:59:53  eric
 * *** empty log message ***
 *
 * Revision 2.50  2000/06/06  12:42:55  phil
 * ajout de Cmp division_xpun (const Cmp&,int)
 *
 * Revision 2.49  2000/05/22  13:33:15  phil
 * ajout des trucs pour poisson avec dzpuis == 3
 *
 * Revision 2.48  2000/04/03  17:01:01  phil
 * ajout de sxpun_1d
 *
 * Revision 2.47  2000/03/16  16:28:30  phil
 * *** empty log message ***
 *
 * Revision 2.46  2000/03/09  13:52:55  phil
 * *** empty log message ***
 *
 * Revision 2.45  2000/03/09  13:42:34  phil
 * vire les trucs relatifs a Poisson compacts
 *
 * Revision 2.44  2000/03/06  10:27:07  eric
 * Ajout des protos som_*_symy et som_*_asymy.
 *
 * Revision 2.43  2000/01/20  14:07:59  phil
 * vire poisson_vect et xksk
 *
 * Revision 2.42  1999/12/15  09:41:52  eric
 * *** empty log message ***
 *
 *
 * $Header: /cvsroot/Lorene/C++/Include/proto.h,v 1.51 2014/10/13 08:52:36 j_novak Exp $
 *
 */

namespace Lorene {
class Param ;
class Tbl ;
class Mtbl ;
class Mtbl_cf ;
class Map_af ;
class Matrice ;
class Valeur ;
class Base_val ;
class Cmp ;
class Tenseur ;
class Scalar ;
class Tensor ;
class Sym_tensor ;
class Vector ;
class Param_elliptic ;


// Routines calcul de coefficients
// -------------------------------
double* cheb_ini(const int) ;
double* chebimp_ini(const int) ;
void four1d(const int, double* ) ;
void chebyf1d(const int, double* ) ; 
void chebyr1d(const int, double* ) ;
void cfpcossin(const int* ,const int* ,  double* ) ;
void cfpcossini(const int* ,const int* ,  double* ) ;
void cftcos(const int*, const int*, double*, const int*, double*) ;
void cftsin(const int*, const int*, double*, const int*, double*) ;
void cftcosp(const int*, const int*, double*, const int*, double*) ;
void cftcosi(const int*, const int*, double*, const int*, double*) ;
void cftsinp(const int*, const int*, double*, const int*, double*) ;
void cftsini(const int*, const int*, double*, const int*, double*) ;
void cftcossincp(const int*, const int*, double*, const int*, double*) ;
void cftcossinsi(const int*, const int*, double*, const int*, double*) ;
void cftcossinsp(const int*, const int*, double*, const int*, double*) ;
void cftcossinci(const int*, const int*, double*, const int*, double*) ;
void cftcossins(const int*, const int*, double*, const int*, double*) ;
void cftcossinc(const int*, const int*, double*, const int*, double*) ;
void cftleg(const int*, const int*, double*, const int*, double*) ;
void cftlegmp(const int*, const int*, double*, const int*, double*) ;
void cftlegmi(const int*, const int*, double*, const int*, double*) ;
void cftlegp(const int*, const int*, double*, const int*, double*) ;
void cftlegpp(const int*, const int*, double*, const int*, double*) ;
void cftlegi(const int*, const int*, double*, const int*, double*) ;
void cftlegip(const int*, const int*, double*, const int*, double*) ;
void cftlegpi(const int*, const int*, double*, const int*, double*) ;
void cftlegii(const int*, const int*, double*, const int*, double*) ;
void cfrcheb(const int*, const int*, double*, const int*, double*) ;
void cfrchebp(const int*, const int*, double*, const int*, double*) ;
void cfrchebi(const int*, const int*, double*, const int*, double*) ;
void cfrchebpimp(const int*, const int*, double*, const int*, double*) ;
void cfrchebpimi(const int*, const int*, double*, const int*, double*) ;
void cfrchebpip(const int*, const int*, double*, const int*, double*) ;
void cfrchebpii(const int*, const int*, double*, const int*, double*) ;
void cipcossin(const int* , const int* , const int* , double* , double* ) ;
void cipcossini(const int* , const int* , const int* , double* , double* ) ;
void citcos(const int*, const int*, double*, const int*, double*) ;
void citcosp(const int*, const int*, double*, const int*, double*) ;
void citcosi(const int*, const int*, double*, const int*, double*) ;
void citsinp(const int*, const int*, double*, const int*, double*) ;
void citsini(const int*, const int*, double*, const int*, double*) ;
void citcossincp(const int*, const int*, double*, const int*, double*) ;
void citcossinsi(const int*, const int*, double*, const int*, double*) ;
void citcossinsp(const int*, const int*, double*, const int*, double*) ;
void citcossinci(const int*, const int*, double*, const int*, double*) ;
void citcossins(const int*, const int*, double*, const int*, double*) ;
void citcossinc(const int*, const int*, double*, const int*, double*) ;
void citleg(const int*, const int*, double*, const int*, double*) ;
void citlegmp(const int*, const int*, double*, const int*, double*) ;
void citlegmi(const int*, const int*, double*, const int*, double*) ;
void citlegp(const int*, const int*, double*, const int*, double*) ;
void citlegpp(const int*, const int*, double*, const int*, double*) ;
void citlegi(const int*, const int*, double*, const int*, double*) ;
void citlegip(const int*, const int*, double*, const int*, double*) ;
void citlegpi(const int*, const int*, double*, const int*, double*) ;
void citlegii(const int*, const int*, double*, const int*, double*) ;
void circheb(const int*, const int*, double*, const int*, double*) ;
void circhebp(const int*, const int*, double*, const int*, double*) ;
void circhebi(const int*, const int*, double*, const int*, double*) ;
void cirleg(const int*, const int*, double*, const int*, double*) ;
void cirlegp(const int*, const int*, double*, const int*, double*) ;
void cirlegi(const int*, const int*, double*, const int*, double*) ;
void circhebpimp(const int*, const int*, double*, const int*, double*) ;
void circhebpimi(const int*, const int*, double*, const int*, double*) ;
void circhebpip(const int*, const int*, double*, const int*, double*) ;
void circhebpii(const int*, const int*, double*, const int*, double*) ;
double* legendre(int , int ) ;
double* legendre_norm(int , int ) ;
double* mat_cossincp_legp(int, int) ;   
double* mat_cossinci_legi(int, int) ; 
double* mat_cossinc_leg(int, int) ;   
double* mat_cosp_legpp(int, int) ;   
double* mat_cosi_legip(int, int) ;   
double* mat_sini_legpi(int, int) ;   
double* mat_sinp_legii(int, int) ;   
double* mat_cos_legmp(int, int) ;   
double* mat_sin_legmi(int, int) ;   
double* mat_legp_cossincp(int,  int) ;
double* mat_legi_cossinci(int,  int) ;
double* mat_leg_cossinc(int,  int) ;
double* mat_legpp_cosp(int, int) ;   
double* mat_legip_cosi(int, int) ;   
double* mat_legpi_sini(int, int) ;   
double* mat_legii_sinp(int, int) ;   
double* mat_legmp_cos(int, int) ;   
double* mat_legmi_sin(int, int) ;   
void chb_cossincp_legp(const int* , const double* , double* ) ;
void chb_legp_cossincp(const int* , const double* , double* ) ;
void chb_cossinc_leg(const int* , const double* , double* ) ;
void chb_leg_cossinc(const int* , const double* , double* ) ;
void chb_cosp_legpp(const int* , const double* , double* ) ;
void chb_legpp_cosp(const int* , const double* , double* ) ;
void chb_cosi_legip(const int* , const double* , double* ) ;
void chb_legip_cosi(const int* , const double* , double* ) ;
void chb_sini_legpi(const int* , const double* , double* ) ;
void chb_legpi_sini(const int* , const double* , double* ) ;
void chb_cossinci_legi(const int* , const double* , double* ) ;
void chb_legi_cossinci(const int* , const double* , double* ) ;
void chb_sinp_legii(const int* , const double* , double* ) ;
void chb_legii_sinp(const int* , const double* , double* ) ;
void chb_cos_legmp(const int* , const double* , double* ) ;
void chb_legmp_cos(const int* , const double* , double* ) ;
void chb_sin_legmi(const int* , const double* , double* ) ;
void chb_legmi_sin(const int* , const double* , double* ) ;

double int1d_chebp(int, const double* ) ;
double int1d_chebi(int, const double* ) ;
double int1d_cheb(int, const double* ) ;

//Routines Legendre en r
void cirleg(const int*, const int*, double*, const int*, double*) ;
void cirlegp(const int*, const int*, double*, const int*, double*) ;
void cirlegi(const int*, const int*, double*, const int*, double*) ;
void cfrleg(const int*, const int*, double*, const int*, double*) ;
void cfrlegp(const int*, const int*, double*, const int*, double*) ;
void cfrlegi(const int*, const int*, double*, const int*, double*) ;
void legendre_collocation_points(int, double*) ;

// Routines Jacobi
double* jacobi(int, double) ;
double* pointsgausslobatto(int) ;
Tbl jacobipointsgl(int) ;
double* coeffjaco(int, double*) ;
void cfrjaco02(const int*, const int*, double*, const int*, double*);

// Routines calcul de coef inverse
void cipcossin(const int* , const int* , const int* , double* , double* ) ;
void citcosp(const int*, const int*, double*, const int*, double*) ;
void citcosi(const int*, const int*, double*, const int*, double*) ;
void citcos(const int*, const int*, double*, const int*, double*) ;
void citsin(const int*, const int*, double*, const int*, double*) ;
void citsinp(const int*, const int*, double*, const int*, double*) ;
void citsini(const int*, const int*, double*, const int*, double*) ;
void citcossincp(const int*, const int*, double*, const int*, double*) ;
void citcossinsi(const int*, const int*, double*, const int*, double*) ;
void citcossinsp(const int*, const int*, double*, const int*, double*) ;
void citcossinci(const int*, const int*, double*, const int*, double*) ;
void citcossins(const int*, const int*, double*, const int*, double*) ;
void citcossinc(const int*, const int*, double*, const int*, double*) ;
void citlegp(const int*, const int*, double*, const int*, double*) ;
void citlegpp(const int*, const int*, double*, const int*, double*) ;
void citlegi(const int*, const int*, double*, const int*, double*) ;
void circheb(const int*, const int*, double*, const int*, double*) ;
void circhebp(const int*, const int*, double*, const int*, double*) ;
void circhebi(const int*, const int*, double*, const int*, double*) ;
void circhebpimp(const int*, const int*, double*, const int*, double*) ;
void circhebpimi(const int*, const int*, double*, const int*, double*) ;
void cirjaco02(const int*, const int*, double* , const int*, double*) ;

// Routines calculant la matrice du laplacien
Matrice _laplacien_mat_pas_prevu(int, int, double, int) ;
Matrice _laplacien_mat_r_chebp(int, int, double, int) ;
Matrice _laplacien_mat_r_chebi(int, int, double, int) ;
Matrice _laplacien_mat_r_chebu(int, int, double, int) ;
Matrice _laplacien_mat_r_chebu_deux(int, int) ;
Matrice _laplacien_mat_r_chebu_trois(int, int) ;
Matrice _laplacien_mat_r_chebu_quatre(int, int) ;
Matrice _laplacien_mat_r_chebu_cinq(int, int) ;
Matrice _laplacien_mat_r_cheb(int, int, double, int) ;
Matrice laplacien_mat(int , int , double , int, int ) ;

//Routines de passage a bande versions Matrice et Tbl
Matrice _cl_pas_prevu (const Matrice&, int, double, int) ;
Matrice _cl_r_cheb (const Matrice&, int, double, int) ;
Matrice _cl_r_chebi (const Matrice&, int, double, int) ;
Matrice _cl_r_chebu (const Matrice&, int, double, int) ;
Matrice _cl_r_chebu_cinq (const Matrice&, int) ;
Matrice _cl_r_chebu_quatre (const Matrice&, int) ;
Matrice _cl_r_chebu_trois (const Matrice&, int) ;
Matrice _cl_r_chebu_deux (const Matrice&, int) ;
Matrice _cl_r_chebp (const Matrice&, int, double, int) ;
Matrice combinaison (const Matrice&, int, double, int, int) ;

Tbl _cl_pas_prevu (const Tbl&, int) ;
Tbl _cl_r_cheb (const Tbl&, int) ;
Tbl _cl_r_chebi (const Tbl&, int) ;
Tbl _cl_r_chebu (const Tbl&, int) ;
Tbl _cl_r_chebu_deux (const Tbl&) ;
Tbl _cl_r_chebu_trois (const Tbl&) ;
Tbl _cl_r_chebu_quatre (const Tbl&) ;
Tbl _cl_r_chebu_cinq (const Tbl&) ;
Tbl _cl_r_chebp (const Tbl&, int) ;
Tbl combinaison (const Tbl&, int, int) ;


// Routines de preparation du laplacien inversible
Matrice _prepa_nondege_pas_prevu(const Matrice &, int , double, int) ;
Matrice _prepa_nondege_r_cheb (const Matrice&, int, double, int) ;
Matrice _prepa_nondege_r_chebp (const Matrice&, int, double, int) ;
Matrice _prepa_nondege_r_chebi (const Matrice&, int, double, int) ;
Matrice _prepa_nondege_r_chebu (const Matrice&, int, double, int) ;
Matrice _prepa_nondege_r_chebu_deux (const Matrice&, int) ;
Matrice _prepa_nondege_r_chebu_trois (const Matrice&, int) ;
Matrice _prepa_nondege_r_chebu_quatre (const Matrice&, int) ;
Matrice _prepa_nondege_r_chebu_cinq (const Matrice&, int) ;
Matrice prepa_nondege (const Matrice&, int, double, int, int) ;

//Routines de calcul de la solution particuliere
Tbl _solp_pas_prevu(const Matrice&, const Matrice&, double, double, const Tbl&, int) ;
Tbl _solp_r_cheb (const Matrice&, const  Matrice&, double, double, const Tbl&, int) ;
Tbl _solp_r_chebp (const Matrice&, const  Matrice&, double, double, const Tbl&, int) ;
Tbl _solp_r_chebi (const Matrice&, const Matrice&, double, double, const Tbl&, int) ;
Tbl _solp_r_chebu (const Matrice&, const Matrice&, double, double, const Tbl&, int) ;
Tbl _solp_r_chebu_deux (const Matrice&, const Matrice&, const Tbl&) ;
Tbl _solp_r_chebu_trois (const Matrice&, const Matrice&, double, const Tbl&) ;
Tbl _solp_r_chebu_quatre (const Matrice&, const Matrice&, double, const Tbl&) ;
Tbl _solp_r_chebu_cinq (const Matrice&, const Matrice&, const Tbl&) ;
Tbl solp (const Matrice&, const Matrice&, double, double, const Tbl&, int, int) ;

//Routines de calcul des solutions homogenes
Tbl _solh_pas_prevu (int, int, double) ;
Tbl _solh_r_cheb (int, int, double) ;
Tbl _solh_r_chebp (int, int, double) ;
Tbl _solh_r_chebi (int, int, double) ;
Tbl _solh_r_chebu (int, int, double) ;
Tbl solh (int, int, double, int) ;

// Routines helmholtz minus :
Matrice helmholtz_minus_mat(int , int, double , double , double, int ) ;
Matrice cl_helmholtz_minus (const Matrice&, int) ;
Tbl cl_helmholtz_minus (const Tbl&, int) ;
Matrice prepa_helmholtz_minus_nondege (const Matrice&, int) ;
Tbl solp_helmholtz_minus (const Matrice&, const Matrice&, const Tbl&, 
			  double, double, int, int) ;
Tbl solh_helmholtz_minus (int, int, double, double, double, int) ;

// Routines helmholtz plus :
Matrice helmholtz_plus_mat(int , int, double , double , double, int ) ;
Matrice cl_helmholtz_plus (const Matrice&, int) ;
Tbl cl_helmholtz_plus (const Tbl&, int) ;
Matrice prepa_helmholtz_plus_nondege (const Matrice&, int) ;
Tbl solp_helmholtz_plus (const Matrice&, const Matrice&, const Tbl&, 
			  double, double, int) ;
Tbl solh_helmholtz_plus (int, int, double, double, double, int) ;


//Routines de calcul des valeurs limites 
Tbl val_solh (int, double, double, int) ;
Tbl val_solp (const Tbl&, double, int) ;

double val1_dern_1d (int, const Tbl&, int) ;
double valm1_dern_1d (int, const Tbl&, int) ;


//Routines de derivations version 1d
void _d2sdx2_1d_pas_prevu(int, double*, double* ) ;
void _d2sdx2_1d_r_chebu(int, double*, double* ) ;
void _d2sdx2_1d_r_cheb(int, double*, double* ) ;
void _d2sdx2_1d_r_chebp(int, double*, double* ) ;
void _d2sdx2_1d_r_chebi(int, double*, double * ) ;
void d2sdx2_1d(int, double** , int) ;

void _dsdx_1d_pas_prevu(int, double*, double* ) ;
void _dsdx_1d_r_chebu(int, double*, double* ) ;
void _dsdx_1d_r_chebp(int, double*, double* ) ;
void _dsdx_1d_r_chebi(int, double*, double* ) ;
void dsdx_1d(int, double** , int) ;

void _multx_1d_pas_prevu(int, double*, double* ) ;
void _multx_1d_r_cheb(int, double*, double* ) ;
void multx_1d(int, double **, int) ;
void multxpun_1d(int, double **, int) ;

void _sx_1d_pas_prevu(int, double*, double* ) ;
void _sx_1d_r_chebi(int, double*, double* ) ;
void _sx_1d_r_chebp(int, double*, double* ) ;
void sx_1d(int, double **, int) ;

void _sx2_1d_pas_prevu(int, double*, double*) ;
void _sx2_1d_identite(int, double*, double*) ;
void _sx2_1d_r_chebp(int, double*, double*) ;
void _sx2_1d_r_chebi(int, double*, double*) ;
void _sxm12_1d_r_chebu(int, double *, double*) ;
void sx2_1d(int, double**, int) ;

void _sxdsdx_1d_pas_prevu(int, double*, double*) ;
void _dsdx_1d_r_cheb(int, double*, double*) ;
void _sxdsdx_1d_r_chebi(int, double*, double*) ;
void _sxdsdx_1d_r_chebp(int, double*, double*) ;
void sxdsdx_1d(int, double** , int) ;

//Routines de derivations (pour sol_dalembert)
void _dsdx_r_chebp(Tbl *, int &) ;
void _dsdx_r_chebi(Tbl *, int &) ;

// Resolution de l'equation de Poisson
int nullite_plm_sym (int, int, int, int) ;
int nullite_plm_nonsym (int, int, int, int) ;
int nullite_plm_nonsym_anti (int, int, int, int) ;
int nullite_plm (int, int, int, int, Base_val) ;

void donne_lm_sym (int, int, int, int, int&, int&, int&) ;
void donne_lm_nonsym (int, int, int, int, int&, int&, int&) ;
void donne_lm_nonsym_anti (int, int, int, int, int&, int&, int&) ;
void donne_lm (int, int, int, int, Base_val, int&, int&, int&) ;


// Les sommations en r :
void som_r_pas_prevu
    (double*, const int, const int, const int, const double, double*) ;
void som_r_cheb
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebi
    (double*, const int, const int, const int, const double, double*) ;    
void som_r_chebp
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebu
    (double*, const int, const int, const int, const double, double*) ; 
void som_r_chebpim_p
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebpim_i
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebpi_p
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebpi_i
    (double*, const int, const int, const int, const double, double*) ;
void som_r_cheb_symy
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebu_symy
    (double*, const int, const int, const int, const double, double*) ; 
void som_r_chebpim_p_symy
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebpim_i_symy
    (double*, const int, const int, const int, const double, double*) ;
void som_r_cheb_asymy
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebu_asymy
    (double*, const int, const int, const int, const double, double*) ; 
void som_r_chebpim_p_asymy
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebpim_i_asymy
    (double*, const int, const int, const int, const double, double*) ;
void som_r_leg
    (double*, const int, const int, const int, const double, double*) ;
void som_r_legi
    (double*, const int, const int, const int, const double, double*) ;    
void som_r_legp
    (double*, const int, const int, const int, const double, double*) ;
void som_r_jaco02
    (double*, const int, const int, const int, const double, double*) ;
    
// Les sommations en theta :
void som_tet_pas_prevu
    (double*, const int, const int, const double, double*) ;
void som_tet_cos
    (double*, const int, const int, const double, double* ) ;
void som_tet_cos_p
    (double*, const int, const int, const double, double* ) ;
void som_tet_cos_i
    (double*, const int, const int, const double, double* ) ;
void som_tet_sin
    (double*, const int, const int, const double, double* ) ;
void som_tet_sin_p
    (double*, const int, const int, const double, double* ) ;
void som_tet_sin_i
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_cp
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_ci
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_c
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_s
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_sp
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_si
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_cp_symy
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_ci_symy
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_cp_asymy
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_ci_asymy
    (double*, const int, const int, const double, double* ) ;

// Les sommations en phi : 
void som_phi_pas_prevu
    (double*, const int, const double, double* ) ;
void som_phi_cossin
    (double*, const int, const double, double* ) ;
void som_phi_cossin_p
    (double*, const int, const double, double* ) ;
void som_phi_cossin_i
    (double*, const int, const double, double* ) ;
void som_phi_cossin_symy
    (double*, const int, const double, double* ) ;
void som_phi_cossin_asymy
    (double*, const int, const double, double* ) ;

// les divisions et multiplications par x-1 :
void sxm1_1d_cheb(int, double*) ;
void mult_xm1_1d_cheb(int, const double*, double*) ;
void mult2_xm1_1d_cheb(int, const double*, double*) ;

// x * dsdx ...
void _xdsdx_1d_pas_prevu (int, double*, double*) ;
void _xdsdx_1d_r_cheb (int, double*, double*) ;
void _xdsdx_1d_r_chebp (int, double*, double*) ;
void _xdsdx_1d_r_chebi (int, double*, double*) ;
void xdsdx_1d(int, double**, int) ;

// Multiplication par x^2
void multx2_1d(int, double **, int)	;
void _multx2_1d_r_cheb(int, double* , double *);
void _multx2_1d_r_chebp(int, double* , double *);
void _multx2_1d_r_chebi(int, double* , double *);
void _multx2_1d_pas_prevu(int, double* , double *);

// division par (x+1)
void sxpun_1d(int, double **, int)	;
void _sxpun_1d_r_cheb(int, double* , double *);
void _sxpun_1d_pas_prevu(int, double* , double *);
Cmp division_xpun (const Cmp&, int) ;

// Fonctions liees a la resolution des l'equation des ondes
void get_operateur_dal(const Param&, const int&, const int&, 
		       int&, Matrice& );
Tbl dal_inverse(const int&, const int&, const Matrice&, const Tbl&, 
		const bool) ;
Mtbl_cf sol_dalembert(Param&, const Map_af&, const Mtbl_cf&) ;

void runge_kutta3_wave_sys(double, const Scalar&, const Scalar&, Scalar& , Scalar&, int dl=0 ) ;
void evolve_outgoing_BC(double, int, const Scalar&, Scalar&, Tbl&, Tbl&, Tbl&, int dl=0) ;
void tilde_laplacian(const Scalar& B_in, Scalar& tilde_lap, int dl=-1) ;
void initialize_outgoing_BC(int, const Scalar& , const Scalar& , Tbl&) ;

// Fonctions liees aux operateurs elliptiques degeneres: obtention d'espaces-temps de type Kerr
void tensorelliptic ( Scalar source, Scalar& resu, double fitd1, double fit2d1, double fit0d2 = 0., double fit1d2 = 0., double fit0d3 = 0., double fit1d3 = 0.);
 void tensorellipticBt ( Scalar source, Scalar& resu, double fitd1, double fit2d1, double fit0d2 = 0., double fit1d2 = 0., double fit0d3 = 0., double fit1d3 = 0.);
void tensorellipticCt ( Scalar source, Scalar& resu, double fitd1, double fit2d1, double fit0d2, double fit1d2, double fit0d3, double fit1d3);

 Sym_tensor secmembre_kerr ( const Sym_tensor& hij, const Sym_tensor& aa,const Scalar& nn,const Scalar& ppsi,const Vector& bb);

 Sym_tensor boundfree_tensBC( Sym_tensor source, Vector Beta, Scalar Psi, Scalar Nn, Sym_tensor hij_guess, double precision , int loopmax = 250) ;

// Trucs utilises pour poisson_compact :
Matrice lap_cpt_mat(int, int, int) ;
Matrice xdsdx_mat(int, int, int) ;
Matrice combinaison_cpt (const Matrice &, int, int) ;
Tbl combinaison_cpt (const Tbl &, int) ;

// Trucs binaires :
void dirichlet_binaire (const Cmp& source_un, const Cmp& source_deux, 
			const Valeur& boundary_un, const Valeur& boundary_deux,
				Cmp& sol_un, Cmp& sol_deux, int num_front, 
				double precision) ;
void dirichlet_binaire (const Cmp& source_un, const Cmp& source_deux, 
			double bound_un, double bound_deux, 
				Cmp& sol_un, Cmp& sol_deux, int num_front, 
				double precision) ;
void dirichlet_binaire (const Scalar& source_un, const Scalar& source_deux, 
			const Valeur& boundary_un, const Valeur& boundary_deux,
			Scalar& sol_un, Scalar& sol_deux, int num_front, 
			double precision) ;

void neumann_binaire (const Cmp& source_un, const Cmp& source_deux, 
		      const Valeur& boundary_un, const Valeur& boundary_deux, 
		      Cmp& sol_un, Cmp& sol_deux, int num_front, 
		      double precision) ;
void neumann_binaire (const Cmp& source_un, const Cmp& source_deux, 
		      double bound_un, double bound_deux, 
		      Cmp& sol_un, Cmp& sol_deux, int num_front, 
		      double precision) ;
void neumann_binaire (const Scalar& source_un, const Scalar& source_deux, 
		      const Valeur& boundary_un, const Valeur& boundary_deux,
		      Scalar& sol_un, Scalar& sol_deux, int num_front, 
		      double precision) ;

void poisson_vect_frontiere (double lambda, const Tenseur& source, Tenseur& shift, 
	    const Valeur& lim_x, const Valeur& lim_y, const Valeur& lim_z, 
	    int num_front, double precision, int itermax) ;
void poisson_vect_boundary (double lambda, const Vector& source, Vector& shift,
 	    const Valeur& lim_x, const Valeur& lim_y, const Valeur& lim_z, 
	    int num_front, double precision, int itermax) ;

void poisson_vect_binaire ( double lambda, 
		const Tenseur& source_un, const Tenseur& source_deux, 
		const Valeur& bound_x_un, const Valeur& bound_y_un, 
		const Valeur& bound_z_un, const Valeur& bound_x_deux, 
		const Valeur& bound_y_deux, const Valeur& bound_z_deux, 
		Tenseur& sol_un, Tenseur& sol_deux, int num_front, double precision) ;
void poisson_vect_binaire ( double lambda, 
		const Vector& source_un, const Vector& source_deux, 
		const Valeur& bound_x_un, const Valeur& bound_y_un, 
		const Valeur& bound_z_un, const Valeur& bound_x_deux, 
		const Valeur& bound_y_deux, const Valeur& bound_z_deux, 
		Vector& sol_un, Vector& sol_deux, int num_front, double precision) ;

// Elliptic solvers :
Mtbl_cf elliptic_solver  (const Param_elliptic&, const Mtbl_cf&) ;

Mtbl_cf elliptic_solver_boundary  (const Param_elliptic& ope_var, const Mtbl_cf& source, const Mtbl_cf& bound, double fact_dir, double fact_neu ) ;

Mtbl_cf elliptic_solver_no_zec  (const Param_elliptic&, const Mtbl_cf&, 
				 double val) ;
Mtbl_cf elliptic_solver_only_zec  (const Param_elliptic&, const Mtbl_cf&, 
				 double val) ;
Mtbl_cf elliptic_solver_sin_zec  (const Param_elliptic&, const Mtbl_cf&, double*, double*) ;
Mtbl_cf elliptic_solver_fixe_der_zero  (double, 
					    const Param_elliptic&, 
					    const Mtbl_cf&) ;

// Integrale 2D pour les etoiles en rotation
double integrale2d(const Scalar&) ;

// Solution de la composante r de Poisson vectoriel, pour l=0 uniquement
Scalar pois_vect_r0(const Scalar& ) ; 

// Regularisation du shift :
double regle (Tenseur& shift_auto, const Tenseur& shift_comp, double omega, double) ;

// Trucs pour la solution de Misner-Lindquist
double serie_lindquist_plus (double rayon, double distance, double xa, double ya, 
	double za, double precision, double itemax) ;

double serie_lindquist_moins (double rayon, double distance, double xa, double ya, 
	double za, double precision, double itemax) ;

double adm_serie (double rayon, double distance, double precision) ;

double bare_serie (double rayon, double distance, double precision) ;

void set_lindquist (Cmp& psi_un, Cmp& psi_deux, double rayon, double precision) ;

void separation (const Cmp& c1, const Cmp& c2, Cmp& res1, Cmp& res2, int decrois, 
    int puiss, int lmax, double precision, const double relax = 0.5, const int itemax = 100, const int flag = 1) ;

// Spectral cutoff used in tensor elliptic solvers, and solving for stationary black hole spacetimes
 
void coupe_l_tous( Sym_tensor& hij,Sym_tensor& aa, Scalar& nn,Scalar& ppsi, Vector& bb, int ntt, int cutoff);
void tensor_coupe_l( Sym_tensor& ten, int ntt, int cutoff);


}
#endif
