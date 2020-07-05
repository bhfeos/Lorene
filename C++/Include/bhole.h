/*
 *  Definition of Lorene classes Bhole
 *				 Bhole_binaire
 *
 */

/*
 *   Copyright (c) 2000-2001 Philippe Grandclement
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


#ifndef __BHOLE_H_ 
#define __BHOLE_H_ 

/*
 * $Id: bhole.h,v 1.18 2014/10/13 08:52:32 j_novak Exp $
 * $Log: bhole.h,v $
 * Revision 1.18  2014/10/13 08:52:32  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.17  2007/05/15 12:44:18  p_grandclement
 * Scalar version of reevaluate
 *
 * Revision 1.16  2007/04/24 20:15:30  f_limousin
 * Implementation of Dirichlet and Neumann BC for the lapse
 *
 * Revision 1.15  2006/04/27 09:12:29  p_grandclement
 * First try at irrotational black holes
 *
 * Revision 1.14  2005/08/29 15:10:12  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.13  2004/12/02 09:33:04  p_grandclement
 * *** empty log message ***
 *
 * Revision 1.12  2004/03/25 12:35:34  j_novak
 * now using namespace Unites
 *
 * Revision 1.11  2004/03/25 09:14:06  j_novak
 * *** empty log message ***
 *
 * Revision 1.10  2004/03/25 08:53:33  j_novak
 * Error in the doc corrected.
 *
 * Revision 1.9  2004/03/24 17:14:04  j_novak
 * Translation of comments to doxygen
 *
 * Revision 1.8  2003/11/25 07:12:58  k_taniguchi
 * Change the argument of update_metric from the class Etoile_bin to Et_bin_nsbh.
 *
 * Revision 1.7  2003/11/13 13:43:53  p_grandclement
 * Addition of things needed for Bhole::update_metric (const Etoile_bin&, double, double)
 *
 * Revision 1.6  2003/10/24 13:05:48  p_grandclement
 * correction of the equations for Bin_ns_bh...
 *
 * Revision 1.5  2003/02/13 16:40:24  p_grandclement
 * Addition of various things for the Bin_ns_bh project, non of them being
 * completely tested
 *
 * Revision 1.4  2003/01/31 16:57:12  p_grandclement
 * addition of the member Cmp decouple used to compute the K_ij auto, once
 * the K_ij total is known
 *
 * Revision 1.3  2002/12/18 10:28:49  e_gourgoulhon
 *
 * Added the set_mp function.
 *
 * Revision 1.2  2002/09/13 09:17:31  j_novak
 * Modif. commentaires
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.45  2001/06/28  15:08:16  eric
 * Ajout de la fonction area().
 *
 * Revision 2.44  2001/05/16  11:32:45  phil
 * ajout de init_bhole_seul
 *
 * Revision 2.43  2001/05/07  12:35:28  phil
 * mmodifi commentaires
 *
 * Revision 2.42  2001/05/07  11:43:43  phil
 * *** empty log message ***
 *
 * Revision 2.41  2001/05/07  11:40:50  phil
 * mise a jour des commentaires
 *
 * Revision 2.40  2001/05/07  09:30:37  phil
 * *** empty log message ***
 *
 * Revision 2.39  2001/05/07  09:28:37  phil
 * *** empty log message ***
 *
 * Revision 2.38  2001/05/07  09:11:31  phil
 * *** empty log message ***
 *
 * Revision 2.37  2001/04/27  09:25:11  phil
 * fait_decouple est devenu public
 *
 * Revision 2.36  2001/04/26  12:04:11  phil
 * *** empty log message ***
 *
 * Revision 2.35  2001/04/05  13:42:40  phil
 * *** empty log message ***
 *
 * Revision 2.34  2001/04/05  13:34:29  phil
 * ajout resolution en utilisant phi
 *
 * Revision 2.33  2001/03/22  10:49:34  phil
 * *** empty log message ***
 *
 * Revision 2.32  2001/03/22  10:41:25  phil
 * pleins de modif
 *
 * Revision 2.31  2001/02/28  13:22:31  phil
 * modif vire kk_auto
 *
 * Revision 2.30  2001/02/12  15:37:42  phil
 * ajout calcul de J a linfini
 *
 * Revision 2.29  2001/01/29  14:29:56  phil
 * ajout type rotation
 *
 * Revision 2.28  2001/01/24  10:57:08  phil
 * ajout de set_rayonm
 *
 * Revision 2.27  2001/01/22  09:29:08  phil
 * vire les procedures qui servent pas
 *
 * Revision 2.26  2001/01/10  09:31:05  phil
 * modification de fait_kk_auto (membre de binaire maintenant)
 *
 * Revision 2.25  2001/01/09  13:38:15  phil
 * ajout memebre kk_auto
 *
 * Revision 2.24  2000/12/21  10:47:07  phil
 * retour eventuel a 2.21
 *
 * Revision 2.21  2000/12/20  09:07:56  phil
 * ajout set_statiques
 *
 * Revision 2.20  2000/12/18  16:40:39  phil
 * *** empty log message ***
 *
 * Revision 2.19  2000/12/18  16:38:04  phil
 * ajout convergence vers une masse donnee
 *
 * Revision 2.18  2000/12/15  16:41:28  phil
 * ajout calcul de la separation
 *
 * Revision 2.17  2000/12/14  10:44:28  phil
 * ATTENTION : PASSAGE DE PHI A PSI
 *
 * Revision 2.16  2000/12/13  15:35:18  phil
 * ajout calcul bare_masse
 *
 * Revision 2.15  2000/12/04  14:29:37  phil
 * ajout de grad_n_tot
 *
 * Revision 2.14  2000/12/01  16:12:52  phil
 * *** empty log message ***
 *
 * Revision 2.13  2000/12/01  14:16:20  phil
 * *** empty log message ***
 *
 * Revision 2.12  2000/11/24  15:14:44  phil
 * ajout de find_horizon
 *
 * Revision 2.11  2000/11/24  09:57:18  phil
 * ajout calcul masse et moment pour systeme
 *
 * Revision 2.10  2000/11/17  10:03:13  phil
 * ajout de coal
 *
 * Revision 2.9  2000/11/15  18:26:29  phil
 * simplifaction resolution du shift : on bosse a omega bloque
 *
 * Revision 2.8  2000/11/15  12:59:50  phil
 * changement solve_shift_omega
 *
 * Revision 2.7  2000/11/15  09:41:03  phil
 * *** empty log message ***
 *
 * Revision 2.6  2000/11/15  09:39:26  phil
 * ajout viriel
 *
 * Revision 2.5  2000/11/03  12:56:41  phil
 * ajout de const
 *
 * Revision 2.4  2000/10/26  08:20:45  phil
 * ajout verifie_shift
 *
 * Revision 2.3  2000/10/23  09:15:12  phil
 * modif commentaires
 *
 * Revision 2.2  2000/10/20  10:51:13  phil
 * Modif commentaires (minimale !)
 *
 * Revision 2.1  2000/10/20  09:27:28  phil
 * *** empty log message ***
 *
 * Revision 2.0  2000/10/20  09:22:04  phil
 * *** empty log message ***
 *
 *
 * $Header: /cvsroot/Lorene/C++/Include/bhole.h,v 1.18 2014/10/13 08:52:32 j_novak Exp $
 *
 */

#include "tenseur.h"

namespace Lorene {

class Tenseur ;
class Tenseur_sym ;
class Map_af ;
class Bhole_binaire ;
class Et_bin_nsbh ;

#define COROT 0
#define IRROT 1

/**
 * Black hole. \ingroup (star)
 * 
 * This class represents a black hole in a comformaly flat approximation. It is 
 * defined with an affine mapping with a nucleus, not used in calculations, and 
 * must have,  at least two shells. The horizon (i.e. the surface where \e N=0 ) 
 * is located at the inner boudary of the innermost shell (i.e. the domain
 * number 1).
 * 
 * It can described either :
 * 
 * \li \b an \b isolated \b black \b hole . The fields of the companion and the total
 * fields being zero.
 * \li \b a \b black \b hole \b in \b a \b binary \b system . The boost having no direct 
 * meaning in this case.
 * 
 * The tensor \f$K^{ij}\f$ is the the extrinsic curavature tensor with the conformal 
 * factor extracted.
 * 
 * The tensor \f$A^{ij}\f$ denotes \f$\nabla^i N^j + \nabla^j N^i-\frac{2}{3}
 * \nabla_k N^k \delta^{ij}\f$,  that is \f$2NK^{ij}\f$.
 *
 */
class Bhole {

    // Data :
    protected:
	
	Map_af& mp ;  ///< Affine mapping.
	double rayon ; ///< Radius of the horizon in LORENE's units.
	double omega ; ///< Angular velocity in LORENE's units.
	double omega_local ; ///< local angular velocity
	int rot_state ; ///< State of rotation
	double lapse_hori ; /// Lapse on the horizon
	
	/**
	 * Cartesian components of the boost in the reference frame.
	 */
	double* boost ;
	double regul ; ///< Intensity of the correction on the shift vector. 
	
	Tenseur n_auto ; ///< Part of \e N  generated by the hole.
	Tenseur n_comp ; ///<Part of \e N  generated by the companion hole.
	Tenseur n_tot ; ///< Total \e N .
	
	Tenseur psi_auto ; ///< Part of \f$\Psi\f$ generated by the hole.
	Tenseur psi_comp ; ///< Part of \f$\Psi\f$ generated by the companion hole.
	Tenseur psi_tot ; ///< Total \f$\Psi\f$.
	
	Tenseur grad_n_tot ; ///< Total gradient of \e N .
	Tenseur grad_psi_tot ; ///< Total gradient of \f$\Psi\f$.
	
	Tenseur shift_auto ; ///< Part of \f$\beta^i\f$ generated by the hole.
	
	Tenseur_sym taij_auto ; ///< Part of \f$A^{ij}\f$ generated by the hole.
	Tenseur_sym taij_comp ;///< Part of \f$A^{ij}\f$ generated by the companion hole.
	/**
	 * Total \f$A^{ij}\f$,  which must be zero on the horizon of the
	 * regularisation on the shift has been done.
	 */
	Tenseur_sym taij_tot ;
	
	Tenseur_sym tkij_auto ; ///< Auto \f$K^{ij}\f$.
	Tenseur_sym tkij_tot ; ///< Total \f$K^{ij}\f$.
	
	/**
	 * Function used to construct the part of \f$K^{ij}\f$ generated by the hole 
	 * from the total \f$K^{ij}\f$. Only used for a binary system.
	 * 
	 * Mainly this \c Cmp  is 1 around the hole and 0 around the companion
	 * and the sum of \c decouple for the hole and his companion is 1 
	 * everywhere.
	 */
	Cmp decouple ;
	
    //Constructors :
    public:
	/**
	 * Standard constructor. All the fields are set to zero, except
	 * \c rayon ,  which value is deducted from the \c mapping 
	 */
	Bhole (Map_af& mapping) ;
	Bhole (const Bhole&) ;	///< Constructor by copy
	Bhole (Map_af&, FILE*, bool old=false) ; ///< Constructor from a \c Map_af  and a file
	~Bhole() ; ///< Destructor
	
    public:
    //sauvegarde
	void sauve (FILE* fich) const ; ///< Write on a file
	
    public:
	void operator= (const Bhole&) ; ///< Affectation
	
	/// Returns the mapping (readonly).
	const Map_af& get_mp() const {return mp;} ; 

	/// Read/write of the mapping
	Map_af& set_mp() {return mp; } ;

	/**
	 * Returns the radius of the horizon.
	 */
	double get_rayon() const {return rayon;} ;
	
	/**
	 * Sets the radius of the horizon to \c ray .
	 */
	void set_rayon(double ray) {rayon = ray ;} ;
	
	/**
	 * Returns the angular velocity.
	 */
	double get_omega() const {return omega;} ;
	/**
	 * Sets the angular velocity to \c ome .
	 */
	void set_omega(double ome) {omega = ome ;} ;
	/**
	 * Returns the local angular velocity.
	 */
	double get_omega_local() const {return omega_local;} ;
	/**
	 * Sets the local angular velocity to \c ome .
	 */
	void set_omega_local(double ome) {omega_local = ome ;} ;
	/**
	 * Returns the state of rotation.
	 */
	double get_rot_state() const {return rot_state;} ;
	/**
	 * Returns the state of rotation.
	 */
	void set_rot_state(int rotation) {rot_state = rotation;} ;
	/**
	 * Returns the cartesian components of the boost with respect to the 
	 * reference frame.
	 */
	double* get_boost () const {return boost;} ;
	void set_boost (double vx, double vy, double vz) {
	    boost[0] = vx ; boost[1] = vy ; boost[2] = vz ;
	}
	
	/**
	 * Returns the intensity of the regularisation on \f$\vec{N}\f$.
	 */
	double get_regul() const {return regul;} ;
	
	/**
	 * Returns the part of \e N  generated by the hole.
	 */
	const Tenseur& get_n_auto() const {return n_auto ;} ;
	/**
	 * Returns the part of \e N  generated by the hole as a cmp
	 */
	Cmp return_n_auto() const {return n_auto() ;} ;
	/**
	 * Returns the part of \e N  generated by the companion hole.
	 */
	const Tenseur& get_n_comp() const {return n_comp ;} ;
	/**
	 * Returns the total \e N .
	 */
	const Tenseur& get_n_tot() const {return n_tot ;} ;
	
	/**
	 * Returns the part of \f$\Psi\f$ generated by the hole.
	 */
	const Tenseur& get_psi_auto() const {return psi_auto ;} ;
	/**
	 * Returns the part of \f$\Psi\f$  generated by the hole as a cmp
	 */
	Cmp return_psi_auto() const {return psi_auto() ;} ;
	/**
	 * Returns the part of \f$\Psi\f$ generated by the companion hole.
	 */
	const Tenseur& get_psi_comp() const {return psi_comp ;} ;
	/**
	 * Returns the total \f$\Psi\f$.
	 */
	const Tenseur& get_psi_tot() const {return psi_tot ;} ;
	
	/**
	 * Returns the gradient of \f$\Psi\f$.
	 */
	const Tenseur& get_grad_psi_tot() const {return grad_psi_tot ;} ;
	
	/**
	 * Returns the gradient of \e N .
	 */
	const Tenseur& get_grad_n_tot() const {return grad_n_tot ;} ;
	
	
	/**
	 * Returns the part of \f$\beta^i\f$ generated by the hole.
	 */
	const Tenseur& get_shift_auto() const {return shift_auto ;} ;
	
	Cmp return_shift_auto(int i) const {return shift_auto(i) ;} ;
	
	/**
	 * Returns the part of \f$A^{ij}\f$ generated by the hole.
	 */
	const Tenseur& get_taij_auto() const {return taij_auto ;} ;
	/**
	 * Returns the part of \f$A^{ij}\f$ generated by the companion hole.
	 */
	const Tenseur& get_taij_comp() const {return taij_comp ;} ;
	/**
	 * Returns the total \f$A^{ij}\f$.
	 */
	const Tenseur& get_taij_tot() const {return taij_tot ;}
	
	/**
	 * Returns the total \f$K^{ij}\f$.
	 */
	const Tenseur& get_tkij_tot() const {return tkij_tot ;}
	/**
	 * Returns the part of \f$K^{ij}\f$ generated by the hole.
	 */
	const Tenseur& get_tkij_auto() const {return tkij_auto ;}
	/**
	 * Returns the function used to construct \c tkij_auto  from \c tkij_tot .
	 */
	const Cmp get_decouple() const {return decouple ;}
	
	Cmp& set_n_auto() {return n_auto.set() ;}
	Cmp& set_psi_auto() {return psi_auto.set() ;}
	
	void set_lapse_hori(double xx) {lapse_hori = xx;} 
    public:
	
	/**
	 * Imports the part of \e N  due to the companion hole \c comp . The 
	 * total \e N  is then calculated.
	 * 
	 * It also imports the gradient of \e N  and construct the total \f$\nabla N\f$.
	 */
	void fait_n_comp (const Bhole& comp) ;
	
	/**
	 * Imports the part of \f$\Psi\f$ due to the companion hole \c comp . The 
	 * total \f$\Psi\f$ is then calculated.
	 * 
	 * It also imports the gradient of \f$\Psi\f$ and construct the total \f$\nabla \Psi\f$.
	 */
	void fait_psi_comp (const Bhole& comp) ;
		/**
	 * Imports the part of \e N  due to the companion ns \c comp . The 
	 * total \e N  is then calculated.
	 * 
	 * It also imports the gradient of \e N  and construct the total \f$\nabla N\f$.
	 */
	void fait_n_comp (const Et_bin_nsbh& comp) ;
	
	/**
	 * Imports the part of \f$\Psi\f$ due to the companion ns \c comp . The 
	 * total \f$\Psi\f$ is then calculated.
	 * 
	 * It also imports the gradient of \f$\Psi\f$ and construct the total \f$\nabla \Psi\f$.
	 */
	void fait_psi_comp (const Et_bin_nsbh& comp) ;

	/**
	 * Calculates the part of \f$A^{ij}\f$ generated by \c shift_auto .
	 */
	void fait_taij_auto () ;
	

	/**
	 * Corrects \c shift_auto  in such a way that the total \f$A^{ij}\f$ is 
	 * equal to zero in the horizon,  which should ensure the regularity 
	 * of \f$K^{ij}\f$,  using the stored values of the boost and the angular 
	 * velocity.
	 * 
	 * \b WARNING :  this should only be used for a black hole in 
	 * a binary system \c Bhole_binaire  or \c Bin_ns_bh .
	 * 
	 * @param comp [input] : the part of \f$\beta^i\f$ generated by the companion 
	 * hole.
	 */
	void regularise_shift (const Tenseur& comp) ;
	
	/**
	 * Computes the viriel error, that is the difference between the ADM and the Komar 
	 * masses,  calculated by the asymptotic behaviours of respectively \f$\Psi\f$ and \e N .
	 * 
	 * \b WARNING  this should only be used for an isolated black hole.
	 */
	 double viriel_seul () const ;
	 
	/**
	 * Sets the values of the fields to :
      	 * \li \c n_auto  \f$= \frac{1}{2}-2\frac{a}{r}\f$
	 * \li \c n_comp   \f$= \frac{1}{2}\f$
	 * \li \c psi_auto  \f$= \frac{1}{2}+\frac{a}{r}\f$
	 * \li \c psi_comp   \f$= \frac{1}{2}\f$
	 * 
	 * \e a  being the radius of the hole, the other fields being set to zero.
	 */
	void init_bhole () ;
	
	/**
	 *  Set the inital values to those of Kerr
	 * @param masse [input] : ADM mass in LORENE's units.
	 * @param moment [input] : \f$\frac{J}{M}=a\f$, in LORENE's units.
	 * 
	 * The radius \f$r_0\f$ of \c *this is supposed to be equal to 
	 * \f$\frac{1}{2}\sqrt{M^2-a^2}\f$.
	 * 
	 * The fields have then the following values~:
	 * 
	 * \li \f$\Omega = \frac{a}{2M\sqrt{M^2-a^2}}\f$
	 * 
	 * \li \f$\Psi^4 = 
	 * 1+\frac{2M}{r}+\frac{3M^2+a^2\cos^2\theta}{2r^2}+\frac{2Mr_0^2}{r^3}+
	 * \frac{r_0^4}{r^4}\f$.
	 * \li \f$N = \left(1-\frac{2MR}{\Sigma}+\frac{4a^2M^2R^2\sin^2\theta}
	 * {\Sigma^2\left(R^2+a^2\right)+2a^2\Sigma R\sin^2\theta}\right)^{\frac{1}{2}}\f$
	 * \li \f$N^\phi = \frac{2aMR}{\Sigma\left(R^2+a^2\right)+2a^2MR\sin^2\theta}\f$
	 * \end {itemize}
	 * where~:
	 * 
	 * \li \f$R = r+\frac{M^2-a^2}{4r}+M\f$.
	 * \li \f$\Sigma = R^2+a^2-2MR\f$
	 * \li \f$r_0\f$ is equal to \c rayon .
	 * 
	 */
	void init_kerr (double masse, double moment) ;
	
	/**
	 * Solves the equation for \e N ~:
	 *    \f[
	 *\Delta N = -\frac{2}{\Psi}\nabla_i \Psi \nabla^i N 
	 * + N \Psi^4 K_{ij} K^{ij}
	 * \f]
	 * with the condition that \c N =0 on the horizon.
	 * @param relax [input] : the relaxation parameter.
	 * 
	 * \b WARNING  this should only be used for an isolated black hole.
	 */
	void solve_lapse_seul (double relax) ;
	
	/**
	 * Solves the equation for \f$\Psi\f$~:
	 *    \f[
	 * \Delta \Psi = - \frac{\Psi^5} {8} K_{ij}K^{ij}
	 *  \f]
	 * with the condition that \f$\partial_r \Psi =-\frac{1}{2 {\tt rayon}} 
	 * f\left(\theta, \phi\right)\f$ on the horizon, where \c f is 
	 * the value of \f$\Psi\f$ on the horizon at the preceeding step.
	 * @param relax [input] : the relaxation parameter.
	 * 
	 * \b WARNING  this should only be used for an isolated black hole.
	 */
	 
	void solve_psi_seul (double relax) ;
	
	/**
	 * Solves the equation for \f$\vec{\beta}\f$~:
	 *    \f[
	 *\Delta \beta^i +\frac{1}{3} \nabla^i \left(\nabla_j \beta^j\right)
	 * = -6A^{ij}\frac{\nabla_j \Psi}{\Psi} + 2 K^{ij}\nabla_j N
	 * \f]
	 * with \f$\vec{\beta} = -\Omega \vec{m} - \vec{V}\f$ on the horizon, 
	 * \f$\vec{V}\f$ being the boost and \f$\vec{m} = \frac{\partial}{\partial \phi}\f$.
	 * The solution is solved using Oohara-scheme and an iteration.
	 * @param precis [input] : parameter for the Oohara-solver which is an iterative 
	 scheme.
	 * @param relax [input] : the relaxation parameter.
	 * 
	 * \b WARNING  this should only be used for an isolated black hole.
	 */
	void solve_shift_seul (double precis, double relax) ;
	
	/**
	 * Corrects the shift in the innermost shell, so that it remains \f$
	 * {\mathcal{C}}^2\f$ and that \f$A^{ij}\f$ equals zero on the horizon.
	 * 
	 * \c regul is then,  the relative difference between the shift before
	 * and after the regularisation.
	 * 
	 * \b WARNING  this should only be used for an isolated black hole.
	 */
	void regularise_seul () ;
	/**
	 * Solves the equation for \e N ~:
	 *    \f[
	 *\Delta N = -\frac{2}{\Psi}\nabla_i \Psi \nabla^i N 
	 * + N \Psi^4 K_{ij} K^{ij}
	 * \f]
	 * with the condition that \c N =0 on the horizon.
	 * @param relax [input] : the relaxation parameter.
	 * 
	 * \b WARNING  this should only be used for BH in a \c Bin_ns_bh .
	 */
	
	void solve_lapse_with_ns (double relax, int bound_nn, 
				  double lim_nn) ;
	
	/**
	 * Solves the equation for \f$\Psi\f$~:
	 *    \f[
	 (\Delta \Psi = - \frac{\Psi^5} {8} K_{ij}K^{ij}
	 *  \f]
	 * with the condition that \f$\partial_r \Psi=-\frac{1}{2 {\tt rayon}} 
	 * f\left(\theta, \phi\right)\f$ on the horizon, where \c f is the value of \f$\Psi\f$ on 
	 * the horizon at the preceeding step.
	 * @param relax [input] : the relaxation parameter.
	 * 
	 * \b WARNING  this should only be used for BH in a \c Bin_ns_bh .
	 */
	 
	void solve_psi_with_ns (double relax) ;
	
	/**
	 * Solves the equation for \f$\vec{\beta}\f$~:
	 *    \f[
	 *\Delta \beta^i +\frac{1}{3} \nabla^i \left(\nabla_j \beta^j\right)
	 * = -6A^{ij}\frac{\nabla_j \Psi}{\Psi} + 2 K^{ij}\nabla_j N
	 * \f]
	 * with \f$\vec{\beta} = -\Omega \vec{m} - \vec{V}\f$ on the horizon, 
	 * \f$\vec{V}\f$ being the boost and \f$\vec{m} = \frac{\partial}{\partial \phi}\f$.
	 * The solution is solved using Oohara-scheme and an iteration.
	 * @param ns [input] : the companion.
	 * @param precis [input] : parameter for the Oohara-solver which is an iterative 
	 scheme.
	 * @param relax [input] : the relaxation parameter.
	 * 
	 * \b WARNING  this should only be used for BH in a \c Bin_ns_bh .
	 */
	void solve_shift_with_ns (const Et_bin_nsbh& ns, 
				  double precis, double relax,
				  int bound_nn, double lim_nn) ;
	
	void equilibrium (const Et_bin_nsbh& ns, double precis, double relax,
			  int bound_nn, double lim_nn) ;
	void update_metric (const Et_bin_nsbh& ns) ;
        
	
	/**
	 * Calculates the total \f$K^{ij}\f$. The regularisation of the shift must be done
	 * before to ensure regularity.
	 */
	void fait_tkij() ;
	
	/// Computes the area of the throat. 
	double area() const ; 
	
	/**
	 *  Calculates the ADM mass of the black hole using :
	 * \f$M = -\frac{1}{2 \pi} \oint_{\infty} \nabla^i \Psi {\mathrm d} S_i\f$.
	 * 
	 * \b WARNING  this should only be used for an isolated black hole.
	 */
	double masse_adm_seul () const ;
	
	/**
	 *  Calculates the Komar mass of the black hole using :
	 * \f$M = \frac{1}{4 \pi} \oint_{\infty} \nabla^i N {\mathrm d} S_i\f$.
	 * 
	 * \b WARNING  this should only be used for an isolated black hole.
	 */
	double masse_komar_seul() const ;
	
	/**
	 * Calculates the angular momentum of the black hole using the formula at infinity :
	 * \f$ J = -\frac{1}{8 \pi} \oint_{\infty} K_j^i m^j {\mathrm d}S_i\f$
	 * where \f$\vec{m} = \frac{\partial}{\partial \phi}\f$.
	 * 
	 *\b WARNING  It supposes that the boost is zero and should only be 
	 * used for an isolated black hole..
	 */
	double moment_seul_inf() const ;
	
	/**
	 * Calculates the angular momentum of the black hole using the formula on the horizon :
	 * \f$ J = -\frac{1}{8 \pi} \oint_{H} \Psi^6 K_j^i m^j {\mathrm d}S_i\f$
	 * where \f$\vec{m} = \frac{\partial}{\partial \phi}\f$ and \e H denotes the horizon.
	 * 
	 *\b WARNING  It supposes that the boost is zero and should only be 
	 * used for an isolated black hole..
	 */
	double moment_seul_hor() const ;
	
	/*
	* local angular momentum
	*/
	double local_momentum() const ;
	
	/**
	 * Initiates the black hole for a resolution with \f$\Phi = \log \Psi\f$.
	 */
	void init_bhole_phi () ;
	
	/**
	 * Initiates for a single the black hole.
	 * 
	 * \b WARNING  It supposes that the boost is zero and should only be 
	 * used for an isolated black hole..
	 */
	void init_bhole_seul () ;
	
	friend class Bhole_binaire ; ///< Binary black hole system.
	friend class Bin_ns_bh ; ///< Binary NS-BH
} ;

/**
 * Binary black holes system. \ingroup (star) 
 * 
 * This class is intended for dealing with binary black holes configurations 
 * in the conformaly flat approximation.
 */
class Bhole_binaire {
    
    // data :
    private:
	// les deux trous noirs.
	Bhole hole1 ;	///< Black hole one
	Bhole hole2 ;	///< Black hole two
	
	// Tableau sur les deux trous.
	Bhole* holes[2] ; ///< Array on the black holes
	
	double pos_axe ; /// Position of the axis of rotation
	double omega ;	///< Angular velocity
	
    public:
	// constructeurs & destructeur
	Bhole_binaire(Map_af& mp1, Map_af& mp2) ;   ///< Standard constructor
	Bhole_binaire(const Bhole_binaire& ) ;	///< Copy
	~Bhole_binaire() ;  ///< Destructor
	
    public:
	// trucs pour modifier
	void operator=(const Bhole_binaire&) ; ///< Affectation operator
	
	/**
	 * Read/write of a component of the system. \c i  must be equal to
	 * 1 or 2.
	 */
	Bhole& set(int i) 
	    { assert( (i==1) || (i==2) ); 
	      return *holes[i-1] ;} ; 
	/**
	 * Sets the orbital velocity to \c ome
	 */
	void set_omega(double ome) {omega = ome ; 
			     hole1.set_omega (ome) ;
			     hole2.set_omega (ome) ;} ;
	
	void set_pos_axe (double) ;
			     
    public:
	// trucs pour lire :
	/**
	 * Read only of a component of the system. \c i must be equal to
	 * 1 or 2.
	 */
	const Bhole& operator()(int i) const 
	    { assert( (i==1) || (i==2) ); 
	      return *holes[i-1] ;} ;
	      
	/// Returns the angular velocity 
	double get_omega() const {return omega; } ; 
    
	/**
	 * Initialisation of the system. Each hole is set close to a Schwarzschild one 
	 * and the parts of the fields generated by 
	 * the companion are calculated.
	 * 
	 * The angular velocity is set to zero.
	 */
	void init_bhole_binaire() ;
	
	/**
	 * Computes the viriel error, that is the difference between the ADM and the Komar 
	 * masses,  calculated by the asymptotic behaviours of respectively \f$\Psi\f$ and \e N .
	 */
	double viriel() const ;
	/**
	 * Solves the equation for the lapse~:
	 *    \f[
	 * \Delta N_a = -\frac{2}{\Psi}\nabla_i \Psi_a \nabla^i N + N \Psi^4 K_{ij}K_a^{ij}
	 *  \f]
	 * The fields are the total values excpet those with subscript \f$_a\f$, which are 
	 * the fields generated by each holes (\e a  = 1, 2). The boundary conditions are 
	 * such that \e N =0 on both horizons.
	 * The companion contributions are then obtained.
	 * @param precis [input] : precision,  for the boudary conditions are
	 * obtained by iteration.
	 * @param relax [input] : relaxation parameter.
	 */
	void solve_lapse (double precis, double relax) ;
	
	/**
	 * Solves the equation for the conformal factor~:
	 *    \f[
	 *\Delta \Psi_a = -\frac {\Psi^5}{8} K_{ij}K_a^{ij}
	 * \f]
	 * The fields are the total values excpet those with subscript \f$_a\f$, which are 
	 * the fields generated by each holes (\e a  = 1, 2). The boundary conditions are such that 
	 * \f$\partial_r \Psi = -\frac{1}{2 r_0} f\left(\theta, \phi\right)\f$ on both horizons,  \f$r_0\f$ being 
	 * the radii of those horizons and \e f the value of \f$\Psi\f$ on the 
	 * horizons at the preceeding step. The companion contributions are then 
	 * obtained.
	 * @param precis [input] : precision,  for the boudary conditions are being 
	 * obtained by iteration.
	 * @param relax [input] : relaxation parameter.
	 */
	void solve_psi (double precis, double relax) ;
   
	/**
	 * Solves the equation for the shift, using the Oohara-Nakarmure scheme~:
	 *    \f[
	 *\Delta \beta_a^i +\frac{1}{3} \nabla^i \left(\nabla_j \beta_a^j\right) = 
	 * -6A^{ij}\frac{\nabla_i \Psi_a}{\Psi} + 2K^{ij}\nabla_j N_a
	 * \f]
	 * using the current \f$\Omega\f$ for the boudary conditions~:
	 * \f$\vec{N} = -\Omega \vec{m}\f$ on both horizons. 
	 * The fields are the total values excpet those with subscript \f$_a\f$, which are 
	 * the fields generated by each holes (\c a = 1, 2). The companion contributions are then 
	 * obtained.
	 * @param precis [input] : precision for the solver, the boundary 
	 * conditions and the inversion of the operator being 
	 * obtained by iteration.
	 * @param relax [input] : relaxation parameter.
	 */
	void solve_shift (double precis, double relax) ;
	
	/**
	 * Calculation af the extrinsic curvature tensor.
	 * 
	 * The regularisation of the shifts must have been done. All the 
	 * contributions of \f$A^{ij}\f$ are then calculated and the total tensor 
	 * must be zero on both horizons. The computation is done to avoid every singularity 
	 * close to the horizons (division by \e N =0) for it is done in the coefficient space 
	 * in the two regions surrounding the holes.
	 */
	void fait_tkij () ;
	
	/**
	 * Calculates {tt decouple} which is used to obtain \c tkij_auto  by the formula : 
	 * \c tkij_auto  = \c decouple  * \c tkij_tot .
	 * (see the membre {tt Cmp decouple  for more precisions about its value).
	 * 
	 */
	void fait_decouple () ;
    
    public:
	 /**
	  * Initialize the systeme to Misner Lindquist solution, that is solving for \e N  and 
	  * \f$Psi\f$ in the case \f$\Omega = 0\f$.
	  * @param precis [input] : precision for the convergence (on \e N ).
	  * @param relax [input] : relaxation parameter.
	  */
	 
	 void set_statiques (double precis, double relax) ;
	 
	 /**
	  * Solves the equation for a particular angular velocity, the systeme being 
	  * initialized to Misner-Lindquist solution.
	  * @param precis [input] : precision for the convergence (on \f$\beta\f$).
	  * @param relax [input] : relaxation parameter.
	  * @param nbre_ome [input] : number of intermediates velocities to go from 0 to 
	  * \c omega , typically 10.
	  * @param sortie [input] : flag for the output on files (0 no output files).
	  * @returns : the virial error.
	  */
	  void coal (double precis, double relax, int nbre_ome, double search_ome, 
	  		double m1, double m2, const int sortie = 0) ;  

			    
	/**
	 *  Calculates the ADM mass of the system using :
	 * \f$M = -\frac{1}{2 \pi} \oint_{\infty} \nabla^i \Psi {\mathrm d} S_i\f$.
	 */
	double adm_systeme() const ;
	
	/**
	 *  Calculates the Komar mass of the system using :
	 * \f$M = \frac{1}{4 \pi} \oint_{\infty} \nabla^i N {\mathrm d} S_i\f$.
	 */
	double komar_systeme() const ;
	
	/**
	 * Calculates the angular momentum of the black hole using the formula at infinity :
	 * \f$ J = -\frac{1}{8 \pi} \oint_{\infty} K_j^i m^j {\mathrm d}S_i\f$
	 * where \f$\vec{m} = \frac{\partial}{\partial \phi}\f$.
	 */
	double moment_systeme_inf() ;
	
	/**
	 * Calculates the angular momentum of the black hole using the formula on the horizon :
	 * \f$ J = -\sum_{a= 1, 2} \frac{1}{8 \pi} \oint_{H_a} \Psi^6 K_j^i m^j {\mathrm d}S_i\f$
	 * where \f$\vec{m} = \frac{\partial}{\partial \phi}\f$ and \f$H_a\f$ denotes the horizon 
	 * of the hole \e a .
	 * 
	 */
	double moment_systeme_hor() const ;
	
	/**
	 * Calculation of the proper distance between the two spheres of inversion, 
	 * along the x axis.
	 * @param nr [input] : number of points used for the calculation.
	 */
	 double distance_propre(const int nr = 65) const ;
	 
	Tbl linear_momentum_systeme_inf() const ;
	
	/**
	 * Solve the equation for the logarithm of \f$\Psi\f$.
	 */
	void solve_phi (double precision, double relax) ;
	/**
	 * Initiates the system for a resolution using the logarithm of \f$\Psi\f$.
	 */
	void init_phi() ;
} ;

}
#endif
