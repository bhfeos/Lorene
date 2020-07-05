/*
 *  Definition of Lorene class Change_var
 *
 */

/*
 *   Copyright (c) 2003 Philippe Grandclement
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

#ifndef __CHANGE_VAR_H_ 
#define __CHANGE_VAR_H_ 

/*
 * $Id: change_var.h,v 1.6 2014/10/13 08:52:32 j_novak Exp $
 * $Log: change_var.h,v $
 * Revision 1.6  2014/10/13 08:52:32  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2004/06/22 08:49:56  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * Revision 1.4  2004/05/14 08:51:00  p_grandclement
 * *** empty log message ***
 *
 * Revision 1.3  2004/03/23 14:54:45  j_novak
 * More documentation
 *
 * Revision 1.2  2004/03/05 09:18:48  p_grandclement
 * Addition of operator sec_order_r2
 *
 * Revision 1.1  2003/12/11 14:57:00  p_grandclement
 * I had forgotten the .h (sorry folks...)
 *
 *
 * $Header: /cvsroot/Lorene/C++/Include/change_var.h,v 1.6 2014/10/13 08:52:32 j_novak Exp $
 *
 */

#include "proto.h"

// Defines the various types of variable changes with the associated 
// functions F, G, F' and G'

#define STD 1
double one (double) ;
double zero (double) ;

#define W_BETA 2
double ide (double) ;

#define W_BETA_INF 3
double part_ln (double) ;
double part_ln_der (double) ;

#define H_BETA 4

#define LAMBDA_RN 5
double plus_log(double) ;
double moins_sur(double) ;

#define NU_RN 6
double moins_log(double) ;
double plus_sur(double) ;

namespace Lorene {
/**
 * This class defines a variable change to be used when solving 
 * elliptic equations. In {\em one particular domain}, 
 * the change is defined as 
 * follows : \ingroup (ellip)
 * 
 * \f$W = F\left(r\right) + G \left(r\right) w\f$
 *
 * where \f$W\f$ is the usual variable, i.e. the one that is \f$\mathcal{C}^1\f$ 
 * and \f$w\f$ is the variable used for solving the elliptic equation.
 * The functions \f$F\f$ and \f$G\f$ are arbitrary functions of \f$r\f$.
 * The type of change must be explicitely implemented in the constructor for 
 * {\tt Change_var}.
 * 
 **/

class Change_var {

 protected:
  double (* func_F) (double) ; /// Pointer on the function \f$F\f$.
  double (* der_F) (double) ; /// Pointer on the derivative of \f$F\f$.
  double (* func_G) (double) ; /// Pointer on the function \f$G\f$.
  double (* der_G) (double) ; /// Pointer on the derivative of \f$G\f$.

  double mult_F ; /// Multiplicative factor for F ## PROVISORY
  double add_F ; /// Additive factor for F ## PROVISORY

 public:

  /**
   * Standard constructor. Case with without parameter.
   * {\tt var} defines explicitely the type of variable to be used.
   * Are currently implemented :
   * \begin{itemize}
   * \item {\tt var} \f$= {\rm STD} \Longrightarrow F=0 \, \& \, G=1\f$.
   * \item {\tt var} \f$= {\rm W\_BETA} \Longrightarrow F=1 \, \& \, G = r\f$.
   * \item {\tt var} \f$= {\rm W\_BETA\_INF} \Longrightarrow F = 1+ \frac{1}{3}r^2\ln r \, \& \, G = r\f$.
   * \item {\tt var} \f$= {\rm H\_BETA} \Longrightarrow F=1 \, \&\, G=1\f$.
   * \item {\tt var} \f$= {\rm LAMBDA_RN} \Longrightarrow F=-\ln r \, \& \, G=1\f$.
   * \item {\tt var} \f$= {\rm NU_RN} \Longrightarrow F=\ln r \, \& \, G=1\f$.
   * \end{itemize} 
   **/
  
  Change_var (int var) ;

  /**
   * Standard constructor. Case with one multiplicative parameter for F
   * {\tt var} defines explicitely the type of variable to be used.
   **/

  Change_var (int var, double) ;

  /**
   * Standard constructor. Case with one multiplicative parameter for F and one additive.
   * {\tt var} defines explicitely the type of variable to be used.
   **/

  Change_var (int var, double, double) ;
 

  Change_var (const Change_var& so) ; /// Constructor by copy.
  ~Change_var() ; /// Standard destructor.

 public:
  double val_F (double x) ; /// Returns the value of \f$F\f$ at {\tt x}.
  double val_der_F (double x) ; /// Returns the value of \f$F'\f$ at {\tt x}.
  double val_G (double x) ; /// Returns the value of \f$G\f$ at {\tt x}.
  double val_der_G (double x) ; /// Returns the value of \f$G'\f$ at {\tt x}.
} ;

}
#endif
