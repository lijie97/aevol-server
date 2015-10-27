// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
// 
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// 
// ****************************************************************************


#ifndef __AE_INFLUENCE_R_H__
#define  __AE_INFLUENCE_R_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>

// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "ae_protein.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================
class ae_influence_R
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_influence_R( ae_rna* rna, ae_protein *protein, double enhancing_coef, double operating_coef );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_influence_R( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline ae_rna* get_rna( void );
    inline double  get_enhancing_coef( void );
    inline double  get_operating_coef( void );

    // =================================================================
    //                            Public Methods
    // =================================================================
    inline double           get_enhancer_activity( void );
    inline double           get_operator_activity( void );
    //inline ae_influence_R*  copy( void );

    // =================================================================
    //                           Public Attributes
    // =================================================================

  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*    ae_influence_R( const ae_influence_R &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
      };*/

    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_protein* protein_;
    ae_rna*     rna_;
    double      enhancing_coef_;
    double      operating_coef_;
};

// =====================================================================
//                          Accessors definitions
// =====================================================================
ae_rna* ae_influence_R::get_rna( void )
{
  return rna_;
}

double ae_influence_R::get_enhancing_coef( void )
{
  return enhancing_coef_;
}

double ae_influence_R::get_operating_coef( void )
{
  return operating_coef_;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
double ae_influence_R::get_enhancer_activity( void )
{
  return protein_->get_concentration() * enhancing_coef_;
}

double ae_influence_R::get_operator_activity( void )
{
  return protein_->get_concentration() * operating_coef_;
}

/*
ae_influence_R* ae_influence_R::copy( void )
{
  return new ae_influence_R( rna_, protein_, enhancing_coef_, operating_coef_ );
}
*/
} // namespace aevol
#endif // __AE_INFLUENCE_R_H__
