//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2004  LIRIS.
// Web: https://liris.cnrs.fr/
// E-mail: carole.knibbe@liris.cnrs.fr
// Original Authors : Guillaume Beslon, Carole Knibbe, Virginie Lefort
//                    David Parsons
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//*****************************************************************************


/** \class
 *  \brief
 */


#ifndef __AE_INFLUENCE_R_H__
#define  __AE_INFLUENCE_R_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>

// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_protein.h>

// =================================================================
//                          Class declarations
// =================================================================
class ae_influence_R : public ae_object
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
    ae_protein* _protein;
    ae_rna*     _rna;
    double      _enhancing_coef;
    double      _operating_coef;
};

// =====================================================================
//                          Accessors definitions
// =====================================================================
ae_rna* ae_influence_R::get_rna( void )
{
  return _rna;
}

double ae_influence_R::get_enhancing_coef( void )
{
  return _enhancing_coef;
}

double ae_influence_R::get_operating_coef( void )
{
  return _operating_coef;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
double ae_influence_R::get_enhancer_activity( void )
{
  return _protein->get_concentration() * _enhancing_coef;
}

double ae_influence_R::get_operator_activity( void )
{
  return _protein->get_concentration() * _operating_coef;
}

/*
ae_influence_R* ae_influence_R::copy( void )
{
  return new ae_influence_R( _rna, _protein, _enhancing_coef, _operating_coef );
}
*/

#endif // __AE_INFLUENCE_R_H__
