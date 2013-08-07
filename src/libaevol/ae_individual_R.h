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


#ifndef __AE_INDIVIDUAL_R_H__
#define  __AE_INDIVIDUAL_R_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>

// =================================================================
//                            Project Files
// =================================================================
#include <ae_individual.h>
#include <ae_rna_R.h>
#include <ae_protein_R.h>

// =================================================================
//                          Class declarations
// =================================================================
class ae_individual_R : public virtual ae_individual
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_individual_R( const ae_individual_R &model );
    ae_individual_R( void );
    ae_individual_R(  ae_individual_R* parent, int32_t id,
                      ae_jumping_mt* mut_prng, ae_jumping_mt* stoch_prng  );
    ae_individual_R( gzFile backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_individual_R( void );

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    virtual void evaluate( ae_environment* envir );
    void    set_influences( void );
    void    update_concentrations( void );
    void    multiply_concentrations( double factor );
    int8_t  get_quadon( ae_genetic_unit* gen_unit, ae_strand strand, int32_t pos );
    void    save( gzFile backup_file );
    // =================================================================
    //                           Public Attributes
    // =================================================================

  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*    ae_individual_R( const ae_individual &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
      };*/

    // =================================================================
    //                           Protected Methods
    // =================================================================
    virtual void make_protein_list( void );
    virtual void make_rna_list( void );

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_list* _inherited_protein_list;
    ae_list* _rna_list_coding;  // Please note that these RNAs are
                                // actually managed via genetic units.

};

// =====================================================================
//                          Accessors definitions
// =====================================================================

#endif // __AE_INDIVIDUAL_R_H__
