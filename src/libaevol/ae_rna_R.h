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
//*****************************************************************************


#ifndef  __AE_RNA_R_H__
#define  __AE_RNA_R_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>

// =================================================================
//                            Project Files
// =================================================================
#include <ae_rna.h>
#include <ae_protein.h>
#include <ae_influence_R.h>

// =================================================================
//                          Class declarations
// =================================================================

class ae_rna_R : public ae_rna
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_rna_R( ae_genetic_unit* gen_unit, const ae_rna_R &model );
    ae_rna_R( void );
    ae_rna_R( ae_genetic_unit* gen_unit );
    ae_rna_R( ae_genetic_unit* gen_unit, ae_strand strand, int32_t index, int8_t diff );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_rna_R( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline  ae_list*  get_influence_list( void );

    // =================================================================
    //                            Public Methods
    // =================================================================
    void    set_influences( ae_list* protein_list );
    double  get_synthesis_rate( void );
    void    remove_influence( ae_influence_R* influence );
    // =================================================================
    //                           Public Attributes
    // =================================================================

  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*    ae_rna( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
      };*/
    ae_rna_R( const ae_rna_R &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };

    // =================================================================
    //                           Protected Methods
    // =================================================================
    //inline  ae_rna_R* copy( void );
            int32_t   get_enhancer_position( void );
            int32_t   get_operator_position( void );
            void      add_influence( ae_protein *protein, double enhancing_coef, double operating_coef );
            double    affinity_with_protein( int32_t index, ae_protein *protein );

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_list*  _influence_list;

};

// =====================================================================
//                          Accessors definitions
// =====================================================================
ae_list* ae_rna_R::get_influence_list( void )
{
  return _influence_list;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
/*
ae_rna_R* ae_rna_R::copy( void )
{
  return new ae_rna_R( this );
}
*/

#endif // __AE_RNA_R_H__
