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
#include "Rna.h"
#include "Protein.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================

class Rna_R : public Rna
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
	Rna_R( GeneticUnit* gen_unit, const Rna_R &model );
	Rna_R( void );
	Rna_R( GeneticUnit* gen_unit );
	Rna_R( GeneticUnit* gen_unit, ae_strand strand, int32_t index, int8_t diff );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Rna_R( void );

    // =================================================================
    //                              Accessors
    // =================================================================
//    inline  ae_list*  get_influence_list( void );
    inline std::vector<ae_protein*> get_protein_list( void );
    // =================================================================
    //                            Public Methods
    // =================================================================
    void    set_influences( std::vector<ae_protein*> protein_list );
    double  get_synthesis_rate( void );
    double  get_affinity_with_protein( int32_t index, ae_protein *protein );
    int32_t get_enhancer_position( void );
    int32_t get_operator_position( void );
    // =================================================================
    //                           Public Attributes
    // =================================================================
    std::vector<ae_protein*> _protein_list;
    std::vector<int> _enhancing_coef_list;
    std::vector<int> _operating_coef_list;

  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*    ae_rna( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
      };*/
    Rna_R( const Rna_R &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };

    // =================================================================
    //                           Protected Methods
    // =================================================================
    //inline  ae_rna_R* copy( void );
    double    affinity_with_protein( int32_t index, ae_protein *protein );

    // =================================================================
    //                          Protected Attributes
    // =================================================================


};

// =====================================================================
//                          Accessors definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================
/*
ae_rna_R* ae_rna_R::copy( void )
{
  return new ae_rna_R( this );
}
*/

} // namespace aevol

#endif // __AE_RNA_R_H__
