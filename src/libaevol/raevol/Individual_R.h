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


#ifndef __AE_INDIVIDUAL_R_H__
#define  __AE_INDIVIDUAL_R_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>

// =================================================================
//                            Project Files
// =================================================================
#include "Individual.h"
#include "ae_rna_R.h"
#include "ae_protein_R.h"
#include "Habitat.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================
class ae_individual_R : public virtual Individual
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
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
    /**
      * Main evaluation method
      */
     virtual void Evaluate();
     /**
      * Evaluate within the provided context
      */
     virtual void EvaluateInContext(const Habitat& habitat);
     virtual void reevaluate();
     virtual void clear_everything_except_dna_and_promoters();
     void do_transcription_translation_folding();
     void do_transcription();
     void do_translation();
     void do_folding();
     void compute_phenotype();
     void compute_distance_to_target(const PhenotypicTarget& target);

    void    set_influences( void );
    void    update_concentrations( void );
    void    multiply_concentrations( double factor );
    int8_t  get_quadon( GeneticUnit* gen_unit, ae_strand strand, int32_t pos );
    void    save( gzFile backup_file );

    inline std::vector<ae_protein_R*> get_inherited_protein_list( void) const;
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
    std::vector<ae_protein_R*> _inherited_protein_list;
    std::vector<ae_rna_R*> _rna_list_coding;   // Please note that these RNAs are
                                // actually managed via genetic units.

};

// =====================================================================
//                          Accessors definitions
// =====================================================================
inline std::vector<ae_protein_R*> ae_individual_R::get_inherited_protein_list( void ) const
{
  assert( _inherited_protein_list );

  return _inherited_protein_list;
}

} // namespace aevol

#endif // __AE_INDIVIDUAL_R_H__
