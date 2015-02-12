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


/*! \class ae_gene_mutation
    \brief Currently used only by post-treatments, on a specific lineage, to monitor the fate of paralogs.
         Each paralog maintains a list of the mutations it underwent. A gene mutation is a mutation, but
         enriched with the generation when it occurred and the position where it occurred in the coding RNA
         (relative to the first bp of the promoter).
*/
 
 
 #ifndef __AE_GENE_MUTATION_H__
#define  __AE_GENE_MUTATION_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>




// =================================================================
//                            Project Files
// =================================================================
#include "ae_mutation.h"
#include "ae_rna.h"
#include "ae_enums.h"


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================



enum ae_gene_loss_type
{
  NOT_LOST_YET = 0,
  LOST_BY_LOCAL_MUTATION  = 1,
  DELETED = 2,
  BROKEN_BY_REAR = 3,
  DUPLICATED = 4
};




enum ae_gene_mutation_region
{
  UPSTREAM = 0,
  CDS = 1,
  BOTH = 2,
};





class ae_gene_mutation : public ae_mutation
{
  friend class ae_gene_tree_node;
  
 public :

  // =================================================================
  //                             Constructors
  // =================================================================
  
  // Creates a copy of the mutation mut, but enriched with the generation when it occured
  // and the position where it occurred in the RNA, relative to the first bp of the promoter
  ae_gene_mutation(ae_mutation const & mut, int32_t gener, int32_t cdsPosBefore, ae_strand strandBefore, ae_gene_mutation_region region );
  
  ae_gene_mutation( const ae_gene_mutation &model );
  
  // =================================================================
  //                             Destructors
  // =================================================================
  
  virtual ~ae_gene_mutation();
  
  // =================================================================
  //                              Accessors
  // =================================================================
  
  inline int32_t get_generation() const;
  inline double get_impact_on_metabolic_error() const;
  inline ae_gene_mutation_region get_region();
  inline void set_impact_on_metabolic_error(double impact);
  
 
  // =================================================================
  //                            Public Methods
  // =================================================================
  void get_description_string_for_gene_mut(char * str); // str must be at least of size 60
  int8_t type_of_event(); // 0 if local mut, 1 if rearrangement, 2 if transfer

     
 protected :
  
  // =================================================================
  //                         Forbidden Constructors
  // =================================================================
  
  ae_gene_mutation( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
  /* ae_gene_mutation( const ae_gene_mutation &model ) */
  /*   { */
  /*     printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ ); */
  /*     exit( EXIT_FAILURE ); */
  /*   }; */
  
  // =================================================================
  //                           Protected Methods
  // =================================================================
  
  // =================================================================
  //                          Protected Attributes
  // =================================================================
  
  int32_t*  _position_relative_to_shine_dal; /* array of positions similar to the _pos array of the ae_mutation class (size 1 for the switch, 2 for an inversion, etc.) */
  int32_t   _generation;
  double    _impact_on_metabolic_error;
  ae_gene_mutation_region _region;
  
};





// =====================================================================
//                         Inline Accessors' definitions
// =====================================================================

inline int32_t ae_gene_mutation::get_generation() const
{
  return _generation;
}

inline double ae_gene_mutation::get_impact_on_metabolic_error() const
{
  return _impact_on_metabolic_error;

}


inline void ae_gene_mutation::set_impact_on_metabolic_error(double impact)
{
  _impact_on_metabolic_error = impact;
}


inline ae_gene_mutation_region ae_gene_mutation::get_region()
{
  return _region;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // __AE_GENE_MUTATION_H__
