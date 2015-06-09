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


/** \class ae_gene_tree_node
 *  \brief Currently used only by post-treatments, on a specific lineage, to monitor the fate of paralogs.
 *         Each node corresponds to a coding RNA. When it is duplicated, two new nodes are added in the tree,
 *         as children of the ancestral version. The left child corresponds to the original DNA segment, while
 *         the right child corresponds to the copy that was reinserted elsewhere, possibly in another genetic 
 *         unit.  
 */
 
 
 #ifndef __AE_GENE_TREE_H__
#define  __AE_GENE_TREE_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>



// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "ae_rna.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================



enum ae_gene_creation_type
{
  INITIALIZATION = 0,
  LOCAL_MUTATION  = 1,
  REARRANGEMENT = 2,
  TRANSFER = 3
};




class ae_gene_tree_node;

 
class ae_gene_tree
{  
  friend class ae_gene_tree_node;

  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_gene_tree( void );
    ae_gene_tree( int32_t nodeCreationDate, Protein * protein, const ae_mutation* mut = NULL ); // Creates a tree with just a root node.
   

  
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_gene_tree( void );
  
    // =================================================================
    //                              Accessors
    // =================================================================
    
    inline ae_gene_creation_type get_creation_type() const;
    inline int32_t get_begin_gener() const;
    inline int32_t get_end_gener() const;
    inline int32_t get_total_nb_nodes() const;
    inline int32_t get_nb_internal_nodes() const;
    inline int32_t get_nb_leaves() const;
    inline int32_t get_nb_active_leaves() const;
  

    // =================================================================
    //                            Public Methods
    // =================================================================

    void set_end_gener_if_active_leaves(int32_t gener);
    void update_pointers_in_tree_leaves(GeneticUnit * unit);
    void anticipate_mutation_effect_on_genes_in_tree_leaves(const ae_mutation* mut, int32_t lengthOfGeneticUnit);
    void register_actual_mutation_effect_on_genes_in_tree_leaves(const ae_mutation* mut, GeneticUnit * unit, int32_t gener, double impact_on_metabolic_error);

    /* void duplicate_this_gene(ae_gene_tree_node * node, int32_t duplicDate, Protein * newProtein); */
    /* void report_gene_mutation(ae_gene_tree_node * node, ae_gene_mutation * geneMut); */
    /* void report_gene_loss(ae_gene_tree_node * node, int32_t geneLossDate, ae_gene_loss_type geneLossType); */

    ae_gene_tree_node * search_in_leaves(const Protein * protein); // Returns NULL if not found
    void print_to_screen(void); // for debug purposes
    void write_to_files(const char * topologyFileName, const char * nodeAttributesFileName, int32_t end_gener); // Newick format for the topology file
    void write_nodes_in_tabular_file(int32_t treeID, FILE * f);  // f must already be open 

  
    // =================================================================
    //                           Public Attributes
    // =================================================================
  
  
  
  
  
  protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================

    ae_gene_tree( const ae_gene_tree &model )
      {
        printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
        exit( EXIT_FAILURE );
      };
    
  
    // =================================================================
    //                           Protected Methods
    // =================================================================

  
    // =================================================================
    //                          Protected Attributes
    // =================================================================
     
    ae_gene_tree_node * _root;
    ae_gene_creation_type _creation_type;
    int32_t _begin_gener;
    int32_t _end_gener;
    int32_t _total_nb_nodes;
    int32_t _nb_internal_nodes;
    int32_t _nb_leaves;
    int32_t _nb_active_leaves;
    
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

inline ae_gene_creation_type ae_gene_tree::get_creation_type() const
{
  return _creation_type;
}

inline int32_t ae_gene_tree::get_begin_gener() const
{
  return _begin_gener;
}

inline int32_t ae_gene_tree::get_end_gener() const
{
  return _end_gener;
}

inline int32_t ae_gene_tree::get_total_nb_nodes() const
{
  return _total_nb_nodes;
}

inline int32_t ae_gene_tree::get_nb_internal_nodes() const
{
  return _nb_internal_nodes;
}

inline int32_t ae_gene_tree::get_nb_leaves() const
{
  return _nb_leaves;
}

inline int32_t ae_gene_tree::get_nb_active_leaves() const
{
  return _nb_active_leaves;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol

#endif // __AE_GENE_TREE_H__
