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


/** \class GeneTreeNode
 *  \brief Currently used only by post-treatments, on a specific lineage, to monitor the fate of paralogs.
 *         Each node corresponds to a coding RNA. When it is duplicated, two new nodes are added in the tree,
 *         as children of the ancestral version. The left child corresponds to the original DNA segment, while
 *         the right child corresponds to the copy that was reinserted elsewhere, possibly in another genetic 
 *         unit.  
 */
 
 
 #ifndef AEVOL_GENE_TREE_H__
#define  AEVOL_GENE_TREE_H__
 
 
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
#include "Rna.h"

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




class GeneTreeNode;

 
class GeneTree
{  
  friend class GeneTreeNode;

  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    GeneTree( void );
    GeneTree( int32_t nodeCreationDate, Protein * protein, const Mutation * mut = NULL ); // Creates a tree with just a root node.
   

  
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~GeneTree( void );
  
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
    void anticipate_mutation_effect_on_genes_in_tree_leaves(const Mutation * mut, int32_t lengthOfGeneticUnit);
    void register_actual_mutation_effect_on_genes_in_tree_leaves(const Mutation * mut, GeneticUnit * unit, int32_t gener, double impact_on_metabolic_error);

    /* void duplicate_this_gene(GeneTreeNode * node, int32_t duplicDate, Protein * newProtein); */
    /* void report_gene_mutation(GeneTreeNode * node, GeneMutation * geneMut); */
    /* void report_gene_loss(GeneTreeNode * node, int32_t geneLossDate, ae_gene_loss_type geneLossType); */

    GeneTreeNode * search_in_leaves(const Protein * protein); // Returns NULL if not found
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

    GeneTree( const GeneTree &model )
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
     
    GeneTreeNode * root_;
    ae_gene_creation_type creation_type_;
    int32_t begin_gener_;
    int32_t end_gener_;
    int32_t total_nb_nodes_;
    int32_t nb_internal_nodes_;
    int32_t nb_leaves_;
    int32_t nb_active_leaves_;
    
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

inline ae_gene_creation_type GeneTree::get_creation_type() const
{
  return creation_type_;
}

inline int32_t GeneTree::get_begin_gener() const
{
  return begin_gener_;
}

inline int32_t GeneTree::get_end_gener() const
{
  return end_gener_;
}

inline int32_t GeneTree::get_total_nb_nodes() const
{
  return total_nb_nodes_;
}

inline int32_t GeneTree::get_nb_internal_nodes() const
{
  return nb_internal_nodes_;
}

inline int32_t GeneTree::get_nb_leaves() const
{
  return nb_leaves_;
}

inline int32_t GeneTree::get_nb_active_leaves() const
{
  return nb_active_leaves_;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol

#endif // AEVOL_GENE_TREE_H__
