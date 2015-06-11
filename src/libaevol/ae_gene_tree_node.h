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
 
 
 #ifndef __AE_GENE_TREE_NODE_H__
#define  __AE_GENE_TREE_NODE_H__
 
 
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
#include "GeneMutation.h"
#include "ae_gene_tree.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================






 
class ae_gene_tree_node
{  
  friend class ae_gene_tree;

  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_gene_tree_node( int32_t nodeCreationDate, Protein * protein);

  
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_gene_tree_node( void );
  
    // =================================================================
    //                            Public Methods
    // =================================================================
    ae_gene_tree_node * search_in_subtree_leaves(const Protein * protein);
    void print_subtree_to_screen(void); // for debug purposes
    void write_subtree_to_files(FILE * topologyFile, FILE * nodeAttributesFile, int32_t end_gener);
    void write_subtree_nodes_in_tabular_file(int32_t treeID, FILE *f); // f must already be open
    void update_pointers_in_subtree_leaves(GeneticUnit * unit);
    void anticipate_mutation_effect_on_genes_in_subtree_leaves(const ae_mutation* mut, int32_t lengthOfGeneticUnit);
    void register_actual_mutation_effect_on_genes_in_subtree_leaves(ae_gene_tree* tree, const ae_mutation* mut, GeneticUnit* unit, int32_t gener, double impact_on_metabolic_error);
  
  protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_gene_tree_node( void )
      {
        printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
        exit( EXIT_FAILURE );
      };
    ae_gene_tree_node( const ae_gene_tree_node &model )
      {
        printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
        exit( EXIT_FAILURE );
      };

    // =================================================================
    //                          Protected Attributes
    // =================================================================
      static int32_t     _nextID;
    
    int32_t             _ID;
    int32_t             _dna_creation_date;   // generation when the DNA of this gene copy was created. For left nodes, equals the creation date of the parent node. 
    int32_t             _node_creation_date;  // generation when this node was created. For right nodes, equals the dna creation date. 
    int32_t             _gene_loss_date;      // generation when the gene was lost: became a pseudogene, or was deleted, or was broken by a rearrangement
    ae_gene_loss_type   _gene_loss_type;      // NOT_LOST_YET or LOST_BY_LOCAL_MUTATION or DELETED or BROKEN_BY_REAR
    ae_strand           _strand;
    int32_t             _shine_dal_position;
    size_t             _nb_promoters;
    int32_t *           _promoter_positions;
    Protein *        _protein_pointer;     // for a leaf (current state of a gene), points to the potein object
                                              // for an internal node (ancestral state of a gene), points to NULL
                                              Rna **           _rna_pointers;         // for a leaf (current state of a gene), points to the RNA object
                                              // for an internal node (ancestral state of a gene), points to NULL
    std::list<GeneMutation *> mutation_list;       // list of ae_gene_mutations since the creation date of the node, i.e. since the last duplication


    ae_gene_tree_node * _left_child;    // NULL until the gene is duplicated, then points to the copy lying on the original DNA segment 
    ae_gene_tree_node * _right_child;   // NULL until the gene is duplicated, then points to the copy lying on the duplicated DNA segment 
                                        // (which was reinserted elsewhere in the genome, possibly on another genetic unit)
    ae_gene_tree_node * _parent_node;   // points to the node that corresponds to the state of the gene before the last duplication (NULL for the root of the gene tree)

    bool _cds_possibly_modified;
    bool _cds_completely_deleted;
    bool _proms_possibly_modified;
    bool _gene_possibly_duplicated;
    int32_t _putative_position_for_the_duplicate;
    
    
};
} // namespace aevol

#endif // __AE_GENE_TREE_NODE_H__
