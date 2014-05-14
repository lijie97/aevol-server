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


// =================================================================
//                              Libraries
// =================================================================
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>



// =================================================================
//                            Project Files
// =================================================================

#include <ae_gene_tree.h>
#include <ae_gene_tree_node.h>





//##############################################################################
//                                                                             #
//                            Class ae_gene_tree                               #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================



// =================================================================
//                             Constructors
// =================================================================

ae_gene_tree::ae_gene_tree( void )
{
  _root = NULL;
  _begin_gener = 0;
  _end_gener = 0;
  _total_nb_nodes = 0;
  _nb_internal_nodes = 0;
  _nb_leaves = 0;
  _nb_active_leaves = 0;
  _creation_type = INITIALIZATION;
}

// Creates a tree with just a root node.
ae_gene_tree::ae_gene_tree( int32_t nodeCreationDate, ae_protein * protein, ae_mutation * mut /* = NULL */ )
{
  _root = new ae_gene_tree_node(nodeCreationDate, protein);
  if (mut == NULL)  _creation_type = INITIALIZATION;
  else if ((mut->get_mut_type() == SWITCH) || (mut->get_mut_type() == S_INS) || (mut->get_mut_type() == S_DEL)) _creation_type = LOCAL_MUTATION;
  else if ((mut->get_mut_type() == DUPL) || (mut->get_mut_type() == DEL) || (mut->get_mut_type() == TRANS) || (mut->get_mut_type() == INV) ) _creation_type = REARRANGEMENT;
  else _creation_type = TRANSFER;
 
    
  _begin_gener = nodeCreationDate;
  _end_gener = nodeCreationDate;
  _total_nb_nodes = 1;
  _nb_internal_nodes = 0;
  _nb_leaves = 1;
  if (protein != NULL) _nb_active_leaves = 1; else _nb_active_leaves = 0;
}




// =================================================================
//                             Destructors
// =================================================================

ae_gene_tree::~ae_gene_tree( void )
{
  delete _root;
}


// =================================================================
//                           Public Methods
// =================================================================

void ae_gene_tree::set_end_gener_if_active_leaves(int32_t gener)
{
  if (_nb_active_leaves > 0) _end_gener = gener;
}

void ae_gene_tree::update_pointers_in_tree_leaves(ae_genetic_unit * unit)
{
  _root->update_pointers_in_subtree_leaves(unit);
}

void ae_gene_tree::anticipate_mutation_effect_on_genes_in_tree_leaves(ae_mutation * mut, int32_t lengthOfGeneticUnit)
{
  _root->anticipate_mutation_effect_on_genes_in_subtree_leaves(mut, lengthOfGeneticUnit);
}

void ae_gene_tree::register_actual_mutation_effect_on_genes_in_tree_leaves( ae_mutation * mut, ae_genetic_unit * unit, int32_t gener, double impact_on_metabolic_error)
{
  _root->register_actual_mutation_effect_on_genes_in_subtree_leaves( this, mut, unit, gener, impact_on_metabolic_error);
}

// void ae_gene_tree::duplicate_this_gene( ae_gene_tree_node * node, int32_t duplicDate, ae_protein * newProtein )
// {
//   if (newProtein == node->_protein_pointer) {fprintf(stderr, "Error, duplication with the same protein\n"); exit(EXIT_FAILURE);}

//   // Create a new node for the "old" DNA segment
//   node->_left_child = new ae_gene_tree_node(duplicDate, node->_protein_pointer);
//   node->_left_child->_node_creation_date = duplicDate;
//   node->_left_child->_dna_creation_date = node->_dna_creation_date;
//   node->_left_child->_parent_node = node;

//   // Create a new node for the "new" DNA segment
//   node->_right_child = new ae_gene_tree_node(duplicDate, newProtein);
//   node->_right_child->_node_creation_date = duplicDate;
//   node->_right_child->_dna_creation_date = duplicDate;
//   node->_right_child->_parent_node = node;
  
//   // This node becomes internal, it represents an ancestral (obsolete) state of the gene
//   node->_protein_pointer = NULL;
//   for (int32_t i = 0; i < node->_nb_promoters; i++) {node->_rna_pointers[i] = NULL;}
//   node->_gene_loss_type = DUPLICATED; 
//   node->_gene_loss_date = duplicDate;

//   // Update tree statistics
//   _total_nb_nodes += 2;
//   _nb_internal_nodes ++;
//   _nb_leaves ++;  // - 1 + 2 (the ex-leaf becomes an internal node, 2 leaves are created)
//   if (newProtein != NULL) _nb_active_leaves ++;
//   if (duplicDate > _end_gener) _end_gener = duplicDate;
// }


// void ae_gene_tree::report_gene_mutation( ae_gene_tree_node * node, ae_gene_mutation * geneMut)
// {
//   node->_mutation_list->add(geneMut);
//   if (geneMut->get_generation() > _end_gener) _end_gener = geneMut->get_generation();
// }
 


// void ae_gene_tree::report_gene_loss( ae_gene_tree_node * node, int32_t geneLossDate, ae_gene_loss_type geneLossType)
// {
//   node->_gene_loss_date = geneLossDate;
//   node->_gene_loss_type = geneLossType;
//   node->_protein_pointer = NULL;
//   for (int32_t i = 0; i < node->_nb_promoters; i++) {node->_rna_pointers[i] = NULL;}
//   if (geneLossDate > _end_gener) _end_gener = geneLossDate;
//   _nb_active_leaves --;
// }
  


ae_gene_tree_node * ae_gene_tree::search_in_leaves(const ae_protein * protein)
{
  return _root->search_in_subtree_leaves(protein);
}
 

void ae_gene_tree::print_to_screen(void)
{
  _root->print_subtree_to_screen();
  printf("\n");
}


void ae_gene_tree::write_to_files(const char * topologyFileName, const char * nodeAttributesFileName, int32_t end_gener)
{
  FILE  * topology_file = fopen(topologyFileName, "w");
  if (topology_file == NULL)
    {
      fprintf(stderr, "Error: cannot create file %s.\n", topologyFileName);
      exit(EXIT_FAILURE);
    }

  FILE  * node_attributes_file = fopen(nodeAttributesFileName, "w");
  if (node_attributes_file == NULL)
    {
      fprintf(stderr, "Error: cannot create file %s.\n", topologyFileName);
      exit(EXIT_FAILURE);
    }  


  _root->write_subtree_to_files(topology_file, node_attributes_file, end_gener);
  fprintf(topology_file, ";");

  fclose(topology_file);
  fclose(node_attributes_file);
}


// f must already be open 
void ae_gene_tree::write_nodes_in_tabular_file(int32_t treeID, FILE * f)
{
  _root->write_subtree_nodes_in_tabular_file(treeID, f);
}


// =================================================================
//                           Protected Methods
// =================================================================

