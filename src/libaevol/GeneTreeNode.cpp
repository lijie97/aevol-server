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
#include <cassert>
#include <list>
#include <algorithm>
#include <memory>

// =================================================================
//                            Project Files
// =================================================================
#include "GeneTreeNode.h"
#include "GeneMutation.h"
#include "GeneticUnit.h"


 using std::list;


namespace aevol {

//##############################################################################
//                                                                             #
//                          Class GeneTreeNode                            #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================
// =================================================================
//                             Constructors
// =================================================================
int32_t GeneTreeNode::_nextID = 0;

GeneTreeNode::GeneTreeNode( int32_t nodeCreationDate, Protein * protein)
{
  _ID = GeneTreeNode::_nextID;
  GeneTreeNode::_nextID ++;

  _dna_creation_date = nodeCreationDate;
  _node_creation_date = nodeCreationDate;
  _gene_loss_date  = -1;
  _gene_loss_type  = NOT_LOST_YET;

  _protein_pointer = protein;
  _shine_dal_position = protein->get_shine_dal_pos();
  _strand = protein->get_strand();

  _nb_promoters = protein->get_rna_list().size();
  //printf("%d promoters at positions", _nb_promoters); // debug
  _promoter_positions = new int32_t[_nb_promoters];
  _rna_pointers = new Rna *[_nb_promoters];
  int32_t i = 0;
  for (const auto& rna: protein->get_rna_list()) {
    _rna_pointers[i] = rna;
    _promoter_positions[i] = rna->get_promoter_pos();
    i++;
  }
  // printf(" \n "); //debug

  _left_child  = NULL;
  _right_child = NULL;
  _parent_node = NULL;

  _cds_possibly_modified = false;
  _proms_possibly_modified = false;
  _gene_possibly_duplicated = false;
  _cds_completely_deleted = false;
  _putative_position_for_the_duplicate = -1;
}


// =================================================================
//                             Destructors
// =================================================================

GeneTreeNode::~GeneTreeNode( void )
{
  if (_left_child != NULL) delete _left_child;
  if (_right_child != NULL) delete _right_child;
  if (_promoter_positions != NULL) delete [] _promoter_positions;
  if (_rna_pointers != NULL) delete [] _rna_pointers;
  std::for_each(mutation_list.begin(), mutation_list.end(), std::default_delete<GeneMutation>());
}



// =================================================================
//                            Public Methods
// =================================================================

GeneTreeNode *GeneTreeNode::search_in_subtree_leaves(const Protein * protein)
{
  GeneTreeNode *result_left = NULL, *result_right = NULL;
  if ((_left_child == NULL) && (_right_child == NULL)) // I am a leaf
    {
      if (_protein_pointer == protein) return this;
      else return NULL;
    }
  else // I am an internal node
    {
      if (_left_child != NULL)  result_left = _left_child->search_in_subtree_leaves(protein);
      if (_right_child != NULL) result_right = _right_child->search_in_subtree_leaves(protein);
      if ((result_left == NULL) && (result_right == NULL)) return NULL;
      else if ((result_left == NULL) && (result_right != NULL)) return result_right;
      else if ((result_left != NULL) && (result_right == NULL)) return result_left;
      else
        {
          fprintf(stderr, "Error, the protein %p should not be found twice in the tree.\n", protein);
          // get the _root address to print the whole tree to screen
          GeneTreeNode * n = _parent_node, *root = this;
          while (n!=NULL) {root = n; n = n->_parent_node; }
          // here, n==NULL and root points on the root of the tree
          root->print_subtree_to_screen();
          exit(EXIT_FAILURE);
        }
    }
}


void GeneTreeNode::print_subtree_to_screen(void)
{
  // Postorder tree traversal

  /* Left subtree */
  if (_left_child != NULL)  _left_child->print_subtree_to_screen();

  /* Right subtree */
  if (_right_child != NULL) _right_child->print_subtree_to_screen();

  /* Current tree node */
  printf("Node ID: %" PRId32 "\n", _ID);
  if (_parent_node != NULL)   printf("Parent ID: %" PRId32 "\n", _parent_node->_ID);
  else                        printf("Parent ID: none\n");
  if (_left_child != NULL)    printf("Left child ID: %" PRId32 "\n", _left_child->_ID);
  else                        printf("Left child ID: none\n");
  if (_right_child != NULL)   printf("Right child ID: %" PRId32 "\n", _right_child->_ID);
  else                        printf("Right child ID: none\n");
  printf("Node creation date: %" PRId32 "\n", _node_creation_date);
  printf("DNA creation date: %" PRId32 "\n", _dna_creation_date);
  switch(_gene_loss_type)
    {
    case NOT_LOST_YET:            printf("Node status: NOT_LOST_YET\n"); break;
    case LOST_BY_LOCAL_MUTATION:  printf("Gene loss type: LOST_BY_LOCAL_MUTATION\n");  printf("Gene loss date: %" PRId32 "\n", _gene_loss_date); break;
    case DELETED :                printf("Gene loss type: DELETED\n");  printf("Gene loss date: %" PRId32 "\n", _gene_loss_date); break;
    case BROKEN_BY_REAR:          printf("Gene loss type: BROKEN_BY_REAR\n");  printf("Gene loss date: %" PRId32 "\n", _gene_loss_date); break;
    case DUPLICATED:              printf("Node status: DUPLICATED\n");  printf("Duplication date: %" PRId32 "\n", _gene_loss_date);break;
    default: break;
    }
  printf("Protein pointer: %p, Shine-Dalgarno position: %" PRId32 "\n", _protein_pointer, _shine_dal_position);
  if( _strand == LEADING) printf("Strand: LEADING\n");
  else                    printf("Strand: LAGGING\n");
  for(size_t i = 0; i < _nb_promoters; i++)
    {
      printf("Promoter at %" PRId32 ", rna pointer %p\n", _promoter_positions[i], _rna_pointers[i]);
    }
  printf("Number of mutations: %zu\n", mutation_list.size());
  for (const auto& mutation: mutation_list) {
    // TODO vld: simplify
    char str[128];
    mutation->get_description_string_for_gene_mut(str);
    printf("  %s\n", str);
  }
  printf("\n\n");
}


void GeneTreeNode::write_subtree_to_files(FILE * topologyFile, FILE * nodeAttributesFile, int32_t end_gener)
{
  // Newick format for the topology file (postorder tree traversal with parentheses and branch lengths)


  if ((_left_child != NULL) || (_right_child != NULL))
    {
      fprintf(topologyFile, "(");
      /* Left subtree */
      if (_left_child != NULL) _left_child->write_subtree_to_files(topologyFile, nodeAttributesFile, end_gener);
     fprintf(topologyFile, ", ");
     /* Right subtree */
     if (_right_child != NULL) _right_child->write_subtree_to_files(topologyFile, nodeAttributesFile, end_gener);
     fprintf(topologyFile, ")");
    }

  /* Current tree node */
  fprintf(topologyFile, "%" PRId32 "", _ID);
  if (_gene_loss_type == NOT_LOST_YET)  fprintf(topologyFile, ":%" PRId32 "", end_gener - _node_creation_date);
  else fprintf(topologyFile, ":%" PRId32 "", _gene_loss_date - _node_creation_date);


  fprintf(nodeAttributesFile, "Node ID: %" PRId32 "\n", _ID);
  if (_parent_node != NULL)   fprintf(nodeAttributesFile, "Parent ID: %" PRId32 "\n", _parent_node->_ID);
  else                        fprintf(nodeAttributesFile, "Parent ID: none\n");
  if (_left_child != NULL)    fprintf(nodeAttributesFile, "Left child ID: %" PRId32 "\n", _left_child->_ID);
  else                        fprintf(nodeAttributesFile, "Left child ID: none\n");
  if (_right_child != NULL)   fprintf(nodeAttributesFile, "Right child ID: %" PRId32 "\n", _right_child->_ID);
  else                        fprintf(nodeAttributesFile, "Right child ID: none\n");
  fprintf(nodeAttributesFile, "Node creation date: %" PRId32 "\n", _node_creation_date);
  fprintf(nodeAttributesFile, "DNA creation date: %" PRId32 "\n", _dna_creation_date);
  switch(_gene_loss_type)
    {
    case NOT_LOST_YET:            fprintf(nodeAttributesFile, "Node status: NOT_LOST_YET\n"); break;
    case LOST_BY_LOCAL_MUTATION:  fprintf(nodeAttributesFile, "Gene loss type: LOST_BY_LOCAL_MUTATION\n"); fprintf(nodeAttributesFile, "Gene loss date: %" PRId32 "\n", _gene_loss_date); break;
    case DELETED :                fprintf(nodeAttributesFile, "Gene loss type: DELETED\n"); fprintf(nodeAttributesFile, "Gene loss date: %" PRId32 "\n", _gene_loss_date); break;
    case BROKEN_BY_REAR:          fprintf(nodeAttributesFile, "Gene loss type: BROKEN_BY_REAR\n"); fprintf(nodeAttributesFile, "Gene loss date: %" PRId32 "\n", _gene_loss_date); break;
    case DUPLICATED:              fprintf(nodeAttributesFile, "Node status:    DUPLICATED\n"); fprintf(nodeAttributesFile, "Duplication date: %" PRId32 "\n", _gene_loss_date); break;
    default: break;
    }
  if( _strand == LEADING) fprintf(nodeAttributesFile, "Strand: LEADING\n");
  else                    fprintf(nodeAttributesFile, "Strand: LAGGING\n");
  fprintf(nodeAttributesFile, "Shine-Dalgarno position: %" PRId32 "\n", _shine_dal_position);
  for (size_t i= 0; i < _nb_promoters; i++)
    {
      fprintf(nodeAttributesFile, "Position of promoter %" PRId32 ": %" PRId32 "\n", static_cast<int32_t>(i+1), _promoter_positions[i]);
    }
  fprintf(nodeAttributesFile, "Number of mutations: %" PRId32 "\n", static_cast<int32_t>(mutation_list.size()));
  for (const auto& mutation: mutation_list) {
    // TODO vld: simplify
    char str[128];
    mutation->get_description_string_for_gene_mut(str);
    fprintf(nodeAttributesFile, "  %s", str);
    fprintf(nodeAttributesFile, "\n");
  }

  if (_gene_loss_type == NOT_LOST_YET)
    {
      assert( _protein_pointer != NULL );
      fprintf(nodeAttributesFile, "  Shine-Dalgarno pos:%" PRId32 ", Stop pos: %" PRId32 ", M: %.8f, W: %.8f, H: %.8f, nb promoters: %" PRId32 ", conc: %.8f \n", \
              _protein_pointer->get_shine_dal_pos(), _protein_pointer->get_last_STOP_base_pos(), \
              _protein_pointer->get_mean(), _protein_pointer->get_width(), _protein_pointer->get_height(),  \
              static_cast<int32_t>(_protein_pointer->get_rna_list().size()), _protein_pointer->get_concentration() );
    }
  fprintf(nodeAttributesFile, "\n\n");
}



// all attributes on a single line
void GeneTreeNode::write_subtree_nodes_in_tabular_file(int32_t treeID, FILE * f)
{
  /* Left subtree */
  if (_left_child != NULL) _left_child->write_subtree_nodes_in_tabular_file(treeID, f);
  /* Right subtree */
  if (_right_child != NULL)  _right_child->write_subtree_nodes_in_tabular_file(treeID, f);


  fprintf(f, "%" PRId32 " ", treeID);
  fprintf(f, "%" PRId32 " ", _ID);
  if (_parent_node != NULL)   fprintf(f, "%" PRId32 " ", _parent_node->_ID);
  else                        fprintf(f, "-1 ");
  if (_left_child != NULL)    fprintf(f, "%" PRId32 " ", _left_child->_ID);
  else                        fprintf(f, "-1 ");
  if (_right_child != NULL)   fprintf(f, "%" PRId32 " ", _right_child->_ID);
  else                        fprintf(f, "-1 ");
  fprintf(f, "%" PRId32 " ", _node_creation_date);
  fprintf(f, "%" PRId32 " ", _dna_creation_date);
  switch(_gene_loss_type)
    {
    case NOT_LOST_YET:            fprintf(f, "NOT_LOST_YET "); break;
    case LOST_BY_LOCAL_MUTATION:  fprintf(f, "LOST_BY_LOCAL_MUTATION "); break;
    case DELETED :                fprintf(f, "DELETED "); break;
    case BROKEN_BY_REAR:          fprintf(f, "BROKEN_BY_REAR "); break;
    case DUPLICATED:              fprintf(f, "DUPLICATED ");  break;
    default: break;
    }
  fprintf(f, "%" PRId32 " ", _gene_loss_date);
  if (_strand == LEADING) fprintf(f, "LEADING ");
  else                    fprintf(f, "LAGGING ");
  fprintf(f, "%" PRId32 " ", _shine_dal_position);
  fprintf(f, "%" PRId32 " ", static_cast<int32_t>(_nb_promoters));
  if (_gene_loss_type == NOT_LOST_YET)
    {
      assert( _protein_pointer != NULL );
      fprintf(f, "%.8f %.8f %.8f %.8f ", \
              _protein_pointer->get_mean(), _protein_pointer->get_width(), _protein_pointer->get_height(),  \
              _protein_pointer->get_concentration() );
    }
  else {fprintf(f, "-1 -1 -1 -1 ");}

  fprintf(f, "%" PRId32 " ", static_cast<int32_t>(mutation_list.size()));
  int32_t nb_localmut_upstream_neutral = 0, nb_localmut_upstream_benef = 0, nb_localmut_upstream_delet = 0;
  int32_t nb_rear_upstream_neutral = 0, nb_rear_upstream_benef = 0, nb_rear_upstream_delet = 0;
  int32_t nb_localmut_cds_neutral = 0, nb_localmut_cds_benef = 0, nb_localmut_cds_delet = 0;
  int32_t nb_rear_cds_neutral = 0, nb_rear_cds_benef = 0, nb_rear_cds_delet = 0;
  for (const auto& mutation: mutation_list) {
    if (mutation == &(*mutation_list.back()))
      // do not count the last event, if it was disruptive
      if ((_gene_loss_type == DELETED) || (_gene_loss_type == LOST_BY_LOCAL_MUTATION) || (_gene_loss_type == BROKEN_BY_REAR))
        break;

    if (mutation->type_of_event() == 0) {
      if (mutation->get_region() == UPSTREAM) {
        if      (mutation->get_impact_on_metabolic_error() == 0.0) nb_localmut_upstream_neutral++;
        else if (mutation->get_impact_on_metabolic_error() < 0.0)  nb_localmut_upstream_benef++;
        else if (mutation->get_impact_on_metabolic_error() > 0.0)  nb_localmut_upstream_delet++;
      }
      else {
        if      (mutation->get_impact_on_metabolic_error() == 0.0) nb_localmut_cds_neutral++;
        else if (mutation->get_impact_on_metabolic_error() < 0.0)  nb_localmut_cds_benef++;
        else if (mutation->get_impact_on_metabolic_error() > 0.0)  nb_localmut_cds_delet++;
      }
    }
    else {
      if (mutation->get_region() == UPSTREAM) {
        if      (mutation->get_impact_on_metabolic_error() == 0.0) nb_rear_upstream_neutral++;
        else if (mutation->get_impact_on_metabolic_error() < 0.0)  nb_rear_upstream_benef++;
        else if (mutation->get_impact_on_metabolic_error() > 0.0)  nb_rear_upstream_delet++;
      }
      else {
        if      (mutation->get_impact_on_metabolic_error() == 0.0) nb_rear_cds_neutral++;
        else if (mutation->get_impact_on_metabolic_error() < 0.0)  nb_rear_cds_benef++;
        else if (mutation->get_impact_on_metabolic_error() > 0.0)  nb_rear_cds_delet++;
      }
    }
  }
  fprintf(f, "%" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " ", \
          nb_localmut_upstream_neutral, nb_localmut_upstream_benef, nb_localmut_upstream_delet,  \
          nb_localmut_cds_neutral, nb_localmut_cds_benef, nb_localmut_cds_delet, \
          nb_rear_upstream_neutral, nb_rear_upstream_benef, nb_rear_upstream_delet, \
          nb_rear_cds_neutral, nb_rear_cds_benef, nb_rear_cds_delet);



  fprintf(f, "\n");

}



// This is an auxiliary function for the method anticipate_mutation_effect_on_genes_in_subtree.
// The segment should go from 'first' to 'last' (included) in the clockwise sense.
// 'first' and 'last' should not be equal.
static bool breakpoint_inside_segment(int32_t pos_brkpt, int32_t first, int32_t last)
{
  if (first < last) // most frequent case
    {
      if( (first <= pos_brkpt) && (pos_brkpt <= last)) {return true;}
      else {return false;}
    }
  else // special case where the segment overlaps ori
    {
      if( (first <= pos_brkpt) || (pos_brkpt <= last) ) {return true;}
      else {return false;}
    }
}


// This is an auxiliary function for the method anticipate_mutation_effect_on_genes_in_subtree.
// It return true if the subsegment [first, last] is totally included in the segment [pos1, pos2].
// The subsegment should go from 'first' to 'last' in the clockwise sense and the segment should
// go from 'pos1' to 'pos2' in the clockwise sense.
static bool subsegment_totally_in_segment(int32_t pos1, int32_t pos2, int32_t first, int32_t last)
{
  if ( (first < last)  && (pos1 <= pos2) )
    {
      if ( ((first >= pos1) && (first <= pos2)) && ((last >= pos1) && (last <= pos2)) ) {return true; }
      else {return false;}
    }
  else if ( (first < last) && (pos1 > pos2) )  // mut seg in 2 pieces but not the gene
    {
      if ( (first >= pos1) || (last <= pos2) )  // the gene is either completely in [pos1, genlen-1] or completely in [0, pos2]
        {
          return true;
        }
      else return false;
    }
  else if ( (first > last) && (pos1 <= pos2) )  // gene in two pieces but not mut seg, the gene cannot be totally included
    {
      return false;
    }
  else // both mut seg and the gene are in 2 pieces
    {
      if ((first >= pos1) && (last <= pos2)) {return true;}
      else {return false;}
    }
}


void GeneTreeNode::update_pointers_in_subtree_leaves(GeneticUnit * unit)
{
 if ((_left_child != NULL) || (_right_child != NULL)) // I am a internal node
    {
      if (_left_child != NULL)   _left_child->update_pointers_in_subtree_leaves(unit);
      if (_right_child != NULL)  _right_child->update_pointers_in_subtree_leaves(unit);
    }
  else // no child => I am a leaf => there is work to do for me !
  {
    if (_gene_loss_type != NOT_LOST_YET)
      // inactive leaf
      return;

    // TODO vld: refactor DUPLICATED CODE (ref dc1)
    auto& pl = unit->get_protein_list(_strand); // shorthand
    auto protein =
        find_if(pl.begin(), pl.end(),
                [this](Protein & p)
                {return p.get_shine_dal_pos() == _shine_dal_position;});
    if (protein != pl.end()) {
      /* The strand and shine dal position are correct */
      /* Update the protein and rna pointers and positions */
      _nb_promoters = protein->get_rna_list().size();
      if (_promoter_positions != NULL) delete [] _promoter_positions;
      if (_rna_pointers != NULL) delete [] _rna_pointers;
      _promoter_positions = new int32_t[_nb_promoters];
      _rna_pointers = new Rna *[_nb_promoters];
      size_t i = 0;
      for (const auto& rna: protein->get_rna_list()) {
        _rna_pointers[i] = rna;
        _promoter_positions[i] = rna->get_promoter_pos();
        i++;
      }
    }
    else {
      fprintf(stderr, "Error: cannot find a protein that should be there.\n");
      exit(EXIT_FAILURE);
    }
  }
}

void GeneTreeNode::anticipate_mutation_effect_on_genes_in_subtree_leaves(const Mutation * mut, int32_t lengthOfGeneticUnit)
{
  if ((_left_child != NULL) || (_right_child != NULL)) // I am a internal node
    {
      if (_left_child != NULL)   _left_child->anticipate_mutation_effect_on_genes_in_subtree_leaves(mut, lengthOfGeneticUnit);
      if (_right_child != NULL)  _right_child->anticipate_mutation_effect_on_genes_in_subtree_leaves(mut, lengthOfGeneticUnit);
    }
  else // no child => I am a leaf => there is work to do for me !
    {
      if (_gene_loss_type != NOT_LOST_YET)
        {
          // inactive leaf
          return;
        }

      int32_t genlen = lengthOfGeneticUnit; // in bp
      int32_t pos0 = -1, pos1 = -1, pos2 = -1, pos2bis = -1, pos3 = -1, mutlength = -1;
      // int32_t pos1donor = -1, pos2donor = -1, pos3donor = -1;  ae_sense sense = DIRECT;  // related to transfer (TO DO)
      bool invert = false;
      MutationType type = mut->get_mut_type();
      switch(type)
        {
        case SWITCH:
          {
            mut->get_infos_point_mutation(&pos0);
            mutlength = 1;
            break;
          }
        case S_INS:
          {
            mut->get_infos_small_insertion(&pos0, &mutlength);
            break;
          }
        case S_DEL:
          {
            mut->get_infos_small_deletion(&pos0, &mutlength);
            break;
          }
        case DUPL:
          {
            mut->get_infos_duplication(&pos1, &pos2, &pos0);
            // pos2 is actually not included in the segment, the real end of the segment is pos2 - 1
            pos2bis = pos2;
            pos2 = ae_utils::mod(pos2 - 1, genlen);
            mutlength = mut->get_length();
            break;
          }
        case DEL:
          {
            mut->get_infos_deletion(&pos1, &pos2);
            pos2bis = pos2;
            pos2 = ae_utils::mod(pos2 - 1, genlen);
            mutlength = mut->get_length();
            break;
          }
        case TRANS:
          {
            mut->get_infos_translocation(&pos1, &pos2, &pos3, &pos0, &invert); // TO DO: check whether the 4 positions are also in absolute coordinates when rearrangements with alignments
            pos2bis = pos2;
            pos2 = ae_utils::mod(pos2 - 1, genlen);
            mutlength = mut->get_length();
            break;
          }
        case INV:
          {
            mut->get_infos_inversion(&pos1, &pos2);
            pos2bis = pos2;
            pos2 = ae_utils::mod(pos2 - 1, genlen);
            mutlength = mut->get_length();
            break;
          }
        case INSERT:
          {
            // TO DO
            break;
          }
        case INS_HT:
          {
            // TO DO
            break;
          }
        case REPL_HT:
          {
            // TO DO
            break;
          }
        default:
          {
            fprintf(stderr, "Error: unknown mutation type in GeneTreeNode::anticipate_mutation_effect_on_genes_in_subtree.\n");
          }
        }


      int32_t first_cds, last_cds;
      int32_t first_upstream, last_upstream; // "upstream region" is the segment between the furthest promoter and the Shine-Dalgarno sequence
      int32_t nbprom = _protein_pointer->get_rna_list().size();
      assert(nbprom != 0);
      assert(nbprom == static_cast<int32_t>(_nb_promoters));
      int32_t position_furthest_prom = -1, currentprompos = -1;
      if (_protein_pointer->get_strand() == LEADING)
        {
          first_cds = _protein_pointer->get_shine_dal_pos();
          last_cds = _protein_pointer->get_last_STOP_base_pos( );
          for (const auto& rna: _protein_pointer->get_rna_list()) {
              currentprompos = rna->get_promoter_pos();
              if (currentprompos > first_cds) currentprompos = currentprompos - genlen; // negative value for promoters on the other side of ori
              if ((position_furthest_prom == -1) || (position_furthest_prom < currentprompos)) // we need the smallest promoter position
                {
                  position_furthest_prom = rna->get_promoter_pos();
                }
          }
          position_furthest_prom = ae_utils::mod(position_furthest_prom, genlen); // restore a positive value if necessary
          first_upstream = position_furthest_prom;
          last_upstream = ae_utils::mod(first_cds - 1, genlen);
        }
      else
        {
          first_cds = _protein_pointer->get_last_STOP_base_pos( );
          last_cds = _protein_pointer->get_shine_dal_pos();
          for (const auto& rna: _protein_pointer->get_rna_list()) {
              currentprompos = rna->get_promoter_pos();
              if (currentprompos < last_cds) currentprompos = currentprompos + genlen; // value larger than genlen for promoters on the other side of ori
              if ((position_furthest_prom == -1) || (position_furthest_prom > currentprompos)) // we need the largest promoter position
                {
                  position_furthest_prom = rna->get_promoter_pos();
                }
          }
          position_furthest_prom = ae_utils::mod(position_furthest_prom, genlen); // restore a value < genlen if necessary
          first_upstream = ae_utils::mod(last_cds - 1, genlen);
          last_upstream = position_furthest_prom;
        }




      switch(type)
        {
        case SWITCH:
          {
            if (breakpoint_inside_segment(pos0, first_cds, last_cds)) _cds_possibly_modified = true;
            if (breakpoint_inside_segment(pos0, first_upstream, last_upstream)) _proms_possibly_modified = true;
            break;
          }
        case S_INS:
          {
            if (breakpoint_inside_segment(pos0, first_cds, last_cds)) _cds_possibly_modified = true;
            if (breakpoint_inside_segment(pos0, first_upstream, last_upstream)) _proms_possibly_modified = true;

            if (_shine_dal_position >= pos0) _shine_dal_position = ae_utils::mod(_shine_dal_position + mutlength, genlen + mutlength);
            for (size_t i = 0; i < _nb_promoters; i++)
              {
                if (_promoter_positions[i] >= pos0){ _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] + mutlength, genlen + mutlength);}
              }

            break;
          }
        case S_DEL:
          {
            // If the Shine-Dalgarno position is in the small deleted segment, mark the cds as deleted
            // (our gene tracking is based on the tracking of the Shine-Dalgarno position,
            // and we cannot predict the position of a bp that was deleted: we lose track of the gene)
            if (mutlength == 1)
              {
                if (_protein_pointer->get_shine_dal_pos() == pos0)  _cds_completely_deleted = true;
              }
            else // mutlength > 1
              {
                if (breakpoint_inside_segment(_protein_pointer->get_shine_dal_pos(), pos0, ae_utils::mod(pos0 + mutlength - 1, genlen))) _cds_completely_deleted = true;
              }

            if (!(_cds_completely_deleted))
              {
                if (breakpoint_inside_segment(pos0, ae_utils::mod(first_cds - mutlength, genlen), last_cds)) _cds_possibly_modified = true;
                if (breakpoint_inside_segment(pos0, ae_utils::mod(first_upstream - mutlength, genlen), last_upstream)) _proms_possibly_modified = true;

                if ( pos0 + mutlength <= genlen ) // the deletion does not contain the replication origin
                  {
                    if (_shine_dal_position >= pos0) _shine_dal_position = ae_utils::mod(_shine_dal_position - mutlength, genlen - mutlength);
                    for (size_t i = 0; i < _nb_promoters; i++)
                      {
                        if (_promoter_positions[i] >= pos0){ _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] - mutlength, genlen - mutlength);}
                      }
                  }
                else // the deletion contains the replication origin
                  {
                    int32_t nb_del_after_ori = mutlength - genlen + pos0;
                    if (_shine_dal_position >= 0) _shine_dal_position = ae_utils::mod(_shine_dal_position - nb_del_after_ori, genlen - mutlength);
                    for (size_t i = 0; i < _nb_promoters; i++)
                      {
                        if (_promoter_positions[i] >= 0){ _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] - nb_del_after_ori, genlen - mutlength);}
                      }
                  }
              }
            break;
          }
        case DUPL:
          {
            if (subsegment_totally_in_segment(pos1, pos2, first_cds, last_cds))
              {
                _gene_possibly_duplicated = true;
                _putative_position_for_the_duplicate = ae_utils::mod(ae_utils::mod(_shine_dal_position - pos1, genlen) + pos0, genlen + mutlength);
              }
            if (breakpoint_inside_segment(pos0, first_cds, last_cds)) _cds_possibly_modified = true;
            if (breakpoint_inside_segment(pos0, first_upstream, last_upstream)) _proms_possibly_modified = true;

            if (_shine_dal_position >= pos0) _shine_dal_position = ae_utils::mod(_shine_dal_position + mutlength, genlen + mutlength);
            for (size_t i = 0; i < _nb_promoters; i++)
              {
                if (_promoter_positions[i] >= pos0){ _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] + mutlength, genlen + mutlength);}
              }

            break;
          }
        case DEL:
          {
            // If the Shine-Dalgarno is in the deleted segment, mark the cds as deleted
            // (our gene tracking is based on the tracking of the Shine-Dalgarno position,
            // and we cannot predict the position of a bp that was deleted: we lose track of the gene)
           if (mutlength == 1)
              {
                if (_protein_pointer->get_shine_dal_pos() == pos1)  _cds_completely_deleted = true;
              }
            else // mutlength > 1
              {
                if (breakpoint_inside_segment(_protein_pointer->get_shine_dal_pos(), pos1, pos2)) _cds_completely_deleted = true;
              }

            if (!(_cds_completely_deleted))
              {
                if (breakpoint_inside_segment(pos1, first_cds, last_cds) || breakpoint_inside_segment(pos2, first_cds, last_cds)) _cds_possibly_modified = true;
                if (breakpoint_inside_segment(pos1, first_upstream, last_upstream) || breakpoint_inside_segment(pos2, first_upstream, last_upstream)) _proms_possibly_modified = true;

                if ( pos1 < pos2bis ) // the deletion does not contain the replication origin
                  {
                    if (_shine_dal_position >= pos1) _shine_dal_position = ae_utils::mod(_shine_dal_position - mutlength, genlen - mutlength);
                    for (size_t i = 0; i < _nb_promoters; i++)
                      {
                        if (_promoter_positions[i] >= pos1){ _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] - mutlength, genlen - mutlength);}
                      }
                  }
                else  // the deletion contains the replication origin
                  {
                    if (_shine_dal_position >= 0) _shine_dal_position = ae_utils::mod(_shine_dal_position - pos2bis, genlen - mutlength);
                    for (size_t i = 0; i < _nb_promoters; i++)
                      {
                        if (_promoter_positions[i] >= 0){ _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] - pos2bis, genlen - mutlength);}
                      }
                  }
              }
            break;
          }
        case TRANS:
          {
            if (breakpoint_inside_segment(pos1, first_cds, last_cds)) _cds_possibly_modified = true;   // beginning of the excised segment
            if (breakpoint_inside_segment(pos2, first_cds, last_cds)) _cds_possibly_modified = true;  // end of the excised segment
            if (breakpoint_inside_segment(pos3, first_cds, last_cds)) _cds_possibly_modified = true;  // breakpoint inside the segment for the reinsertion
            if (breakpoint_inside_segment(pos0, first_cds, last_cds)) _cds_possibly_modified = true;  // reinsertion point in the genetic unit
            if (breakpoint_inside_segment(pos1, first_upstream, last_upstream)) _proms_possibly_modified = true;   // beginning of the excised segment
            if (breakpoint_inside_segment(pos2, first_upstream, last_upstream)) _proms_possibly_modified = true;  // end of the excised segment
            if (breakpoint_inside_segment(pos3, first_upstream, last_upstream)) _proms_possibly_modified = true;  // breakpoint inside the segment for the reinsertion
            if (breakpoint_inside_segment(pos0, first_upstream, last_upstream)) _proms_possibly_modified = true;  // reinsertion point in the genetic unit

            int32_t pos_min = ae_utils::min( pos1, ae_utils::min( pos2bis, ae_utils::min( pos3, pos0 ) ) );
            int32_t pos_B, pos_C, pos_D, pos_E;
            int32_t len_B, len_C, len_D;
            if ( ! invert )
              {
                if ( pos_min == pos1 )          { pos_B = pos1; pos_C = pos3; pos_D = pos2bis; pos_E = pos0; }
                else if ( pos_min == pos2bis )  { pos_B = pos2bis; pos_C = pos0; pos_D = pos1; pos_E = pos3; }
                else if ( pos_min == pos3 )     { pos_B = pos3; pos_C = pos2bis; pos_D = pos0; pos_E = pos1; }
                else                            { pos_B = pos0; pos_C = pos1; pos_D = pos3; pos_E = pos2bis; } // if ( pos_min == pos0 )
                len_B = pos_C - pos_B;
                len_C = pos_D - pos_C;
                len_D = pos_E - pos_D;
                if      ((_shine_dal_position >= pos_B) && (_shine_dal_position < pos_C))   _shine_dal_position = ae_utils::mod(_shine_dal_position + len_D + len_C, genlen);
                else if ((_shine_dal_position >= pos_C) && (_shine_dal_position < pos_D))   _shine_dal_position = ae_utils::mod(_shine_dal_position + len_D - len_B, genlen);
                else if ((_shine_dal_position >= pos_D) && (_shine_dal_position < pos_E))   _shine_dal_position = ae_utils::mod(_shine_dal_position - len_B - len_C, genlen);
                for (size_t i = 0; i < _nb_promoters; i++)
                  {
                    if      ((_promoter_positions[i] >= pos_B) && (_promoter_positions[i] < pos_C))  _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] + len_D + len_C, genlen);
                    else if ((_promoter_positions[i] >= pos_C) && (_promoter_positions[i] < pos_D))  _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] + len_D - len_B, genlen);
                    else if ((_promoter_positions[i] >= pos_D) && (_promoter_positions[i] < pos_E))  _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] - len_B - len_C, genlen);
                  }
              }
            else // invert
              {
                if ( pos_min == pos1 )
                  {
                    pos_B = pos1; pos_C = pos3; pos_D = pos2bis; pos_E = pos0;
                    len_B = pos_C - pos_B;
                    len_C = pos_D - pos_C;
                    len_D = pos_E - pos_D;
                    if  ((_shine_dal_position >= pos_B) && (_shine_dal_position < pos_C))
                      {
                        if (_strand == LEADING) {_strand = LAGGING;} else {_strand = LEADING;}
                        _shine_dal_position = pos_B + pos_C - _shine_dal_position -1 ;
                        _shine_dal_position = ae_utils::mod(_shine_dal_position + len_D, genlen);
                      }
                    else if ((_shine_dal_position >= pos_C) && (_shine_dal_position < pos_D))
                      {
                        if (_strand == LEADING) {_strand = LAGGING;} else {_strand = LEADING;}
                        _shine_dal_position = pos_C + pos_D - _shine_dal_position - 1;
                        _shine_dal_position = ae_utils::mod(_shine_dal_position + len_D, genlen);
                      }
                    else if ((_shine_dal_position >= pos_D) && (_shine_dal_position < pos_E))
                      {
                        _shine_dal_position = ae_utils::mod(_shine_dal_position - len_B - len_C, genlen);
                      }
                    for (size_t i = 0; i < _nb_promoters; i++)
                      {
                        if  ((_promoter_positions[i] >= pos_B) && (_promoter_positions[i] < pos_C))
                          {
                            _promoter_positions[i] = pos_B + pos_C - _promoter_positions[i] - 1;
                            _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] + len_D, genlen);
                          }
                        else if ((_promoter_positions[i] >= pos_C) && (_promoter_positions[i] < pos_D))
                          {
                            _promoter_positions[i] = pos_C + pos_D - _promoter_positions[i] - 1;
                            _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] + len_D, genlen);
                          }
                        else if ((_promoter_positions[i] >= pos_D) && (_promoter_positions[i] < pos_E))
                          {
                            _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] - len_B - len_C, genlen);
                          }
                      }
                  }
                else if ( pos_min == pos2bis )
                  {
                    pos_B = pos2bis; pos_C = pos0; pos_D = pos1; pos_E = pos3;
                    len_B = pos_C - pos_B;
                    len_C = pos_D - pos_C;
                    len_D = pos_E - pos_D;
                    if  ((_shine_dal_position >= pos_B) && (_shine_dal_position < pos_C))
                      {
                        if (_strand == LEADING) {_strand = LAGGING;} else {_strand = LEADING;}
                        _shine_dal_position = pos_B + pos_C - _shine_dal_position - 1;
                        _shine_dal_position = ae_utils::mod(_shine_dal_position + len_D, genlen);
                      }
                    else if ((_shine_dal_position >= pos_C) && (_shine_dal_position < pos_D))
                      {
                        if (_strand == LEADING) {_strand = LAGGING;} else {_strand = LEADING;}
                        _shine_dal_position = pos_C + pos_D - _shine_dal_position - 1;
                        _shine_dal_position = ae_utils::mod(_shine_dal_position + len_D, genlen);
                      }
                    else if ((_shine_dal_position >= pos_D) && (_shine_dal_position < pos_E))
                      {
                        _shine_dal_position = ae_utils::mod(_shine_dal_position - len_B - len_C, genlen);
                      }
                    for (size_t i = 0; i < _nb_promoters; i++)
                      {
                        if  ((_promoter_positions[i] >= pos_B) && (_promoter_positions[i] < pos_C))
                          {
                            _promoter_positions[i] = pos_B + pos_C - _promoter_positions[i] - 1;
                            _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] + len_D, genlen);
                          }
                        else if ((_promoter_positions[i] >= pos_C) && (_promoter_positions[i] < pos_D))
                          {
                            _promoter_positions[i] = pos_C + pos_D - _promoter_positions[i] - 1;
                            _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] + len_D, genlen);
                          }
                        else if ((_promoter_positions[i] >= pos_D) && (_promoter_positions[i] < pos_E))
                          {
                            _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] - len_B - len_C, genlen);
                          }
                      }
                  }
                else if ( pos_min == pos3 )
                  {
                    pos_B = pos3; pos_C = pos2bis; pos_D = pos0; pos_E = pos1;
                    len_B = pos_C - pos_B;
                    len_C = pos_D - pos_C;
                    len_D = pos_E - pos_D;
                    if  ((_shine_dal_position >= pos_B) && (_shine_dal_position < pos_C))
                      {
                        _shine_dal_position = ae_utils::mod(_shine_dal_position + len_C + len_D, genlen);
                      }
                    else if ((_shine_dal_position >= pos_C) && (_shine_dal_position < pos_D))
                      {
                        if (_strand == LEADING) {_strand = LAGGING;} else {_strand = LEADING;}
                        _shine_dal_position = pos_C + pos_D - _shine_dal_position - 1;
                        _shine_dal_position = ae_utils::mod(_shine_dal_position - len_B, genlen);
                      }
                    else if ((_shine_dal_position >= pos_D) && (_shine_dal_position < pos_E))
                      {
                        if (_strand == LEADING) {_strand = LAGGING;} else {_strand = LEADING;}
                        _shine_dal_position = pos_D + pos_E - _shine_dal_position - 1;
                        _shine_dal_position = ae_utils::mod(_shine_dal_position - len_B, genlen);
                      }
                    for (size_t i = 0; i < _nb_promoters; i++)
                      {
                        if  ((_promoter_positions[i] >= pos_B) && (_promoter_positions[i] < pos_C))
                          {
                            _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] + len_C + len_D, genlen);
                          }
                        else if ((_promoter_positions[i] >= pos_C) && (_promoter_positions[i] < pos_D))
                          {
                            _promoter_positions[i] = pos_C + pos_D - _promoter_positions[i] - 1;
                            _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] - len_B, genlen);
                          }
                        else if ((_promoter_positions[i] >= pos_D) && (_promoter_positions[i] < pos_E))
                          {
                            _promoter_positions[i] = pos_D + pos_E - _promoter_positions[i] - 1;
                            _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] - len_B, genlen);
                          }
                      }

                  }
                else // if ( pos_min == pos0 )
                  {
                    pos_B = pos0; pos_C = pos1; pos_D = pos3; pos_E = pos2bis;
                    len_B = pos_C - pos_B;
                    len_C = pos_D - pos_C;
                    len_D = pos_E - pos_D;
                    if  ((_shine_dal_position >= pos_B) && (_shine_dal_position < pos_C))
                      {
                        _shine_dal_position = ae_utils::mod(_shine_dal_position + len_C + len_D, genlen);
                      }
                    else if ((_shine_dal_position >= pos_C) && (_shine_dal_position < pos_D))
                      {
                        if (_strand == LEADING) {_strand = LAGGING;} else {_strand = LEADING;}
                        _shine_dal_position = pos_C + pos_D - _shine_dal_position - 1;
                        _shine_dal_position = ae_utils::mod(_shine_dal_position - len_B, genlen);
                      }
                    else if ((_shine_dal_position >= pos_D) && (_shine_dal_position < pos_E))
                      {
                        if (_strand == LEADING) {_strand = LAGGING;} else {_strand = LEADING;}
                        _shine_dal_position = pos_D + pos_E - _shine_dal_position - 1;
                        _shine_dal_position = ae_utils::mod(_shine_dal_position - len_B, genlen);
                      }
                    for (size_t i = 0; i < _nb_promoters; i++)
                      {
                        if  ((_promoter_positions[i] >= pos_B) && (_promoter_positions[i] < pos_C))
                          {
                            _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] + len_C + len_D, genlen);
                          }
                        else if ((_promoter_positions[i] >= pos_C) && (_promoter_positions[i] < pos_D))
                          {
                            _promoter_positions[i] = pos_C + pos_D - _promoter_positions[i] - 1;
                            _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] - len_B, genlen);
                          }
                        else if ((_promoter_positions[i] >= pos_D) && (_promoter_positions[i] < pos_E))
                          {
                            _promoter_positions[i] = pos_D + pos_E - _promoter_positions[i] - 1;
                            _promoter_positions[i] = ae_utils::mod(_promoter_positions[i] - len_B, genlen);
                          }
                      }
                  }
              }

            break;
          }
        case INV:
          {
            if (breakpoint_inside_segment(pos1, first_cds, last_cds)) _cds_possibly_modified = true;
            if (breakpoint_inside_segment(pos2, first_cds, last_cds)) _cds_possibly_modified = true;
            if (breakpoint_inside_segment(pos1, first_upstream, last_upstream)) _proms_possibly_modified = true;   // beginning of the excised segment
            if (breakpoint_inside_segment(pos2, first_upstream, last_upstream)) _proms_possibly_modified = true;  // end of the excised segment

            if  ((_shine_dal_position >= pos1) && (_shine_dal_position < pos2bis))
              {
                if (_strand == LEADING) {_strand = LAGGING;} else {_strand = LEADING;}
                _shine_dal_position = ae_utils::mod(pos1 + pos2bis - _shine_dal_position - 1, genlen);
              }
            break;
          }
        case INSERT:
          {
            // TO DO
            break;
          }
        case INS_HT:
          {
            // TO DO
            break;
          }
        case REPL_HT:
          {
            // TO DO
            break;
          }
        default:
          // Only simple mutation types are considered.
          break;
        }
    }
}




void GeneTreeNode::register_actual_mutation_effect_on_genes_in_subtree_leaves(
    GeneTree * tree, const Mutation * mut, GeneticUnit* unit, int32_t gener, double impact_on_metabolic_error)
{
  if ((_left_child != NULL) || (_right_child != NULL)) // I am a internal node, just delegate work to others
    {
      if (_left_child != NULL)   _left_child->register_actual_mutation_effect_on_genes_in_subtree_leaves( tree, mut, unit, gener, impact_on_metabolic_error);
      if (_right_child != NULL)  _right_child->register_actual_mutation_effect_on_genes_in_subtree_leaves( tree, mut, unit, gener, impact_on_metabolic_error);
    }
  else // no child => I am a leaf => there is work to do for me !
    {

      if (_gene_loss_type != NOT_LOST_YET)
        {
          // inactive leaf
          return;
        }

      GeneMutation * genemut = NULL;

      if (_cds_completely_deleted)
        {
          genemut = new GeneMutation(*mut, gener, _shine_dal_position, _strand, CDS);
          genemut->set_impact_on_metabolic_error(impact_on_metabolic_error);
          mutation_list.push_back(genemut);
          if (gener > tree->_end_gener) tree->_end_gener = gener;

          _gene_loss_date = gener;
          _gene_loss_type = DELETED;
          _protein_pointer = NULL;
          for (size_t i = 0; i < _nb_promoters; i++) {_rna_pointers[i] = NULL;}
          if (gener > tree->_end_gener) tree->_end_gener = gener;
          (tree->_nb_active_leaves) --;
          return;
        }


      if ( (!_cds_completely_deleted) && (!_cds_possibly_modified) && (!_proms_possibly_modified) )
        {
          // This CDS was not affected by the mutation (it could have be moved or duplicated however).
          // Just make sure that we have correctly predicted the positions of the SD sequence and of the promoters.

          // TODO vld: refactor DUPLICATED CODE (ref dc1)
          auto& pl = unit->get_protein_list(_strand);
          auto protein =
              find_if(pl.begin(), pl.end(),
                      [this](Protein & p)
                      { return p.get_shine_dal_pos() == _shine_dal_position; });
          if (protein != pl.end()) {
            /* The strand and shine dal position are correct */
            /* Update the protein and rna pointers and positions */
            _nb_promoters = protein->get_rna_list().size();
            if (_promoter_positions != NULL) delete [] _promoter_positions;
            if (_rna_pointers != NULL) delete [] _rna_pointers;
            _promoter_positions = new int32_t[_nb_promoters];
            _rna_pointers = new Rna *[_nb_promoters];
            size_t i = 0;
            for (const auto& rna: protein->get_rna_list()) {
              _rna_pointers[i] = rna;
              _promoter_positions[i] = rna->get_promoter_pos();
              i++;
            }
          }
          else
            {
              fprintf(stderr, "Error: cannot find a protein that should have survived.\n");
              char str[100];
              mut->get_generic_description_string(str);
              printf("Mutation : %s\n", str);
              printf("CDS should be at %d ", _shine_dal_position);
              if (_strand == LEADING) printf("LEADING\n");
              else printf("LAGGING\n");
              unit->print_proteins();
              exit(EXIT_FAILURE);
            }
        }


      if (_cds_possibly_modified || _proms_possibly_modified)
        {
          /* Record the impact of the mutation on the metabolic error */
          if      ((_cds_possibly_modified) && (_proms_possibly_modified))   genemut = new GeneMutation(*mut, gener, _shine_dal_position, _strand, BOTH);
          else if ((_cds_possibly_modified) && (!_proms_possibly_modified))  genemut = new GeneMutation(*mut, gener, _shine_dal_position, _strand, CDS);
          else if ((!_cds_possibly_modified) && (_proms_possibly_modified))  genemut = new GeneMutation(*mut, gener, _shine_dal_position, _strand, UPSTREAM);
          genemut->set_impact_on_metabolic_error(impact_on_metabolic_error);
          mutation_list.push_back(genemut);
          if (gener > tree->_end_gener) tree->_end_gener = gener;

          /* Check whether the protein survived the event */
          // TODO vld: refactor DUPLICATED CODE (ref dc1)
          auto& pl = unit->get_protein_list(_strand); // shorthand
          auto protein =
              find_if(pl.begin(), pl.end(),
                      [this](Protein & p)
                      { return p.get_shine_dal_pos() == _shine_dal_position; });
          if (protein != pl.end()) {
            /* The strand and shine dal position are correct */
            /* Update the protein and rna pointers and positions */
            _nb_promoters = protein->get_rna_list().size();
            if (_promoter_positions != NULL) delete [] _promoter_positions;
            if (_rna_pointers != NULL) delete [] _rna_pointers;
            _promoter_positions = new int32_t[_nb_promoters];
            _rna_pointers = new Rna *[_nb_promoters];
            size_t i = 0;
            for (const auto& rna: protein->get_rna_list()) {
              _rna_pointers[i] = rna;
              _promoter_positions[i] = rna->get_promoter_pos();
              i++;
            }
          }
          else
            {
              /* The protein does not exist anymore, the gene was killed by the event */
              _gene_loss_date = gener;
              if ((mut->get_mut_type() == SWITCH) || (mut->get_mut_type() == S_INS) || (mut->get_mut_type() == S_DEL)) _gene_loss_type = LOST_BY_LOCAL_MUTATION;
              else if ((mut->get_mut_type() == DUPL) || (mut->get_mut_type() == DEL) || (mut->get_mut_type() == TRANS) || (mut->get_mut_type() == INV)) _gene_loss_type = BROKEN_BY_REAR;
              _protein_pointer = NULL;
              for (size_t i = 0; i < _nb_promoters; i++) {_rna_pointers[i] = NULL;}
              if (gener > tree->_end_gener) (tree->_end_gener) = gener;
              (tree->_nb_active_leaves) --;

            }

        }


      if (_gene_possibly_duplicated) {
        /* Check whether the duplicated CDS found a promoter */
        /* It should be on the same strand as myself, at the _putative_position_for_the_duplicate */

        auto& pl = unit->get_protein_list(_strand); // shorthand
        auto protein =
            find_if(pl.begin(), pl.end(),
                    [this](Protein & p)
                    { return p.get_shine_dal_pos() == _putative_position_for_the_duplicate; });

        if (protein != pl.end()) {
          if (_protein_pointer != NULL) {
            // Create a new node for the "old" DNA segment
            _left_child = new GeneTreeNode(gener, _protein_pointer);
            _left_child->_node_creation_date = gener;
            _left_child->_dna_creation_date = _dna_creation_date;
            _left_child->_parent_node = this;
          }
          else {
            // Poor old gene was killed by the insertion of the duplicated segment
            // We do not create a left child
            _left_child = NULL;
          }

          // Create a new node for the "new" DNA segment
          _right_child = new GeneTreeNode(gener, &*protein);
          _right_child->_node_creation_date = gener;
          _right_child->_dna_creation_date = gener;
          _right_child->_parent_node = this;

          // This node becomes internal, it represents an ancestral (obsolete) state of the gene
          _protein_pointer = NULL;
          for (size_t i = 0; i < _nb_promoters; i++)
            _rna_pointers[i] = NULL;
          _gene_loss_type = DUPLICATED;
          _gene_loss_date = gener;

          // Update tree statistics
          if (_left_child != NULL) {
            (tree->_total_nb_nodes) += 2;
            (tree->_nb_internal_nodes) ++;
            (tree->_nb_leaves) ++;  // - 1 + 2 (the ex-leaf becomes an internal node, 2 leaves are created)

          }
          else {
            (tree->_total_nb_nodes) += 1;
            (tree->_nb_internal_nodes) ++;
            // (tree->_nb_leaves) remains unchanged  <==  - 1 + 1 (the ex-leaf becomes an internal node, 1 leave is created)
          }
          if (protein != pl.end())
            tree->_nb_active_leaves++;
          if (gener > tree->_end_gener)
            tree->_end_gener = gener;
        }
        // else nothing to do, the duplication was only partial, not a complete gene duplication
      }


      /* Get ready for the next mutation */
      _cds_possibly_modified = false;
      _proms_possibly_modified = false;
      _gene_possibly_duplicated = false;
      _cds_completely_deleted = false;
      _putative_position_for_the_duplicate = -1;
    }

}



// =================================================================
//                           Protected Methods
// =================================================================


} // namespace aevol