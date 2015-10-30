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

#ifndef AEVOL_PROTEIN_R_H_
#define AEVOL_PROTEIN_R_H_

// =================================================================
//                              Libraries
// =================================================================

// =================================================================
//                            Project Files
// =================================================================
#include "ae_protein_R.h"
#include "ae_influence_R.h"
#include "ae_codon.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                           Class ae_protein_R                                #
//                                                                             #
//##############################################################################

// =================================================================
//                             Constructors
// =================================================================
ae_protein_R::ae_protein_R(ae_genetic_unit* gen_unit, const ae_protein_R &model) : ae_protein(gen_unit, model)
{
  delta_concentration_  = model.delta_concentration_;
  influence_list_       = new ae_list();
  inherited_            = false;
}

ae_protein_R::ae_protein_R(ae_genetic_unit* gen_unit, ae_list* codon_list, ae_strand strand, int32_t shine_dal_pos,
                            ae_rna* rna)  :
  ae_protein(gen_unit, codon_list, strand, shine_dal_pos, rna)
{
  influence_list_       = new ae_list();
  delta_concentration_  = 0;
  inherited_            = false;
}

/*
ae_protein_R::ae_protein_R(ae_protein_R* parent) :
ae_protein(parent)
{
  //delta_concentration_  = parent->delta_concentration_;
  codon_list_           = parent->codon_list_->copy();
//  rna_list_             = parent->rna_list_->copy();
  rna_list_             = new ae_list();
  //influence_list_       = parent->influence_list_->copy();
  influence_list_       = new ae_list();
}
*/

ae_protein_R::ae_protein_R(gzFile backup_file) : ae_protein::ae_protein(backup_file)
{
  // the Influence list is re-calculate afterward, and then is not saved, nor use in this consctructor.
  gzread(backup_file, &delta_concentration_,   			sizeof(delta_concentration_));
  gzread(backup_file, &inherited_,   			sizeof(inherited_));

  influence_list_       = new ae_list();
}

// =================================================================
//                             Destructors
// =================================================================
ae_protein_R::~ae_protein_R()
{
  remove_influences();
  influence_list_->erase(NO_DELETE);
  delete influence_list_;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_protein_R::compute_delta_concentration()
{
  delta_concentration_ = 0;

  ae_list_node* rna_node  = rna_list_->first();
  ae_rna_R*     rna       = NULL;

  while (rna_node != NULL)
  {
    assert(inherited_ == false);
    rna = (ae_rna_R*)rna_node->obj();

    delta_concentration_ += rna->synthesis_rate();

    rna_node = rna_node->next();
  }

  //printf("degradation rate : %f \n", ae_common::degradation_rate);
  delta_concentration_ -= ae_common::degradation_rate * concentration_;
  delta_concentration_ *= ae_common::degradation_step;
}

int8_t ae_protein_R::codon(int32_t index)
{
  return dynamic_cast< ae_codon* >(AA_list_->object(index))->value();
}

void ae_protein_R::add_influence(ae_influence_R *influence)
{
  influence_list_->add(influence);
}

void ae_protein_R::save(gzFile backup_file)
{
  ae_protein::save(backup_file);

  // the Influence list is re-calculate afterward, and then is not saved.
  gzwrite(backup_file, &delta_concentration_,   	sizeof(delta_concentration_));
  gzwrite(backup_file, &inherited_,   			sizeof(inherited_));
}
// =================================================================
//                           Protected Methods
// =================================================================
void ae_protein_R::remove_influences()
{
  ae_list_node*   influence_node;
  ae_influence_R* influence;
  ae_rna_R*       rna;

  influence_node = influence_list_->first();
  while (influence_node != NULL)
  {
    influence = (ae_influence_R*)influence_node->obj();
    rna       = (ae_rna_R*)influence->rna();
    rna->remove_influence(influence);

    influence_node = influence_node->next();
  }

}

} // namespace aevol

#endif // AEVOL_PROTEIN_R_H_
