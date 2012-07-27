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

// =================================================================
//                            Project Files
// =================================================================
#include <ae_protein_R.h>
#include <ae_influence_R.h>
#include <ae_codon.h>

//##############################################################################
//                                                                             #
//                           Class ae_protein_R                                #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_protein_R::ae_protein_R( ae_genetic_unit* gen_unit, const ae_protein_R &model ) : ae_protein( gen_unit, model )
{
  _delta_concentration  = model._delta_concentration;
  _influence_list       = new ae_list();
  _inherited            = false;
}

ae_protein_R::ae_protein_R( ae_genetic_unit* gen_unit, ae_list* codon_list, ae_strand strand, int32_t shine_dal_pos,
                            ae_rna* rna )  :
  ae_protein( gen_unit, codon_list, strand, shine_dal_pos, rna )
{
  _influence_list       = new ae_list();
  _delta_concentration  = 0;
  _inherited            = false;
}

/*
ae_protein_R::ae_protein_R( ae_protein_R* parent ) :
ae_protein( parent )
{
  //_delta_concentration  = parent->_delta_concentration;
  _codon_list           = parent->_codon_list->copy();
//  _rna_list             = parent->_rna_list->copy();
  _rna_list             = new ae_list();
  //_influence_list       = parent->_influence_list->copy();
  _influence_list       = new ae_list();
}
*/

ae_protein_R::ae_protein_R( gzFile* backup_file ) : ae_protein::ae_protein( backup_file )
{
  // the Influence list is re-calculate afterward, and then is not saved, nor use in this consctructor.
  gzread( backup_file, &_delta_concentration,   			sizeof(_delta_concentration) );
  gzread( backup_file, &_inherited,   			sizeof(_inherited) );

  _influence_list       = new ae_list();
}
 
// =================================================================
//                             Destructors
// =================================================================
ae_protein_R::~ae_protein_R( void )
{
  remove_influences();
  _influence_list->erase( NO_DELETE );
  delete _influence_list;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_protein_R::compute_delta_concentration( void )
{
  _delta_concentration = 0;

  ae_list_node* rna_node  = _rna_list->get_first();
  ae_rna_R*     rna       = NULL;

  while ( rna_node != NULL )
  {
    assert( _inherited == false );
    rna = (ae_rna_R*)rna_node->get_obj();

    _delta_concentration += rna->get_synthesis_rate();

    rna_node = rna_node->get_next();
  }

  //printf("degradation rate : %f \n", ae_common::degradation_rate);
  _delta_concentration -= ae_common::degradation_rate * _concentration;
  _delta_concentration *= ae_common::degradation_step;
}

int8_t ae_protein_R::get_codon( int32_t index )
{
  return dynamic_cast< ae_codon* >( _AA_list->get_object( index ) )->get_value();
}

void ae_protein_R::add_influence( ae_influence_R *influence )
{
  _influence_list->add( influence );
}

void ae_protein_R::save( gzFile* backup_file ) 
{
  ae_protein::save( backup_file );

  // the Influence list is re-calculate afterward, and then is not saved.
  gzwrite( backup_file, &_delta_concentration,   	sizeof(_delta_concentration) );
  gzwrite( backup_file, &_inherited,   			sizeof(_inherited) );
}
// =================================================================
//                           Protected Methods
// =================================================================
void ae_protein_R::remove_influences( void )
{
  ae_list_node*   influence_node;
  ae_influence_R* influence;
  ae_rna_R*       rna;

  influence_node = _influence_list->get_first();
  while ( influence_node != NULL )
  {
    influence = (ae_influence_R*)influence_node->get_obj();
    rna       = (ae_rna_R*)influence->get_rna();
    rna->remove_influence( influence );

    influence_node = influence_node->get_next();
  }

}
