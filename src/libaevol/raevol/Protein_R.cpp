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




// =================================================================
//                              Libraries
// =================================================================

// =================================================================
//                            Project Files
// =================================================================
#include "Protein_R.h"
#include "Codon.h"
#include <algorithm>

namespace aevol {

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

Protein_R::Protein_R( GeneticUnit* gen_unit, const Protein_R &model ) : Protein::Protein( gen_unit, model )
{
  _concentration         = model._concentration;
  _initial_concentration = model._concentration;
  _delta_concentration   = model._delta_concentration;
  _signal                = model._signal;
  _inherited             = model._inherited;
//  _influence_list        = new ae_list();
  not_pure_TF			 = false;
}


Protein_R::Protein_R( GeneticUnit* gen_unit, const std::list<Codon*> codon_list,
							Strand strand, int32_t shine_dal_pos,
                            ae_rna* rna )  :
		Protein::Protein( gen_unit, codon_list, strand, shine_dal_pos, rna )
{
	_initial_concentration = 0;
//  _influence_list       = new ae_list();
  _delta_concentration  = 0;
  _inherited            = false;
  _signal               = false;
  not_pure_TF			 = false;
}

//used to build the signal protein
Protein_R::Protein_R( const std::list<ae_codon*> codon_list, double concentration)  :
		Protein::Protein( codon_list, concentration )
{
//  _influence_list       = new ae_list();
  _initial_concentration = 0;
  _delta_concentration  = 0;
  _inherited            = false;
  _signal               = true;
  not_pure_TF			 = false;
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

Protein_R::Protein_R( gzFile backup_file ) : Protein::Protein( backup_file )
{
	_initial_concentration = 0;
  // the Influence list is re-calculate afterward, and then is not saved, nor use in this consctructor.
  gzread( backup_file, &_delta_concentration,   	sizeof(_delta_concentration) );
  gzread( backup_file, &_inherited,   			sizeof(_inherited) );
  gzread( backup_file, &_signal,   			sizeof(_signal) );
  not_pure_TF			 = false;

//  _influence_list       = new ae_list();
}
 
// =================================================================
//                             Destructors
// =================================================================
Protein_R::~Protein_R( void )
{
	_rna_R_list.clear();
}

// =================================================================
//                            Public Methods
// =================================================================
void Protein_R::compute_delta_concentration( void )
{
  _delta_concentration = 0;
  if( _signal == 0 )
  {
	for (const auto& rna: rna_list)
    {
      assert( _inherited == false);

      _delta_concentration += ((ae_rna_R)rna).get_synthesis_rate();
    }
    _delta_concentration -= _gen_unit->get_exp_m()->get_exp_s()->get_degradation_rate() * _concentration;
    _delta_concentration *= _gen_unit->get_exp_m()->get_exp_s()->get_degradation_step();
  }
}

int8_t Protein_R::get_codon( int32_t index )
{
  return _AA_list[index];
}

void Protein_R::save( gzFile backup_file )
{
  Protein::save( backup_file );

  // the Influence list is re-calculate afterward, and then is not saved.
  gzwrite( backup_file, &_delta_concentration,   	sizeof(_delta_concentration) );
  gzwrite( backup_file, &_inherited,   			sizeof(_inherited) );
  gzwrite( backup_file, &_signal,   			sizeof(_signal) );
}
// =================================================================
//                           Protected Methods
// =================================================================
void Protein_R::remove_influences( void )
{
  printf("ALERTE la proteine veut détruire une influence !!!\n");

  _rna_R_list.clear();
}

} // namespace aevol
