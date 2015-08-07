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
#include "GeneticUnit.h"
#include "ExpManager.h"

#include <algorithm>

namespace aevol {

//##############################################################################
//                                                                             #
//                           Class Protein_R                                #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================

Protein_R::Protein_R( GeneticUnit* gen_unit, const Protein_R &model ) : Protein::Protein( gen_unit, model ) {
  _concentration = model._concentration;
  _initial_concentration = model._concentration;
  _delta_concentration = model._delta_concentration;
  _signal = model._signal;
  _inherited = model._inherited;
  not_pure_TF = false;

  if (!_AA_list.empty()) {
    _cod_tab = new int8_t[_AA_list.size()];
    int i = 0;
    for (auto cod : _AA_list) {
      _cod_tab[i] = cod->get_value();
      i++;
    }
  }
}


Protein_R::Protein_R( GeneticUnit* gen_unit, const std::vector<Codon*> codon_list,
							Strand strand, int32_t shine_dal_pos,
                            Rna* rna, double w_max )  :
		Protein::Protein( gen_unit, codon_list, strand, shine_dal_pos, rna, w_max )
{
  _rna_R_list.push_back((Rna_R*)rna);

	_initial_concentration = 0;
  _delta_concentration   = 0;
  _inherited             = false;
  _signal                = false;
  not_pure_TF			       = false;

  if (!_AA_list.empty()) {
    _cod_tab = new int8_t[_AA_list.size()];
    int i = 0;
    for (auto cod : _AA_list) {
      _cod_tab[i] = cod->get_value();
      i++;
    }
  }
}

//used to build the signal protein
Protein_R::Protein_R( const std::vector<Codon*> codon_list, double concentration)  :
		Protein::Protein( codon_list, concentration )
{
  _initial_concentration = 0;
  _delta_concentration  = 0;
  _inherited            = false;
  _signal               = true;
  not_pure_TF			 = false;

  if (!_AA_list.empty()) {
    _cod_tab = new int8_t[_AA_list.size()];
    int i = 0;
    for (auto cod : _AA_list) {
      _cod_tab[i] = cod->get_value();
      i++;
    }
  }
}

Protein_R::Protein_R( gzFile backup_file ) : Protein::Protein( backup_file )
{
  // the Influence list is re-calculate afterward, and then is not saved, nor use in this consctructor.
  gzread( backup_file, &_delta_concentration,   	sizeof(_delta_concentration) );
  gzread( backup_file, &_inherited,   			sizeof(_inherited) );
  gzread( backup_file, &_signal,   			sizeof(_signal) );
  not_pure_TF			 = false;

  if (!_AA_list.empty()) {
    _cod_tab = new int8_t[_AA_list.size()];
    int i = 0;
    for (auto cod : _AA_list) {
      _cod_tab[i] = cod->get_value();
      i++;
    }
  }
}
 
// =================================================================
//                             Destructors
// =================================================================
Protein_R::~Protein_R( void )
{
	_rna_R_list.clear();

  delete [] _cod_tab;
}

// =================================================================
//                            Public Methods
// =================================================================
void Protein_R::compute_delta_concentration( void )
{
  _delta_concentration = 0;

  if( _signal == false )
  {
	for (auto& rna: _rna_R_list)
    {
      assert( _inherited == false);

      _delta_concentration += rna->get_synthesis_rate();

    }

    _delta_concentration -= _gen_unit->get_exp_m()->get_exp_s()->get_degradation_rate() * _concentration;
    _delta_concentration *= 1/_gen_unit->get_exp_m()->get_exp_s()->get_degradation_step();

  }
}

int8_t Protein_R::get_codon( int32_t index )
{
  return _cod_tab[index];
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


void Protein_R::add_RNA( Rna * rna )
{
  Protein::add_RNA(rna);
  _rna_R_list.push_back((Rna_R*)rna);
}

} // namespace aevol
