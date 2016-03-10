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
long Protein_R::id = 0;

// =================================================================
//                             Constructors
// =================================================================
Protein_R::Protein_R( GeneticUnit* gen_unit, const Protein_R &model ) : Protein::Protein( gen_unit, model ) {
  concentration_ = model.concentration_;
  _initial_concentration = model.concentration_;
  _delta_concentration = model._delta_concentration;
  _signal = model._signal;
  _inherited = model._inherited;
  is_TF_ = model.is_TF_;
  _id = id++;

  if (!AA_list_.empty()) {
    _cod_tab = new int8_t[AA_list_.size()];
    int i = 0;
    for (auto cod : AA_list_) {
      _cod_tab[i] = cod->value();
      i++;
    }
  }
}


Protein_R::Protein_R( GeneticUnit* gen_unit, const std::list<Codon*> codon_list,
							Strand strand, int32_t shine_dal_pos,
                            Rna* rna, double w_max )  :
		Protein::Protein( gen_unit, codon_list, strand, shine_dal_pos, rna, w_max )
{
  _rna_R_list.push_back((Rna_R*)rna);

	_initial_concentration = concentration_;
  _delta_concentration   = 0;
  _inherited             = false;
  _signal                = false;
  is_TF_			       = false;
  _id = id++;

  if (!AA_list_.empty()) {
    _cod_tab = new int8_t[AA_list_.size()];
    int i = 0;
    for (auto cod : AA_list_) {
      _cod_tab[i] = cod->value();
      i++;
    }
  }
}

//used to build the signal protein
Protein_R::Protein_R( const std::list<Codon*> codon_list, double concentration, double w_max)  :
		Protein::Protein( codon_list, concentration, w_max )
{
  _initial_concentration = 0;
  _delta_concentration  = 0;
  _inherited            = false;
  _signal               = true;
  is_TF_			 = false;
  _id = id++;

  if (!AA_list_.empty()) {
    _cod_tab = new int8_t[AA_list_.size()];
    int i = 0;
    for (auto cod : AA_list_) {
      _cod_tab[i] = cod->value();
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
  is_TF_			 = false;
  _id = id++;

  if (!AA_list_.empty()) {
    _cod_tab = new int8_t[AA_list_.size()];
    int i = 0;
    for (auto cod : AA_list_) {
      _cod_tab[i] = cod->value();
      i++;
    }
  }
}

// =================================================================
//                             Destructors
// =================================================================
Protein_R::~Protein_R( void )
{
  for (unsigned int i = 0; i < _rna_R_list.size(); i++)
    _rna_R_list[i] = NULL;

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
    //printf("Protein %ld is generated by ",_id);
	for (auto& rna: _rna_R_list)
    {
      //if (_id == 34483) printf("%ld (influenced by %ld) at %f  - ",rna->get_id(),rna->_operating_coef_list.size(),rna->get_synthesis_rate());
      assert( _inherited == false);

      _delta_concentration += rna->get_synthesis_rate();
      ///if (_id == 34483) printf("Prot %ld synthesis by %ld at rate %e\n",_id,rna->get_id(),rna->get_synthesis_rate());
    }
    //if (_id == 34483) printf("\n");
    //if (_id == 34483)  printf("Prot %ld BEFORE DEGRADATION concentration %f %f\n",_id,concentration_,_delta_concentration);

    _delta_concentration -= gen_unit_->exp_m()->exp_s()->get_degradation_rate() * concentration_;
    _delta_concentration *= 1/((double)gen_unit_->exp_m()->exp_s()->get_nb_degradation_step());

    //if (_id == 34483)  printf("Prot %ld AFTER degradation concentration %f %f\n",_id,concentration_,_delta_concentration);
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
  _initial_concentration += rna->basal_level();
  //printf("Add RNA %ld to protein %ld (influence by %ld)\n",((Rna_R*)rna)->get_id(),_id,((Rna_R*)rna)->_operating_coef_list.size());
  _rna_R_list.push_back((Rna_R*)rna);
}

} // namespace aevol
