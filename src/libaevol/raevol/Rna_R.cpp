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
#include <math.h>

// =================================================================
//                            Project Files
// =================================================================
#include "Rna_R.h"
#include "Individual_R.h"
#include "ExpManager.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                           Class ae_rna_R                                    #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
Rna_R::Rna_R( GeneticUnit* gen_unit, const Rna_R &model ) : Rna( gen_unit, model )
{
  _protein_list = model._protein_list;
}

Rna_R::Rna_R( GeneticUnit* gen_unit, Strand strand, int32_t index, int8_t diff ) :
		Rna( gen_unit, strand, index, diff )
{
}

/*
ae_rna_R::ae_rna_R( ae_rna_R* parent ) :
ae_rna( parent )
{
  _influence_list = parent->_influence_list->copy();
}
*/

// =================================================================
//                             Destructors
// =================================================================
Rna_R::~Rna_R( void )
{
	_protein_list.clear();
	_enhancing_coef_list.clear();
	_operating_coef_list.clear();
}

// =================================================================
//                            Public Methods
// =================================================================
void Rna_R::set_influences( const std::list<Protein*> protein_list )
{
    int32_t enhancer_position = get_enhancer_position();
	  int32_t operator_position = get_operator_position();

	  _protein_list.clear();
    _protein_list.resize(protein_list.size());
    _enhancing_coef_list.clear();
	  _enhancing_coef_list.resize(protein_list.size());
    _operating_coef_list.clear();
	  _operating_coef_list.resize(protein_list.size());

    int i = 0;
	  for (auto& prot : protein_list) {
      _protein_list[i] = (Protein_R*) prot;

		  _enhancing_coef_list[i] = affinity_with_protein( enhancer_position, prot );
		  _operating_coef_list[i] = affinity_with_protein( operator_position, prot );

	    if ( _enhancing_coef_list[i] != 0.0 || _operating_coef_list[i] != 0.0 )
	    	((Protein_R*)prot)->not_pure_TF = true;

      i++;
	  }

}

double Rna_R::get_synthesis_rate( void )
{
  double enhancer_activity  = 0;
  double operator_activity  = 0;

  for (unsigned int i = 0; i < _enhancing_coef_list.size(); i++) {
  	enhancer_activity  += _enhancing_coef_list[i];
    operator_activity  += _operating_coef_list[i];
  }

  double enhancer_activity_pow_n  = pow( enhancer_activity, _gen_unit->get_exp_m()->get_exp_s()->get_hill_shape_n() );
  double operator_activity_pow_n  = pow( operator_activity, _gen_unit->get_exp_m()->get_exp_s()->get_hill_shape_n() );

  return   _basal_level
           * (_gen_unit->get_exp_m()->get_exp_s()->get_hill_shape()
              / (operator_activity_pow_n + _gen_unit->get_exp_m()->get_exp_s()->get_hill_shape()))
           * (1 + ((1 / _basal_level) - 1)
                  * (enhancer_activity_pow_n /
                     (enhancer_activity_pow_n + _gen_unit->get_exp_m()->get_exp_s()->get_hill_shape())));
}

// =================================================================
//                           Protected Methods
// =================================================================
int32_t Rna_R::get_enhancer_position( void )
{
  if(_strand == LEADING)
  {
    return (_pos - 20)  % ( _gen_unit->get_dna()->get_length() );
  }
  else  // _strand = LAGGING
  {
    return (_pos + 20)  % ( _gen_unit->get_dna()->get_length() );
  }
}

int32_t Rna_R::get_operator_position( void )
{
  if(_strand == LEADING)
  {
    return (_pos + PROM_SIZE)  % ( _gen_unit->get_dna()->get_length() );
  }
  else  // _strand = LAGGING
  {
    return (_pos - PROM_SIZE)  % ( _gen_unit->get_dna()->get_length() );
  }
}

double Rna_R::affinity_with_protein( int32_t index, Protein *protein )
{
	  double  max = 0;
	  double  temp;
	  int32_t len = protein->get_length();
	  int32_t quadon_tab[5];
	//  int32_t* codon_tab;
	//  codon_tab = new int32_t[len];

	  Individual_R*  indiv = NULL;
	  Protein_R*     prot  = NULL;

	  // printf("affinity_with_protein - len = %d\n",len);

	  // Putting the quadons and the codons on local tab
	  indiv = dynamic_cast< Individual_R* >( _gen_unit->get_indiv() );
	  prot  = ( Protein_R* )( protein );
	  for ( int32_t i = 0 ; i < 5; i++ )
	  {
	    quadon_tab[i] = indiv->get_quadon( _gen_unit, _strand, (index+i) );
	  }
	//  for (int32_t i = 0 ; i < len ; i++ )
	//  {
	//    codon_tab[i] = prot->get_codon(i);
	//  }
	//
	//  t2 = high_resolution_clock::now();
	//    	  duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	//    	  ae_logger::addLog(CODON,duration);
	//    	  t1 = t2;

	  // Calculate the affinity
	  for ( int32_t i = 0 ; i < len - 4; i++ )
	  {
	    temp  = 1 *
	    		_gen_unit->get_exp_m()->get_exp_s()->get_binding_matrix(quadon_tab[0],prot->get_cod_tab(i)) *
	    		_gen_unit->get_exp_m()->get_exp_s()->get_binding_matrix(quadon_tab[1],prot->get_cod_tab(i+1)) *
	    		_gen_unit->get_exp_m()->get_exp_s()->get_binding_matrix(quadon_tab[2],prot->get_cod_tab(i+2)) *
	    		_gen_unit->get_exp_m()->get_exp_s()->get_binding_matrix(quadon_tab[3],prot->get_cod_tab(i+3)) *
	    		_gen_unit->get_exp_m()->get_exp_s()->get_binding_matrix(quadon_tab[4],prot->get_cod_tab(i+4));

	//    for ( int32_t j = 0 ; j < 5 ; j++ )
	//    {
	//      temp *= ae_common::_binding_matrix[quadon_tab[j]][prot->_AA_list[i+j]];
	//
	//      if ( temp == 0.0 ) break;
	//    }

	    max = (max < temp) ? temp : max;
	  }
	//  delete [] codon_tab;
	  return max;
}
} // namespace aevol