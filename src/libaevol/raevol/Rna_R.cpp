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
#ifdef __BLAS__
#include <cblas.h>
#endif
// =================================================================
//                            Project Files
// =================================================================
#include "Rna_R.h"
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

long Rna_R::id = 0;
// =================================================================
//                             Constructors
// =================================================================
Rna_R::Rna_R( GeneticUnit* gen_unit, const Rna_R &model ) : Rna( gen_unit, model )
{
  _protein_list = model._protein_list;
  _enhancing_coef_list = model._enhancing_coef_list;
  _operating_coef_list = model._operating_coef_list;
  _id = model._id;
  _nb_influences = model._nb_influences;
}

Rna_R::Rna_R( GeneticUnit* gen_unit, Strand strand, int32_t index, int8_t diff ) :
		Rna( gen_unit, strand, index, diff )
{
  _id = id++;
  _nb_influences = 0;
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
void Rna_R::set_influences( std::list<Protein*>& protein_list )
{
	int32_t enhancer_position = get_enhancer_position();
	int32_t operator_position = get_operator_position();

	_protein_list.clear();
  std::vector<Protein_R*>::iterator itprot = _protein_list.begin();
  _protein_list.resize(0);
 _enhancing_coef_list.clear();
  std::vector<ProteinConcentration>::iterator itenh = _enhancing_coef_list.begin();
	_enhancing_coef_list.resize(0);
  _operating_coef_list.clear();
  std::vector<ProteinConcentration>::iterator itope = _operating_coef_list.begin();
	_operating_coef_list.resize(0);


  int i = 0;
  ProteinConcentration enhance=0,operate=0;

  //#pragma omp simd
	for (auto& prot : protein_list) {
    enhance = affinity_with_protein( enhancer_position, prot );
    operate = affinity_with_protein( operator_position, prot );

    if (enhance != 0.0 || operate != 0.0) {

      _protein_list.insert(itprot,(Protein_R*) prot);

      _enhancing_coef_list.insert(itenh, enhance);
      _operating_coef_list.insert(itope, operate);

      _protein_list[i]->is_TF_ = true;
     // _protein_concentration_list[i] = prot->concentration();
      i++;
      itprot = _protein_list.begin()+i;
      itenh = _enhancing_coef_list.begin()+i;
      itope = _operating_coef_list.begin()+i;
    }
    //else
    // _protein_list[i] = nullptr;
	}
  _nb_influences = i;

  //TODO NOT USEFUL ??? _nb_influences = i==0 ? 0 : i-1;

  /*if (gen_unit_->indiv()->id() == 12608)
    printf("12608 RNA %d is influenced by %d proteins\n",_id,_nb_influences);*/

  /*if (protein_list.size() > 0)
    printf("Set Influences of RNA %ld with %ld %ld %ld\n",_id,_enhancing_coef_list.size(),_operating_coef_list.size(),
         _protein_list.size());*/

}

ProteinConcentration Rna_R::get_synthesis_rate( void )
{
  ProteinConcentration enhancer_activity  = 0;
  ProteinConcentration operator_activity  = 0;

//#ifndef __BLAS__
  for (int i = 0; i < _nb_influences; i++) {
    /*if (gen_unit_->indiv()->id() == 12608)
      printf("12608 RNA %d  due to Protein %d is %f %f %f\n",
             _id,_protein_list[i]->get_id(),
             _enhancing_coef_list[i],_operating_coef_list[i],_protein_list[i]->concentration_);*/

  	enhancer_activity  += _enhancing_coef_list[i] * _protein_list[i]->concentration_;
    operator_activity  += _operating_coef_list[i] * _protein_list[i]->concentration_;
  }
/*#else
  ProteinConcentration enhancer_tab[_nb_influences];
  ProteinConcentration operator_tab[_nb_influences];

#ifdef __SIMD
  #pragma omp simd
#endif
  for (int i = 0; i < _nb_influences; i++) {
  	enhancer_tab[i] = _enhancing_coef_list[i] * _protein_list[i]->concentration_;
  	}

#ifdef __SIMD
  #pragma omp simd
#endif
  for (int i = 0; i < _nb_influences; i++) {
    operator_tab[i] = _operating_coef_list[i] * _protein_list[i]->concentration_;
  }

#ifdef __FLOAT_CONCENTRATION
  ProteinConcentration enhancer_activity_blas = cblas_sasum(_nb_influences,enhancer_tab,1);
  ProteinConcentration operator_activity_blas = cblas_sasum(_nb_influences,operator_tab,1);
#else
  ProteinConcentration enhancer_activity_blas = cblas_dasum(_nb_influences,enhancer_tab,1);
  ProteinConcentration operator_activity_blas = cblas_dasum(_nb_influences,operator_tab,1);
#endif
#endif
*/
  ProteinConcentration enhancer_activity_pow_n  = pow( enhancer_activity, gen_unit_->exp_m()->exp_s()->get_hill_shape_n() );
  ProteinConcentration operator_activity_pow_n  = pow( operator_activity, gen_unit_->exp_m()->exp_s()->get_hill_shape_n() );

  //if (enhancer_activity != 0.0 || operator_activity != 0.0)
  /*if (_id == 132073) printf("Synthesis of RNA %ld : E %f O %f EP %f OP %f SN %f S %f B %f\n",_id,enhancer_activity,operator_activity,enhancer_activity_pow_n,
                                                operator_activity_pow_n,gen_unit_->exp_m()->exp_s()->get_hill_shape_n(),
         gen_unit_->exp_m()->exp_s()->get_hill_shape(),basal_level_);*/
  /*if (gen_unit_->indiv()->id() == 12608)
    printf("12608 RNA %d is %f %f %f\n",_id,enhancer_activity,operator_activity,basal_level_);*/

  return   basal_level_
           * (gen_unit_->exp_m()->exp_s()->get_hill_shape()
              / (operator_activity_pow_n + gen_unit_->exp_m()->exp_s()->get_hill_shape()))
           * (1 + ((1 / basal_level_) - 1)
                  * (enhancer_activity_pow_n /
                     (enhancer_activity_pow_n + gen_unit_->exp_m()->exp_s()->get_hill_shape())));
}

// =================================================================
//                           Protected Methods
// =================================================================
int32_t Rna_R::get_enhancer_position( void )
{
  int32_t length = gen_unit_->dna()->length();

  if(strand_ == LEADING)
  {
    return (pos_ - 20)  % ( length ) < 0 ?
           ((pos_ - 20)  % ( length )) + ( length ) :
           (pos_ - 20)  % ( length );
  }
  else  // strand_ = LAGGING
  {
    return (pos_ + 20)  % ( length ) < 0 ?
           ((pos_ + 20)  % ( length )) + ( length ) :
           (pos_ + 20)  % ( length );
  }
}

int32_t Rna_R::get_operator_position( void )
{
  int32_t length = gen_unit_->dna()->length();

  if(strand_ == LEADING)
  {
    return (pos_ + PROM_SIZE)  % ( length ) < 0 ?
           (pos_ + PROM_SIZE)  % ( length ) + (length) :
           (pos_ + PROM_SIZE)  % ( length );
  }
  else  // strand_ = LAGGING
  {
    return (pos_ - PROM_SIZE)  % ( length ) < 0 ?
           (pos_ - PROM_SIZE)  % ( length ) + (length) :
           (pos_ - PROM_SIZE)  % ( length );
  }
}

ProteinConcentration Rna_R::affinity_with_protein( int32_t index, Protein *protein )
{
  int32_t len = protein->length();

  if (len > 5) {

#ifndef __BLAS__
    ProteinConcentration max = 0;
    ProteinConcentration temp = 1;
#else
    ProteinConcentration  tab_temp[len-4];
#endif


    int32_t quadon_tab[5];
    Protein_R* prot = NULL;
    prot = (Protein_R*) (protein);
    for (int32_t i = 0; i < 5; i++) {
      quadon_tab[i] = gen_unit_->indiv_r_->get_quadon(gen_unit_, strand_, (index + i));
    }


    ProteinConcentration (* binding_matrix)[MAX_QUADON][MAX_CODON] = &(gen_unit_->exp_m()->exp_s()->_binding_matrix);
#ifdef __SIMD
#pragma omp simd
#endif
    for (int32_t i = 0; i < len - 4; i++) {

#ifndef __BLAS__
      temp = 1;
#else
      tab_temp[i]  = (*binding_matrix)[quadon_tab[0]][prot->_cod_tab[i]];
      tab_temp[i]  *= (*binding_matrix)[quadon_tab[1]][prot->_cod_tab[i+1]];
      tab_temp[i]  *= (*binding_matrix)[quadon_tab[2]][prot->_cod_tab[i+2]];
      tab_temp[i]  *= (*binding_matrix)[quadon_tab[3]][prot->_cod_tab[i+3]];
      tab_temp[i]  *= (*binding_matrix)[quadon_tab[4]][prot->_cod_tab[i+4]];
#endif

#ifndef __BLAS__
      for (int8_t j = 0; j < 5; j++) {
        temp *= gen_unit_->exp_m()->exp_s()->_binding_matrix[quadon_tab[j]][prot->_cod_tab[
            i + j]];
      }

      max = (max < temp) ? temp : max;
#endif
    }

#ifdef __BLAS__
#ifdef __FLOAT_CONCENTRATION
    return tab_temp[cblas_isamax(len-4,tab_temp,1)];
#else
    return tab_temp[cblas_idamax(len-4,tab_temp,1)];
#endif
#else
    return max;
#endif
  } else {
    return 0.0;
  }
}
} // namespace aevol
