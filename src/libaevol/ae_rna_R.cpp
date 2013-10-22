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
#include <ae_rna_R.h>
#include <ae_influence_R.h>
#include <ae_individual_R.h>

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
ae_rna_R::ae_rna_R( void ) : ae_rna()
{
  _influence_list = new ae_list();
}

ae_rna_R::ae_rna_R( ae_genetic_unit* gen_unit, const ae_rna_R &model ) : ae_rna( gen_unit, model )
{
  _influence_list = new ae_list();
}

ae_rna_R::ae_rna_R( ae_genetic_unit* gen_unit ) : ae_rna( gen_unit )
{
  _influence_list = new ae_list();
}

ae_rna_R::ae_rna_R( ae_genetic_unit* gen_unit, ae_strand strand, int32_t index, int8_t diff ) :
  ae_rna( gen_unit, strand, index, diff )
{
  _influence_list = new ae_list();
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
ae_rna_R::~ae_rna_R( void )
{
  //printf("nb_el = %ld\n", _influence_list->get_nb_elts());
  _influence_list->erase( DELETE_OBJ );
  delete _influence_list;

}

// =================================================================
//                            Public Methods
// =================================================================
void ae_rna_R::set_influences( ae_list* protein_list )
{
  int32_t enhancer_position = get_enhancer_position();
  int32_t operator_position = get_operator_position();

  ae_list_node<ae_protein*>* prot_node;
  ae_protein* prot;

  prot_node = protein_list->get_first();
  while ( prot_node != NULL )
  {
    prot = prot_node->get_obj();
    
    double activation = affinity_with_protein( enhancer_position, prot );
    double inhibition = affinity_with_protein( operator_position, prot );

    if ( activation != 0.0 || inhibition != 0.0 )
    {
      add_influence(  prot, activation, inhibition );
    }

    prot_node = prot_node->get_next();
  }
}

double ae_rna_R::get_synthesis_rate( void )
{
  double enhancer_activity  = 0;
  double operator_activity  = 0;

  ae_list_node<ae_influence_R*>*   influence_node;
  ae_influence_R* influence;

  influence_node = _influence_list->get_first();
  while ( influence_node != NULL )
  {
    influence = influence_node->get_obj();

    enhancer_activity  += influence->get_enhancer_activity();
    operator_activity  += influence->get_operator_activity();

    influence_node = influence_node->get_next();
  }

  double enhancer_activity_pow_n  = pow( enhancer_activity, ae_common::hill_shape_n );
  double operator_activity_pow_n  = pow( operator_activity, ae_common::hill_shape_n );

  return   _basal_level
        * (ae_common::hill_shape / (operator_activity_pow_n + ae_common::hill_shape))
        * (1 + ((1 / _basal_level) - 1) * (enhancer_activity_pow_n / (enhancer_activity_pow_n + ae_common::hill_shape)));
}

void ae_rna_R::remove_influence( ae_influence_R* influence )
{
  //printf("remove_influence, rna = %p, nb_el = %ld\n", this, _influence_list->get_nb_elts());
  _influence_list->remove( influence, DELETE_OBJ /*delete_node*/, DELETE_OBJ /*delete_obj*/ );
}

// =================================================================
//                           Protected Methods
// =================================================================
int32_t ae_rna_R::get_enhancer_position( void )
{
  if(_strand == LEADING)
  {
    return utils::mod( (_pos - 20) , _gen_unit->get_dna()->get_length() );
  }
  else  // _strand = LAGGING
  {
    return utils::mod( (_pos + 20) , _gen_unit->get_dna()->get_length() );
  }
}

int32_t ae_rna_R::get_operator_position( void )
{
  if(_strand == LEADING)
  {
    return utils::mod( (_pos + PROM_SIZE) , _gen_unit->get_dna()->get_length() );
  }
  else  // _strand = LAGGING
  {
    return utils::mod( (_pos - PROM_SIZE) , _gen_unit->get_dna()->get_length() );
  }
}

void ae_rna_R::add_influence( ae_protein *protein, double enhancing_coef, double operating_coef )
{
  ae_influence_R* influence = new ae_influence_R( this, protein, enhancing_coef, operating_coef );
  //printf("add_influence, rna = %p, influence = %p\n", this, influence);
  //printf("nb_el = %ld\n", _influence_list->get_nb_elts());
  _influence_list->add( influence );
  //printf("/ nb_el = %ld\n", _influence_list->get_nb_elts());
  dynamic_cast< ae_protein_R* >( protein )->add_influence( influence );
}

double ae_rna_R::affinity_with_protein( int32_t index, ae_protein *protein )
{
  double  max = 0;
  double  temp;
  int32_t len = protein->get_length();
  int32_t quadon_tab[5];
  int32_t* codon_tab;
  codon_tab = new int32_t[len];

  ae_individual_R*  indiv = NULL;
  ae_protein_R*     prot  = NULL;

  // Putting the quadons and the codons on local tab
  indiv = dynamic_cast< ae_individual_R* >( _gen_unit->get_indiv() );
  prot  = ( ae_protein_R* )( protein );

  for ( int32_t i = 0 ; i < 5; i++ )
  {
    quadon_tab[i] = indiv->get_quadon( _gen_unit, _strand, (index+i) );
  }

  for (int32_t i = 0 ; i < len ; i++ )
  {
    codon_tab[i] = prot->get_codon(i);
  }

  // Calculate the affinity
  for ( int32_t i = 0 ; i < len - 4; i++ )
  {
    temp  = 1;
    
    for ( int32_t j = 0 ; j < 5 ; j++ )
    {
      temp *= ae_common::get_binding( quadon_tab[j] , codon_tab[i+j] );

      if ( temp == 0.0 ) break;
    }
    
    max = (max < temp) ? temp : max;
  }

  delete [] codon_tab;
  return max;

/*
double  max = 0;
  double  temp;
  int32_t len = protein->get_length();

  ae_individual_R*  indiv =           NULL;
  ae_protein_R*     prot  =           NULL;

  for ( int32_t i = 0 ; i < len - 4 ; i++ )
  {
    temp  = 1;
    
    for ( int32_t j = 0 ; (temp != 0.0) && (j < 5) ; j++ )
    {
      indiv = dynamic_cast< ae_individual_R* >( _gen_unit->get_indiv() );
      prot  = ( ae_protein_R* )( protein );
      temp *= ae_common::get_binding( indiv->get_quadon( _gen_unit, _strand, (index+j) ), prot->get_codon(i+j));
    }
    
    max = (max < temp) ? temp : max;
  }
  
  return max;
*/
}
