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
#include <ae_individual_R.h>

//##############################################################################
//                                                                             #
//                           Class ae_individual_R                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_individual_R::ae_individual_R( const ae_individual_R &model ) : ae_individual( model ) // TODO : When is this used?
{
  //printf("ae_individual_R( model )");
  
  _rna_list_coding = new ae_list();
  
  if( ae_common::with_heredity )
  {
    _inherited_protein_list = new ae_list();
    
    // We copy all the proteins from model
    for ( int8_t strand = LEADING ; strand <= LAGGING ; strand++ )
    {
      ae_list_node<ae_protein_R*>* prot_node = model._protein_list->get_first();
      ae_protein_R* prot;
      
      while ( prot_node != NULL )
      {
        prot = prot_node->get_obj();
        
        ae_protein_R* inherited_prot = new ae_protein_R( NULL, *prot );
        inherited_prot->set_inherited( true );
        
        _inherited_protein_list->add( inherited_prot );

        prot_node = prot_node->get_next();
      }
    }
  }
  else
  {
    _inherited_protein_list = NULL;
  }
}

/*
 * Used at initialization
*/
ae_individual_R::ae_individual_R( void ) : ae_individual()
{
  _rna_list_coding = new ae_list();
  
  if( ae_common::with_heredity )
  {
    _inherited_protein_list = new ae_list();
  }
  else
  {
    _inherited_protein_list = NULL;
  }
}

ae_individual_R::ae_individual_R( ae_individual_R* parent, int32_t id,
                                  ae_jumping_mt* mut_prng, ae_jumping_mt* stoch_prng )
        : ae_individual( parent, id, mut_prng, stoch_prng )
{
  //~ printf( "ae_individual_R( parent ) : I have %d inherited proteins\n", parent->get_protein_list()->get_nb_elts() );
  
  _rna_list_coding = new ae_list();
  
  if( ae_common::with_heredity )
  {
    _inherited_protein_list = new ae_list();
    
    // We copy all the proteins from parent
    for ( int8_t strand = LEADING ; strand <= LAGGING ; strand++ )
    {
      ae_list_node<ae_protein_R*>* prot_node = parent->_protein_list->get_first();
      ae_protein_R* prot;
      
      while ( prot_node != NULL )
      {
        prot = prot_node->get_obj();
        
        ae_protein_R* inherited_prot = new ae_protein_R( NULL, *prot );
        inherited_prot->set_inherited( true );
        
        _inherited_protein_list->add( inherited_prot );

        prot_node = prot_node->get_next();
      }
    }
  }
  else
  {
    _inherited_protein_list = NULL;
  }
}

ae_individual_R::ae_individual_R( gzFile backup_file ) : ae_individual( backup_file )
{
  _rna_list_coding = new ae_list();
  
  _inherited_protein_list = NULL;
  if( ae_common::with_heredity)
  {
    // Retreive inherited proteins
    _inherited_protein_list = new ae_list();
    int16_t nb_inherited_proteins = 0;
    gzread( backup_file, &nb_inherited_proteins,  sizeof(nb_inherited_proteins) );
  
    for ( int16_t i = 0 ; i < nb_inherited_proteins ; i++ )
    {
	  _inherited_protein_list->add( new ae_protein_R( backup_file ) );
    }
  }  
}

// =================================================================
//                             Destructors
// =================================================================
ae_individual_R::~ae_individual_R( void )
{
  assert( !ae_common::with_heredity || _inherited_protein_list != NULL );
  
  if ( _inherited_protein_list != NULL )
  {
    _inherited_protein_list->erase( DELETE_OBJ );
    delete _inherited_protein_list;
    _inherited_protein_list = NULL;
  }
  
  if ( _rna_list_coding != NULL )
  {
    _rna_list_coding->erase( NO_DELETE );
    delete _rna_list_coding;
    _rna_list_coding = NULL;
  }
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_individual_R::evaluate( ae_environment* envir )
{
  // ---------------------------------------------------------------------------
  // 1) Transcription - Translation - Folding
  // ---------------------------------------------------------------------------
  ae_list_node<ae_genetic_unit*>*     gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit = NULL;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = gen_unit_node->get_obj();
    
    gen_unit->do_transcription();
    gen_unit->do_translation();
    gen_unit->compute_phenotypic_contribution(); // Compute basal levels
    
    gen_unit_node = gen_unit_node->get_next();
  }
  
  //----------------------------------------------------------------------------
  // 2) Make a list of all the proteins present in the individual
  //----------------------------------------------------------------------------
  make_protein_list();
  make_rna_list();

  //----------------------------------------------------------------------------
  // 3) Create influence graph
  //----------------------------------------------------------------------------
  set_influences();

  // Computes some statistical data about the individual
  //~ _genome->compute_statistical_data();

/*
  printf("number of protein : inherited : %d \n",_inherited_protein_list->get_nb_elts());
  printf("number of protein : leading : %d \n",_genetic_unit_list->get_first()->get_obj()->get_protein_list()[LEADING]->get_nb_elts());
  printf("number of protein : lagging : %d \n",_genetic_unit_list->get_first()->get_obj()->get_protein_list()[LAGGING]->get_nb_elts());
  printf("number of protein : total : %d \n",_protein_list->get_nb_elts());
*/

  //----------------------------------------------------------------------------
  // 4) Make the individual "live its life" and compute partial phenotypes and 
  //    fitnesses
  //----------------------------------------------------------------------------
  int16_t indiv_age   = 0;
  double fitness_temp = 0;
  
  // Go from an evaluation date to the next
  for( int16_t evaluation_index = 0 ; evaluation_index < ae_common::individual_evaluation_nbr ; evaluation_index++ )
  {
    // Let the individual evolve until the evaluation date
    while( ( indiv_age < ae_common::individual_evaluation_dates->get_value( evaluation_index ) ) )
    {
	  //Updating the concentrations in order to respect the degradation step.
	  for( int i = 0; i < 1/ae_common::degradation_step; i++ )
	  {
        update_concentrations();
        indiv_age++;
	  }
    }

    // Evaluate the individual's phenotype
    compute_phenotype();
    compute_distance_to_target( envir );
    compute_fitness( envir );
    fitness_temp += _fitness;
  }

  // The individual may have some years left to live
  // TODO : Useless if no mutation during lifetime => commented
  //~ while( indiv_age <= ae_common::individual_life_time )
  //~ {
    //~ update_concentrations();
    //~ indiv_age++;
  //~ }

  //----------------------------------------------------------------------------
  // 5) Compute final fitness
  //----------------------------------------------------------------------------
  _fitness = fitness_temp / (double)ae_common::individual_evaluation_nbr;
}

void ae_individual_R::set_influences( void )
// Compute the influence of each protein over each coding RNA
// As non-coding RNAs are completely inert, we don't care about their concentration
// so we don't care if proteins activate or inhibit their transcription.
{
  ae_list_node<ae_rna_R*>* rna_node  = NULL;
  ae_rna_R*     rna       = NULL;
  
  //
  rna_node = _rna_list_coding->get_first();
  while ( rna_node != NULL )
  {
    rna = rna_node->get_obj();

    //~ printf( "%d proteins\n", _protein_list->get_nb_elts() );
    rna->set_influences( _protein_list );

    rna_node = rna_node->get_next();
  }
}

void ae_individual_R::update_concentrations( void )
{
  //_phenotype->print_points();

  ae_list_node<ae_protein_R*>* prot_node       = NULL;
  ae_list_node>ae_protein_R*>* next_prot_node  = NULL;
  ae_protein_R* prot            = NULL;

  // Compute all the changes that will be applied to the concentrations
  // Concentrations must not be changed at this stage
  prot_node = _protein_list->get_first();
  while ( prot_node != NULL )
  {
    prot = prot_node->get_obj();

    prot->compute_delta_concentration();

    prot_node = prot_node->get_next();
  }

  // Apply the changes in concentrations we have just computed
  prot_node = _protein_list->get_first();
  while ( prot_node != NULL )
  {
    prot = prot_node->get_obj();

    prot->update_concentration();

    // Keep track of the next node in the list, in case we need to remove the current one
    next_prot_node = prot_node->get_next();

    if( ae_common::with_heredity )
    {
      if( prot->is_inherited() && prot->get_concentration() < ae_common::protein_presence_limit )
      {
        // The protein has to be removed from both the individual's _protein_list and _inherited_proteins_list
        _protein_list->remove( prot_node, DELETE_OBJ /*delete_node*/, NO_DELETE /*delete_obj*/ );
        _inherited_protein_list->remove( prot, DELETE_OBJ /*delete_node*/, DELETE_OBJ /*delete_obj*/ );
      }
    }

    // ready to study the next node
    prot_node = next_prot_node;
  }
}

// Multiply the concentration of each protein by <factor>
void ae_individual_R::multiply_concentrations( double factor )
{
  ae_list_node<ae_protein_R*>* prot_node = _protein_list->get_first();
  ae_protein_R* prot      = NULL;
  
  while ( prot_node != NULL )
  {
    prot = prot_node->get_obj();

    prot->multiply_concentration( factor );

    prot_node = prot_node->get_next();
  }
}

int8_t ae_individual_R::get_quadon( ae_genetic_unit* gen_unit, ae_strand strand, int32_t pos )
{
  const char* dna = gen_unit->get_dna()->get_data();
  int32_t  len    = gen_unit->get_dna()->get_length();
  int8_t quadon   = 0;

  if ( strand == LEADING )
  {
    for ( int8_t i = 0 ; i < QUADON_SIZE ; i++ )
    {
      if ( dna[utils::mod((pos+i),len)] == '1' )
      {
        quadon += 1 << (QUADON_SIZE - i - 1);  //pow( 2, QUADON_SIZE - i - 1 );
      }
    }
  }
  else  // ( strand == LAGGING )
  {
    for ( int8_t i = 0 ; i < QUADON_SIZE ; i++ )
    {
      if ( dna[utils::mod((pos-i),len)] != '1' ) // == and not != because we are on the complementary strand...
      {
        quadon += 1 << (QUADON_SIZE - i - 1);  //pow( 2, QUADON_SIZE - i - 1 );
      }
    }
  }

  return quadon;
}

void ae_individual_R::save( gzFile backup_file )
{
  ae_individual::save( backup_file );
  // Test if there is heredity, and if the generation is the first one (no inherited protein list).
  if (ae_common::with_heredity && _inherited_protein_list != NULL )
  {
    // Write inherited proteins
    int16_t nb_inherited_proteins = _inherited_protein_list->get_nb_elts();
    gzwrite( backup_file, &nb_inherited_proteins,  sizeof(nb_inherited_proteins) );

    ae_list_node<ae_protein_R*>*  inherited_protein_node = _inherited_protein_list->get_first();
    ae_protein_R*  inherited_protein;

    for ( int16_t i = 0 ; i < nb_inherited_proteins ; i++ )
    {
    inherited_protein = inherited_protein_node->get_obj();
    
    inherited_protein->save( backup_file );
    
    inherited_protein_node = inherited_protein_node->get_next();
    }
  }
}
// =================================================================
//                           Protected Methods
// =================================================================
void ae_individual_R::make_protein_list( void )
{
  ae_individual::make_protein_list();
  _protein_list->add_list( _inherited_protein_list );
}

void ae_individual_R::make_rna_list( void )
{
  ae_individual::make_rna_list();
  
  // Parse the newly created RNA list and copy the coding RNAs in _rna_list_coding.
  ae_list_node<ae_rna*>* rna_node  = _rna_list->get_first();
  ae_rna*       rna       = NULL;
  
  while ( rna_node != NULL )
  {
    rna = rna_node->get_obj();
    
    if ( rna->is_coding() == true )
    {
      _rna_list_coding->add( rna );
    }
    
    rna_node = rna_node->get_next();
  }
}
