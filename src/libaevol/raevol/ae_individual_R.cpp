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
#include "ae_individual_R.h"
namespace aevol {

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

/*
 * Used at initialization
*/
ae_individual_R::ae_individual_R() : ae_individual()
{
  rna_list_coding_ = new ae_list();
  
  if( ae_common::with_heredity )
  {
    inherited_protein_list_ = new ae_list();
  }
  else
  {
    inherited_protein_list_ = NULL;
  }
}

ae_individual_R::ae_individual_R( ae_individual_R* parent, int32_t id,
                                  ae_jumping_mt* mut_prng, ae_jumping_mt* stoch_prng )
        : ae_individual( parent, id, mut_prng, stoch_prng )
{
  //~ printf( "ae_individual_R( parent ) : I have %d inherited proteins\n", parent->protein_list()->nb_elts() );
  
  rna_list_coding_ = new ae_list();
  
  if( ae_common::with_heredity )
  {
    inherited_protein_list_ = new ae_list();
    
    // We copy all the proteins from parent
    for ( int8_t strand = LEADING ; strand <= LAGGING ; strand++ )
    {
      ae_list_node<ae_protein_R*>* prot_node = parent->protein_list_->first();
      ae_protein_R* prot;
      
      while ( prot_node != NULL )
      {
        prot = prot_node->obj();
        
        ae_protein_R* inherited_prot = new ae_protein_R( NULL, *prot );
        inherited_prot->set_inherited( true );
        
        inherited_protein_list_->add( inherited_prot );

        prot_node = prot_node->next();
      }
    }
  }
  else
  {
    inherited_protein_list_ = NULL;
  }
}

ae_individual_R::ae_individual_R( gzFile backup_file ) : ae_individual( backup_file )
{
  rna_list_coding_ = new ae_list();
  
  inherited_protein_list_ = NULL;
  if( ae_common::with_heredity)
  {
    // Retreive inherited proteins
    inherited_protein_list_ = new ae_list();
    int16_t nb_inherited_proteins = 0;
    gzread( backup_file, &nb_inherited_proteins,  sizeof(nb_inherited_proteins) );
  
    for ( int16_t i = 0 ; i < nb_inherited_proteins ; i++ )
    {
	  inherited_protein_list_->add( new ae_protein_R( backup_file ) );
    }
  }  
}

// =================================================================
//                             Destructors
// =================================================================
ae_individual_R::~ae_individual_R()
{
  assert( !ae_common::with_heredity || inherited_protein_list_ != NULL );
  
  if ( inherited_protein_list_ != NULL )
  {
    inherited_protein_list_->erase( DELETE_OBJ );
    delete inherited_protein_list_;
    inherited_protein_list_ = NULL;
  }
  
  if ( rna_list_coding_ != NULL )
  {
    rna_list_coding_->erase( NO_DELETE );
    delete rna_list_coding_;
    rna_list_coding_ = NULL;
  }
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_individual_R::evaluate( Environment* envir )
{
  // ---------------------------------------------------------------------------
  // 1) Transcription - Translation - Folding
  // ---------------------------------------------------------------------------
  ae_list_node<GeneticUnit*>*     gen_unit_node = genetic_unit_list_->first();
  GeneticUnit*  gen_unit = NULL;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = gen_unit_node->obj();
    
    gen_unit->do_transcription();
    gen_unit->do_translation();
    gen_unit->compute_phenotypic_contribution(); // Compute basal levels
    
    gen_unit_node = gen_unit_node->next();
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
  //~ genome_->compute_statistical_data();

/*
  printf("number of protein : inherited : %d \n",inherited_protein_list_->nb_elts());
  printf("number of protein : leading : %d \n",genetic_unit_list_->first()->obj()->protein_list()[LEADING]->nb_elts());
  printf("number of protein : lagging : %d \n",genetic_unit_list_->first()->obj()->protein_list()[LAGGING]->nb_elts());
  printf("number of protein : total : %d \n",protein_list_->nb_elts());
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
    while( ( indiv_age < ae_common::individual_evaluation_dates->value( evaluation_index ) ) )
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
    fitness_temp += fitness_;
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
  fitness_ = fitness_temp / (double)ae_common::individual_evaluation_nbr;
}

void ae_individual_R::set_influences()
// Compute the influence of each protein over each coding RNA
// As non-coding RNAs are completely inert, we don't care about their concentration
// so we don't care if proteins activate or inhibit their transcription.
{
  ae_list_node<ae_rna_R*>* rna_node  = NULL;
  ae_rna_R*     rna       = NULL;
  
  //
  rna_node = rna_list_coding_->first();
  while ( rna_node != NULL )
  {
    rna = rna_node->obj();

    //~ printf( "%d proteins\n", protein_list_->nb_elts() );
    rna->set_influences( protein_list_ );

    rna_node = rna_node->next();
  }
}

void ae_individual_R::update_concentrations()
{
  //phenotype_->print_points();

  ae_list_node<ae_protein_R*>* prot_node       = NULL;
  ae_list_node>ae_protein_R*>* next_prot_node  = NULL;
  ae_protein_R* prot            = NULL;

  // Compute all the changes that will be applied to the concentrations
  // Concentrations must not be changed at this stage
  prot_node = protein_list_->first();
  while ( prot_node != NULL )
  {
    prot = prot_node->obj();

    prot->compute_delta_concentration();

    prot_node = prot_node->next();
  }

  // Apply the changes in concentrations we have just computed
  prot_node = protein_list_->first();
  while ( prot_node != NULL )
  {
    prot = prot_node->obj();

    prot->update_concentration();

    // Keep track of the next node in the list, in case we need to remove the current one
    next_prot_node = prot_node->next();

    if( ae_common::with_heredity )
    {
      if( prot->is_inherited() && prot->concentration() < ae_common::protein_presence_limit )
      {
        // The protein has to be removed from both the individual's protein_list_ and inherited_proteins_list_
        protein_list_->remove( prot_node, DELETE_OBJ /*delete_node*/, NO_DELETE /*delete_obj*/ );
        inherited_protein_list_->remove( prot, DELETE_OBJ /*delete_node*/, DELETE_OBJ /*delete_obj*/ );
      }
    }

    // ready to study the next node
    prot_node = next_prot_node;
  }
}

// Multiply the concentration of each protein by <factor>
void ae_individual_R::multiply_concentrations( double factor )
{
  ae_list_node<ae_protein_R*>* prot_node = protein_list_->first();
  ae_protein_R* prot      = NULL;
  
  while ( prot_node != NULL )
  {
    prot = prot_node->obj();

    prot->multiply_concentration( factor );

    prot_node = prot_node->next();
  }
}

int8_t ae_individual_R::quadon( GeneticUnit* gen_unit, ae_strand strand, int32_t pos )
{
  const char* dna = gen_unit->dna()->data();
  int32_t  len    = gen_unit->dna()->length();
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
  if (ae_common::with_heredity && inherited_protein_list_ != NULL )
  {
    // Write inherited proteins
    int16_t nb_inherited_proteins = inherited_protein_list_->nb_elts();
    gzwrite( backup_file, &nb_inherited_proteins,  sizeof(nb_inherited_proteins) );

    ae_list_node<ae_protein_R*>*  inherited_protein_node = inherited_protein_list_->first();
    ae_protein_R*  inherited_protein;

    for ( int16_t i = 0 ; i < nb_inherited_proteins ; i++ )
    {
    inherited_protein = inherited_protein_node->obj();
    
    inherited_protein->save( backup_file );
    
    inherited_protein_node = inherited_protein_node->next();
    }
  }
}
// =================================================================
//                           Protected Methods
// =================================================================
void ae_individual_R::make_protein_list()
{
  ae_individual::make_protein_list();
  protein_list_->add_list( inherited_protein_list_ );
}

void ae_individual_R::make_rna_list()
{
  ae_individual::make_rna_list();
  
  // Parse the newly created RNA list and copy the coding RNAs in rna_list_coding_.
  ae_list_node<ae_rna*>* rna_node  = rna_list_->first();
  ae_rna*       rna       = NULL;
  
  while ( rna_node != NULL )
  {
    rna = rna_node->obj();
    
    if ( rna->is_coding() == true )
    {
      rna_list_coding_->add( rna );
    }
    
    rna_node = rna_node->next();
  }
}
} // namespace aevol
