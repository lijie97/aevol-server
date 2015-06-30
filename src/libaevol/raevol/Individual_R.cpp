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
#include "Individual_R.h"
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
Individual_R::Individual_R( void ) : Individual()
{

}

Individual_R::Individual_R( Individual_R* parent, int32_t id,
                                  ae_jumping_mt* mut_prng, ae_jumping_mt* stoch_prng )
        : Individual( parent, id, mut_prng, stoch_prng )
{
  //~ printf( "ae_individual_R( parent ) : I have %d inherited proteins\n", parent->get_protein_list()->get_nb_elts() );
  
    for (int i = 0; i < parent->_protein_list.size(); i++)
    {
    	if( parent->_protein_list[i]->get_concentration() > parent->get_exp_m()->get_exp_s()->get_protein_presence_limit() )
    	{
    		Protein_R* inherited_prot = new Protein_R( (Protein_R*)parent->_protein_list[i] );
    		inherited_prot->set_inherited( true );
    		_inherited_protein_list.push_back( inherited_prot );
    	}
    }
}

Individual_R::Individual_R( gzFile backup_file ) : Individual( backup_file )
{
  if( get_exp_m()->get_exp_s()->get_with_heredity() )
  {
    // Retreive inherited proteins
//    _inherited_protein_list = new ae_list();
    int16_t nb_inherited_proteins = 0;
    gzread( backup_file, &nb_inherited_proteins,  sizeof(nb_inherited_proteins) );
  
    for ( int16_t i = 0 ; i < nb_inherited_proteins ; i++ )
    {
	  _inherited_protein_list.push_back( new Protein_R( backup_file ) );
    }
  }  
}

// =================================================================
//                             Destructors
// =================================================================
Individual_R::~Individual_R( void )
{
  assert( !get_exp_m()->get_exp_s()->get_with_heredity()  );
  
  for (int i = 0; i < _inherited_protein_list.size(); i++) {
	  Protein_R* dp = _inherited_protein_list[i];
	  delete dp;
  }

  _inherited_protein_list.clear();

  _rna_list_coding.clear();
}

// =================================================================
//                            Public Methods
// =================================================================
void Individual_R::evaluate( Environment* envir )
{
	  //protections
	//  if ( _nb_env_list == 0 )
	//  {
	//    printf( "ERROR in ae_individual_R::evaluate env_list is empty\n " );
	//    exit( EXIT_FAILURE );
	//  }

	  // la protection est génante dans le cas de modification environnementale au cours de la simulation (overload)
	  /*
	  if (env_list->get_nb_elts() != ae_common::individual_environment_nbr)
	  {
	    printf( "ERROR in ae_individual_R::evaluate the number of elements in env_list does not match wich ae_common::individual_environment_nbr\n " );
	    exit( EXIT_FAILURE );
	  }
	  */
	  //debug
	  //printf("évaluation de l'individu %d\n", _index_in_population);

	  // ---------------------------------------------------------------------------
	  // 1) Transcription - Translation - Folding - make_protein_list
	  // ---------------------------------------------------------------------------

	  _transcribed = false;
	  _translated = false;
	  _folded = false;


	  //  printf("avant : do_transcription_translation_folding() \n");
	  //contain make protein list
	  do_transcription_translation_folding();

	  // printf("taille liste _rna_list_coding : %d\n",_rna_list_coding->get_nb_elts());
	  //printf("taille liste _protein_list : %d\n",_protein_list->get_nb_elts());


	  if(_phenotype != NULL)
	  {
	    delete _phenotype;
	    _phenotype = NULL;
	  }
	  _phenotype = new Phenotype();


	  //----------------------------------------------------------------------------
	  // 2) Make a list of all the rna present in the individual
	  //    and initialise the concentrations of the proteins
	  //----------------------------------------------------------------------------
	#ifdef __TRACING__
	  high_resolution_clock::time_point t1 = high_resolution_clock::now();
	#endif
	  make_rna_list();


	//  ae_list_node* prot_node       = NULL;
	//  ae_protein_R* prot            = NULL;

	//  prot_node = _protein_list->get_first();
	//  while ( prot_node != NULL )
	//  {
	//    prot = (ae_protein_R*)prot_node->get_obj();
	  for (int i = 0; i < _protein_list.size(); i++) {
		  ((ae_protein_R*)_protein_list[i])->set_initial_concentration();
	//    prot_node = prot_node->get_next();
	  }
	#ifdef __TRACING__
	  high_resolution_clock::time_point t2 = high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	  ae_logger::addLog(MAKERNALIST,duration);
	  t1 = t2;
	#endif
	  //----------------------------------------------------------------------------
	  // 3) Create influence graph (including the signals)
	  //----------------------------------------------------------------------------
	//	ae_list* init_prot_list = new ae_list;
	//  init_prot_list->add_list(_protein_list);
	  ae_environment* envir = NULL;
	  std::vector<ae_protein*> init_protein = _protein_list;
	  std::vector<int> env_switch(env_list.size());
	  std::map<int,std::vector<ae_protein_R*>*> _cloned_signals;
	  ae_protein_R* cloned_signal = NULL;

	  for(int i = 0; i < env_list.size(); i++)
	  {
	     envir = env_list[i];
	     env_switch[i] = envir->get_id();
	//     printf("Env at age %d : %d (%d)\n",i,envir->get_id(),env_list.size());

	//     if (i > 0) printf("switch at %d\n",ae_common::individual_environment_dates->get_value(i-1));

	     if (_cloned_signals.find(envir->get_id()) == _cloned_signals.end()) {
	    	 std::vector<ae_protein_R*>* loc_signals = new std::vector<ae_protein_R*>();
	    	 for ( int8_t j = 0; j < envir->get_signals().size(); j++)
	    	 {
	    			 cloned_signal = new ae_protein_R(NULL,*(envir->_signals[j]));
	    			 cloned_signal->set_concentration(0.0);
	    			 _protein_list.push_back(cloned_signal);
	    			 loc_signals->push_back(cloned_signal);
	    	 }
	    	 _cloned_signals.insert(std::make_pair(envir->get_id(),loc_signals));
	//    	 printf("Cloned signals for %d is %ld (or %ld)\n",
	//    			 envir->get_id(),
	//    			 loc_signals->size(),
	//    			 _cloned_signals[envir->get_id()]->size());
	     }
	  }

	  cloned_signal = NULL;
	#ifdef __TRACING__
	  t2 = high_resolution_clock::now();
	  duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	  ae_logger::addLog(CLONESIGNAL,duration);
	  t1 = t2;
	#endif

	  set_influences();

	#ifdef __TRACING__
	  t2 = high_resolution_clock::now();
	  duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	  ae_logger::addLog(SETINFLUS,duration);
	  t1 = t2;
	#endif
	  //----------------------------------------------------------------------------
	  // 4) Make the individual "live its life" and compute partial phenotypes and
	  //    fitnesses
	  //----------------------------------------------------------------------------

	  int16_t indiv_age     = 0;
	  double dist_temp      = 0;
	  int8_t compteur_env   = 0;
	  //set envir to the first environment
	  envir = env_list[compteur_env];
	//  for ( auto it = _cloned_signals.begin(); it != _cloned_signals.end(); ++it ) {
	//  	for ( int8_t i = 0; i < it->second->size(); i++)
	//  	      {
	//  	    	  printf("-- Concentration %d %d %d : %f\n",indiv_age,
	//  	    			  it->first,i,it->second->at(i)->get_concentration());
	//  	      }
	//  }
	  //Set the concentrations of the signals proteins
	//  _signals = envir->get_signals();
	  for ( int8_t i = 0; i < _cloned_signals[env_switch[compteur_env]]->size(); i++)
	  {
	    ((ae_protein_R*) _cloned_signals[env_switch[compteur_env]]->at(i))->set_concentration(0.9);
	  }


	  while( indiv_age < ae_common::individual_evaluation_dates->get_value( ae_common::individual_evaluation_nbr - 1) )
	  {
	#ifdef __TRACING__
		  t1 = high_resolution_clock::now();
	#endif
	//	printf("ENV %d : %d / %d %d\n",indiv_age,env_switch[compteur_env],compteur_env,env_list.size());
	    //Updating the concentrations in order to respect the degradation step.
	    for( int16_t i = 0; i < 1/ae_common::degradation_step; i++ )
	    {
	      update_concentrations();
	    }
	    indiv_age++;

	//    for ( auto it = _cloned_signals.begin(); it != _cloned_signals.end(); ++it ) {
	//    	for ( int8_t i = 0; i < it->second->size(); i++)
	//    	      {
	//    	    	  printf("Concentration %d %d %d : %f\n",indiv_age,
	//    	    			  it->first,i,it->second->at(i)->get_concentration());
	//    	      }
	//    }
	#ifdef __TRACING__
	    t2 = high_resolution_clock::now();
	    duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	    ae_logger::addLog(UPDATECONCENT,duration);
	    t1 = t2;
	#endif
	    if( ae_common::individual_evaluation_dates->search(indiv_age) != -1)
	    {
	      // Evaluate the individual's phenotype
	//    	printf("START Update Phenotype %d\n",this->_index_in_population);
	      update_phenotype();
	//      printf("END Update Phenotype %d\n",this->_index_in_population);
	      _distance_to_target_computed = false;
	      _phenotype_computed = true;
	//      printf("START Compute DIST %d\n",this->_index_in_population);
	      compute_distance_to_target( envir ); /// <<<<STILL ISSUE WITH ENVIR
	//      printf("END Compute DIST %d\n",this->_index_in_population);
	      dist_temp += _dist_to_target_by_feature[METABOLISM];
	      //  printf("indiv_age : %d dist_to_target : %lf\n",indiv_age,_dist_to_target_by_feature[METABOLISM]);
	    }

	#ifdef __TRACING__
	    t2 = high_resolution_clock::now();
	    duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	    ae_logger::addLog(UPDATEPHENO,duration);
	    t1 = t2;
	#endif

	    if( ae_common::individual_environment_dates->search(indiv_age) != -1)
	    {
	      //Remove the signals of this environment
	      for ( int8_t i = 0; i < _cloned_signals[env_switch[compteur_env]]->size(); i++)
	      {
	    	  _cloned_signals[env_switch[compteur_env]]->at(i)->set_concentration(0.);
	      }

	      // Change the environment at is next value
	      compteur_env+=1;
	      envir = env_list[compteur_env];

	      // Add the signals of this new environment
	//      _signals = envir->get_signals();
	      for ( int8_t i = 0; i < _cloned_signals[env_switch[compteur_env]]->size(); i++)
	      {
	    	  _cloned_signals[env_switch[compteur_env]]->at(i)->set_concentration(0.9);
	      }
	    }

	#ifdef __TRACING__
	    t2 = high_resolution_clock::now();
	    duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	    ae_logger::addLog(CHANGEENV,duration);
	    t1 = t2;
	#endif
	  }

	  //Remove the signals of the last environment
	  for ( int8_t i = 0; i < _cloned_signals[env_switch[compteur_env]]->size(); i++)
	  {
		  _cloned_signals[env_switch[compteur_env]]->at(i)->set_concentration(0.0);
	  }


	  //----------------------------------------------------------------------------
	  // 5) Compute final fitness and final dist to target
	  //----------------------------------------------------------------------------
	#ifdef __TRACING__
	  t1 = high_resolution_clock::now();
	#endif
	  // On devrait faire la somme du carré des erreurs afin d'éviter qu'elles puissent se compenser
	  _dist_to_target_by_feature[METABOLISM] = dist_temp / (double)ae_common::individual_evaluation_nbr;
	  _fitness_computed=false;
	  // yoram attention il peut y avoir des soucis si on utilise des environnements segmentés ici
	  compute_fitness( envir );
	  // set _protein list to its standard value in order to not disturb...

	  // TODO Delete all the signals !
	  for ( auto it = _cloned_signals.begin(); it != _cloned_signals.end(); ++it )
	//  for(int i = 0; i < env_list.size(); i++)
	  {
	//	  printf("Cloned signals for %d is %ld\n",it->first,((std::vector<ae_protein_R*>) it->second).size());
		  for (int j = 0; j < ((std::vector<ae_protein_R*>*) it->second)->size(); j++) {
			  ae_protein_R *to_delete = ((std::vector<ae_protein_R*>*) it->second)->at(j);
	//		  printf("Cloned DELETE signal %p (%d)\n",to_delete,it->first);
			  delete to_delete;
		  }
		  ((std::vector<ae_protein_R*>) it->first).clear();
	  }
	  _cloned_signals.clear();

	  _protein_list.clear();
	  _protein_list = init_protein;
	  init_protein.clear();
	//  _protein_list->add_list(init_prot_list);
	//  delete init_prot_list;
	//  init_prot_list = NULL;
	  _phenotype_computed = true;

	//  printf("Fitness %f\n",this->get_fitness());

	#ifdef __TRACING__
	  t2 = high_resolution_clock::now();
	    duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	    ae_logger::addLog(FITNESS,duration);
	    t1 = t2;
	#endif
}

void Individual_R::set_influences( void )
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

void Individual_R::update_concentrations( void )
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
void Individual_R::multiply_concentrations( double factor )
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

int8_t Individual_R::get_quadon( GeneticUnit* gen_unit, Strand strand, int32_t pos )
{
  const char* dna = gen_unit->get_dna()->get_data();
  int32_t  len    = gen_unit->get_dna()->get_length();
  int8_t quadon   = 0;

  if ( strand == LEADING )
  {
    for ( int8_t i = 0 ; i < QUADON_SIZE ; i++ )
    {
      if ( dna[(pos+i) % len] == '1' )
      {
        quadon += 1 << (QUADON_SIZE - i - 1);  //pow( 2, QUADON_SIZE - i - 1 );
      }
    }
  }
  else  // ( strand == LAGGING )
  {
    for ( int8_t i = 0 ; i < QUADON_SIZE ; i++ )
    {
      if ( dna[(pos-i) % len] != '1' ) // == and not != because we are on the complementary strand...
      {
        quadon += 1 << (QUADON_SIZE - i - 1);  //pow( 2, QUADON_SIZE - i - 1 );
      }
    }
  }

  return quadon;
}

void Individual_R::save( gzFile backup_file )
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
void Individual_R::make_protein_list( void )
{
  Individual::make_protein_list();
  _protein_list->add_list( _inherited_protein_list );
}

void Individual_R::make_rna_list( void )
{
  Individual::make_rna_list();
  
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
} // namespace aevol
