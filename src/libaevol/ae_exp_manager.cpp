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


/*! \class ae_exp_manager
    \brief This class allows to manage an experiment.
    
    An experiment manager allows to... manage an experiment.
    It owns a population and an experimental_setup that can be loaded from a pair of aevol binary files (pop and exp_setup)
*/
 
 
// =================================================================
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================
#include <ae_exp_manager.h>




//##############################################################################
//                                                                             #
//                             Class ae_exp_manager                            #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_exp_manager::ae_exp_manager( void )
{
  _num_gener    = 0;
  _first_gener  = 0;
  _last_gener   = 0;
  
  _pop      = new ae_population( this );
  _exp_s    = new ae_exp_setup( this );
  _output_m = new ae_output_manager( this );
  
  _quit_signal_received = false;
}

// =================================================================
//                             Destructors
// =================================================================
ae_exp_manager::~ae_exp_manager( void )
{
  delete _pop;
  delete _exp_s;
  delete _output_m;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_exp_manager::save_experiment( void ) const
{
  _output_m->save_experiment();
}

void ae_exp_manager::load_experiment( char* pop_file_name, char* exp_setup_file_name, char* out_man_file_name, bool verbose )
{
  _output_m->load_experiment( pop_file_name, exp_setup_file_name, out_man_file_name, verbose );

  // Evaluate individuals
  _pop->evaluate_individuals( _exp_s->get_env() );

  // If the population is spatially structured, then the individuals are saved and loaded in the order of the grid and not in increasing order of fitness
  // so we have to sort the individuals
  if ( is_spatially_structured() )
  {
    _pop->sort_individuals();
  }


  // If the simulation is being continued (not just post-processed), prepare output data accordingly
  //~ if ( to_be_run )
  //~ {
    //~ // Prepare stat files
    //~ _stats  = new ae_stats( _first_gener );
    
    //~ // Prepare tree
    //~ if ( ae_common::rec_params->get_record_tree() == true )
    //~ { 
      //~ mkdir( "tree", 0755 );
      //~ _tree  = new ae_tree(); 
    //~ }
    //~ else
    //~ {
      //~ _tree = NULL;
    //~ }
    
    //~ // Prepare dump
    //~ if ( ae_common::rec_params->get_dump_period() > 0 )
    //~ {
      //~ mkdir( "dump", 0755 );
      //~ _dump = new ae_dump();
    //~ }
    //~ else
    //~ {
      //~ _dump = NULL;
    //~ }
    
    //~ mkdir( "backup", 0755 );
    
    //~ if ( ae_common::rec_params->is_logged( LOG_LOADS ) == true )
    //~ {
      //~ // Write an entry in the LOADS log file
      //~ fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "GENERATION_OVERLOAD %"PRId32"\n", _num_gener );
      //~ if ( param_overloader->get_nb_overloaded() > 0 )
      //~ {
        //~ //fprintf( _logs->get_log( LOG_LOADS ), "  Overloaded parameters:\n" );
        //~ param_overloader->write_log( ae_common::rec_params->get_log( LOG_LOADS ) );
        //~ fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "\n" );
      //~ }
      //~ else
      //~ {
        //~ fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "  No overloaded parameters\n\n" );
      //~ }
    //~ }
  //~ }
  //~ else
  //~ {
    //~ // We just want to inspect the state of the simulation at this moment
    //~ ae_common::rec_params->init_logs( 0 );
  //~ }
  
  // Initialize display
  #ifdef __X11
    //~ ((ae_exp_setup_X11*) this)->initialize( ae_common::pop_structure, ae_common::params->get_allow_plasmids() );
  #endif // def __X11
}

void ae_exp_manager::run_evolution( void )
{
  // dump the initial state of the population; useful for restarts
  _output_m->write_current_generation_outputs();
  
  while ( _num_gener < _last_gener )
  {
    printf( "============================== %"PRId32" ==============================\n", _num_gener );
    printf( "  distance to target (metabolic) : %f\n", ((ae_individual *) _pop->get_indivs()->get_last()->get_obj())->get_dist_to_target_by_feature( METABOLISM ) );

    if ( quit_signal_received() ) break;
    
    #ifdef __X11
      display();
    #endif
    
    // Take one step in the evolutionary loop
    _exp_s->step_to_next_generation();
    _num_gener++;

    // Write statistical data and store phylogenetic data (tree)
    _output_m->write_current_generation_outputs();
  }
  
  _output_m->flush();
  
  printf( "============================== %"PRId32" ==============================\n", _num_gener );
  printf( "  distance to target (metabolic) : %f\n", ((ae_individual *) _pop->get_indivs()->get_last()->get_obj())->get_dist_to_target_by_feature( METABOLISM ) );
  printf( "===================================================================\n");
  printf ("  The run is finished. \n"); 
  printf ("  Printing the final best individual into best_last_org.txt \n"); 
  const char* out_file_name = "best_last_org.txt"; 
  FILE* org_file = fopen( out_file_name, "w" );
  fputs( ((ae_individual *) _pop->get_indivs()->get_last()->get_obj())->get_genetic_unit(0)->get_dna()->get_data(), org_file); 
  fclose ( org_file ); 
}

// =================================================================
//                           Protected Methods
// =================================================================

// =================================================================
//                          Non inline accessors
// =================================================================
