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
#include <ae_tree.h>

#include <ae_macros.h>
#include <ae_exp_manager.h>
#include <ae_exp_setup.h>
#include <ae_population.h>
#include <ae_individual.h>
#include <ae_utils.h>




//##############################################################################
//                                                                             #
//                                Class ae_tree                                #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================
const int32_t ae_tree::NO_PARENT = -1;

// =================================================================
//                             Constructors
// =================================================================
ae_tree::ae_tree( ae_exp_manager* exp_m, ae_tree_mode tree_mode, int32_t tree_step )
{
  _exp_m = exp_m;
    
  _tree_mode = tree_mode;
  _tree_step = tree_step;

  switch ( _tree_mode )
  {
    case NORMAL :
    {
      _nb_indivs    = new int32_t [_tree_step];
      _replics      = new ae_replication_report** [_tree_step];
      
      // All pointers in the _replics table must be set to NULL, otherwise
      // the destructor won't work properly if called before the matrix was 
      // filled 
      for ( int32_t i = 0 ; i < tree_step ; i++ )
      {
        _replics[i] = NULL;
      }
      
      break;
    }
    case LIGHT :
    {
      // TO DO

      //~ // Creates the big arrays to store number of individuals and child->parent relations
      //~ _nb_indivs  = new int32_t[ae_common::nb_generations];
      //~ _parent     = new int32_t*[ae_common::nb_generations];
      
      //~ _parent[0]    = new int32_t[ae_common::init_params->get_init_pop_size()];
      
      //~ _nb_indivs[0] = ae_common::init_params->get_init_pop_size();
      
      //~ // Individuals from the initial population don't have parents
      //~ for ( int32_t i = 0 ; i < _nb_indivs[0] ; i++ )
      //~ {
        //~ _parent[0][i] = NO_PARENT;
      //~ }
      
      break;
    }
  }
}



ae_tree::ae_tree( ae_exp_manager* exp_m, char* tree_file_name )
{
  _exp_m = exp_m;

  _tree_mode = _exp_m->get_tree_mode();
  _tree_step = _exp_m->get_tree_step();
  switch ( _tree_mode )
  {
    case NORMAL :
    {
      gzFile tree_file = gzopen( tree_file_name, "r" );
      if ( tree_file == Z_NULL )
      {
        printf( "ERROR : Could not read tree file %s\n", tree_file_name );
        exit( EXIT_FAILURE );
      }
      
      ae_replication_report * replic_report = NULL;
            
      _nb_indivs    = new int32_t[_tree_step];
      _replics      = new ae_replication_report**[_tree_step];      
      
      gzread( tree_file, _nb_indivs, _tree_step * sizeof(_nb_indivs[0]) );

      for ( int32_t gener_i = 0 ; gener_i < _tree_step ; gener_i++ )
      {
        _replics[gener_i] = new ae_replication_report*[_nb_indivs[gener_i]];
        for ( int32_t indiv_i = 0 ; indiv_i < _nb_indivs[gener_i] ; indiv_i++ )
        {
          // Retreive a replication report
          replic_report = new ae_replication_report( tree_file, NULL );
          
          // Put it at its rightful position
          _replics[gener_i][replic_report->get_id()] = replic_report;
        }
      }
      gzclose( tree_file );
       
      break;
    }
    case LIGHT :
    {
      // TODO
      break;
    }
  }
}




// =================================================================
//                             Destructors
// =================================================================
ae_tree::~ae_tree( void )
{
  switch ( _tree_mode )
  {
    case NORMAL :
    {
      if ( _replics != NULL )
      {
        for ( int32_t i = 0 ; i < _tree_step ; i++ )
        {
          if ( _replics[i] != NULL )
          {
            for ( int32_t j = 0 ; j < _nb_indivs[i] ; j++ )
            {
              delete _replics[i][j];
            }
            
            delete [] _replics[i];
          }
        }
        delete [] _replics;
      }

      break;
    }
    case LIGHT :
    {
      //~ for( int32_t n = 0 ; n < ae_common::nb_generations ; n++ )
      //~ {
        //~ delete _parent[n];
      //~ }
      
      //~ delete [] _parent;
      
      break;
    }
  }

  delete [] _nb_indivs;
}

// =================================================================
//                            Public Methods
// =================================================================



int32_t ae_tree::get_nb_indivs( int32_t generation ) const
{
  return _nb_indivs[ae_utils::mod(generation - 1, _tree_step)];
}


ae_replication_report * ae_tree::get_report_by_index( int32_t generation, int32_t index ) const
{
  assert( _tree_mode == NORMAL );
  
  return _replics[ae_utils::mod(generation - 1, _tree_step)][index];
}


ae_replication_report * ae_tree::get_report_by_rank( int32_t generation, int32_t rank ) const
{
  assert( _tree_mode == NORMAL );
  int32_t nb_indivs = get_nb_indivs( generation );
  assert( rank <= nb_indivs );
  
  for ( int32_t i = 0 ; i < nb_indivs ; i++ )
  {
    if ( _replics[ae_utils::mod(generation - 1, _tree_step)][i]->get_rank() == rank )
    {
      return _replics[ae_utils::mod(generation - 1, _tree_step)][i];
    }
  }
  
  fprintf( stderr, "ERROR: Couldn't find indiv with rank %" PRId32 " in file %s:%d\n", rank, __FILE__, __LINE__ );
  return NULL;
}


void ae_tree::set_nb_indivs (int32_t nb_indivs, int32_t generation)
{
  _nb_indivs[ae_utils::mod(generation - 1, _tree_step)] = nb_indivs;
}


void ae_tree::fill_tree_with_cur_gener( void )
{
  assert( _exp_m != NULL && _exp_m->get_num_gener() > 0 );
  
  switch ( _tree_mode )
  {
    case NORMAL :
    {
      //  -1 because the tree's arrays contain information on 
      // generations n*TREE_STEP+1 --> (n+1)*_tree_step
      // (for _tree_step == 100, information on generations
      // 1 to 100, or 101 to 200, or 201 to 300, etc)
      int32_t gener_i     = ae_utils::mod( _exp_m->get_num_gener() - 1, _tree_step );
      _nb_indivs[gener_i] = _exp_m->get_pop()->get_nb_indivs();
      _replics[gener_i]   = new ae_replication_report* [_nb_indivs[gener_i]];
      
      
      ae_list_node<ae_individual*>* indiv_node = _exp_m->get_indivs()->get_first();
      ae_individual*  indiv       = NULL;
      int32_t         num_indiv   = 0;
      
      while ( indiv_node != NULL )
      {
        indiv = indiv_node->get_obj();
        
        _replics[gener_i][num_indiv++] = indiv->get_replic_report();

        indiv_node = indiv_node->get_next();
      }
      
      break;
    }
    case LIGHT :
    {
      // TO DO

      // TO CHECK !!
      // not sure that gener_i should be used in this block...

      // int32_t gener_i     = ae_utils::mod( _exp_m->get_num_gener() - 1, _tree_step );
      // _nb_indivs[gener_i] = _exp_m->get_nb_indivs();
      // _parent[gener_i] = new int32_t [_nb_indivs[gener_i]];
      
      // ae_list_node<ae_individual*>* indiv_node = _exp_m->get_indivs()->get_first();
      // ae_individual*  indiv       = NULL;
      // int32_t         num_indiv   = 0;
      
      // while( indiv_node != NULL )
      // {
      //   indiv = indiv_node->get_obj();
        
      //   _parent[gener_i][num_indiv++] = indiv->get_replic_report()->get_parent_id();
        
      //   indiv_node = indiv_node->get_next();
      // }
      
      break;
    }
  }
}

void ae_tree::write_to_tree_file( gzFile tree_file )
{
  switch ( _tree_mode )
  {
    case NORMAL :
    {
      // Write the tree in the backup
      gzwrite( tree_file, &_nb_indivs[0], _tree_step * sizeof(_nb_indivs[0]) );
      
      for ( int32_t gener_i = 0 ; gener_i < _tree_step ; gener_i++ )
      {
        for ( int32_t indiv_i = 0 ; indiv_i < _nb_indivs[gener_i] ; indiv_i++ )
        {
          assert(_replics[gener_i][indiv_i] != NULL);
          _replics[gener_i][indiv_i]->write_to_tree_file( tree_file );
        }
      }
      
      // Reinitialize the tree
      for ( int32_t gener_i = 0 ; gener_i < _tree_step ; gener_i++ )
      {
        for ( int32_t indiv_i = 0 ; indiv_i < _nb_indivs[gener_i] ; indiv_i++ )
        {
          delete _replics[gener_i][indiv_i];
          _replics[gener_i][indiv_i] = NULL;
        }
      
        delete [] _replics[gener_i];
        _replics[gener_i] = NULL;
      }
      
      break;
    }
    case LIGHT :
    {
      // TODO : ?
      break;
    }
  }
}



// =================================================================
//                  Non-inline accessors' definition
// =================================================================
void ae_tree::set_replic_report( int32_t id, ae_replication_report* replic_report )
{
  assert( _tree_mode == NORMAL );
  
  int32_t gener_i = ae_utils::mod(_exp_m->get_num_gener() - 1, _tree_step); // CK: BUGFIX. Previous expression was: _exp_m->get_num_gener() % _tree_step;
  
  if ( _replics[gener_i] == NULL )
  {
    _replics[gener_i] = new ae_replication_report* [_exp_m->get_nb_indivs()];
    
    memset( _replics[gener_i], 0, _exp_m->get_nb_indivs() * sizeof( *_replics ) );  
  }

  assert( _replics[gener_i][id] == NULL );
  _replics[gener_i][id] = new ae_replication_report(*replic_report);
}




// CK: Added for aevol_modify
void ae_tree::set_replic_report( int32_t generation, int32_t id, ae_replication_report* replic_report )
{
  assert( _tree_mode == NORMAL );

  int32_t g = ae_utils::mod(generation - 1, _tree_step);

  if ( _replics[g] == NULL )
  {
    _replics[g] = new ae_replication_report* [_exp_m->get_nb_indivs()];
    memset( _replics[g], 0, _exp_m->get_nb_indivs() * sizeof( *_replics ) );  
  }

  if ( _replics[g][id] != NULL )
    {
      printf("Erased previous replication report for indiv %d of generation %d (%d)\n", id, generation, g);
      delete _replics[g][id];
    }
  _replics[g][id] = new ae_replication_report(*replic_report);
  _replics[g][id]->set_id(id);

  // debug
  printf("Added replication report for indiv %d of generation %d (%d) :\n", id, generation, g);
  printf("  ID            %d \n", _replics[g][id]->get_id() );
  printf("  Rank          %d \n", _replics[g][id]->get_rank() );
  printf("  Genome size   %d \n", _replics[g][id]->get_genome_size() );
  printf("  P. ID         %d \n", _replics[g][id]->get_parent_id() );
  printf("  P. met.err.   %f \n", _replics[g][id]->get_parent_metabolic_error() );
  printf("  P. size       %d \n", _replics[g][id]->get_parent_genome_size() );
}




// =================================================================
//                           Protected Methods
// =================================================================
