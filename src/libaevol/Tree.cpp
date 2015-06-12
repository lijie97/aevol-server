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
#include "Tree.h"

#include "macros.h"
#include "ExpManager.h"
#include "ExpSetup.h"
#include "Individual.h"
#include "ae_utils.h"


namespace aevol {


//##############################################################################
//                                                                             #
//                                Class Tree                                #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================
const int32_t Tree::NO_PARENT = -1;

// =================================================================
//                             Constructors
// =================================================================
Tree::Tree(ExpManager * exp_m, ae_tree_mode tree_mode, int64_t tree_step)
{
  _exp_m = exp_m;

  _tree_mode = tree_mode;
  _tree_step = tree_step;

  switch ( _tree_mode )
  {
    case NORMAL :
    {
      _nb_indivs    = new int32_t [_tree_step];
      _replics      = new ReplicationReport ** [_tree_step];

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



Tree::Tree( ExpManager * exp_m, char* tree_file_name )
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

      ReplicationReport * replic_report = NULL;

      _nb_indivs    = new int32_t[_tree_step];
      _replics      = new ReplicationReport **[_tree_step];

      gzread( tree_file, _nb_indivs, _tree_step * sizeof(_nb_indivs[0]) );

      for ( int64_t t = 0 ; t < _tree_step ; t++ )
      {
        _replics[t] = new ReplicationReport *[_nb_indivs[t]];
        for ( int32_t indiv_i = 0 ; indiv_i < _nb_indivs[t] ; indiv_i++ )
        {
          // Retreive a replication report
          replic_report = new ReplicationReport( tree_file, NULL );

          // Put it at its rightful position
          _replics[t][replic_report->get_id()] = replic_report;
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
Tree::~Tree( void )
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



int32_t Tree::get_nb_indivs(int64_t t) const
{
  return _nb_indivs[ae_utils::mod(t - 1, _tree_step)];
}


ReplicationReport *Tree::get_report_by_index(int64_t t, int32_t index) const
{
  assert( _tree_mode == NORMAL );

  return _replics[ae_utils::mod(t - 1, _tree_step)][index];
}


ReplicationReport *Tree::get_report_by_rank(int64_t t, int32_t rank) const
{
  assert( _tree_mode == NORMAL );
  int32_t nb_indivs = get_nb_indivs(t);
  assert( rank <= nb_indivs );

  for ( int32_t i = 0 ; i < nb_indivs ; i++ )
  {
    if ( _replics[ae_utils::mod(t - 1, _tree_step)][i]->get_rank() == rank )
    {
      return _replics[ae_utils::mod(t - 1, _tree_step)][i];
    }
  }

  fprintf( stderr, "ERROR: Couldn't find indiv with rank %" PRId32 " in file %s:%d\n", rank, __FILE__, __LINE__ );
  return NULL;
}


void Tree::set_nb_indivs (int32_t nb_indivs, int64_t t)
{
  _nb_indivs[ae_utils::mod(t - 1, _tree_step)] = nb_indivs;
}


void Tree::fill_tree_with_cur_gener(void)
{
  assert(_exp_m != NULL && Time::get_time() > 0);

  switch (_tree_mode)
  {
    case NORMAL :
    {
      //  -1 because the tree's arrays contain information on
      // generations n*TREE_STEP+1 --> (n+1)*_tree_step
      // (for _tree_step == 100, information on generations
      // 1 to 100, or 101 to 200, or 201 to 300, etc)
      int64_t time_i = ae_utils::mod(Time::get_time() - 1, _tree_step);
      _nb_indivs[time_i] = _exp_m->get_nb_indivs();
      _replics[time_i]   = new ReplicationReport * [_nb_indivs[time_i]];


      for (const auto& indiv: _exp_m->get_indivs()) {
        size_t num_indiv = 0;
        _replics[time_i][num_indiv++] = indiv->get_replic_report();
      }

      break;
    }
    case LIGHT :
    {
      // TO DO

      // TO CHECK !!
      // not sure that time_i should be used in this block...

      // int32_t gener_i     = ae_utils::mod( _exp_m->get_num_gener() - 1, _tree_step );
      // _nb_indivs[gener_i] = _exp_m->get_nb_indivs();
      // _parent[gener_i] = new int32_t [_nb_indivs[gener_i]];

      // int32_t         num_indiv   = 0;
      // for (const auto& indiv: _exp_m->get_indivs_std())
      //   _parent[gener_i][num_indiv++] = indiv->get_replic_report()->get_parent_id();
      
      break;
    }
  }
}

void Tree::write_to_tree_file( gzFile tree_file )
{
  switch ( _tree_mode )
  {
    case NORMAL :
    {
      // Write the tree in the backup
      gzwrite( tree_file, &_nb_indivs[0], _tree_step * sizeof(_nb_indivs[0]) );

      for ( int64_t t = 0 ; t < _tree_step ; t++ )
      {
        for ( int32_t indiv_i = 0 ; indiv_i < _nb_indivs[t] ; indiv_i++ )
        {
          assert(_replics[t][indiv_i] != NULL);
          _replics[t][indiv_i]->write_to_tree_file( tree_file );
        }
      }

      // Reinitialize the tree
      for ( int32_t t = 0 ; t < _tree_step ; t++ )
      {
        for ( int32_t indiv_i = 0 ; indiv_i < _nb_indivs[t] ; indiv_i++ )
        {
          delete _replics[t][indiv_i];
          _replics[t][indiv_i] = NULL;
        }

        delete [] _replics[t];
        _replics[t] = NULL;
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
void Tree::set_replic_report(int32_t id, ReplicationReport * replic_report)
{
  assert(_tree_mode == NORMAL);

  int64_t t = ae_utils::mod(Time::get_time() - 1, _tree_step);

  if (_replics[t] == NULL)
  {
    _replics[t] = new ReplicationReport * [_exp_m->get_nb_indivs()];

    memset(_replics[t], 0, _exp_m->get_nb_indivs() * sizeof(*_replics));
  }

  assert(_replics[t][id] == NULL);
  _replics[t][id] = new ReplicationReport(*replic_report);
}




// CK: Added for aevol_modify
void Tree::set_replic_report(int64_t t, int32_t id, ReplicationReport * replic_report)
{
  assert( _tree_mode == NORMAL );

  t = ae_utils::mod(t - 1, _tree_step);

  if (_replics[t] == NULL)
  {
    _replics[t] = new ReplicationReport * [_exp_m->get_nb_indivs()];
    memset( _replics[t], 0, _exp_m->get_nb_indivs() * sizeof( *_replics ) );
  }

  if ( _replics[t][id] != NULL )
    {
      printf("Erased previous replication report for indiv %" PRId32 " at relative time %" PRId64 "\n", id, t);
      delete _replics[t][id];
    }
  _replics[t][id] = new ReplicationReport(*replic_report);
  _replics[t][id]->set_id(id);

  // debug
  // printf("Added replication report for indiv %d at relative time %d :\n", id, t);
  // printf("  ID            %d \n", _replics[t][id]->get_id() );
  // printf("  Rank          %d \n", _replics[t][id]->get_rank() );
  // printf("  Genome size   %d \n", _replics[t][id]->get_genome_size() );
  // printf("  P. ID         %d \n", _replics[t][id]->get_parent_id() );
  // printf("  P. met.err.   %f \n", _replics[t][id]->get_parent_metabolic_error() );
  // printf("  P. size       %d \n", _replics[t][id]->get_parent_genome_size() );
}




// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol