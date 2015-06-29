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
#include "Utils.h"


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
Tree::Tree(ExpManager* exp_m, TreeMode tree_mode, int64_t tree_step)
{
  _exp_m = exp_m;

  _tree_mode = tree_mode;
  _tree_step = tree_step;

  switch ( _tree_mode )
  {
    case NORMAL :
    {
      _replics = new ReplicationReport** [_tree_step];

      for (int32_t time = 0 ; time < tree_step ; time++) {
        _replics[time] = new ReplicationReport* [exp_m->get_nb_indivs()];
        for (int32_t num_indiv = 0 ;
             num_indiv < exp_m->get_nb_indivs() ;
             num_indiv++) {
          _replics[time][num_indiv] = new ReplicationReport();
        }
      }

      break;
    }
    case LIGHT :
    {
      // TODO

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


/**
 *
 */
Tree::Tree(ExpManager * exp_m, char* tree_file_name)
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

      _replics = new ReplicationReport** [_tree_step];

      for ( int64_t t = 0 ; t < _tree_step ; t++ )
      {
        _replics[t] = new ReplicationReport* [_exp_m->get_nb_indivs()];
        for (int32_t indiv_i = 0 ;
             indiv_i < _exp_m->get_nb_indivs() ;
             indiv_i++) {
          // Retrieve a replication report
          ReplicationReport* replic_report =
              new ReplicationReport(tree_file, nullptr);

          // Put it at its rightful position
          _replics[t][replic_report->id()] = replic_report;
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
            for ( int32_t j = 0 ; j < _exp_m->get_nb_indivs() ; j++ )
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
}

// =================================================================
//                            Public Methods
// =================================================================
ReplicationReport** Tree::get_reports(int64_t t) const {
  assert(_tree_mode == NORMAL);
  return _replics[Utils::mod(t - 1, _tree_step)];
}

ReplicationReport* Tree::get_report_by_index(int64_t t, int32_t index) const {
  assert(_tree_mode == NORMAL);
  return _replics[Utils::mod(t - 1, _tree_step)][index];
}


ReplicationReport* Tree::get_report_by_rank(int64_t t, int32_t rank) const {
  assert(_tree_mode == NORMAL);
  int32_t nb_indivs = _exp_m->get_nb_indivs();
  assert(rank <= nb_indivs);

  for (int32_t i = 0 ; i < nb_indivs ; i++) {
    if (_replics[Utils::mod(t - 1, _tree_step)][i]->rank() == rank) {
      return _replics[Utils::mod(t - 1, _tree_step)][i];
    }
  }

  fprintf(stderr,
          "ERROR: Couldn't find indiv with rank %" PRId32 " in file %s:%d\n",
          rank, __FILE__, __LINE__);
  return NULL;
}

void Tree::signal_end_of_generation() {
  auto cur_reports = get_reports(Time::get_time());
  for (int32_t i = 0; i < _exp_m->get_nb_indivs(); i++) {
    cur_reports[i]->signal_end_of_generation();
  }
}

void Tree::write_to_tree_file(gzFile tree_file)
{
  switch ( _tree_mode )
  {
    case NORMAL :
    {
      // Write the tree in the backup
      for ( int64_t t = 0 ; t < _tree_step ; t++ )
      {
        for ( int32_t indiv_i = 0 ; indiv_i < _exp_m->get_nb_indivs() ; indiv_i++ )
        {
          assert(_replics[t][indiv_i] != NULL);
          _replics[t][indiv_i]->write_to_tree_file( tree_file );
        }
      }

      // Reinitialize the tree
      for ( int32_t t = 0 ; t < _tree_step ; t++ )
      {
        for ( int32_t indiv_i = 0 ; indiv_i < _exp_m->get_nb_indivs() ; indiv_i++ )
        {
          delete _replics[t][indiv_i];
          _replics[t][indiv_i] = new ReplicationReport();
        }
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

void Tree::update(Observable& o, void* arg) {
  auto indivs = reinterpret_cast<Individual**>(arg);
  get_report_by_index(Time::get_time(), indivs[0]->get_id())->init(indivs[0],
                                                                   indivs[1]);
}



// =================================================================
//                  Non-inline accessors' definition
// =================================================================

// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
