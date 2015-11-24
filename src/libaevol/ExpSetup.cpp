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
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <set>

// =================================================================
//                            Project Files
// =================================================================
#include "ExpSetup.h"
#include "JumpingMT.h"


#include "GzHelpers.h"

namespace aevol {
//##############################################################################
//                                                                             #
//                              Class ExpSetup                                 #
//                                                                             #
//##############################################################################

// ===========================================================================
//                         Definition of static attributes
// ===========================================================================

// ===========================================================================
//                                  Constructors
// ===========================================================================
ExpSetup::ExpSetup( ExpManager * exp_m )
{
  _exp_m  = exp_m;
  
  // -------------------------------------------------------------- Selection
  _sel = new Selection( exp_m );
  
  // --------------------------------------------------------------- Transfer
  _with_HT                    = false;
  _repl_HT_with_close_points  = false;
  _HT_ins_rate                = 0.0;
  _HT_repl_rate               = 0.0;
  _repl_HT_detach_rate         = 0.0;
  
  // --------------------------------------------------------------- Plasmids
  _with_plasmids    = false;
  _prob_plasmid_HT  = 0.0;
  _tune_donor_ability     = 0.0;
  _tune_recipient_ability = 0.0;
  _donor_cost       = 0.0;
  _recipient_cost   = 0.0;
  _swap_GUs         = false;
  
  // -------------------------------------------------------------- Secretion
  _with_secretion = false;
  _secretion_contrib_to_fitness = 0.0;
  _secretion_cost               = 0.0;

  _fuzzy_flavor                 = 0;

#ifdef __REGUL
  _protein_presence_limit = 1e-2;
  _degradation_rate  = 1;
  _degradation_step  = 10;
  _with_heredity          = false;
  _nb_indiv_age      = 20*_degradation_step;
  _eval_step         = 5;

  _hill_shape_n      = 4;
  _hill_shape_theta  = 0.5;
  _hill_shape        = std::pow( _hill_shape_theta, _hill_shape_n );

  _list_eval_step    = new std::set<int>();
#endif
}
  

// ===========================================================================
//                                 Destructor
// ===========================================================================

/*!
*/
void ExpSetup::write_setup_file(gzFile exp_setup_file) const {
  // --------------------------------------------------------------- Transfer
  gzwrite(exp_setup_file,
          static_cast<int8_t>(with_HT_),
          static_cast<int8_t>(repl_HT_with_close_points_));
  if (with_HT_)
    gzwrite(exp_setup_file,
            HT_ins_rate_,
            HT_repl_rate_);
  if(repl_HT_with_close_points_)
    gzwrite(exp_setup_file,
            repl_HT_detach_rate_);

  // --------------------------------------------------------------- Plasmids
  int8_t tmp_with_plasmids = with_plasmids();
  gzwrite(exp_setup_file, tmp_with_plasmids);
  if (tmp_with_plasmids)
    gzwrite(exp_setup_file,
            prob_plasmid_HT_,
            tune_donor_ability_,
            tune_recipient_ability_,
            donor_cost_,
            recipient_cost_,
            static_cast<int8_t>(swap_GUs_));

  // -------------------------------------------------------------- Secretion
  gzwrite(exp_setup_file,
          static_cast<int8_t>(with_secretion_),
          secretion_contrib_to_fitness_,
          secretion_cost_);

  sel()->write_setup_file(exp_setup_file);

#ifdef __REGUL
  gzwrite( exp_setup_file, &_hill_shape,  sizeof(_hill_shape) );
  gzwrite( exp_setup_file, &_hill_shape_n,  sizeof(_hill_shape_n) );
  gzwrite( exp_setup_file, &_hill_shape_theta,  sizeof(_hill_shape_theta) );

  gzwrite( exp_setup_file, &_degradation_rate,  sizeof(_degradation_rate) );
  gzwrite( exp_setup_file, &_degradation_step,  sizeof(_degradation_step) );

  gzwrite( exp_setup_file, &_nb_indiv_age,  sizeof(_nb_indiv_age) );

  gzwrite( exp_setup_file, &_with_heredity,  sizeof(_with_heredity) );
  gzwrite( exp_setup_file, &_protein_presence_limit,  sizeof(_protein_presence_limit) );

  gzwrite( exp_setup_file, &_eval_step,  sizeof(_eval_step) );

  char* binding_matrix_file_name = new char[100];

  sprintf( binding_matrix_file_name, "binding_matrix.rae" );

  gzFile binding_matrix_file = gzopen( binding_matrix_file_name, "w" );

  if ( binding_matrix_file == Z_NULL )
  {
    printf( "ERROR : Could not write binding matrix file %s\n", binding_matrix_file_name );
    exit( EXIT_FAILURE );
  }

  write_binding_matrix_to_backup( binding_matrix_file );
  gzclose( binding_matrix_file );

  delete[] binding_matrix_file_name;

  unsigned int eval_step_size = _list_eval_step->size();
  gzwrite(exp_setup_file, &eval_step_size,  sizeof(eval_step_size));

  for(auto eval_step : *_list_eval_step) {
    gzwrite(exp_setup_file, &eval_step,  sizeof(eval_step));
  }
#endif

}

void ExpSetup::load(gzFile setup_file, gzFile backup_file, bool verbose) {
  // -------------------------------------------- Retrieve transfer parameters
  int8_t tmp_with_HT;
  int8_t tmp_repl_HT_with_close_points;
  gzread(setup_file,
         tmp_with_HT,
         tmp_repl_HT_with_close_points);
  with_HT_ = static_cast<bool>(tmp_with_HT);
  repl_HT_with_close_points_ = static_cast<bool>(tmp_repl_HT_with_close_points);
  if (with_HT_)
  {
    gzread(setup_file,
           HT_ins_rate_,
           HT_repl_rate_);
  }
   if(repl_HT_with_close_points_)
    gzread(setup_file, repl_HT_detach_rate_);


  // -------------------------------------------- Retrieve plasmid parameters
  int8_t tmp_with_plasmids;
  gzread(setup_file, tmp_with_plasmids);
  with_plasmids_ = static_cast<bool>(tmp_with_plasmids);
  if (with_plasmids_)
  {
    int8_t tmp_swap_GUs;
    gzread(setup_file,
           prob_plasmid_HT_,
           tune_donor_ability_,
           tune_recipient_ability_,
           donor_cost_,
           recipient_cost_,
           tmp_swap_GUs);
    swap_GUs_ = static_cast<bool>(tmp_swap_GUs);
  }

  // ------------------------------------------ Retrieve secretion parameters
  int8_t tmp_with_secretion;
  gzread(setup_file, tmp_with_secretion,
         secretion_contrib_to_fitness_,
         secretion_cost_);
  with_secretion_ = static_cast<bool>(tmp_with_secretion);

  // ---------------------------------------------- Retrieve selection context
  sel()->load(setup_file, backup_file, verbose);

#ifdef __REGUL
  gzread( setup_file, &_hill_shape,  sizeof(_hill_shape) );
  gzread( setup_file, &_hill_shape_n,  sizeof(_hill_shape_n) );
  gzread( setup_file, &_hill_shape_theta,  sizeof(_hill_shape_theta) );

  gzread( setup_file, &_degradation_rate,  sizeof(_degradation_rate) );
  gzread( setup_file, &_degradation_step,  sizeof(_degradation_step) );

  gzread( setup_file, &_nb_indiv_age,  sizeof(_nb_indiv_age) );

  gzread( setup_file, &_with_heredity,  sizeof(_with_heredity) );
  gzread( setup_file, &_protein_presence_limit,  sizeof(_protein_presence_limit) );

  gzread( setup_file, &_eval_step,  sizeof(_eval_step) );

  char* binding_matrix_file_name = new char[100];
//    _binding_matrix = new double[MAX_QUADON][MAX_CODON];

  sprintf( binding_matrix_file_name, "binding_matrix.rae" );

  gzFile binding_matrix_file = gzopen( binding_matrix_file_name, "r" );

  if ( binding_matrix_file == Z_NULL )
  {
    printf( "ERROR : Could not read binding matrix file %s\n", binding_matrix_file_name );
    exit( EXIT_FAILURE );
  }

  read_binding_matrix_from_backup( binding_matrix_file );
  gzclose( binding_matrix_file );

  delete[] binding_matrix_file_name;

  unsigned int eval_step_size;
  gzread(setup_file, &eval_step_size,  sizeof(eval_step_size));

  int eval_val;
  for(unsigned int i = 0; i < eval_step_size; i++) {
    gzread(setup_file, &eval_val,  sizeof(eval_val));
    _list_eval_step->insert(eval_val);
  }

#endif

}

#ifdef __REGUL
void ExpSetup::init_binding_matrix( bool random_binding_matrix, double binding_zeros_percentage,
		std::shared_ptr<JumpingMT> prng)
{
  if(random_binding_matrix==1)
  {
    for( int8_t i = 0; i < MAX_QUADON; i++ )  // i for the quadons
    {
      for( int8_t j = 0; j < MAX_CODON; j++ )  // j for the codons
      {
        if( prng->random() > binding_zeros_percentage)
        {
        	_binding_matrix[i][j] = prng->random();
        }
        else
        {
        	_binding_matrix[i][j] = 0;
        }
//  	  printf("m[%d][%d] = %f\n",MAX_QUADON,MAX_CODON, _binding_matrix[i][j]);
      }
    }
  }
  else // random_binding_matrix == 0
  {
    char* binding_matrix_file_name = new char[100];
//    _binding_matrix = new double[MAX_QUADON][MAX_CODON];

    sprintf( binding_matrix_file_name, "binding_matrix.rae" );

    gzFile binding_matrix_file = gzopen( binding_matrix_file_name, "r" );

    if ( binding_matrix_file == Z_NULL )
    {
      printf( "ERROR : Could not read binding matrix file %s\n", binding_matrix_file_name );
      exit( EXIT_FAILURE );
    }

    read_binding_matrix_from_backup( binding_matrix_file );
    gzclose( binding_matrix_file );

    delete[] binding_matrix_file_name;
  }




}

void ExpSetup::read_binding_matrix_from_backup(gzFile binding_matrix_file) {
	for (int i=0; i < MAX_QUADON; i++)
		for (int j=0; j < MAX_CODON; j++) {
			gzread( binding_matrix_file, &(_binding_matrix[i][j]), sizeof(double) );
		}
}

void ExpSetup::write_binding_matrix_to_backup(gzFile binding_matrix_file) const {
	double value;
	for (int i=0; i < MAX_QUADON; i++)
		for (int j=0; j < MAX_CODON; j++) {
			value = _binding_matrix[i][j];
			gzwrite( binding_matrix_file, &value, sizeof(double) );
		}
}

void ExpSetup::write_binding_matrix_to_file( FILE* file ) const
{
  for( int16_t row = 0 ; row < MAX_QUADON ; row++ )
  {
    for( int16_t column = 0 ; column < MAX_CODON ; column++ )
    {
      fprintf( file, "\t%e", _binding_matrix[row][column] );
    }
    fprintf( file, "\n");
  }
}

void ExpSetup::print_binding_matrix()
{
  for( int16_t row = 0 ; row < MAX_QUADON ; row++ )
  {
    for( int16_t column = 0 ; column < MAX_CODON ; column++ )
    {
      printf("M[%d][%d] = %e\n", row, column, _binding_matrix[row][column] );
    }
  }
}
#endif

// ===========================================================================
//                                Protected Methods
// ===========================================================================
} // namespace aevol
