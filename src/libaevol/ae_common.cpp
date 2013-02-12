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
#include <stdlib.h>
#include <math.h>

// =================================================================
//                            Project Files
// =================================================================
#include <ae_common.h>
#include <ae_exp_setup.h>
#include <ae_population.h>
#include <ae_list.h>
#include <ae_gaussian.h>
#include <ae_point_2d.h>

//##############################################################################
//                                                                             #
//                               Class ae_common                               #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

ae_exp_setup*     ae_common::sim          = NULL;
ae_population*    ae_common::pop          = NULL;
ae_params_init*   ae_common::init_params  = NULL;
ae_params*        ae_common::params       = NULL;
ae_params_record* ae_common::rec_params   = NULL;


#ifdef __X11
  ae_exp_setup_X11*  ae_common::sim_display = NULL;
  ae_population_X11* ae_common::pop_display = NULL;
#endif


// Number of generations to be simulated
int32_t ae_common::nb_generations         = 1000000;

// Population properties
bool    ae_common::pop_structure     = false;
int16_t ae_common::grid_x            = 25; 
int16_t ae_common::grid_y            = 40; 
int32_t ae_common::migration_number  = 0;



#ifdef __REGUL
  // Regulation factors
  double ae_common::hill_shape_n      = 4;
  double ae_common::hill_shape_theta  = 0.5;
  double ae_common::hill_shape        = pow( hill_shape_theta, hill_shape_n );

  // Degradation equation
  double ae_common::degradation_rate  = 1;
  double ae_common::degradation_step  = 0.1;

  // Individual life
  int16_t         ae_common::individual_evaluation_nbr   = 1;
  ae_array_short* ae_common::individual_evaluation_dates = NULL;

  // Binding matrix
  ae_matrix_double* ae_common::_binding_matrix = NULL;

  // Heredity
  bool    ae_common::with_heredity          = false;
  double  ae_common::protein_presence_limit = 1e-5;
#endif

// =================================================================
//                             Constructors
// =================================================================

// =================================================================
//                             Destructors
// =================================================================

// =================================================================
//                            Public Methods
// =================================================================
void ae_common::write_to_backup( gzFile backup_file )
{
  // Version info
  int8_t aevol_version_len = strlen(PACKAGE_VERSION);
  gzwrite( backup_file, &aevol_version_len, sizeof(aevol_version_len) );
  char* aevol_version = new char[aevol_version_len + 1];
  strcpy( aevol_version, PACKAGE_VERSION );
  gzwrite( backup_file, aevol_version, (aevol_version_len+1) * sizeof(aevol_version[0])  );
  delete aevol_version;
  
  // Initialization parameters
  init_params->write_to_backup( backup_file ); // TODO : get rid of !
  
  // Permanent parameters
  //~ params->write_to_backup( backup_file );
  
  // Recording parameters
  rec_params->write_to_backup( backup_file );
  
  // Population properties 
  int8_t tmp_pop_structure = pop_structure? 1 : 0;
  gzwrite( backup_file, &tmp_pop_structure,           sizeof(tmp_pop_structure)         );
  gzwrite( backup_file, &grid_x,                      sizeof(grid_x)                    );
  gzwrite( backup_file, &grid_y,                      sizeof(grid_y)                    );
  gzwrite( backup_file, &migration_number,            sizeof(migration_number)          );
  

  // R-AEVOL specific
  #ifdef __REGUL
    // Regulation factors
    gzwrite( backup_file, &hill_shape_n,                sizeof(hill_shape_n)      );
    gzwrite( backup_file, &hill_shape_theta,            sizeof(hill_shape_theta)  );
    gzwrite( backup_file, &hill_shape,                  sizeof(hill_shape)        );

    // Degradation rate
    gzwrite( backup_file, &degradation_rate,            sizeof(degradation_rate)  );
    gzwrite( backup_file, &degradation_step,            sizeof(degradation_step)  );

    // Individual life
    gzwrite( backup_file, &individual_evaluation_nbr,   sizeof(individual_evaluation_nbr)  );
    individual_evaluation_dates->write_to_backup( backup_file );

    // Binding matrix
    _binding_matrix->write_to_backup( backup_file );
    
    // Heredity
    int8_t tmp_with_heredity = with_heredity? 1 : 0;
    gzwrite( backup_file, &tmp_with_heredity,           sizeof(tmp_with_heredity)           );
    gzwrite( backup_file, &protein_presence_limit,      sizeof(protein_presence_limit)      );
  #endif
}

void ae_common::read_from_backup( gzFile backup_file, bool verbose )
{
  // Version info
  int8_t aevol_version_len;
  gzread( backup_file, &aevol_version_len, sizeof(aevol_version_len) );
  char* aevol_version = new char[aevol_version_len + 1];
  gzread( backup_file, aevol_version, (aevol_version_len+1) * sizeof(aevol_version[0]) );
  int major, minor, cur_major, cur_minor;
  sscanf( aevol_version, "%d.%d.", &major, &minor );
  sscanf( PACKAGE_VERSION, "%d.%d.", &cur_major, &cur_minor );
  if ( cur_major != major || cur_minor != minor )
  {
    printf( "ERROR in %s:%d :\n  This backup was generated with version %d.%d, current version is %d.%d\n", __FILE__, __LINE__, major, minor, cur_major, cur_minor );
    exit( EXIT_FAILURE );
  }
  delete aevol_version;
  
  // Initialization parameters
  init_params = new ae_params_init(); // TODO : get rid of !
  init_params->read_from_backup( backup_file );
  
  // Permanent parameters
  params = new ae_params();
  //~ params->read_from_backup( backup_file );
  
  // Recording parameters
  rec_params = new ae_params_record();
  rec_params->read_from_backup( backup_file );
  
  // Population properties 
  int8_t tmp_pop_structure;
  gzread( backup_file, &tmp_pop_structure,            sizeof(tmp_pop_structure)         );
  pop_structure = (tmp_pop_structure!=0);
  gzread( backup_file, &grid_x,                       sizeof(grid_x)                    );
  gzread( backup_file, &grid_y,                       sizeof(grid_y)                    );
  gzread( backup_file, &migration_number,             sizeof(migration_number)          );
  

#ifdef __REGUL
  // Regulation factors
  gzread( backup_file, &hill_shape_n,                 sizeof(hill_shape_n)      );
  gzread( backup_file, &hill_shape_theta,             sizeof(hill_shape_theta)  );
  gzread( backup_file, &hill_shape,                   sizeof(hill_shape)        );
  
  // Degradation rate
  gzread( backup_file, &degradation_rate,             sizeof(degradation_rate)  );
  gzread( backup_file, &degradation_step,             sizeof(degradation_step)  );
  
  // Individual life
  gzread( backup_file, &individual_evaluation_nbr,    sizeof(individual_evaluation_nbr)    );
  individual_evaluation_dates  = new ae_array_short( individual_evaluation_nbr );
  individual_evaluation_dates->read_from_backup( backup_file );


  // Binding matrix
  _binding_matrix = new ae_matrix_double( MAX_QUADON, MAX_CODON );
  _binding_matrix->read_from_backup( backup_file );

      
  // Heredity
  int8_t tmp_with_heredity;
  gzread( backup_file, &tmp_with_heredity,           sizeof(tmp_with_heredity)           );
  with_heredity = (tmp_with_heredity!=0);
  gzread( backup_file, &protein_presence_limit,      sizeof(protein_presence_limit)      );
#endif
}

void ae_common::print_to_file( void )
{
  FILE* param_out = fopen( "param.out", "w" );
  assert( param_out );
  
  // Initialization parameters
  fprintf( param_out, "\n========================== INITIALIZATION PARAMETERS ==========================\n" );
  init_params->print_to_file( param_out );
  
  fprintf( param_out, "\n============================ SIMULATION PARAMETERS =============================\n" );
  params->print_to_file( param_out );
  
  fprintf( param_out, "\n============================= RECORDING PARAMETERS =============================\n" );
  rec_params->print_to_file( param_out );
  
  fprintf( param_out, "\n=============================== OTHER PARAMETERS ===============================\n" );
  
  // Number of generations to compute
  fprintf( param_out, "\nNumber of generations to compute ------------------------\n" );
  fprintf( param_out, "nb_generations               %"PRId32"\n", nb_generations          );

  // Population properties
  fprintf( param_out, "\nPopulation properties -----------------------------------\n" );
  fprintf( param_out, "pop_structure :              %s\n",  pop_structure? "true" : "false" );
  fprintf( param_out, "grid_x :                     %"PRId16"\n", grid_x                    ); 
  fprintf( param_out, "grid_y :                     %"PRId16"\n", grid_y                    ); 
  fprintf( param_out, "migration_number :           %"PRId32"\n", migration_number          );

  #ifdef __REGUL
    fprintf( param_out, "\nRAEVOL SPECIFIC PARAMETERS ==============================\n" );
  
    // Regulation factors
    fprintf( param_out, "\nRegulation factors --------------------------------------\n" );
    fprintf( param_out, "hill_shape_n :                %e\n",                 hill_shape_n  );
    fprintf( param_out, "hill_shape_theta :            %e\n",             hill_shape_theta  );
    fprintf( param_out, "hill_shape :                  %e\n",                   hill_shape  );

    // Degradation rate
    fprintf( param_out, "\nDegradation rate ----------------------------------------\n" );
    fprintf( param_out, "degradation_rate :            %e\n",             degradation_rate  );
    fprintf( param_out, "degradation_step :            %e\n",             degradation_step  );

    // Individual life
    fprintf( param_out, "\nIndividual life -----------------------------------------\n" );
    fprintf( param_out, "individual_evaluation_nbr :   %"PRId16"\n",    individual_evaluation_nbr  );
    fprintf( param_out, "individual_evaluation_dates : ");
    individual_evaluation_dates->print_to_file( param_out );

    // Binding matrix
    fprintf( param_out, "\nBinding matrix ------------------------------------------\n" );
    fprintf( param_out, "binding_matrix :\n");
    _binding_matrix->print_to_file( param_out );

    // Heredity
    fprintf( param_out, "\nHeredity ------------------------------------------------\n" );
    fprintf( param_out, "with_heredity :               %s\n", with_heredity? "true" : "false"  );
    fprintf( param_out, "protein_presence_limit :      %e\n", protein_presence_limit           );
  #endif
  
  fclose( param_out );
}

void ae_common::clean( void )
{
  if (sim != NULL)
  {
    delete sim;
    sim = NULL;
  }
  
  #ifdef __REGUL
    if ( individual_evaluation_dates != NULL )
    {
      delete individual_evaluation_dates;
    }
    if ( _binding_matrix != NULL )
    {
      delete _binding_matrix;
      _binding_matrix = NULL;
    }
  #endif
}

#ifdef __REGUL
void ae_common::init_binding_matrix( void )
{
  _binding_matrix = new ae_matrix_double( MAX_QUADON, MAX_CODON );

  for( int8_t i = 0; i < MAX_QUADON; i++ )  // i for the quadons
  {
    for( int8_t j = 0; j < MAX_CODON; j++ )  // j for the codons
    {
      if( sim->alea->random() > binding_zeros_percentage)
      {
        _binding_matrix->set_value( i, j, sim->alea->random() );
      }
      else
      {
        _binding_matrix->set_value( i, j, 0 );
      }
    }
  }
}
#endif

// =================================================================
//                           Protected Methods
// =================================================================
