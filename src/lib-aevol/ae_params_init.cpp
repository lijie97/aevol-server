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


/*! \class ae_params_init
    \brief Contains all the parameters needed only at startup e.g. init_genome_size
*/
 
 
// =================================================================
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================
#include <ae_params_init.h>
#include <ae_gaussian.h>
#include <ae_point_2d.h>




//##############################################################################
//                                                                             #
//                             Class ae_params_init                            #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_params_init::ae_params_init( void )
{
  // ==================== Initialize with default values ====================
  // PseudoRandom Number Generator
  _seed     = 0;
  _env_seed = 0;
  
  // Initial conditions
  _initial_genome_length  = 5000;
  _init_method            = ONE_GOOD_GENE | CLONE;
  _init_pop_size          = 1000;
  
  // Whether to delete the existing statistics file (otherwise kept with the suffix ".old")
  _delete_old_stats       = false;
  
  // Environment
  _env_gaussians      = new ae_list();
  _env_custom_points  = new ae_list();
  _env_sampling           = 300;
  
  // Environment x-axis segmentation
  _env_axis_is_segmented        = false;
  _env_axis_nb_segments         = 1;
  _env_axis_segment_boundaries  = NULL;
  _env_axis_features            = NULL;
  _env_axis_separate_segments   = false;

  // Environment variation
  _env_var_method = NONE;
  _env_var_sigma  = 0.01;
  _env_var_tau    = 1000;
  
  // Environment noise TODO
  
  // Secretion
  // starting configuration of secretion grid; 0, all are 0; 1, point source of secreted compund
  _secretion_init         = 0;
  
  // Plasmids
  _plasmid_initial_length = 1000;
  _plasmid_initial_gene   = 0;
  
  #ifdef __REGUL
    // Binding matrix
    _binding_zeros_percentage = 0.75;
  #endif
}

// =================================================================
//                             Destructors
// =================================================================
ae_params_init::~ae_params_init( void )
{
  _env_gaussians->erase( DELETE_OBJ );
  delete _env_gaussians;
  
  _env_custom_points->erase( DELETE_OBJ );
  delete _env_custom_points;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_params_init::write_to_backup( gzFile* backup_file ) // Usefull?
{
  // PseudoRandom Number Generator
  gzwrite( backup_file, &_seed,                       sizeof(_seed)                     );
  gzwrite( backup_file, &_env_seed,                   sizeof(_env_seed)                 );

  // Initial conditions
  gzwrite( backup_file, &_initial_genome_length,      sizeof(_initial_genome_length)    );
  gzwrite( backup_file, &_init_method,                sizeof(_init_method)              );
  gzwrite( backup_file, &_init_pop_size,              sizeof(_init_pop_size)            );

  // Whether to delete the existing statistics file (otherwise kept with the suffix ".old")
  int8_t tmp_delete_old_stats = _delete_old_stats? 1 : 0;
  gzwrite( backup_file, &tmp_delete_old_stats,        sizeof(tmp_delete_old_stats)      );
  
  // Environment gaussians
  int16_t nb_gaussians = _env_gaussians->get_nb_elts();
  gzwrite( backup_file, &nb_gaussians, sizeof(nb_gaussians) );
  
  ae_list_node* gaussian_node = _env_gaussians->get_first();
  ae_gaussian*  gaussian;
  for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
  {
    gaussian = ( ae_gaussian* ) gaussian_node->get_obj();

    gaussian->write_to_backup( backup_file );

    gaussian_node = gaussian_node->get_next();
  }
  
  // Environment custom points
  int16_t nb_points = _env_custom_points->get_nb_elts();
  gzwrite( backup_file, &nb_points, sizeof(nb_points) );
  ae_list_node* point_node = _env_custom_points->get_first();
  ae_point_2d*  point;
  for ( int16_t i = 0 ; i < nb_points ; i++ )
  {
    point = ( ae_point_2d* ) point_node->get_obj();

    point->write_to_backup( backup_file );

    point_node = point_node->get_next();
  }
  
  // Environment sampling
  gzwrite( backup_file, &_env_sampling, sizeof(_env_sampling) );
  
  
  // Environment x-axis segmentation
  int8_t tmp_env_axis_is_segmented = _env_axis_is_segmented? 1 : 0;
  gzwrite( backup_file, &tmp_env_axis_is_segmented, sizeof(tmp_env_axis_is_segmented) );
  
  if ( _env_axis_is_segmented )
  {
    gzwrite( backup_file, &_env_axis_nb_segments, sizeof(_env_axis_nb_segments) );
    
    // Write segment boundaries
    for ( int16_t i = 0 ; i < _env_axis_nb_segments + 1 ; i++ )
    {
      gzwrite( backup_file, &_env_axis_segment_boundaries[i], sizeof(_env_axis_segment_boundaries[i]) );
    }
    
    // Write segment features
    for ( int16_t i = 0 ; i < _env_axis_nb_segments ; i++ )
    {
      int8_t tmp_env_axis_features = _env_axis_features[i];
      gzwrite( backup_file, &tmp_env_axis_features, sizeof(tmp_env_axis_features) );
    }
    
    int8_t tmp_env_axis_separate_segments = _env_axis_separate_segments? 1 : 0;
    gzwrite( backup_file, &tmp_env_axis_separate_segments, sizeof(tmp_env_axis_separate_segments) );
  }
  
  
  // Environment variation
  int8_t tmp_env_var_method = _env_var_method;
  gzwrite( backup_file, &tmp_env_var_method,  sizeof(tmp_env_var_method)  );
  gzwrite( backup_file, &_env_var_sigma,      sizeof(_env_var_sigma)      );
  gzwrite( backup_file, &_env_var_tau,        sizeof(_env_var_tau)        );
  
  // Environment noise TODO
  
  
  // Secretion
  gzwrite( backup_file, &_secretion_init, sizeof(_secretion_init) );

  // Plasmids
  gzwrite( backup_file, &_plasmid_initial_length, sizeof(_plasmid_initial_length) );
  gzwrite( backup_file, &_plasmid_initial_gene,   sizeof(_plasmid_initial_gene)   );  

  // R-AEVOL specific
  #ifdef __REGUL
    // Binding matrix
    gzwrite( backup_file, &_binding_zeros_percentage, sizeof(_binding_zeros_percentage) );
  #endif
}

void ae_params_init::read_from_backup( gzFile* backup_file, bool verbose ) // Usefull?
{
  // PseudoRandom Number Generator
  gzread( backup_file, &_seed,                      sizeof(_seed)                     );
  gzread( backup_file, &_env_seed,                  sizeof(_env_seed)                 );

  // Initial conditions
  gzread( backup_file, &_initial_genome_length,     sizeof(_initial_genome_length)    );
  gzread( backup_file, &_init_method,               sizeof(_init_method)              );
  gzread( backup_file, &_init_pop_size,             sizeof(_init_pop_size)            );

  // Statistics collection
  int8_t tmp_delete_old_stats;
  gzread( backup_file, &tmp_delete_old_stats,       sizeof(tmp_delete_old_stats)      );
  _delete_old_stats = (tmp_delete_old_stats != 0);
  
  // Environment gaussians
  if ( verbose )
  {
    printf( "    Loading initial environment\n" );
  }
  int16_t nb_gaussians;
  gzread( backup_file, &nb_gaussians, sizeof(nb_gaussians) );
  
  if ( _env_gaussians->is_empty() == false )
  {
    _env_gaussians->erase( DELETE_OBJ );
  }
  
  for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
  {
    _env_gaussians->add( new ae_gaussian( backup_file ) );
  }
  
  // Environment custom points
  int16_t nb_points;
  gzread( backup_file, &nb_points, sizeof(nb_points) );

  if ( _env_custom_points->is_empty() == false )
  {
    _env_custom_points->erase( DELETE_OBJ );
  }

  for ( int16_t i = 0 ; i < nb_points ; i++ )
  {
    _env_custom_points->add( new ae_point_2d( backup_file ) );
  }
  
  // Environment sampling
  gzread( backup_file, &_env_sampling,    sizeof(_env_sampling) );
  
  
  // Environment x-axis segmentation
  if ( verbose )
  {
    printf( "    Loading x-axis segmentation parameters\n" );
  }
  int8_t tmp_env_axis_is_segmented;
  gzread( backup_file, &tmp_env_axis_is_segmented, sizeof(tmp_env_axis_is_segmented) );
  _env_axis_is_segmented = (tmp_env_axis_is_segmented != 0);
  
  if ( _env_axis_is_segmented )
  {
    gzread( backup_file, &_env_axis_nb_segments, sizeof(_env_axis_nb_segments) );
    _env_axis_segment_boundaries = new double [_env_axis_nb_segments + 1];
    for ( int16_t i = 0 ; i < _env_axis_nb_segments + 1 ; i++ )
    {
      gzread( backup_file, &_env_axis_segment_boundaries[i], sizeof(_env_axis_segment_boundaries[i]) );
    }
    _env_axis_features = new ae_env_axis_feature [_env_axis_nb_segments];
    for ( int16_t i = 0 ; i < _env_axis_nb_segments ; i++ )
    {
      int8_t tmp_env_axis_features;
      gzread( backup_file, &tmp_env_axis_features, sizeof(tmp_env_axis_features) );
      _env_axis_features[i] = (ae_env_axis_feature) tmp_env_axis_features;
    }
    int8_t tmp_env_axis_separate_segments;
    gzread( backup_file, &tmp_env_axis_separate_segments,   sizeof(tmp_env_axis_separate_segments) );
    _env_axis_separate_segments = (tmp_env_axis_separate_segments!=0);
  }
  
  
  // Environment variation
  int8_t tmp_env_var_method;
  gzread( backup_file, &tmp_env_var_method,          sizeof(tmp_env_var_method) );
  _env_var_method = (ae_env_var) tmp_env_var_method;
  gzread( backup_file, &_env_var_sigma,              sizeof(_env_var_sigma)     );
  gzread( backup_file, &_env_var_tau,                sizeof(_env_var_tau)       );
  
  // Environment noise TODO
  
  
  // Secretion
  gzread( backup_file, &_secretion_init,  sizeof(_secretion_init) );

  // Plasmids
  gzread( backup_file, &_plasmid_initial_length, sizeof(_plasmid_initial_length) );
  gzread( backup_file, &_plasmid_initial_gene,   sizeof(_plasmid_initial_gene)   );

  
  // R-AEVOL specific
  #ifdef __REGUL
    // Binding matrix
    gzread( backup_file, &_binding_zeros_percentage, sizeof(_binding_zeros_percentage) );
  #endif
}

void ae_params_init::print_to_file( FILE* file )
{
  // PseudoRandom Number Generator
  fprintf( file, "\nPseudoRandom Number Generator ---------------------------\n" );
  fprintf( file, "seed :                       %"PRId32"\n", _seed     );
  fprintf( file, "env_seed :                   %"PRId32"\n", _env_seed );
  
  // Initial conditions
  fprintf( file, "\nInitial conditions --------------------------------------\n" );
  fprintf( file, "initial_genome_length :      %"PRId32"\n", _initial_genome_length );
  fprintf( file, "init_method :               " );
  if ( _init_method & ONE_GOOD_GENE )
  {
    fprintf( file, " ONE_GOOD_GENE" );
  }
  if ( _init_method & CLONE )
  {
    fprintf( file, " CLONE" );
  }
  if ( _init_method & WITH_INS_SEQ )
  {
    fprintf( file, " WITH_INS_SEQ" );
  }
  fprintf( file, "\n" );
  fprintf( file, "init_pop_size :              %"PRId32"\n", _init_pop_size );
  
  // Environment gaussians
  fprintf( file, "\nEnvironment gaussians -----------------------------------\n" );
  ae_list_node * node = _env_gaussians->get_first();
  while ( node != NULL )
  {
    ae_gaussian * gauss = (ae_gaussian *) node->get_obj();
    fprintf( file, "env_add_gaussian :           %f %f %f \n",gauss->get_height(),gauss->get_mean(),gauss->get_width());
    node = node->get_next();
  }
  
  // Environment custom_points
  fprintf( file, "\nEnvironment custom_points -------------------------------\n" );
  node = _env_custom_points->get_first();
  while ( node != NULL )
  {
    ae_point_2d * point = (ae_point_2d *)node->get_obj();
    fprintf( file, "env_add_point:    %f %f \n",point->x,point->y);
    node = node->get_next();	
  }
  
  // Environment sampling
  fprintf( file, "\nEnvironment sampling ------------------------------------\n" );
  fprintf( file, "env_sampling :               %"PRId16"\n", _env_sampling );
  
  // Environment segmentation
  fprintf( file, "\nEnvironment segmentation --------------------------------\n" );
  fprintf( file, "env_axis_is_segmented  :     %s\n", _env_axis_is_segmented? "true" : "false");
  fprintf( file, "env_axis_nb_segments   :     %"PRId16"\n", _env_axis_nb_segments);
  
  // Environment axis segment boundaries and features
  fprintf( file, "\nEnvironment axis segment boundaries and features --------\n" );
  fprintf( file, "env_axis_segment_boundaries :");
  if ( _env_axis_nb_segments > 1 )
  {
    for ( int k = 0 ; k < _env_axis_nb_segments + 1 ; k++ )
    {
      fprintf( file, " %f ", _env_axis_segment_boundaries[k] );
    }

    fprintf( file, "\n");
    fprintf( file, "env_axis_features :");
    for ( int k = 0 ; k < _env_axis_nb_segments ; k++ )
    {
      switch ( _env_axis_features[k] )
	    {
        case NEUTRAL :
	      {
          fprintf( file, " NEUTRAL " );
          break;
	      }
        case METABOLISM :
	      {
          fprintf( file, " METABOLISM ");
          break;
	      }
        case SECRETION :
	      {
          fprintf( file, " SECRETION ");
          break;
	      }
        default :
	      {
          fprintf( file, " UNKNOWN " );
          break;
	      }
	    }
    }
  }
  fprintf( file, "\n");
  fprintf( file, "env_axis_separate_segments :      %s\n", _env_axis_separate_segments? "true" : "false");
  
  
  // Environment variation
  fprintf( file, "\nEnvironment variation -----------------------------------\n" );
  switch ( _env_var_method )
  {
    case NONE :
    {
      fprintf( file, "env_var_method :             NONE\n" );
      break;
    }
    case AUTOREGRESSIVE_MEAN_VAR :
    {
      fprintf( file, "env_var_method :             AUTOREGRESSIVE_MEAN_VAR\n" );
      break;
    }
    case LOCAL_GAUSSIANS_VAR :
    {
      fprintf( file, "env_var_method :             LOCAL_GAUSSIANS_VAR\n" );
      break;
    }

    default :
    {
      fprintf( file, "env_var_method :             UNKNOWN\n" );
      break;
    }
  }
  fprintf ( file, "env_var_sigma :                  %f \n",       _env_var_sigma );
  fprintf ( file, "env_var_tau   :                  %"PRId32"\n", _env_var_tau   );
  
  
  // Environment noise
  fprintf( file, "\nEnvironment noise ---------------------------------------\n" );
  
  
  // Secretion
  fprintf( file, "\nSecretion -----------------------------------------------\n" );
  fprintf( file, "secretion_init :             %e\n", _secretion_init );

  
  // Plasmids
  fprintf( file, "\nPlasmids ------------------------------------------------\n" );
  fprintf( file, "plasmid_initial_length :     %"PRId32"\n", _plasmid_initial_length                 );
  fprintf( file, "plasmid_initial_gene :       %"PRId32"\n", _plasmid_initial_gene                   );
}

void ae_params_init::clean( void )
{
  _env_gaussians->erase( DELETE_OBJ );
  delete _env_gaussians;
  _env_gaussians = NULL;
  _env_custom_points->erase( DELETE_OBJ );
  delete _env_custom_points;
  _env_custom_points = NULL;
  
  if ( _env_axis_features != NULL )
  {
    delete [] _env_axis_features;
    _env_axis_features = NULL;
  }

  if ( _env_axis_segment_boundaries != NULL )
  {
    delete [] _env_axis_segment_boundaries;
    _env_axis_segment_boundaries = NULL;
  }
}

// =================================================================
//                           Protected Methods
// =================================================================

// =================================================================
//                          Non inline accessors
// =================================================================
