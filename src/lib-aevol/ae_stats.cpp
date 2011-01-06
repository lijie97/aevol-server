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
#include <stdio.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_stats.h>
#include <ae_stat_record.h>
#include <ae_simulation.h>
#include <ae_population.h>
#include <ae_individual.h>
#include <ae_genetic_unit.h>
#ifdef __REGUL
  #include <ae_influence_R.h>
  #include <ae_protein_R.h>
#endif




//##############################################################################
//                                                                             #
//                                Class ae_stats                               #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_stats::ae_stats( void )
{
  _stat_files_best_names = new char*[12];
  _stat_files_glob_names = new char*[12];

  _stat_files_best = new FILE*[12];
  _stat_files_glob = new FILE*[12];

  for ( int8_t i = 0 ; i < 12 ; i++ )
  {
    _stat_files_glob[i] = NULL;
    _stat_files_best[i] = NULL;
    _stat_files_best_names[i] = new char[80];
    _stat_files_glob_names[i] = new char[80];
    _stat_files_best_names[i][0] = '\0';
    _stat_files_glob_names[i][0] = '\0';
  }
    
  strcpy( _stat_files_best_names[0], "stat_fitness_best.out" );
  strcpy( _stat_files_best_names[1], "stat_mutation_best.out" );
  strcpy( _stat_files_best_names[2], "stat_genes_best.out" );
  strcpy( _stat_files_best_names[3], "stat_bp_best.out" );

  strcpy( _stat_files_glob_names[0], "stat_fitness_glob.out" );
  strcpy( _stat_files_glob_names[1], "stat_mutation_glob.out" );
  strcpy( _stat_files_glob_names[2], "stat_genes_glob.out" );
  strcpy( _stat_files_glob_names[3], "stat_bp_glob.out" );
  
  for (int32_t i = 0; i < 4; i++) 
  {
    _stat_files_best[i] = fopen(_stat_files_best_names[i], "w"); 
    _stat_files_glob[i] = fopen(_stat_files_glob_names[i], "w"); 
    assert (_stat_files_glob[i]); 
    assert (_stat_files_best[i]); 
  }
  
  if ( ae_common::allow_plasmids )
  {
  
    strcpy( _stat_files_best_names[4], "stat_fitness_chromosome_best.out" );
    strcpy( _stat_files_best_names[5], "stat_mutation_chromosome_best.out" );
    strcpy( _stat_files_best_names[6], "stat_genes_chromosome_best.out" );
    strcpy( _stat_files_best_names[7], "stat_bp_chromosome_best.out" );
    
    strcpy( _stat_files_best_names[8], "stat_fitness_plasmid_best.out" );
    strcpy( _stat_files_best_names[9], "stat_mutation_plasmid_best.out" );
    strcpy( _stat_files_best_names[10], "stat_genes_plasmid_best.out" );
    strcpy( _stat_files_best_names[11], "stat_bp_plasmid_best.out" );
    
    strcpy( _stat_files_glob_names[4], "stat_fitness_chromosome_glob.out" );
    strcpy( _stat_files_glob_names[5], "stat_mutation_chromosome_glob.out" );
    strcpy( _stat_files_glob_names[6], "stat_genes_chromosome_glob.out" );
    strcpy( _stat_files_glob_names[7], "stat_bp_chromosome_glob.out" );
    
    strcpy( _stat_files_glob_names[8], "stat_fitness_plasmid_glob.out" );
    strcpy( _stat_files_glob_names[9], "stat_mutation_plasmid_glob.out" );
    strcpy( _stat_files_glob_names[10], "stat_genes_plasmid_glob.out" );
    strcpy( _stat_files_glob_names[11], "stat_bp_plasmid_glob.out" );
  
    for (int32_t i = 4; i < 12; i++) 
    {
      _stat_files_best[i] = fopen(_stat_files_best_names[i], "w"); 
      _stat_files_glob[i] = fopen(_stat_files_glob_names[i], "w"); 
      assert (_stat_files_glob[i]); 
      assert (_stat_files_best[i]); 
    } 
  }
}



// Give the filename without the .out extension
// NB: This is used only in ancstats.cpp, for post-treatment analysis
ae_stats::ae_stats( const char * ancstat_file_name )
{
  _stat_files_best_names = NULL;
  _stat_files_glob_names = NULL;
  _stat_files_glob = NULL;
  _stat_files_best = new FILE*[12];
  
  char filename[100];

  char** _stat_data_type_names;
  _stat_data_type_names = new char*[4];

  for ( int8_t i = 0 ; i < 4 ; i++ )
  {
    _stat_data_type_names[i]  = new char[30];
  }
  
  strcpy( _stat_data_type_names[0], "_fitness" ); 
  strcpy( _stat_data_type_names[1], "_mutation" ); 
  strcpy( _stat_data_type_names[2], "_genes" ); 
  strcpy( _stat_data_type_names[3], "_bp" ); 
  
  for (int32_t i = 0; i < 4; i++) 
  {
    strncpy(filename, ancstat_file_name, 100);
    filename[99] = '\0';
    strcat(filename, _stat_data_type_names[i]);
    strcat(filename, ".out");
    _stat_files_best[i] = fopen(filename, "w"); 
    assert (_stat_files_best[i]);
  }
  
  delete [] _stat_data_type_names; 
  
  if ( ae_common::allow_plasmids )
  {
    char** _stat_GU_names;
    _stat_GU_names   = new char*[2];
    _stat_GU_names[0]  = new char[30];
    _stat_GU_names[1]  = new char[30];
    strcpy( _stat_GU_names[0], "_chromosome" );
    strcpy( _stat_GU_names[1], "_plasmid" );

    
    for (int32_t i = 0; i < 2; i++) 
    {
      for (int32_t j = 0; j < 4; i++) 
      {
        strncpy(filename, ancstat_file_name, 100);
        filename[99] = '\0';
        strcat(filename, _stat_data_type_names[i]);
        strcat(filename, _stat_GU_names[j]);
        strcat(filename, ".out");
        _stat_files_best[4+i*4+j] = fopen(filename, "w"); 
        assert(_stat_files_best[4+i*4+j]);
      }
    }
    
    delete [] _stat_GU_names; 
  }
  else
    {
      for (int32_t i = 4; i < 12; i++) _stat_files_best[i] = NULL;
    }
}


ae_stats::ae_stats( int32_t num_gener )
{
  _stat_files_best_names = new char*[12];
  _stat_files_glob_names = new char*[12];

  _stat_files_glob = new FILE*[12];
  _stat_files_best = new FILE*[12];

  for ( int8_t i = 0 ; i < 12 ; i++ )
  {
    _stat_files_glob[i] = NULL;
    _stat_files_best[i] = NULL;
    _stat_files_best_names[i] = new char[80];
    _stat_files_glob_names[i] = new char[80];
    _stat_files_best_names[i][0] = '\0';
    _stat_files_glob_names[i][0] = '\0';
  }
  
  strcpy( _stat_files_best_names[0], "stat_fitness_best.out" );
  strcpy( _stat_files_best_names[1],"stat_mutation_best.out" );
  strcpy( _stat_files_best_names[2], "stat_genes_best.out" );
  strcpy( _stat_files_best_names[3], "stat_bp_best.out" );

  strcpy( _stat_files_glob_names[0], "stat_fitness_glob.out" );
  strcpy( _stat_files_glob_names[1], "stat_mutation_glob.out" );
  strcpy( _stat_files_glob_names[2], "stat_genes_glob.out" );
  strcpy( _stat_files_glob_names[3], "stat_bp_glob.out" );
   
  if ( ae_common::allow_plasmids )
  {
    strcpy( _stat_files_best_names[4], "stat_fitness_chromosome_best.out" );
    strcpy( _stat_files_best_names[5], "stat_mutation_chromosome_best.out" );
    strcpy( _stat_files_best_names[6], "stat_genes_chromosome_best.out" );
    strcpy( _stat_files_best_names[7], "stat_bp_chromosome_best.out" );
    
    strcpy( _stat_files_best_names[8], "stat_fitness_plasmid_best.out" );
    strcpy( _stat_files_best_names[9], "stat_mutation_plasmid_best.out" );
    strcpy( _stat_files_best_names[10], "stat_genes_plasmid_best.out" );
    strcpy( _stat_files_best_names[11], "stat_bp_plasmid_best.out" );
    
    strcpy( _stat_files_glob_names[4], "stat_fitness_chromosome_glob.out" );
    strcpy( _stat_files_glob_names[5], "stat_mutation_chromosome_glob.out" );
    strcpy( _stat_files_glob_names[6], "stat_genes_chromosome_glob.out" );
    strcpy( _stat_files_glob_names[7], "stat_bp_chromosome_glob.out" );
    
    strcpy( _stat_files_glob_names[8], "stat_fitness_plasmid_glob.out" );
    strcpy( _stat_files_glob_names[9], "stat_mutation_plasmid_glob.out" );
    strcpy( _stat_files_glob_names[10], "stat_genes_plasmid_glob.out" );
    strcpy( _stat_files_glob_names[11], "stat_bp_plasmid_glob.out" );
   }
  
  char old_stat_file_glob_name[60];
  char old_stat_file_best_name[60];

  FILE* old_stat_file_best = NULL; 
  FILE* old_stat_file_glob = NULL; 
  
  char  line[500];
  char* ret;
  
  int16_t num_files; 
  if ( ae_common::allow_plasmids )
  {
    num_files = 12; 
  }
  else 
  {
    num_files = 4; 
  }
  
  char* tmp_file_name; 
  for ( int8_t i = 0 ; i < num_files ; i++ )
  {
    strncpy(old_stat_file_glob_name, _stat_files_glob_names[i], 60);
    old_stat_file_glob_name[59] = '\0';
    strcat(old_stat_file_glob_name, ".old");
    
    // Backup stat files
    rename ( _stat_files_glob_names[i], old_stat_file_glob_name );
    // Open backed up files and new stat files
    old_stat_file_glob = fopen( old_stat_file_glob_name, "r" );
    assert( old_stat_file_glob );
    _stat_files_glob[i] = fopen( _stat_files_glob_names[i], "w" );
    assert( _stat_files_glob[i] ); 

    strncpy( old_stat_file_best_name, _stat_files_best_names[i], 60 );
    old_stat_file_best_name[59] = '\0';
    strcat( old_stat_file_best_name, ".old" );
    
    // Backup stat files
    rename ( _stat_files_best_names[i], old_stat_file_best_name );
    // Open backed up files and new stat files
    old_stat_file_best = fopen( old_stat_file_best_name, "r" );
    assert( old_stat_file_best );
    _stat_files_best[i] = fopen( _stat_files_best_names[i], "w" );
    assert( _stat_files_best[i] );
    
    // Copy header for each BEST file
    ret = fgets( line, 500, old_stat_file_best );
    while ( !feof( old_stat_file_best ) && line[0] == '#' )
    {
      fputs( line, _stat_files_best[i] );
      ret = fgets( line, 500, old_stat_file_best );
    }
    // This is the empty line between the header and the values
    fputs( line, _stat_files_best[i] );

    // Copy header for each GLOB file
    ret = fgets( line, 500, old_stat_file_glob );
    while ( !feof( old_stat_file_glob ) && line[0] == '#' )
    {
      fputs( line, _stat_files_glob[i] );
      ret = fgets( line, 500, old_stat_file_glob );
    }
    // This is the empty line between the header and the values
    fputs( line, _stat_files_glob[i] );
    
    // Copy stats until num_gener (included) in BEST files
    ret = fgets( line, 500, old_stat_file_best );
    while ( (int32_t)atol(line) <= ae_common::sim->get_first_gener() && !feof(old_stat_file_best) )
    {
      fputs( line, _stat_files_best[i] );
      ret = fgets( line, 500, old_stat_file_best );
    }
    
    // Copy stats until num_gener (included) in GLOB files
    ret = fgets( line, 500, old_stat_file_glob );
    while ( (int32_t)atol(line) <= ae_common::sim->get_first_gener() && !feof(old_stat_file_glob) )
    {
      fputs( line, _stat_files_glob[i] );
      ret = fgets( line, 500, old_stat_file_glob );
    }
    
    fclose( old_stat_file_best ); 
    fclose( old_stat_file_glob ); 
  }  

  // Flush the new stat files
  flush();
}

// =================================================================
//                             Destructors
// =================================================================
ae_stats::~ae_stats( void )
{
  if(_stat_files_best != NULL)
    {
      for (int8_t i = 0; i < 4; i++) 
        {
          if ( _stat_files_best[i] != NULL ) fclose ( _stat_files_best[i] ); 
          
        }  
      if ( ae_common::allow_plasmids )
        { 
          for (int8_t i = 4; i < 12; i++) 
            {
              if ( _stat_files_best[i] != NULL ) fclose ( _stat_files_best[i] ); 
            }  
        }
      delete [] _stat_files_best;
      _stat_files_best = NULL;
    }

  if(_stat_files_glob != NULL)
    {
      for (int8_t i = 0; i < 4; i++) 
        {
          if ( _stat_files_glob[i] != NULL ) fclose ( _stat_files_glob[i] ); 
          
        }  
      if ( ae_common::allow_plasmids )
        { 
          for (int8_t i = 4; i < 12; i++) 
            {
              if ( _stat_files_glob[i] != NULL ) fclose ( _stat_files_glob[i] ); 
            }  
        }
      delete [] _stat_files_glob;
      _stat_files_glob = NULL;
    }


  if (_stat_files_best_names != NULL)
    {
      for (int8_t i = 0; i < 12; i++)  
        {
          if (_stat_files_best_names[i] != NULL) 
            {
              delete [] _stat_files_best_names[i];
            }
        }
      delete [] _stat_files_best_names;
      _stat_files_best_names = NULL;
    }

  if (_stat_files_glob_names != NULL)
    {
      for (int8_t i = 0; i < 12; i++)  
        {
          if (_stat_files_glob_names[i] != NULL) 
            {
              delete [] _stat_files_glob_names[i];
            }
        }
      delete [] _stat_files_glob_names;
      _stat_files_glob_names = NULL;
    }

}   


// =================================================================
//                            Public Methods
// =================================================================

inline double sqr(double x) {
  return x*x;
}

inline double rsqr(double x) {
  return (x < 0.00000001) ? 0. : sqrt(x);
}


void ae_stats::write_headers( void )
{
  int8_t i;
  
  for ( int8_t j = 0; j < 12; j++) 
  {
    if ( j == 0 || ( ae_common::allow_plasmids && (j == 4 || j == 8) ) )
    {
      if (_stat_files_best != NULL)
      {
        i = 1; 
        write_header( _stat_files_best[j], "-------------------------------------" );
        write_header( _stat_files_best[j], "Fittest individual fitness statistics" );
        write_header( _stat_files_best[j], "-------------------------------------" );
        write_header( _stat_files_best[j], "" );
        write_header( _stat_files_best[j], "Generation", i++ );
        write_header( _stat_files_best[j], "Population size", i++ );
        write_header( _stat_files_best[j], "Fitness", i++ );
        write_header( _stat_files_best[j], "Genome size (amount of DNA)", i++ );
        write_header( _stat_files_best[j], "Metabolic error", i++ );
        write_header( _stat_files_best[j], "Parent's metabolic error", i++ );
        write_header( _stat_files_best[j], "Secretion error", i++ );
        write_header( _stat_files_best[j], "Parent's secretion error", i++ );
        write_header( _stat_files_best[j], "Amount of compound secreted by an individual", i++ );
        write_header( _stat_files_best[j], "Amount of compound present in the grid-cell", i++ );

#ifdef __REGUL
        write_header( _stat_files_best[j], "Number of links in the regulation graph", i++ );
        write_header( _stat_files_best[j], "Number of positive links in the regulation graph", i++ );
        write_header( _stat_files_best[j], "Number of negative links in the regulation graph", i++ );
        write_header( _stat_files_best[j], "Average value of links in the regulation graph", i++ );
        write_header( _stat_files_best[j], "Average value of positive links in the regulation graph", i++ );
        write_header( _stat_files_best[j], "Average value of negative links in the regulation graph", i++ );
#endif
        write_header( _stat_files_best[j], ""); 
      }


      if ( _stat_files_glob != NULL ) // remember, no glob files for ancstats
      {
        i = 1; 
        write_header( _stat_files_glob[j], "--------------------------" );
        write_header( _stat_files_glob[j], "Average fitness statistics" );
        write_header( _stat_files_glob[j], "--------------------------" );
        write_header( _stat_files_glob[j], "" );
        write_header( _stat_files_glob[j], "Generation", i++ );
        write_header( _stat_files_glob[j], "Population size", i++ );
        write_header( _stat_files_glob[j], "Fitness", i++ );
        write_header( _stat_files_glob[j], "Genome size, (amount of DNA)", i++ );
        write_header( _stat_files_glob[j], "Metabolic error", i++ );
        write_header( _stat_files_glob[j], "Parent's metabolic error", i++ );
        write_header( _stat_files_glob[j], "Secretion error", i++ );
        write_header( _stat_files_glob[j], "Parent's secretion error", i++ );
        write_header( _stat_files_glob[j], "Amount of compound secreted by an individual", i++ );
        write_header( _stat_files_glob[j], "Amount of compound present in a grid-cell", i++ );
        
#ifdef __REGUL
        write_header( _stat_files_glob[j], "Number of links in the regulation graph", i++ );
        write_header( _stat_files_glob[j], "Number of positive links in the regulation graph", i++ );
        write_header( _stat_files_glob[j], "Number of negative links in the regulation graph", i++ );
        write_header( _stat_files_glob[j], "Average value of links in the regulation graph", i++ );
        write_header( _stat_files_glob[j], "Average value of positive links in the regulation graph", i++ );
        write_header( _stat_files_glob[j], "Average value of negative links in the regulation graph", i++ );
#endif
        write_header( _stat_files_glob[j], "" );
      }
    }

    if ( j == 1 || ( ae_common::allow_plasmids && (j == 5 || j == 9) ) )
    {
      if ( _stat_files_glob != NULL ) // remember, no glob files for ancstats
      {
        i = 1; 
        write_header( _stat_files_glob[j], "---------------------------" );
        write_header( _stat_files_glob[j], "Average mutation statistics" );
        write_header( _stat_files_glob[j], "---------------------------" );
        write_header( _stat_files_glob[j], "" );
        write_header( _stat_files_glob[j], "Generation", i++ );
        write_header( _stat_files_glob[j], "Number of local mutations undergone", i++ );
        write_header( _stat_files_glob[j], "Number of chromosomic rearrangements undergone", i++ );
        write_header( _stat_files_glob[j], "Number of switch undergone", i++ );
        write_header( _stat_files_glob[j], "Number of indels undergone", i++ );
        write_header( _stat_files_glob[j], "Number of duplications undergone", i++ );
        write_header( _stat_files_glob[j], "Number of deletions undergone", i++ );
        write_header( _stat_files_glob[j], "Number of translocations undergone", i++ );
        write_header( _stat_files_glob[j], "Number of inversions undergone", i++ );
        write_header( _stat_files_glob[j], "" );
      }

      if ( _stat_files_best != NULL )
      {
        i = 1; 
        write_header( _stat_files_best[j], "--------------------------------------" );
        write_header( _stat_files_best[j], "Fittest individual mutation statistics" );
        write_header( _stat_files_best[j], "--------------------------------------" );
        write_header( _stat_files_best[j], "" );
        write_header( _stat_files_best[j], "Generation", i++ );
        write_header( _stat_files_best[j], "Number of local mutations undergone", i++ );
        write_header( _stat_files_best[j], "Number of chromosomic rearrangements undergone", i++ );
        write_header( _stat_files_best[j], "Number of switch undergone", i++ );
        write_header( _stat_files_best[j], "Number of indels undergone", i++ );
        write_header( _stat_files_best[j], "Number of duplications undergone", i++ );
        write_header( _stat_files_best[j], "Number of deletions undergone", i++ );
        write_header( _stat_files_best[j], "Number of translocations undergone", i++ );
        write_header( _stat_files_best[j], "Number of inversions undergone", i++ );
        write_header( _stat_files_best[j], "" );
      }
    }
    
    if ( j == 2 || ( ae_common::allow_plasmids && (j == 6 || j == 10) ) )
    {
      if ( _stat_files_glob != NULL ) // remember, no glob files for ancstats
      {
        i = 1; 
        write_header( _stat_files_glob[j], "" );
        write_header( _stat_files_glob[j], "-----------------------" );
        write_header( _stat_files_glob[j], "Average gene statistics" );
        write_header( _stat_files_glob[j], "-----------------------" );
        write_header( _stat_files_glob[j], "" );
        write_header( _stat_files_glob[j], "Generation", i++ );
        write_header( _stat_files_glob[j], "Number of coding RNAs (at least one gene on RNA)", i++ );
        write_header( _stat_files_glob[j], "Number of non-coding RNAs", i++ );
        write_header( _stat_files_glob[j], "Average size of coding RNAs (at least one gene on RNA)", i++ );
        write_header( _stat_files_glob[j], "Average size of non-coding RNAs", i++ );
        write_header( _stat_files_glob[j], "Number of functional genes", i++ );
        // Non functional genes are those with _width == 0 or _height == 0 or those that lack one kind of codons (M, W or H)
        write_header( _stat_files_glob[j], "Nb of non functional genes", i++ );
        write_header( _stat_files_glob[j], "Average size of functional genes", i++ );
        write_header( _stat_files_glob[j], "Average size of non functional genes (WARNING : bias towards 0)", i++ );
        write_header( _stat_files_glob[j], "" );
      }

      if ( _stat_files_best != NULL )
      {
        i = 1; 
        write_header( _stat_files_best[j], "----------------------------------" );
        write_header( _stat_files_best[j], "Fittest individual gene statistics" );
        write_header( _stat_files_best[j], "----------------------------------" );
        write_header( _stat_files_best[j], "" );      
        write_header( _stat_files_best[j], "Generation", i++ );
        write_header( _stat_files_best[j], "Number of coding RNAs (at least one gene on RNA)", i++ );
        write_header( _stat_files_best[j], "Number of non-coding RNAs", i++ );
        write_header( _stat_files_best[j], "Average size of coding RNAs (at least one gene on RNA)", i++ );
        write_header( _stat_files_best[j], "Average size of non-coding RNAs", i++ );
        write_header( _stat_files_best[j], "Number of functional genes", i++ );
        // Non functional genes are those with _width == 0 or _height == 0 or those that lack one kind of codons (M, W or H)
        write_header( _stat_files_best[j], "Nb of non functional genes", i++ );
        write_header( _stat_files_best[j], "Average size of functional genes", i++ );
        write_header( _stat_files_best[j], "Average size of non functional genes (WARNING : bias towards 0)", i++ );
        write_header( _stat_files_best[j], "" );
      }
    }



    if ( j == 3 || ( ae_common::allow_plasmids && (j == 7 || j == 11) ) )
    {
      if ( _stat_files_glob != NULL ) // remember, no glob files for ancstats
      {
        // TO DO (if needed) . This file is not yet implemented 
        i = 1; 
        write_header( _stat_files_glob[j], "----------------------------" );
        write_header( _stat_files_glob[j], "Average base pair statistics" );
        write_header( _stat_files_glob[j], "----------------------------" );
        write_header( _stat_files_glob[j], "" );
        write_header( _stat_files_glob[j], " This data is not available"); 
        write_header( _stat_files_glob[j], " Computing bp stats for all individuals is extremely costly computationaly" );
        write_header( _stat_files_glob[j], "" ); 
        // write_header( _stat_files_glob[j], "Generation", i++ );
        // write_header( _stat_files_glob[j], "Number of bp not included in any CDS", i++ );
        // write_header( _stat_files_glob[j], "Number of bp not included in any functional CDS", i++ );
        // write_header( _stat_files_glob[j], "Number of bp not included in any non functional CDS", i++ );
        // write_header( _stat_files_glob[j], "Number of bp not included in any RNA", i++ );
        // write_header( _stat_files_glob[j], "Number of bp not included in any coding RNA", i++ );
        // write_header( _stat_files_glob[j], "Number of bp not included in any non coding RNA", i++ );
        // write_header( _stat_files_glob[j], "" );
      }

      if ( _stat_files_best!= NULL )
      {          
        i = 1; 
        write_header( _stat_files_best[j], "---------------------------------------" );
        write_header( _stat_files_best[j], "Fittest individual base pair statistics" );
        write_header( _stat_files_best[j], "---------------------------------------" );
        write_header( _stat_files_best[j], "" ); 
        write_header( _stat_files_best[j], "Generation", i++ );
        write_header( _stat_files_best[j], "Number of bp not included in any CDS", i++ );
        write_header( _stat_files_best[j], "Number of bp not included in any functional CDS", i++ );
        write_header( _stat_files_best[j], "Number of bp not included in any non functional CDS", i++ );
        write_header( _stat_files_best[j], "Number of bp not included in any RNA", i++ );
        write_header( _stat_files_best[j], "Number of bp not included in any coding RNA", i++ );
        write_header( _stat_files_best[j], "Number of bp not included in any non coding RNA", i++ );
        write_header( _stat_files_best[j], "" );
      }
    }
  }    
  flush();
}




// For post-treatments
void ae_stats::write_statistics_of_this_indiv( ae_individual * indiv, int32_t t )
{
  if ( ! ae_common::allow_plasmids )
  {
    // ----------------------------------------------------------
    // Compute and write statistical data for the best individual
    // ----------------------------------------------------------
    ae_stat_record* stat_record = new ae_stat_record( indiv, CHROM, true, t );
    stat_record->write_to_file( _stat_files_best[0], FITNESS_STATS  );
    stat_record->write_to_file( _stat_files_best[1], MUTATION_STATS );
    stat_record->write_to_file( _stat_files_best[2], GENES_STATS    );
    stat_record->write_to_file( _stat_files_best[3], BP_STATS       );
    delete stat_record;
  }
  else
  {
    ae_stat_record* stat_record;

    // i: 0 = ALL GUs; 1 = CHROM; 2 = PLASMIDS
    for ( int8_t i = 0 ; i < 3 ; i++ )
    {
      stat_record = new ae_stat_record( ae_common::sim->get_pop()->get_best(), (chrom_or_gen_unit) i, true, t ); 
      // j : 0 = FITNESS_STATS; 1 = MUTATION_STATS; 2 = GENES_STATS; 3 = BP_STATS
      for ( int8_t j = 0 ; j < 4 ; j++ )
      {
        stat_record->write_to_file( _stat_files_best[i*4+j], (stats_type) j );
      }
    }

    delete stat_record;
  }
}


void ae_stats::write_current_generation_statistics( void )
{
  ae_stat_record* stat_record;
  
  if ( ! ae_common::allow_plasmids )
  {
    // ----------------------------------------------------------
    // Compute and write statistical data for the best individual
    // ----------------------------------------------------------
    stat_record = new ae_stat_record( ae_common::sim->get_pop()->get_best() );
    stat_record->write_to_file( _stat_files_best[0], FITNESS_STATS  );
    stat_record->write_to_file( _stat_files_best[1], MUTATION_STATS );
    stat_record->write_to_file( _stat_files_best[2], GENES_STATS    );
    stat_record->write_to_file( _stat_files_best[3], BP_STATS       );
    
    delete stat_record;
 
    // -----------------------------------------------------
    // Compute and write statistical data for the population
    // -----------------------------------------------------
    stat_record = new ae_stat_record( ae_common::sim->get_pop() );
    stat_record->write_to_file( _stat_files_glob[0], FITNESS_STATS  );
    stat_record->write_to_file( _stat_files_glob[1], MUTATION_STATS );
    stat_record->write_to_file( _stat_files_glob[2], GENES_STATS    );
    // NB: It's too expensive to compute bp stats for every indivs in the pop!
    //~ stat_record->write_to_file( _stat_files_glob[3], BP_STATS       );

    delete stat_record;
  }
  else
  {
    // i: 0 = ALL GUs; 1 = CHROM; 2 = PLASMIDS
    for ( int8_t i = 0 ; i < 3 ; i++ )
    {
      stat_record = new ae_stat_record( ae_common::sim->get_pop()->get_best(), (chrom_or_gen_unit) i ); 

      // j : 0 = FITNESS_STATS; 1 = MUTATION_STATS; 2 = GENES_STATS; 3 = BP_STATS
      for ( int8_t j = 0 ; j < 4 ; j++ )
      {
        stat_record->write_to_file( _stat_files_best[i*4+j], (stats_type) j );
      }
      delete stat_record;

      stat_record = new ae_stat_record( ae_common::sim->get_pop(), (chrom_or_gen_unit) i );    
      // j : 0 = FITNESS_STATS; 1 = MUTATION_STATS; 2 = GENES_STATS; 3 = BP_STATS
      // TO DO (if needed): add global base-pair stats; for now, avoid it
      for ( int8_t j = 0 ; j < 3 ; j++ )
      {
        stat_record->write_to_file( _stat_files_glob[i*4+j], (stats_type) j );
      }
      delete stat_record;
    }  
  }
}


void ae_stats::flush( void )
{

  int16_t num_files;

  if ( ae_common::allow_plasmids ) 
  {
    num_files = 12; 
  } 
  else
  {
    num_files = 4; 
  }
    

  if (_stat_files_glob != NULL) // remember, no glob files for ancstats
    {
      for ( int16_t i = 0 ; i < num_files ; i++ )
        {
          if ( _stat_files_glob[i] != NULL)   fflush( _stat_files_glob[i] );
        }
    }

  if (_stat_files_best != NULL) // remember, no glob files for ancstats
    {
      for ( int16_t i = 0 ; i < num_files ; i++ )
        {
          if ( _stat_files_best[i] != NULL)   fflush( _stat_files_best[i] );
        }
    }
    
}

// =================================================================
//                           Protected Methods
// =================================================================
