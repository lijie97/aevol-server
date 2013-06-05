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
#include <err.h>
#include <errno.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_stats.h>
#include <ae_stat_record.h>
#include <ae_exp_manager.h>
#include <ae_exp_setup.h>
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
/*!
  Create a NEW stat manager
 */
ae_stats::ae_stats( ae_exp_manager* exp_m,
                    const char * prefix /* = "stat" */,
                    bool best_indiv_only /* = false */ )
{
  _exp_m = exp_m;
  init_data();
  set_file_names( prefix, best_indiv_only );
  open_files();
  write_headers();
}

/*!
  Create a stat manager to append existing stats
 */
ae_stats::ae_stats( ae_exp_manager* exp_m,
                    int32_t num_gener,
                    const char * prefix /* = "stat" */,
                    bool best_indiv_only /* = false */,
                    bool addition_old_stats /* = true */,
                    bool delete_old_stats /* = true */)
{
  _exp_m = exp_m;
  init_data();
  set_file_names( prefix, best_indiv_only );
  
  // ---------------------------------------------------------------------------
  //  Make a backup copy (named <original_name>.old) of each file
  //  and copy its content into the new stat file untill <num_gener> is reached
  // ---------------------------------------------------------------------------
  if(addition_old_stats)
  {
    char* old_file_name = new char[100];
    FILE* old_file;
    char* cur_file_name;  // Syntaxic sugar for _stat_files_names[][][]
    FILE* cur_file;       // Syntaxic sugar for _stat_files[][][]
    char  line[500];
    char* trash;
    
    for ( int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++ )
    { 
      for ( int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++ )
      {
        for ( int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++ )
        {
          cur_file_name = _stat_files_names[chrom_or_GU][best_or_glob][stat_type];
          if ( cur_file_name != NULL )
          {
            sprintf( old_file_name, "%s.old", cur_file_name );
            int8_t exist_file = rename( cur_file_name, old_file_name );
            if ( exist_file != 0 )
            {
              printf( "ERROR : Could not rename %s as %s.\n", cur_file_name, old_file_name );
              exit( EXIT_FAILURE );
            }
            
            old_file = fopen( old_file_name, "r" );
            cur_file = fopen( cur_file_name, "w" );
            
            // Copy file header
            trash = fgets( line, 500, old_file );
            while ( !feof( old_file ) && line[0] == '#' )
            {
              fputs( line, cur_file );
              trash = fgets( line, 500, old_file );
            }
            
            // Copy the empty line between the header and the values
            fputs( line, cur_file );
            
            // Copy stats until num_gener (included)
            trash = fgets( line, 500, old_file );
            while ( (int32_t)atol(line) < _exp_m->get_first_gener() && !feof(old_file) )
            {
              fputs( line, cur_file );
              trash = fgets( line, 500, old_file );
            }
            
            fclose( old_file );
            
            _stat_files[chrom_or_GU][best_or_glob][stat_type] = cur_file;
            
            if ( delete_old_stats )
            {
              remove( old_file_name );
            }
          }
        }
      }
    }
    delete old_file_name;
  }
  else // ancstat case
  {
    open_files();
    write_headers(true);
  }

  // Flush the new stat files
  flush();
}

// =================================================================
//                             Destructors
// =================================================================
ae_stats::~ae_stats( void )
{
  for ( int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++ )
  {
    for ( int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++ )
    {
      for ( int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++ )
      {
        if( _stat_files_names[chrom_or_GU][best_or_glob][stat_type] != NULL )
        {
          assert( _stat_files[chrom_or_GU][best_or_glob][stat_type] != NULL );
          
          fclose( _stat_files[chrom_or_GU][best_or_glob][stat_type] );
          _stat_files[chrom_or_GU][best_or_glob][stat_type] = NULL;
          
          delete [] _stat_files_names[chrom_or_GU][best_or_glob][stat_type];
          _stat_files_names[chrom_or_GU][best_or_glob][stat_type] = NULL;
        }
      }
      
      delete [] _stat_files[chrom_or_GU][best_or_glob];
      _stat_files[chrom_or_GU][best_or_glob] = NULL;
      
      delete [] _stat_files_names[chrom_or_GU][best_or_glob];
      _stat_files_names[chrom_or_GU][best_or_glob] = NULL;
    }
    
    delete [] _stat_files[chrom_or_GU];
    _stat_files[chrom_or_GU] = NULL;
    
    delete [] _stat_files_names[chrom_or_GU];
    _stat_files_names[chrom_or_GU] = NULL;
  }
  
  delete [] _stat_files;
  _stat_files = NULL;
  
  delete [] _stat_files_names;
  _stat_files_names = NULL;
}   


// =================================================================
//                            Public Methods
// =================================================================

inline double sqr( double x )
{
  return x*x;
}

inline double rsqr( double x )
{
  return (x < 0.00000001) ? 0. : sqrt(x);
}


void ae_stats::write_headers( bool ancstats_stats /* = false */ )
{
  // Column key in the stat files
  int8_t key;
  
  // --------------------------------------
  //  Write headers in FITNESS_STATS files
  // --------------------------------------
  for ( int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++ )
  {
    if( ancstats_stats)
    {
      write_header( _stat_files[chrom_or_GU][BEST][FITNESS_STATS], "----------------------------------------" );
      write_header( _stat_files[chrom_or_GU][BEST][FITNESS_STATS], " Lineage individuals fitness statistics " );
      write_header( _stat_files[chrom_or_GU][BEST][FITNESS_STATS], "----------------------------------------" );
      write_header( _stat_files[chrom_or_GU][BEST][FITNESS_STATS], "" );
    }
    else
    {
      if ( _stat_files_names[chrom_or_GU][BEST][FITNESS_STATS] != NULL )
      {
        write_header( _stat_files[chrom_or_GU][BEST][FITNESS_STATS], "---------------------------------------" );
        write_header( _stat_files[chrom_or_GU][BEST][FITNESS_STATS], " Fittest individual fitness statistics " );
        write_header( _stat_files[chrom_or_GU][BEST][FITNESS_STATS], "---------------------------------------" );
        write_header( _stat_files[chrom_or_GU][BEST][FITNESS_STATS], "" );
      }
      if ( _stat_files_names[chrom_or_GU][GLOB][FITNESS_STATS] != NULL )
      {
        write_header( _stat_files[chrom_or_GU][GLOB][FITNESS_STATS], "------------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][GLOB][FITNESS_STATS], " Average fitness statistics over the population " );
        write_header( _stat_files[chrom_or_GU][GLOB][FITNESS_STATS], "------------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][GLOB][FITNESS_STATS], "" );
      }
      
      if ( _stat_files_names[chrom_or_GU][SDEV][FITNESS_STATS] != NULL )
      {
        write_header( _stat_files[chrom_or_GU][SDEV][FITNESS_STATS], "------------------------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][SDEV][FITNESS_STATS], " Standard deviation, fitness statistics over the population " );
        write_header( _stat_files[chrom_or_GU][SDEV][FITNESS_STATS], "------------------------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][SDEV][FITNESS_STATS], "" );
      }
      
      if ( _stat_files_names[chrom_or_GU][SKEW][FITNESS_STATS] != NULL )
      {
        write_header( _stat_files[chrom_or_GU][SKEW][FITNESS_STATS], "--------------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][SKEW][FITNESS_STATS], " Skewness statistics, fitness over the population " );
        write_header( _stat_files[chrom_or_GU][SKEW][FITNESS_STATS], "--------------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][SKEW][FITNESS_STATS], "" );
      }
    }
    
    
    for ( int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++ )
    {
      if ( _stat_files_names[chrom_or_GU][best_or_glob][FITNESS_STATS] != NULL )
      {
        assert( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS] != NULL );
        key = 1;
        
        write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Generation", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Population size", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Fitness", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Genome size (amount of DNA)", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Metabolic error", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Parent's metabolic error", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Metabolic fitness", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Secretion error", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Parent's secretion error", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Secretion fitness", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Amount of compound present in the grid-cell", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Int probe", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Double probe", key++ );

        #ifdef __REGUL
          write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Number of links in the regulation graph", key++ );
          write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Number of positive links in the regulation graph", key++ );
          write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Number of negative links in the regulation graph", key++ );
          write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Average value of links in the regulation graph", key++ );
          write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Average value of positive links in the regulation graph", key++ );
          write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], "Average value of negative links in the regulation graph", key++ );
        #endif
        
        write_header( _stat_files[chrom_or_GU][best_or_glob][FITNESS_STATS], ""); 
      }
    }
  }
  
  // ---------------------------------------
  //  Write headers in MUTATION_STATS files
  // ---------------------------------------
  for ( int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++ )
  {
    if( ancstats_stats)
    {
      write_header( _stat_files[chrom_or_GU][BEST][MUTATION_STATS], "-----------------------------------------" );
      write_header( _stat_files[chrom_or_GU][BEST][MUTATION_STATS], " Lineage individuals mutation statistics " );
      write_header( _stat_files[chrom_or_GU][BEST][MUTATION_STATS], "-----------------------------------------" );
      write_header( _stat_files[chrom_or_GU][BEST][MUTATION_STATS], "" );
    }
    else
    {
      if ( _stat_files_names[chrom_or_GU][BEST][MUTATION_STATS] != NULL )
      {
        write_header( _stat_files[chrom_or_GU][BEST][MUTATION_STATS], "----------------------------------------" );
        write_header( _stat_files[chrom_or_GU][BEST][MUTATION_STATS], " Fittest individual mutation statistics " );
        write_header( _stat_files[chrom_or_GU][BEST][MUTATION_STATS], "----------------------------------------" );
        write_header( _stat_files[chrom_or_GU][BEST][MUTATION_STATS], "" );
      }
      if ( _stat_files_names[chrom_or_GU][GLOB][MUTATION_STATS] != NULL )
      {
        write_header( _stat_files[chrom_or_GU][GLOB][MUTATION_STATS], "-------------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][GLOB][MUTATION_STATS], " Average mutation statistics over the population " );
        write_header( _stat_files[chrom_or_GU][GLOB][MUTATION_STATS], "-------------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][GLOB][MUTATION_STATS], "" );
      }
    }
    
    for ( int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++ )
    {
      if ( _stat_files_names[chrom_or_GU][best_or_glob][MUTATION_STATS] != NULL )
      {
        assert( _stat_files[chrom_or_GU][best_or_glob][MUTATION_STATS] != NULL );
        key = 1;
        
        write_header( _stat_files[chrom_or_GU][best_or_glob][MUTATION_STATS], "Generation", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of local mutations undergone", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of chromosomic rearrangements undergone", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of switch undergone", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of indels undergone", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of duplications undergone", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of deletions undergone", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of translocations undergone", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of inversions undergone", key++ );
        
        write_header( _stat_files[chrom_or_GU][best_or_glob][MUTATION_STATS], ""); 
      }
    }
  }
  
  // ---------------------------------------
  //  Write headers in GENES_STATS files
  // ---------------------------------------
  for ( int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++ )
  {
    if( ancstats_stats)
    {
      write_header( _stat_files[chrom_or_GU][BEST][GENES_STATS], "-------------------------------------" );
      write_header( _stat_files[chrom_or_GU][BEST][GENES_STATS], " Lineage individuals gene statistics " );
      write_header( _stat_files[chrom_or_GU][BEST][GENES_STATS], "-------------------------------------" );
      write_header( _stat_files[chrom_or_GU][BEST][GENES_STATS], "" );
    }
    else
    {
      if ( _stat_files_names[chrom_or_GU][BEST][GENES_STATS] != NULL )
      {
        write_header( _stat_files[chrom_or_GU][BEST][GENES_STATS], "------------------------------------" );
        write_header( _stat_files[chrom_or_GU][BEST][GENES_STATS], " Fittest individual gene statistics " );
        write_header( _stat_files[chrom_or_GU][BEST][GENES_STATS], "------------------------------------" );
        write_header( _stat_files[chrom_or_GU][BEST][GENES_STATS], "" );
      }
      if ( _stat_files_names[chrom_or_GU][GLOB][GENES_STATS] != NULL )
      {
        write_header( _stat_files[chrom_or_GU][GLOB][GENES_STATS], "---------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][GLOB][GENES_STATS], " Average gene statistics over the population " );
        write_header( _stat_files[chrom_or_GU][GLOB][GENES_STATS], "---------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][GLOB][GENES_STATS], "" );
      }
    }
    
    for ( int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++ )
    {
      if ( _stat_files_names[chrom_or_GU][best_or_glob][GENES_STATS] != NULL )
      {
        assert( _stat_files[chrom_or_GU][best_or_glob][GENES_STATS] != NULL );
        key = 1;
        
        write_header( _stat_files[chrom_or_GU][best_or_glob][GENES_STATS], "Generation", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][GENES_STATS], "Number of coding RNAs (at least one gene on RNA)", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][GENES_STATS], "Number of non-coding RNAs", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][GENES_STATS], "Average size of coding RNAs (at least one gene on RNA)", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][GENES_STATS], "Average size of non-coding RNAs", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][GENES_STATS], "Number of functional genes", key++ );
        // Non functional genes are those with _width == 0 or _height == 0 or those that lack one kind of codons (M, W or H)
        write_header( _stat_files[chrom_or_GU][best_or_glob][GENES_STATS], "Nb of non functional genes", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][GENES_STATS], "Average size of functional genes", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][GENES_STATS], "Average size of non functional genes (WARNING : bias towards 0)", key++ );
        
        write_header( _stat_files[chrom_or_GU][best_or_glob][GENES_STATS], ""); 
      }
    }
  }
  
  // ---------------------------------------
  //  Write headers in BP_STATS files
  // ---------------------------------------
  for ( int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++ )
  {
    if( ancstats_stats)
    {
      write_header( _stat_files[chrom_or_GU][BEST][BP_STATS], "-------------------------------------------" );
      write_header( _stat_files[chrom_or_GU][BEST][BP_STATS], " Lineage individuals non-coding statistics " );
      write_header( _stat_files[chrom_or_GU][BEST][BP_STATS], "-------------------------------------------" );
      write_header( _stat_files[chrom_or_GU][BEST][BP_STATS], "" );
    }
    else
    {
      if ( _stat_files_names[chrom_or_GU][BEST][BP_STATS] != NULL )
      {
        write_header( _stat_files[chrom_or_GU][BEST][BP_STATS], "------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][BEST][BP_STATS], " Fittest individual non-coding statistics " );
        write_header( _stat_files[chrom_or_GU][BEST][BP_STATS], "------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][BEST][BP_STATS], "" );
      }
      if ( _stat_files_names[chrom_or_GU][GLOB][BP_STATS] != NULL )
      {
        write_header( _stat_files[chrom_or_GU][GLOB][BP_STATS], "---------------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][GLOB][BP_STATS], " Average non-coding statistics over the population " );
        write_header( _stat_files[chrom_or_GU][GLOB][BP_STATS], "---------------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][GLOB][BP_STATS], "" );
        write_header( _stat_files[chrom_or_GU][GLOB][BP_STATS], " This data is not available"); 
        write_header( _stat_files[chrom_or_GU][GLOB][BP_STATS], " Computing bp stats for all individuals is extremely costly computationaly" );
        write_header( _stat_files[chrom_or_GU][GLOB][BP_STATS], "" );
        
        // Mark file as "not to be written into" and close it
        delete [] _stat_files_names[chrom_or_GU][GLOB][BP_STATS];
        _stat_files_names[chrom_or_GU][GLOB][BP_STATS] = NULL;
        fclose( _stat_files[chrom_or_GU][GLOB][BP_STATS] );
        _stat_files[chrom_or_GU][GLOB][BP_STATS] = NULL;
      }
    }
    
    for ( int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++ )
    {
      if ( _stat_files_names[chrom_or_GU][best_or_glob][BP_STATS] != NULL )
      {
        assert( _stat_files[chrom_or_GU][best_or_glob][BP_STATS] != NULL );
        key = 1;
        
        write_header( _stat_files[chrom_or_GU][best_or_glob][BP_STATS], "Generation", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][BP_STATS], "Number of bp not included in any CDS", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][BP_STATS], "Number of bp not included in any functional CDS", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][BP_STATS], "Number of bp not included in any non functional CDS", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][BP_STATS], "Number of bp not included in any RNA", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][BP_STATS], "Number of bp not included in any coding RNA", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][BP_STATS], "Number of bp not included in any non coding RNA", key++ );
        
        write_header( _stat_files[chrom_or_GU][best_or_glob][BP_STATS], ""); 
      }
    }
  }
  
  // ---------------------------------------
  //  Write headers in REAR_STATS files
  // ---------------------------------------
  for ( int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++ )
  {
    if( ancstats_stats)
    {
      write_header( _stat_files[chrom_or_GU][BEST][REAR_STATS], "----------------------------------------------" );
      write_header( _stat_files[chrom_or_GU][BEST][REAR_STATS], " Lineage individuals rearrangement statistics " );
      write_header( _stat_files[chrom_or_GU][BEST][REAR_STATS], "----------------------------------------------" );
      write_header( _stat_files[chrom_or_GU][BEST][REAR_STATS], "" );
    }
    else
    {
      if ( _stat_files_names[chrom_or_GU][BEST][REAR_STATS] != NULL )
      {
        write_header( _stat_files[chrom_or_GU][BEST][REAR_STATS], "---------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][BEST][REAR_STATS], " Fittest individual rearrangement statistics " );
        write_header( _stat_files[chrom_or_GU][BEST][REAR_STATS], "---------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][BEST][REAR_STATS], "" );
      }
      if ( _stat_files_names[chrom_or_GU][GLOB][REAR_STATS] != NULL )
      {
        write_header( _stat_files[chrom_or_GU][GLOB][REAR_STATS], "------------------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][GLOB][REAR_STATS], " Average rearrangement statistics over the population " );
        write_header( _stat_files[chrom_or_GU][GLOB][REAR_STATS], "------------------------------------------------------" );
        write_header( _stat_files[chrom_or_GU][GLOB][REAR_STATS], "" );
      }
    }
    
    for ( int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++ )
    {
      if ( _stat_files_names[chrom_or_GU][best_or_glob][REAR_STATS] != NULL )
      {
        assert( _stat_files[chrom_or_GU][best_or_glob][REAR_STATS] != NULL );
        key = 1;
        
        write_header( _stat_files[chrom_or_GU][best_or_glob][REAR_STATS], "Generation", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][REAR_STATS], "Actual duplication rate", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][REAR_STATS], "Actual deletion rate", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][REAR_STATS], "Actual translocation rate", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][REAR_STATS], "Actual inversion rate", key++ );
        write_header( _stat_files[chrom_or_GU][best_or_glob][REAR_STATS], "Average alignment score (needed score)", key++ );
        
        write_header( _stat_files[chrom_or_GU][best_or_glob][REAR_STATS], ""); 
      }
    }
  }
  
  flush();
}

void ae_stats::write_current_generation_statistics( void )
{
  ae_stat_record** stat_records;
  
  for ( int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++ )
  {
    stat_records = new ae_stat_record* [NB_BEST_OR_GLOB];
    
    stat_records[BEST] = new ae_stat_record( _exp_m, _exp_m->get_best_indiv(), (chrom_or_gen_unit) chrom_or_GU );
    stat_records[GLOB] = new ae_stat_record( _exp_m, _exp_m->get_pop(), (chrom_or_gen_unit) chrom_or_GU );
    stat_records[SDEV] = new ae_stat_record( _exp_m, _exp_m->get_pop(), stat_records[GLOB], (chrom_or_gen_unit) chrom_or_GU );
    stat_records[SKEW] = new ae_stat_record( _exp_m, _exp_m->get_pop(), stat_records[GLOB], stat_records[SDEV], (chrom_or_gen_unit) chrom_or_GU );
    
    for ( int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++ )
    {
      for ( int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++ )
      {
        if ( _stat_files_names[chrom_or_GU][best_or_glob][stat_type] != NULL )
        {
          stat_records[best_or_glob]->write_to_file( _stat_files[chrom_or_GU][best_or_glob][stat_type], (stats_type) stat_type );
        }
      }
      
      delete stat_records[best_or_glob];
    }
    
    delete [] stat_records;
  }
}

void ae_stats::write_statistics_of_this_indiv( ae_individual * indiv, int32_t num_gener )
{
  ae_stat_record* stat_record;
  
  for ( int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++ )
  {
    stat_record = new ae_stat_record( _exp_m, indiv, (chrom_or_gen_unit) chrom_or_GU, true, num_gener );
    
    for ( int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++ )
    {
      if ( _stat_files_names[chrom_or_GU][BEST][stat_type] != NULL )
      {
        assert( _stat_files[chrom_or_GU][BEST][stat_type] != NULL );
        
        stat_record->write_to_file( _stat_files[chrom_or_GU][BEST][stat_type], (stats_type) stat_type );
      }
    }
    
    delete stat_record;
  }
}

void ae_stats::flush( void )
{
  for ( int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++ )
  {
    for ( int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++ )
    {
      for ( int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++ )
      {
        if ( _stat_files_names[chrom_or_GU][best_or_glob][stat_type] != NULL )
        {
          assert( _stat_files[chrom_or_GU][best_or_glob][stat_type] != NULL );
          fflush( _stat_files[chrom_or_GU][best_or_glob][stat_type] );
        }
      }
    }
  }
}

// =================================================================
//                           Protected Methods
// =================================================================
/**
 * Allocate memory and initialize file handlers and file names to NULL
 */
void ae_stats::init_data( void )
{
  _stat_files       = new FILE***[NB_CHROM_OR_GU];
  _stat_files_names = new char***[NB_CHROM_OR_GU];
  
  for ( int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++ )
  {
    _stat_files[chrom_or_GU]        = new FILE**[NB_BEST_OR_GLOB];
    _stat_files_names[chrom_or_GU]  = new char**[NB_BEST_OR_GLOB];
    
    for ( int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++ )
    {
      _stat_files[chrom_or_GU][best_or_glob]        = new FILE*[NB_STATS_TYPES];
      _stat_files_names[chrom_or_GU][best_or_glob]  = new char*[NB_STATS_TYPES];
      
      for ( int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++ )
      {
        _stat_files[chrom_or_GU][best_or_glob][stat_type]       = NULL;
        _stat_files_names[chrom_or_GU][best_or_glob][stat_type] = NULL;
      }
    }
  }
}

/**
 * Construct file names
 * NB: Here is where we "choose" which files we will write.
 *     Files that are not wanted must be left with a NULL name (don't new char[] them)
 *     There is an exception though: for the non-coding file for the population, we will
 *       give it a name temporarily so that we can write the warning headers. Once this
 *       is done, the name will be deleted to mark the file as "not to be written into"
 */
void ae_stats::set_file_names( const char * prefix, bool one_lambda_indiv_only )
{
  // 1) Create stats directory
  int status;
  status = mkdir( STATS_DIR, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, STATS_DIR, errno );
  }
  
  const char* chrom_or_gu_name[NB_CHROM_OR_GU]    = { "", "_chromosome", "_plasmids" };
  const char* best_or_glob_name[NB_BEST_OR_GLOB]  = { "_best", "_glob", "_sdev", "_skew" };
  const char* stat_type_name[NB_STATS_TYPES]      = { "_fitness", "_mutation", "_genes", "_bp", "_rear"};
  
  for ( int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++ )
  {
    for ( int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++ )
    {
      if ( one_lambda_indiv_only && best_or_glob != BEST ) continue;
      
      for ( int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++ )
      {
        //~ // We don't want REAR_STATS when rearrangements are done without alignments
        //~ if ( stat_type == REAR_STATS && ! _exp_m->get_with_alignments() ) continue;
        
        // For now, we only want sdev and skew for fitness data
        if ( best_or_glob > GLOB && stat_type > FITNESS_STATS) continue; 
        if ( (chrom_or_GU != ALL_GU || best_or_glob != GLOB) && stat_type > REAR_STATS) continue;
        
        _stat_files_names[chrom_or_GU][best_or_glob][stat_type] = new char[255];
        
        // Construct the correct name
        if ( one_lambda_indiv_only )
        {
          sprintf(  _stat_files_names[chrom_or_GU][best_or_glob][stat_type],
                    STATS_DIR"/%s%s%s.out",
                    prefix,
                    stat_type_name[stat_type],
                    chrom_or_gu_name[chrom_or_GU] );
        }
        else
        {
          sprintf(  _stat_files_names[chrom_or_GU][best_or_glob][stat_type],
                    STATS_DIR"/%s%s%s%s.out",
                    prefix,
                    stat_type_name[stat_type],
                    chrom_or_gu_name[chrom_or_GU],
                    best_or_glob_name[best_or_glob] );
        }
        
      }
    }
    
  }
}

/**
 * Open files that have a non NULL name
 */
void ae_stats::open_files( void )
{
  for ( int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++ )
  {
    for ( int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++ )
    {
      for ( int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++ )
      {
        if ( _stat_files_names[chrom_or_GU][best_or_glob][stat_type] != NULL )
        {
          _stat_files[chrom_or_GU][best_or_glob][stat_type] = fopen( _stat_files_names[chrom_or_GU][best_or_glob][stat_type], "w" );
        }
      }
    }
  }
}
