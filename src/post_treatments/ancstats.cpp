//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2009  LIRIS.
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
#include <inttypes.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>




// =================================================================
//                            Project Files
// =================================================================
#include <ae_utils.h>
#include <ae_simulation.h>
#include <ae_individual.h>
#include <ae_genetic_unit.h>
#include <ae_list.h>
#include <ae_tree.h>
#include <ae_replication_report.h>
#include <ae_dna_replic_report.h>
#include <ae_mutation.h>
#include <ae_param_loader.h>

//debug
#include <ae_gaussian.h>



void open_environment_stat_file( void );
void write_environment_stats( int32_t num_gener, ae_environment * env );

void open_terminators_stat_file( void );
void write_terminators_stats( int32_t num_gener, ae_individual * indiv );

void open_zones_stat_file( void );
void write_zones_stats( int32_t num_gener, ae_individual * indiv, ae_environment * env );

void open_operons_stat_file( void );
void write_operons_stats( int32_t num_gener, ae_individual * indiv );

void print_help( void );



FILE* env_output_file     = NULL;
FILE* term_output_file    = NULL;
FILE* zones_output_file   = NULL;
FILE* operons_output_file = NULL;

int32_t begin_gener = 0;
int32_t end_gener   = 0;
int32_t final_index = 0;




int main(int argc, char** argv) 
{
  // The input file (lineage.ae or lineage.rae) must contain the following informations:
  //
  // - common data                                                (ae_common::write_to_backup)
  // - begin gener                                                (int32_t)
  // - end gener                                                  (int32_t)
  // - final individual index                                     (int32_t)
  // - initial genome size                                        (int32_t)
  // - initial ancestor (nb genetic units + sequences)            (ae_individual::write_to_backup)
  // - replication report of ancestor at generation begin_gener+1 (ae_replic_report::write_to_backup)    
  // - replication report of ancestor at generation begin_gener+2 (ae_replic_report::write_to_backup)    
  // - replication report of ancestor at generation begin_gener+3 (ae_replic_report::write_to_backup)    
  // - ...
  // - replication report of ancestor at generation end_gener     (ae_replic_report::write_to_backup)    


  
  // ***** The parameters are now loaded from the lineage file *****
   

 
  // =====================
  //  Parse command line
  // =====================

  // Default values
  char*   lineage_file_name   = NULL;
  bool    verbose             = false;
  bool    check               = true;

  const char * short_options = "hvnf:"; 
  static struct option long_options[] =
  {
    {"help",    no_argument,       NULL, 'h'},
    {"verbose", no_argument,       NULL, 'v'},
    {"nocheck", no_argument,       NULL, 'n'},
    {"file",    required_argument, NULL, 'f'},
    {0, 0, 0, 0}
  };

  int option;
  while ( (option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1 ) 
  {
    switch( option ) 
    {
      case 'h' : print_help(); exit(EXIT_SUCCESS);  break;
      case 'v' : verbose = true;                    break;
      case 'n' : check = false;                     break;
      case 'f' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          fprintf( stderr, "ERROR : Option -f or --file : missing argument.\n" );
          exit( EXIT_FAILURE );
        }
        
        lineage_file_name = new char[strlen(optarg) + 1];
        sprintf( lineage_file_name, "%s", optarg );
        break;      
      }
      default :
      {
        fprintf( stderr, "ERROR : Unknown option, check your syntax.\n" );
        print_help();
        exit( EXIT_FAILURE );
      }
    }
  }
  
  if ( lineage_file_name == NULL )
  {
    fprintf( stderr, "ERROR : Option -f or --file missing. \n" );
    exit( EXIT_FAILURE );
  }


  // =======================
  //  Open the lineage file
  // =======================

  gzFile * lineage_file = (gzFile*) gzopen( lineage_file_name, "r" );
  if ( lineage_file == Z_NULL )
  {
    fprintf( stderr, "ERROR : Could not read the lineage file %s\n", lineage_file_name );
    exit( EXIT_FAILURE );
  }
  
  ae_common::read_from_backup( lineage_file );

  gzread( lineage_file, &begin_gener, sizeof(begin_gener) );
  gzread( lineage_file, &end_gener,   sizeof(end_gener)   );
  gzread( lineage_file, &final_index, sizeof(final_index) );

  if ( verbose )
  {
    printf("\n\n");
    printf("================================================================================\n");
    printf(" Statistics of the ancestors of indiv. #%"PRId32" (t=%"PRId32" to %"PRId32")\n", final_index, begin_gener, end_gener);
    printf("================================================================================\n");
  }


  // =========================
  //  Open the output file(s)
  // =========================
  char output_file_name[60];
  snprintf( output_file_name, 60, "ancstats-b%06"PRId32"-e%06"PRId32"-i%"PRId32, begin_gener, end_gener, final_index );
  ae_stats * mystats = new ae_stats( output_file_name );
  mystats->write_headers();
  
  // Optional outputs
  //~ open_environment_stat_file();
  //~ open_terminators_stat_file();
  //~ open_zones_stat_file();
  //~ open_operons_stat_file();
  


  // =========================
  //  Prepare the environment
  // =========================

  if ( verbose )
  {
    printf( "Preparing the environment... " );  
    fflush( stdout );
  }


  // Initialize the environment according the ae_common::... values
  // This instruction also creates the random generator and sets its
  // seed, according to the ae_common::env_seed value.
  ae_environment * env = new ae_environment();

  int32_t t = 0;
  if ( ae_common::env_var_method != NONE )
  {
    for ( t = 0 ; t < begin_gener ; t++ )
    {
      env->apply_variation();
    }
  }

  if ( verbose ) printf("OK\n");

  char backup_file_name[50];
  if ( check )
  {
    // check that the environment is now identical to the one stored
    // in the backup file of generation begin_gener
    
#ifdef __REGUL
    sprintf(backup_file_name,"backup/gen_%06"PRId32".rae", begin_gener);
#else
    sprintf(backup_file_name,"backup/gen_%06"PRId32".ae",  begin_gener);
#endif
    if ( verbose )
    {
      printf("Comparing the environment with the one in %s... ", backup_file_name);  
      fflush(NULL);
    }
    ae_simulation * sim = new ae_simulation ( backup_file_name, false );
    
    ae_environment * stored_env = sim->get_env();

    if ( ! (env->is_identical_to(stored_env)) )
    {
      fprintf(stderr, "ERROR: The replayed environment is not the same\n");
      fprintf(stderr, "       as the one in %s\n", backup_file_name);
      exit(EXIT_FAILURE);
    }
    
    if ( verbose ) printf("OK\n");
    delete sim; 
  }

 


  // ==================================================
  //  Prepare the initial ancestor and write its stats
  // ==================================================

  ae_individual * indiv = new ae_individual( lineage_file );
  indiv->evaluate( env );
  indiv->compute_statistical_data();
  indiv->compute_non_coding();
  
  mystats->write_statistics_of_this_indiv(indiv, t);
  
  
  // Optional outputs
  //~ write_environment_stats( 0, env );
  //~ write_terminators_stats( 0, indiv );
  //~ write_zones_stats( 0, indiv, env );
  //~ write_operons_stats( 0, indiv );
  
  
  if ( verbose )
  {
    printf("Initial fitness     = %f\n", indiv->get_fitness());
    printf("Initial genome size = %"PRId32"\n", indiv->get_total_genome_size());
  }



  // ===============================================================================
  //  Replay the mutations to get the successive ancestors and analyze them
  //  (and, optionally, check that the rebuilt envir is correct each time a backup
  //  is available)
  // ===============================================================================


  ae_replication_report * rep         = NULL;
  ae_list_node *          dnarepnode  = NULL;
  ae_dna_replic_report *  dnarep      = NULL;

  ae_list_node * mnode  = NULL;
  ae_mutation *  mut    = NULL;

  ae_list_node*     unitnode  = NULL;
  ae_genetic_unit * unit      = NULL;

  ae_individual *     stored_indiv    = NULL;
  ae_list_node*       storedunitnode  = NULL;
  ae_genetic_unit *   storedunit      = NULL;

  int32_t index;
  int32_t nb_gener = end_gener - begin_gener;
  
  for ( int32_t i = 0 ; i < nb_gener ; i++ )
  {
    t = begin_gener + i + 1;  // where we are in time...
    rep = new ae_replication_report( lineage_file );
    index = rep->get_index(); // who we are building...
    indiv->set_replication_report( rep );

    if ( verbose ) printf("Rebuilding ancestor at generation %"PRId32" (index %"PRId32")...", t, index); 

    // 1) Rebuild environment
    env->apply_variation();
    
    if ( check && utils::mod(t, ae_common::backup_step) == 0 )
    {
      // check that the environment is now identical to the one stored
      // in the backup file of generation begin_gener
      
      #ifdef __REGUL
        sprintf(backup_file_name,"backup/gen_%06"PRId32".rae", t);
      #else
        sprintf(backup_file_name,"backup/gen_%06"PRId32".ae",  t);
      #endif
      
      if ( verbose )
      {
        printf("Comparing the environment with the one in %s... ", backup_file_name);  
        fflush(NULL);
      }
      
      ae_simulation * sim = new ae_simulation ( backup_file_name, false );
      stored_indiv = new ae_individual( * (ae_individual *)sim->get_pop()->get_indivs()->get_object(index) );
      

      ae_environment * stored_env = sim->get_env();


      if ( ! (env->is_identical_to(stored_env)) )
      {
        fprintf(stderr, "ERROR: The replayed environment is not the same\n");
        fprintf(stderr, "       as the one in %s\n", backup_file_name);
        exit(EXIT_FAILURE);
      }
      
      if ( verbose ) printf("OK\n");
      delete sim;
    }

    // Warning: this portion of code won't work if the number of units changes
    // during the evolution
    
    // 2) Replay replication (create current individual's child)
    dnarepnode  = rep->get_dna_replic_reports()->get_first();
    unitnode    = indiv->get_genetic_unit_list()->get_first();
    
    if ( check && utils::mod(t, ae_common::backup_step) == 0 )
    {
      storedunitnode = stored_indiv->get_genetic_unit_list()->get_first();
    }
    
    // For each genetic unit, replay the replication (undergo all mutations)
    while ( dnarepnode != NULL )
    {
      assert( unitnode != NULL );

      dnarep  = (ae_dna_replic_report *)  dnarepnode->get_obj();
      unit    = (ae_genetic_unit *)       unitnode->get_obj();
      
      unit->get_dna()->set_replic_report( dnarep );

      mnode = dnarep->get_rearrangements()->get_first();              
      while ( mnode != NULL )
      {
        mut = (ae_mutation *) mnode->get_obj();
        unit->get_dna()->undergo_this_mutation(mut); 
        mnode = mnode->get_next();
      }
      
      mnode = dnarep->get_mutations()->get_first();              
      while ( mnode != NULL )
      {
        mut = (ae_mutation *) mnode->get_obj();
        unit->get_dna()->undergo_this_mutation( mut ); 
        mnode = mnode->get_next();
      }

      if ( check && utils::mod(t, ae_common::backup_step) == 0 )
      {
        if ( verbose )
        {
          printf("Checking the sequence of the unit...");
          fflush(NULL);
        }
        
        assert( storedunitnode != NULL );
        storedunit = (ae_genetic_unit *) storedunitnode->get_obj();
        
        char * str1 = new char[unit->get_dna()->get_length() + 1];
        memcpy(str1, unit->get_dna()->get_data(), \
               unit->get_dna()->get_length()*sizeof(char));
        str1[unit->get_dna()->get_length()] = '\0';
        
        char * str2 = new char[(storedunit->get_dna())->get_length() + 1];
        memcpy(str2, (storedunit->get_dna())->get_data(), (storedunit->get_dna())->get_length()*sizeof(char));
        str2[(storedunit->get_dna())->get_length()] = '\0';
        
        if ( strncmp( str1, str2, storedunit->get_dna()->get_length() ) == 0 )
        {
          if ( verbose ) printf(" OK\n");
        }
        else
        {
          if ( verbose ) printf( " ERROR !\n" );
          fprintf( stderr, "Error: the rebuilt unit is not the same as \n");
          fprintf( stderr, "the one stored in backup file %s\n", backup_file_name);
          fprintf( stderr, "Rebuilt unit : %"PRId32" bp\n %s\n", (int32_t)strlen(str1), str1 );
          fprintf( stderr, "Stored unit  : %"PRId32" bp\n %s\n", (int32_t)strlen(str2), str2 );
          delete [] str1;
          delete [] str2;
          gzclose(lineage_file);
          delete indiv;
          delete stored_indiv;
          delete env;
          ae_common::clean();
          exit(EXIT_FAILURE);
        }
        
        delete [] str1;
        delete [] str2;
        
        storedunitnode = storedunitnode->get_next();
      }
      
      
      dnarepnode = dnarepnode->get_next();
      unitnode = unitnode->get_next();
    }
        
    assert( unitnode == NULL );

    // 3) All the mutations have been replayed, we can now evaluate the new individual
    indiv->reevaluate( env );
    indiv->compute_statistical_data();
    indiv->compute_non_coding();

    mystats->write_statistics_of_this_indiv( indiv, t );


    // Optional outputs
    //~ write_environment_stats( t, env );
    //~ write_terminators_stats( t, indiv );
    //~ write_zones_stats( t, indiv, env );
    //~ write_operons_stats( t, indiv );
    

    if ( verbose ) printf(" OK\n");

    delete rep;
    if ( utils::mod(t, ae_common::backup_step) == 0 )
    {
      assert(storedunitnode == NULL);
      delete stored_indiv;
    }
  }

  
  ae_common::clean();

  gzclose(lineage_file);
  delete mystats;
  delete indiv;
  delete env;

  fclose( env_output_file );
  fclose( term_output_file );
  fclose( zones_output_file );
  fclose( operons_output_file );


  exit(EXIT_SUCCESS);
  
}




void open_environment_stat_file( void )
{
  // Open file
  char env_output_file_name[60];
  snprintf( env_output_file_name, 60, "envir-b%06"PRId32"-e%06"PRId32"-i%"PRId32".out", begin_gener, end_gener, final_index );
  env_output_file = fopen( env_output_file_name, "w" );
  
  // Write headers
  fprintf( env_output_file, "# Each line contains : Generation, and then, for each gaussian: M W H.\n" );
  fprintf( env_output_file, "#\n" );
}

void write_environment_stats( int32_t num_gener, ae_environment * env )
{
  // Num gener
  fprintf( env_output_file, "%"PRId32, num_gener );

  // For each gaussian : M W H
  ae_list_node * gaussnode  = env->get_gaussians()->get_first();
  ae_gaussian *  gauss      = NULL;
  while ( gaussnode != NULL )
  {
    gauss = (ae_gaussian *) gaussnode->get_obj();
    fprintf( env_output_file, "     %.16f %.16f %.16f", gauss->get_mean(), gauss->get_width(), gauss->get_height() );
    gaussnode = gaussnode->get_next();
  }
  
  fprintf( env_output_file, "\n" );
}



void open_terminators_stat_file( void )
{
  char term_output_file_name[60];
  snprintf( term_output_file_name, 60, "nb_term-b%06"PRId32"-e%06"PRId32"-i%"PRId32, begin_gener, end_gener, final_index );
  term_output_file = fopen( term_output_file_name, "w" );
}

void write_terminators_stats( int32_t num_gener, ae_individual * indiv )
{
  fprintf(  term_output_file, "%"PRId32" %"PRId32" %"PRId32"\n",
            num_gener,
            indiv->get_total_genome_size(),
            indiv->get_nb_terminators() );
}



void open_zones_stat_file( void )
{
  // Open file
  char zones_output_file_name[60];
  snprintf( zones_output_file_name, 60, "zones-b%06"PRId32"-e%06"PRId32"-i%"PRId32, begin_gener, end_gener, final_index );
  zones_output_file = fopen( zones_output_file_name, "w" );
  
  // Write headers
  fprintf( zones_output_file, "# Each line contains : Generation, and then, for each zone:\n" );
  fprintf( zones_output_file, "#   * Number of activation genes\n" );
  fprintf( zones_output_file, "#   * Number of inhibition genes\n" );
  fprintf( zones_output_file, "#   * Geometric area of the activation genes\n" );
  fprintf( zones_output_file, "#   * Geometric area of the inhibition genes\n" );
  fprintf( zones_output_file, "#   * Geometric area of the resulting phenotype\n" );
  fprintf( zones_output_file, "#\n" );
}

void write_zones_stats( int32_t num_gener, ae_individual * indiv, ae_environment * env )
{
  int16_t           nb_segments = env->get_nb_segments();
  int16_t           num_segment = 0;
  ae_env_segment ** segments    = env->get_segments();
  
  ae_list*          prot_list = indiv->get_protein_list();
  ae_list_node*     prot_node = NULL;
  ae_protein*       prot      = NULL;
    
  // Tables : index 0 for the 0 segment
  //                1 for the neutral segment
  int32_t nb_genes_activ[nb_segments];
  int32_t nb_genes_inhib[nb_segments];
  double  geom_area_activ[nb_segments];
  double  geom_area_inhib[nb_segments];
  double  geom_area_phen[nb_segments];
  
  for ( num_segment = 0 ; num_segment < nb_segments ; num_segment++ )
  {
    nb_genes_activ[num_segment]   = 0;
    nb_genes_inhib[num_segment]   = 0;
    geom_area_activ[num_segment]  = 0.0;
    geom_area_inhib[num_segment]  = 0.0;
    geom_area_phen[num_segment]   = 0.0;
  }
  
  
  ae_fuzzy_set* activ = NULL;
  ae_fuzzy_set* inhib = NULL;
  ae_phenotype* phen  = NULL;
  
  
  
  // Compute number of genes in each segment
  prot_node = prot_list->get_first();
  prot      = NULL;
  
  while ( prot_node != NULL )
  {
    prot = (ae_protein*) prot_node->get_obj();
    
    // Go to the corresponding segment
    num_segment = 0;
    while ( prot->get_mean() > segments[num_segment]->stop )
    {
      num_segment++;
    }
    
    // Add a genes (activ or inhib)
    if ( ! prot->is_degenerated() )
    {
      if ( prot->get_height() > 0 )
      {
        nb_genes_activ[num_segment]++;
      }
      else if ( prot->get_height() < 0 )
      {
        nb_genes_inhib[num_segment]++;
      }
      
      // It the gene is exactly at the frontier between 2 zones, mark it in both
      if ( prot->get_mean() == segments[num_segment]->stop && num_segment < nb_segments - 1 )
      {
        if ( prot->get_height() > 0 )
        {
          nb_genes_activ[num_segment+1]++;
        }
        else if ( prot->get_height() < 0 )
        {
          nb_genes_inhib[num_segment+1]++;
        }
      }
    }
    
    prot_node = prot_node->get_next();
  }
  
  // Compute the geometric areas
  activ = indiv->get_phenotype_activ();
  inhib = indiv->get_phenotype_inhib();
  phen  = indiv->get_phenotype();
  
  for ( num_segment = 0 ; num_segment < nb_segments ; num_segment++ )
  {
    geom_area_activ[num_segment]  = activ->get_geometric_area( segments[num_segment]->start, segments[num_segment]->stop );
    geom_area_inhib[num_segment]  = inhib->get_geometric_area( segments[num_segment]->start, segments[num_segment]->stop );
    geom_area_phen[num_segment]   = phen->get_geometric_area(  segments[num_segment]->start, segments[num_segment]->stop );
  }
  
  
  // Print stats to file
  fprintf(  zones_output_file, "%"PRId32, num_gener );
  
  for ( num_segment = 0 ; num_segment < nb_segments ; num_segment++ )
  {
    fprintf(  zones_output_file, "     %"PRId32" %"PRId32" %lf %lf %lf",
              nb_genes_activ[num_segment],
              nb_genes_inhib[num_segment],
              geom_area_activ[num_segment],
              geom_area_inhib[num_segment],
              geom_area_phen[num_segment] );
  }
  
  fprintf(  zones_output_file, "\n" );
}



void open_operons_stat_file( void )
{
  char operons_output_file_name[60];
  snprintf( operons_output_file_name, 60, "operons-b%06"PRId32"-e%06"PRId32"-i%"PRId32, begin_gener, end_gener, final_index );
  operons_output_file = fopen( operons_output_file_name, "w" );
}

void write_operons_stats( int32_t num_gener, ae_individual * indiv )
{
  int32_t nb_genes_per_rna[20];
  for ( int i = 0 ; i < 20 ; i++ )
  {
    nb_genes_per_rna[i] = 0;
  }
  
  ae_list_node* rna_node = indiv->get_rna_list()->get_first();
  ae_rna*       rna      = NULL;
  
  while ( rna_node != NULL )
  {
    rna = (ae_rna*) rna_node->get_obj();
    
    if ( rna->get_transcribed_proteins()->get_nb_elts() >= 20 )
    {
      printf( "Found operon with 20 genes or more : %"PRId16"\n", rna->get_transcribed_proteins()->get_nb_elts() );
      getchar();
    }
    
    nb_genes_per_rna[rna->get_transcribed_proteins()->get_nb_elts()]++;
    
    rna_node = rna_node->get_next();
  }
  
  fprintf(  operons_output_file, "%"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32"\n",
            num_gener,
            nb_genes_per_rna[0],
            nb_genes_per_rna[1],
            nb_genes_per_rna[2],
            nb_genes_per_rna[3],
            nb_genes_per_rna[4],
            nb_genes_per_rna[5],
            nb_genes_per_rna[6],
            nb_genes_per_rna[7],
            nb_genes_per_rna[8],
            nb_genes_per_rna[9],
            nb_genes_per_rna[10],
            nb_genes_per_rna[11],
            nb_genes_per_rna[12],
            nb_genes_per_rna[13],
            nb_genes_per_rna[14],
            nb_genes_per_rna[15],
            nb_genes_per_rna[16],
            nb_genes_per_rna[17],
            nb_genes_per_rna[18],
            nb_genes_per_rna[19] );
}



void print_help( void )
{
  // TODO
}
