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
#include <inttypes.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <err.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>




// =================================================================
//                            Project Files
// =================================================================
#include <ae_utils.h>
#ifndef __NO_X
  #include <ae_exp_manager_X11.h>
#else
  #include <ae_exp_manager.h>
#endif
#include <ae_individual.h>
#include <ae_genetic_unit.h>
#include <ae_list.h>
#include <ae_tree.h>
#include <ae_replication_report.h>
#include <ae_dna_replic_report.h>
#include <ae_mutation.h>
//#include <ae_param_loader.h>
#include <ae_environment.h>
#include <ae_enums.h>
//#include <ae_common.h>

//debug
#include <ae_gaussian.h>




enum check_type
{
  FULL_CHECK  = 0,
  LIGHT_CHECK = 1,
  ENV_CHECK   = 2,
  NO_CHECK    = 3
};



FILE* open_environment_stat_file( const char * prefix );
void write_environment_stats( int32_t num_gener, ae_environment * env, FILE* env_file);

FILE* open_terminators_stat_file( const char * prefix );
void write_terminators_stats( int32_t num_gener, ae_individual * indiv, FILE* terminator_file );

FILE* open_zones_stat_file( const char * prefix );
void write_zones_stats( int32_t num_gener, ae_individual * indiv, ae_environment * env, FILE* zone_file );

FILE* open_operons_stat_file( const char * prefix );
void write_operons_stats( int32_t num_gener, ae_individual * indiv, FILE* operon_file );

void print_help( void );


double* dist_to_target_segment;


int main(int argc, char** argv) 
{
  // The input file (lineage.ae or lineage.rae) must contain the following information:
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
  char*       lineage_file_name   = NULL;
  bool        verbose             = false;
  check_type  check               = LIGHT_CHECK;   
  double      tolerance           = 0;
  
  const char * short_options = "hvncf:lt:"; 
  static struct option long_options[] =
  {
    {"help",        no_argument,       NULL, 'h'},
    {"verbose",     no_argument,       NULL, 'v'},
    {"nocheck",     no_argument,       NULL, 'n'},
    {"fullcheck",   no_argument,       NULL, 'c'},
    {"file",        required_argument, NULL, 'f'},
    {"tolerance",   required_argument, NULL, 't'},
    {0, 0, 0, 0}
  };

  int option;
  while ( (option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1 ) 
  {
    switch( option ) 
    {
      case 'h' : print_help(); exit(EXIT_SUCCESS);  break;
      case 'v' : verbose = true;                    break;
      case 'n' : check = NO_CHECK;                  break;
      case 'c' : check = FULL_CHECK;                break;
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
      case 't' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          fprintf( stderr, "ERROR : Option -t or --tolerance : missing argument.\n" );
          exit( EXIT_FAILURE );
        }
        check = ENV_CHECK;
        tolerance = atof(optarg);
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
  
  
  printf("\n");
  printf( "WARNING : Parameter change during simulation is not managed (consider -l option)\n" );
  printf("\n");

  // =======================
  //  Open the lineage file
  // =======================
  gzFile lineage_file = gzopen( lineage_file_name, "r" );
  if ( lineage_file == Z_NULL )
  {
    fprintf( stderr, "ERROR : Could not read the lineage file %s\n", lineage_file_name );
    exit( EXIT_FAILURE );
  }

	int32_t begin_gener       = 0;
	int32_t end_gener         = 0;
	int32_t final_indiv_index = 0;
	int32_t final_indiv_rank  = 0;
  

  gzread( lineage_file, &begin_gener,       sizeof(begin_gener)       );
  gzread( lineage_file, &end_gener,         sizeof(end_gener)         );
  gzread( lineage_file, &final_indiv_index, sizeof(final_indiv_index) );
  gzread( lineage_file, &final_indiv_rank,  sizeof(final_indiv_rank)  );

  if ( verbose )
  {
    printf("\n\n");
    printf( "===============================================================================\n" );
    printf( " Statistics of the ancestors of indiv. %"PRId32" (rank %"PRId32") from generation %"PRId32" to %"PRId32"\n",
            final_indiv_index, final_indiv_rank, begin_gener, end_gener );
    printf("================================================================================\n");
  }



  // =========================
  //  Open the experience manager
  // =========================
  
  #ifndef __NO_X
    ae_exp_manager* exp_manager = new ae_exp_manager_X11();
  #else
    ae_exp_manager* exp_manager = new ae_exp_manager();
  #endif

  exp_manager->load( begin_gener, false, true, false);
  ae_environment* env = exp_manager->get_env();
  
  int32_t backup_step = exp_manager->get_backup_step();
  
  // =========================
  //  Open the output file(s)
  // =========================
  // Create missing directories
  int status;
  status = mkdir( "stats/ancstats/", 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, "stats/ancstats/" );
  }
  
  char prefix[50];
  snprintf( prefix, 50, "ancstats/ancstats-b%06"PRId32"-e%06"PRId32"-i%"PRId32"-r%"PRId32,begin_gener, end_gener, final_indiv_index , final_indiv_rank);
  bool best_indiv_only = true;
  bool addition_old_stats = false;
  bool delete_old_stats = true;
  ae_stats * mystats = new ae_stats(exp_manager, begin_gener, prefix, best_indiv_only, addition_old_stats, delete_old_stats);
  //mystats->write_headers();
  
  // Optional outputs
  FILE* env_output_file = open_environment_stat_file(prefix);
  FILE* term_output_file = open_terminators_stat_file(prefix);
  FILE* zones_output_file = NULL;
  if(env->get_nb_segments() > 1)
  {
  	zones_output_file = open_zones_stat_file(prefix);
  }
  FILE* operons_output_file = open_operons_stat_file(prefix);


  // ==================================================
  //  Prepare the initial ancestor and write its stats
  // ==================================================
  ae_individual * indiv = new ae_individual(exp_manager, lineage_file );
  indiv->evaluate( env );
  indiv->compute_statistical_data();
  indiv->compute_non_coding();
  
  mystats->write_statistics_of_this_indiv( indiv, begin_gener );
  
  
  // Optional outputs
  write_environment_stats( begin_gener, env, env_output_file );
  write_terminators_stats( begin_gener, indiv, term_output_file );
  if(env->get_nb_segments() > 1)
  {
  	write_zones_stats( begin_gener, indiv, env, zones_output_file );
  }
  write_operons_stats( begin_gener, indiv, operons_output_file );
  
  
  if ( verbose )
  {
    printf("Initial fitness     = %f\n", indiv->get_fitness());
    printf("Initial genome size = %"PRId32"\n", indiv->get_total_genome_size());
  }

  //delete exp_manager;

  // ===============================================================================
  //  Replay the mutations to get the successive ancestors and analyze them
  //  (and, optionally, check that the rebuilt envir is correct each time a backup
  //  is available)
  // ===============================================================================

  int32_t num_gener = 0;
  
  ae_replication_report* rep = NULL;
  ae_list_node<ae_dna_replic_report*>* dnarepnode  = NULL;
  ae_dna_replic_report* dnarep = NULL;

  ae_list_node<ae_mutation*>* mnode  = NULL;
  ae_mutation* mut = NULL;

  ae_list_node<ae_genetic_unit*>* unitnode  = NULL;
  ae_genetic_unit* unit = NULL;

  ae_individual* stored_indiv = NULL;
  ae_list_node<ae_genetic_unit*>* storedunitnode  = NULL;
  ae_genetic_unit* storedunit = NULL;

  int32_t index;
  int32_t nb_gener = end_gener - begin_gener;
  
  ae_exp_manager* exp_manager_backup = NULL;
  ae_environment* backup_env = NULL;
  
  bool check_now = false;
  
  for ( int32_t i = 0 ; i < nb_gener ; i++ )
  {
    num_gener = begin_gener + i + 1;  // where we are in time..
    
    env->build();
    rep = new ae_replication_report( lineage_file, indiv );
    index = rep->get_id(); // who we are building...
    indiv->set_replication_report( rep );
    
    // Check now?
    check_now = ( ( check == FULL_CHECK && ae_utils::mod( num_gener, backup_step ) == 0 ) || 
                  ( check == ENV_CHECK && ae_utils::mod( num_gener, backup_step ) == 0 ) ||
                  ( check == LIGHT_CHECK && num_gener == end_gener ) );

    if ( verbose ) printf("Rebuilding ancestor at generation %"PRId32" (index %"PRId32")...", num_gener, index); 

    // 1) Rebuild environment
    env->apply_variation();
    
    if ( check_now )
    {
      // check that the environment is now identical to the one stored
      // in the backup file of generation begin_gener
      
      // Load the simulation
      #ifndef __NO_X
      	exp_manager_backup = new ae_exp_manager_X11();
      #else
      	exp_manager_backup = new ae_exp_manager();
      #endif

      exp_manager_backup->load( num_gener, false, true, false );
      backup_env = exp_manager_backup->get_env();
      stored_indiv = new ae_individual( * (ae_individual *)exp_manager_backup->get_indiv_by_id( index ) );
      //delete exp_manager;
  
      if ( verbose )
      {
        printf( "Comparing the environment with the one saved at generation %"PRId32"... ", num_gener );  
        fflush(NULL);
      }

      if ( ! env->is_identical_to(backup_env, tolerance) )
      {
        fprintf(stderr, "ERROR: The replayed environment is not the same\n");
        fprintf(stderr, "       as the one saved at generation %"PRId32"... \n", num_gener );  
        fprintf(stderr, "       with tolerance of %lg\n", tolerance);
        exit( EXIT_FAILURE );
      }
      
      if ( verbose ) printf("OK\n");
    }

    // Warning: this portion of code won't work if the number of units changes
    // during the evolution
    
    // 2) Replay replication (create current individual's child)
    dnarepnode  = rep->get_dna_replic_reports()->get_first();
    unitnode    = indiv->get_genetic_unit_list()->get_first();
    
    if ( check_now )
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
      
      mnode = dnarep->get_HT()->get_first();              
      while ( mnode != NULL )
      {
        mut = (ae_mutation *) mnode->get_obj();
        unit->get_dna()->undergo_this_mutation(mut); 
        mnode = mnode->get_next();
      }

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

      if ( check_now )
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
          fprintf( stderr, "the one saved at generation %"PRId32"... ", num_gener );
          fprintf( stderr, "Rebuilt unit : %"PRId32" bp\n %s\n", (int32_t)strlen(str1), str1 );
          fprintf( stderr, "Stored unit  : %"PRId32" bp\n %s\n", (int32_t)strlen(str2), str2 );
          
          delete [] str1;
          delete [] str2;
          gzclose(lineage_file);
          delete indiv;
          delete stored_indiv;
          delete exp_manager_backup;
          delete exp_manager; 
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

    mystats->write_statistics_of_this_indiv( indiv, num_gener );

    // Optional outputs
    write_environment_stats( num_gener, env, env_output_file );
  	write_terminators_stats( num_gener, indiv, term_output_file );
  	if(env->get_nb_segments() > 1)
  	{
  		write_zones_stats( num_gener, indiv, env, zones_output_file );
  	}
  	write_operons_stats( num_gener, indiv, operons_output_file );  

    if ( verbose ) printf(" OK\n");

    delete rep;
    if ( check_now )
    {
      assert( storedunitnode == NULL );
      delete stored_indiv;
      delete exp_manager_backup;
    }
  }

  gzclose(lineage_file);
  delete exp_manager;
  delete mystats;
  delete indiv;
  //delete env; // already done in ae_exp_setup destructor called by exp_manager destructor

  // Optional outputs
  fclose( env_output_file );
  fclose( term_output_file );
  if(env->get_nb_segments() > 1)
  {
  	fclose( zones_output_file );
  }
  fclose( operons_output_file );

  exit(EXIT_SUCCESS);
}




FILE* open_environment_stat_file( const char * prefix  )
{
  // Open file
  char* env_output_file_name = new char[80];
  sprintf( env_output_file_name, "stats/%s_envir.out",prefix );
  FILE* env_output_file = fopen( env_output_file_name, "w" );
  delete env_output_file_name;
  
  // Write headers
  fprintf( env_output_file, "# Each line contains : Generation, and then, for each gaussian: M W H.\n" );
  fprintf( env_output_file, "#\n" );
  
  return env_output_file;
}

void write_environment_stats( int32_t num_gener, ae_environment * env, FILE*  env_output_file)
{
  // Num gener
  fprintf( env_output_file, "%"PRId32, num_gener );

  // For each gaussian : M W H
  ae_list_node<ae_gaussian*>* gaussnode  = env->get_gaussians()->get_first();
  ae_gaussian*  gauss      = NULL;
  while ( gaussnode != NULL )
  {
    gauss = gaussnode->get_obj();
    fprintf( env_output_file, "     %.16f %.16f %.16f", gauss->get_mean(), gauss->get_width(), gauss->get_height() );
    gaussnode = gaussnode->get_next();
  }
  
  fprintf( env_output_file, "\n" );
}



FILE* open_terminators_stat_file( const char * prefix )
{
  char* term_output_file_name = new char[80];
  sprintf( term_output_file_name, "stats/%s_nb_term.out",prefix );
  FILE* term_output_file = fopen( term_output_file_name, "w" );
  delete [] term_output_file_name;
  
  // Write headers
  fprintf( term_output_file, "# Each line contains : \n" );
  fprintf( term_output_file, "#   * Generation\n" );
  fprintf( term_output_file, "#   * Genome size\n" );
  fprintf( term_output_file, "#   * Terminator number\n");
  fprintf( term_output_file, "#\n" );
  
  return term_output_file;
}

void write_terminators_stats( int32_t num_gener, ae_individual * indiv, FILE* term_output_file  )
{
  fprintf(  term_output_file, "%"PRId32" %"PRId32" %"PRId32"\n",
            num_gener,
            indiv->get_total_genome_size(),
            indiv->get_nb_terminators() );
}



FILE* open_zones_stat_file( const char * prefix  )
{
  // Open file
  char* zones_output_file_name = new char[80];
  sprintf( zones_output_file_name, "stats/%s_zones.out",prefix );
  FILE* zones_output_file = fopen( zones_output_file_name, "w" );
  delete [] zones_output_file_name;
  
  // Write headers
  fprintf( zones_output_file, "# Each line contains : Generation, and then, for each zone:\n" );
  fprintf( zones_output_file, "#   * Number of activation genes\n" );
  fprintf( zones_output_file, "#   * Number of inhibition genes\n" );
  fprintf( zones_output_file, "#   * Geometric area of the activation genes\n" );
  fprintf( zones_output_file, "#   * Geometric area of the inhibition genes\n" );
  fprintf( zones_output_file, "#   * Geometric area of the resulting phenotype\n" );
  fprintf( zones_output_file, "#\n" );
  
  return zones_output_file;
}

void write_zones_stats( int32_t num_gener, ae_individual * indiv, ae_environment * env, FILE* zones_output_file )
{
  assert( env->get_nb_segments() > 1 );
  
  int16_t nb_segments = env->get_nb_segments();
  int16_t num_segment = 0;
  ae_env_segment** segments = env->get_segments();
  
  ae_list<ae_protein*>* prot_list = indiv->get_protein_list();
  ae_list_node<ae_protein*>* prot_node = NULL;
  ae_protein* prot = NULL;
    
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
    prot = prot_node->get_obj();
    
    // Go to the corresponding segment
    num_segment = 0;
    while ( prot->get_mean() > segments[num_segment]->stop )
    {
      num_segment++;
    }
    
    // Add a genes (activ or inhib)
    if ( prot->get_is_functional() )
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



FILE* open_operons_stat_file( const char * prefix  )
{
  char* operons_output_file_name = new char[80];
  sprintf( operons_output_file_name, "stats/%s_operons.out",prefix);
  FILE* operons_output_file = fopen( operons_output_file_name, "w" );
  delete [] operons_output_file_name,
  
  // Write headers
  fprintf( operons_output_file, "# Each line contains : Generation, and then, for 20 RNA, the number of genes inside the RNA\n" );
  return operons_output_file;
}

void write_operons_stats( int32_t num_gener, ae_individual * indiv, FILE*  operons_output_file)
{
  int32_t nb_genes_per_rna[20];
  for ( int i = 0 ; i < 20 ; i++ )
  {
    nb_genes_per_rna[i] = 0;
  }
  
  ae_list_node<ae_rna*>* rna_node = indiv->get_rna_list()->get_first();
  ae_rna* rna = NULL;
  
  while ( rna_node != NULL )
  {
    rna = rna_node->get_obj();
    
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
  printf( "\n" ); 
  printf( "*********************** aevol - Artificial Evolution ******************* \n" );
  printf( "*                                                                      * \n" );
  printf( "*                      Ancstats post-treatment program                 * \n" );
  printf( "*                                                                      * \n" );
  printf( "************************************************************************ \n" );
  printf( "\n\n" ); 
  printf( "This program is Free Software. No Warranty.\n" );
  printf( "Copyright (C) 2009  LIRIS.\n" );
  printf( "\n" ); 
#ifdef __REGUL
  printf( "Usage : rancstats -h\n");
  printf( "or :    rancstats [-vn] -f lineage_file \n" );
#else
  printf( "Usage : ancstats -h\n");
  printf( "or :    ancstats [-vn] -f lineage_file \n" );
#endif
  printf( "\n" ); 
  printf( "This program compute some statistics for the individuals within lineage_file.\n" );
  printf( "\n" ); 
  printf( "WARNING: This program should not be used for simulations run with lateral\n" ); 
  printf( "transfer. When an individual has more than one parent, the notion of lineage\n" ); 
  printf( "used here is not relevant.\n" );
  printf( "\n" );  
  printf( "\t-h or --help       : Display this help.\n" );
  printf( "\n" ); 
  printf( "\t-v or --verbose    : Be verbose, listing generations as they are \n" );
  printf( "\t                       treated.\n" );
  printf( "\n" );
  printf( "\t-n or --nocheck    : Disable genome sequence checking. Makes the \n"); 
  printf( "\t                       program faster, but it is not recommended. \n");
  printf( "\t                       It is better to let the program check that \n");
  printf( "\t                       when we rebuild the genomes of the ancestors\n");
  printf( "\t                       from the lineage file, we get the same sequences\n");
  printf( "\t                       as those stored in the backup files.\n" );
  printf( "\n" );
  printf( "\t-c or --fullcheck  : Will perform the genome and environment checks every\n" );
  printf( "\t                       <BACKUP_STEP> generations. Default behaviour is\n" );
  printf( "\t                       lighter as it only perform sthese checks at the\n" );
  printf( "\t                       ending generation.\n" );
  printf( "\n" ); 
  printf( "\t-f lineage_file or --file lineage_file : \n" );
  printf( "\t                     	Compute the statistics for the individuals within lineage_file.\n" );
  printf( "\t-t tolerance or --tolerance tolerance : \n");
  printf( "\t                       Tolerance used to compare the replayed environment to environment in backup\n");
  printf( "\n" );
}
