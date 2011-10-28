//
// This program extracts data about the individuals and write
// them into text files easy to parse with matlab.
//
// Options are :
//     * -f / --file   backup_file
//            backup file where individuals should be extracted
//     * -o / --output   main_output_file
//            output directory where statistics should be printed
//     * -n / --neutral  output_file
//            print position of neutral regions in 'ouput_file'
//     * -b / --best
//            process only the best individual
// 
// Data concerning genomes is printed in 'main_output_file'. A space delimits two pieces
// of information, a new line two individuals. Format is 
//       "nc1 nc2 nc3 nc4 nc5 nc6 nc7 total\n"
// where:
//     * nc1: number of bases in neutral regions
//     * nc2: number of bases outside CDS
//     * nc3: number of bases outside functional CDS
//     * nc4: number of bases outside non functional CDS
//     * nc5: number of bases outside RNAs
//     * nc6: number of bases outside functional RNAs
//     * nc7: number of bases outside non functional RNAs
//     * total: total size of genome
//
// It is also possible to print neutral regions with the '-n output_file' option.
// In this case, "output_file" contains:
//       "# chromosome length: cl, nb neutral bases, nr neutral regions\n
//       "bnr_1 bnr_2 ... bnr_n\n"
//       "enr_1 enr_2 ... enr_n\n"
// where:
//     * cl:   number of bases in chromosome
//     * nb:   number of bases in neutral regions
//     * nr:   number of neutral regions
//     * bnr_i: beginning of ith neutral region
//     * enr_i: end of ith neutral region
//
// Examples :
//
// For generation 20000, write info about the genomes of all the
// individuals in out_020000 :
//
//    genome_stats -f backup/gen_020000.ae -o out_020000
//
// For generation 20000, write the best individual's info in
// out_020000_best :
//
//    genome_stats -b -f backup/gen_020000.ae -o out_020000_best
//
// For generation 20000, write beginning and ends of neutral regions
// in neutral.out for all the individuals
//
//    genome_stats -f backup/gen_020000.ae -n neutral.out
//


/** \class
 *  \brief
 */

// =================================================================
//                              Libraries
// =================================================================
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>


// =================================================================
//                            Project Files
// =================================================================
#include <ae_common.h>
#include <ae_population.h>
#include <ae_individual.h>
#include <ae_list.h>
#include <ae_simulation.h>



// =================================================================
//                        Secondary Functions
// =================================================================

void print_help( char* prog_name );

// print information about the indivdual's genome to file
void print_genome_info( ae_individual* indiv, FILE* output_file );

// print information about the indivdual's neutral regions to file
void print_neutral_regions( ae_individual* indiv, FILE* output_file);


// =================================================================
//                         Main Function
// =================================================================

int main( int argc, char* argv[] )
{
  // ---------------------------------------
  //      command-line option parsing
  // ---------------------------------------
  // Initialize command-line option variables with default values
  char* backup_file_name           = NULL;
  char* main_output_name           = NULL;
  char* neutral_region_output_name = NULL;
  bool best_only = false;
  
  // Define allowed options
  const char * options_list = "hf:o:bn:";
  static struct option long_options_list[] =
  {
    { "help", 0, NULL, 'h' },
    { "file", 1, NULL, 'f' },
    { "output", 1, NULL, 'o' },
    { "best", 0, NULL, 'b' },
    { "neutral", 0, NULL, 'n' },
    { 0, 0, 0, 0 }
  };

  // Get actual values of the command-line options
  int option;
  while ( ( option = getopt_long(argc, argv, options_list, long_options_list, NULL) ) != -1 ) 
  {
    switch ( option ) 
    {
    case 'h' :
      print_help( argv[0] );
      exit( EXIT_SUCCESS );
      break;
    case 'f' :
      backup_file_name = new char[strlen(optarg) + 1];
      sprintf( backup_file_name, "%s", optarg );
      break;
    case 'o' :
      main_output_name = new char[strlen(optarg) + 1];
      sprintf( main_output_name, "%s", optarg );
      break;
    case 'n' :
      neutral_region_output_name = new char[strlen(optarg) + 1];
      sprintf( neutral_region_output_name, "%s", optarg );
      break;
    case 'b' :
      best_only = true;
      break;
    }
  }
  
  // -------------------------------
  //          Initialize
  // -------------------------------
  FILE* main_output           = NULL;
  FILE* neutral_region_output = NULL;
  
  if ( backup_file_name == NULL )
  {
    printf("You must specify a backup file. Please use the option -f or --file.\n");
    exit(EXIT_FAILURE);
  }
  if ( main_output_name != NULL )
  {
    main_output = fopen(main_output_name,"w");
    if ( main_output == NULL)
    {
      fprintf( stderr, "Warning: Could not open file %s.\n", main_output_name );
    }
  }
  if ( neutral_region_output_name != NULL )
  {
    neutral_region_output = fopen(neutral_region_output_name,"w");
    if ( neutral_region_output == NULL)
    {
      fprintf( stderr, "Warning: Could not open file %s.\n", neutral_region_output_name );
    }
  }
  fflush( stderr );

  printf( "Reading backup file <%s>... \n", backup_file_name );
  fflush( stdout );

  // Load the simulation from backup
  ae_common::sim = new ae_simulation( backup_file_name, false );
  printf("done\n");
  delete [] backup_file_name;

  printf( "Computing phenotypes... \n" );
  fflush( stdout );

  // Evaluate the individuals
  (ae_common::sim->get_pop())->evaluate_individuals(ae_common::sim->get_env());
  
  int i = 0;
  int nb_indiv = (ae_common::sim->get_pop())->get_nb_indivs();

  // --------------------------------
  //         Parse individuals
  // --------------------------------
  if (best_only)
  {
    ae_individual* best = ae_common::sim->get_pop()->get_best();
    if ( main_output != NULL)           { print_genome_info(best, main_output); }
    if ( neutral_region_output != NULL) { print_neutral_regions(best, neutral_region_output); }
  }
  else
  {
    if (ae_common::pop_structure)
    {
      ae_grid_cell*** _pop_grid = (ae_common::sim->get_pop())->get_pop_grid();
      for ( int16_t x = 0 ; x < ae_common::grid_x ; x++ )
      {
        for ( int16_t y = 0 ; y < ae_common::grid_y ; y++ )
        {
          ae_individual* indiv = (_pop_grid[x][y]->get_individual());
	  if ( main_output != NULL)           { print_genome_info(indiv, main_output); }
	  if ( neutral_region_output != NULL) { print_neutral_regions(indiv, neutral_region_output); }
          i++;
        }  
      }
    }
    else
    {
      ae_list_node  * indiv_node = ((ae_common::sim->get_pop())->get_indivs())->get_first();
      ae_individual * indiv      = NULL;

      while( indiv_node != NULL )
      {
        indiv = (ae_individual *) indiv_node->get_obj();
        if ( main_output != NULL)           { print_genome_info(indiv, main_output); }
	if ( neutral_region_output != NULL) { print_neutral_regions(indiv, neutral_region_output); }
        indiv_node = indiv_node->get_next();
	i++;
      }
    }
  }
  
  if ( main_output != NULL ) { fclose(main_output); }
  if ( neutral_region_output != NULL ) { fclose(neutral_region_output); }

  if ( main_output_name != NULL )           { delete [] main_output_name; }
  if ( neutral_region_output_name != NULL ) { delete [] neutral_region_output_name; }
  
  ae_common::clean();
  
  return EXIT_SUCCESS;
}


// =================================================================
//              Implementation of Secondary Functions
// =================================================================


// The export fonction
inline void print_genome_info( ae_individual* indiv, FILE* output_file )
{
  int32_t nb_bases = indiv->get_total_genome_size();
  int32_t nb_bases_in_neutral_regions = indiv->get_nb_bases_in_neutral_regions();
  int32_t nb_bases_in_0_CDS = indiv->get_nb_bases_in_0_CDS();
  int32_t nb_bases_in_0_functional_CDS = indiv->get_nb_bases_in_0_functional_CDS();
  int32_t nb_bases_in_0_non_functional_CDS = indiv->get_nb_bases_in_0_non_functional_CDS();
  int32_t nb_bases_in_0_RNA = indiv->get_nb_bases_in_0_RNA();
  int32_t nb_bases_in_0_coding_RNA = indiv->get_nb_bases_in_0_coding_RNA();
  int32_t nb_bases_in_0_non_coding_RNA = indiv->get_nb_bases_in_0_non_coding_RNA();
  if ( output_file != NULL )
  {
    fprintf( output_file, "%"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32"\n",
	     nb_bases_in_neutral_regions, nb_bases_in_0_CDS, nb_bases_in_0_functional_CDS,
	     nb_bases_in_0_non_functional_CDS, nb_bases_in_0_RNA, nb_bases_in_0_coding_RNA,
	     nb_bases_in_0_non_coding_RNA, nb_bases);
  }
}

inline void print_neutral_regions( ae_individual* indiv, FILE* output_file )
{
  if ( output_file == NULL ) return;
    
  //header
  ae_genetic_unit* chromosome = (ae_genetic_unit*) indiv->get_genetic_unit_list()->get_first()->get_obj();
  int32_t nb_neutral_regions  = chromosome->get_nb_neutral_regions();
  fprintf( output_file, "# chromosome length: %"PRId32", %"PRId32" neutral bases, %"PRId32" neutral regions\n",
	   chromosome->get_dna()->get_length(), chromosome->get_nb_bases_in_neutral_regions(), nb_neutral_regions);
  
  //neutral regions
  if ( nb_neutral_regions > 0)
  {
    int32_t* beginning_nr = chromosome->get_beginning_neutral_regions();
    int32_t* end_nr = chromosome->get_end_neutral_regions();
    
    for (int32_t i=0; i<nb_neutral_regions; i++) fprintf(output_file, "%"PRId32"\t", beginning_nr[i]);
    fprintf(output_file, "\n");
    for (int32_t i=0; i<nb_neutral_regions; i++) fprintf(output_file, "%"PRId32"\t", end_nr[i]);
    fprintf(output_file, "\n");
  }   
}

// TODO: update
void print_help( char* prog_name ) 
{
  printf( "\n\
Usage : genome_stats -h\n\
   or : genome_stats -f source [-o output_file] [-b] \n\n\
\t-h : display this screen\n\
\t--file source : read from the backup file source\n\
\t--output of : extract and save some infos about the genomes of the individuals to file of\
\t--best : only treat the best individual\n");
}  
