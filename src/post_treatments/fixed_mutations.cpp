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

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include "ae_utils.h"
#include "ae_exp_manager.h"
#include "ae_individual.h"
#include "ae_genetic_unit.h"
#include "ae_tree.h"
#include "ae_replication_report.h"
#include "dna_replic_report.h"
#include "ae_mutation.h"

using namespace aevol;

enum check_type
{
  FULL_CHECK  = 0,
  LIGHT_CHECK = 1,
  ENV_CHECK   = 2,
  NO_CHECK    = 3
};



// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);
void print_version( void );

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




  // =====================
  //  Parse command line
  // =====================

  // Default values
  check_type  check               = LIGHT_CHECK;
  char*       lineage_file_name   = NULL;
  bool        verbose             = false;
  double      tolerance           = 0;


  const char * short_options = "hVvncf:t:";
  static struct option long_options[] =
  {
    {"help",      no_argument,       NULL, 'h'},
    {"version",   no_argument,       NULL, 'V' },
    {"verbose",   no_argument,       NULL, 'v'},
    {"nocheck",   no_argument,       NULL, 'n'},
    {"fullcheck", no_argument,       NULL, 'c'},
    {"file",      required_argument, NULL, 'f'},
    {"tolerance",   required_argument, NULL, 't'},
    {0, 0, 0, 0}
  };

  int option;
  while( (option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1 )
  {
    switch( option )
    {
      case 'h' :
      {
        print_help(argv[0]);
        exit( EXIT_SUCCESS );
      }
      case 'V' :
      {
        print_version();
        exit( EXIT_SUCCESS );
      }
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
        print_help(argv[0]);
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
  printf( "WARNING : Parameter change during simulation is not managed in general.\n" );
  printf( "          Only changes in environmental target done with aevol_modify are handled.\n" );
  printf("\n");


  // =======================
  //  Open the lineage file
  // =======================

  gzFile lineage_file = gzopen( lineage_file_name, "r" );
  if (lineage_file == Z_NULL)
  {
    fprintf( stderr, "ERROR : Could not read the lineage file %s\n", lineage_file_name );
    exit( EXIT_FAILURE );
  }

  int32_t begin_gener, end_gener, final_index, num_gener, final_indiv_rank;
  gzread(lineage_file, &begin_gener, sizeof(begin_gener));
  gzread(lineage_file, &end_gener, sizeof(end_gener));
  gzread(lineage_file, &final_index, sizeof(final_index) );
  gzread(lineage_file, &final_indiv_rank,   sizeof(final_indiv_rank) );

  if ( verbose )
  {
    printf("\n\n");
    printf("================================================================================\n");
    printf(" Statistics of the ancestors of indiv. #%" PRId32 " (t=%" PRId32 " to %" PRId32 ")\n",
            final_index, begin_gener, end_gener);
    printf("================================================================================\n");
  }


  // =========================
  //  Open the experience manager
  // =========================

  // Open the experiment manager
  ae_exp_manager* exp_manager = new ae_exp_manager();
  exp_manager->load( begin_gener, false, true, false );
  Environment* env = new Environment( *(exp_manager->get_env()) ); // independent copy

  int32_t backup_step = exp_manager->get_backup_step();


  // =========================
  //  Open the output file(s)
  // =========================

  char output_file_name[60];
  snprintf( output_file_name, 60, "stats/fixedmut-b%06" PRId32 "-e%06" PRId32 "-i%" PRId32 "-r%" PRId32 ".out",begin_gener, end_gener, final_index, final_indiv_rank );

  FILE * output = fopen( output_file_name, "w" );
  if ( output == NULL )
  {
    fprintf( stderr, "ERROR : Could not create the output file %s\n", output_file_name );
    exit( EXIT_FAILURE );
  }


  // Write the header
  fprintf( output, "# #################################################################\n" );
  fprintf( output, "#  Mutations in the lineage of the best indiv at generation %" PRId32 "\n", end_gener );
  fprintf( output, "# #################################################################\n" );
  fprintf( output, "#  1.  Generation       (mut. occurred when producing the indiv. of this generation)\n" );
  fprintf( output, "#  2.  Genetic unit     (which underwent the mutation, 0 = chromosome) \n" );
  fprintf( output, "#  3.  Mutation type    (0: switch, 1: smallins, 2: smalldel, 3:dupl, 4: del, 5:trans, 6:inv, 7:insert, 8:ins_HT, 9:repl_HT) \n" );
  fprintf( output, "#  4.  pos_0            (position for the small events, begin_segment for the rearrangements, begin_segment of the inserted segment for ins_HT, begin_segment of replaced segment for repl_HT) \n" );
  fprintf( output, "#  5.  pos_1            (-1 for the small events, end_segment for the rearrangements, end_segment of the inserted segment for ins_HT, begin_segment of donor segment for repl_HT) \n" );
  fprintf( output, "#  6.  pos_2            (reinsertion point for duplic., cutting point in segment for transloc., insertion point in the receiver for ins_HT, end_segment of the replaced segment for repl_HT, -1 for other events)\n" );
  fprintf( output, "#  7.  pos_3            (reinsertion point for transloc., breakpoint in the donor for ins_HT, end_segment of the donor segment for repl_HT, -1 for other events)\n" );
  fprintf( output, "#  8.  invert           (transloc, was the segment inverted (0/1)?, sense of insertion for ins_HT (0=DIRECT, 1=INDIRECT), sense of the donor segment for repl_HT (0=DIRECT, 1=INDIRECT),-1 for other events)\n" );
  fprintf( output, "#  9.  align_score      (score that was needed for the rearrangement to occur, score of the first alignment for ins_HT and repl_HT)\n" );
  fprintf( output, "#  10. align_score2     (score for the reinsertion for transloc, score of the second alignment for ins_HT and repl_HT)\n" );
  fprintf( output, "#  11. seg_len          (segment length for rearrangement, donor segment length for ins_HT and repl_HT)\n" );
  fprintf( output, "#  12. repl_seg_len     (replaced segment length for repl_HT, -1 for the others)\n" );
  fprintf( output, "#  13. GU_length        (before the event)\n" );
  fprintf( output, "#  14. Impact of the mutation on the metabolic error (negative value = smaller gap after = beneficial mutation) \n" );
  fprintf( output, "#  15. Number of coding RNAs possibly disrupted by the breakpoints \n" );
  fprintf( output, "#  16. Number of coding RNAs completely included in the segment (donor segment in the case of a transfer) \n" );
  fprintf( output, "#  17. Number of coding RNAs that were completely included in the replaced segment (meaningful only for repl_HT) \n" );
  fprintf( output, "####################################################################################################################\n" );
  fprintf( output, "#\n" );
  fprintf( output, "# Header for R\n" );
  fprintf( output, "gener gen_unit mut_type pos_0 pos_1 pos_2 pos_3 invert align_score align_score_2 seg_len repl_seg_len GU_len impact nbgenesatbreak nbgenesinseg nbgenesinreplseg\n" );



  // ==============================
  //  Prepare the initial ancestor
  // ==============================

  ae_individual * indiv = new ae_individual(exp_manager, lineage_file );
  indiv->evaluate( env );
  indiv->compute_statistical_data();

  if ( verbose )
  {
    printf("Initial fitness     = %f\n", indiv->get_fitness());
    printf("Initial genome size = %" PRId32 "\n", indiv->get_total_genome_size());
  }



  // ===============================================================================
  //  Replay the mutation to get the successive ancestors and analyze them
  //  (and, optionally, check that the rebuilt envir and genome are correct each
  //   time a backup is available)
  // ===============================================================================

  ae_individual* stored_indiv = NULL;
  std::list<ae_genetic_unit>::const_iterator stored_unit;

  int32_t i, index, genetic_unit_number, unitlen_before;
  int32_t nb_genes_at_breakpoints, nb_genes_in_segment, nb_genes_in_replaced_segment;
  double metabolic_error_before, metabolic_error_after, impact_on_metabolic_error;
  char mut_descr_string[80];

  ae_exp_manager* exp_manager_backup = NULL;
  Environment* backup_env = NULL;

  bool check_now = false;

  for ( i = 0; i < end_gener - begin_gener; i++ )
  {
    num_gener = begin_gener + i + 1;  // where are we in time...

    ae_replication_report* rep = new ae_replication_report( lineage_file, indiv );
    index = rep->get_id(); // who are we building...
    indiv->set_replication_report(rep);

    // Check now?
    check_now = ( ( check == FULL_CHECK && ae_utils::mod( num_gener, backup_step ) == 0 ) ||
                  ( check == ENV_CHECK && ae_utils::mod( num_gener, backup_step ) == 0 ) ||
                  ( check == LIGHT_CHECK && num_gener == end_gener ) );


    if ( verbose ) printf("Rebuilding ancestor at generation %" PRId32 " (index %" PRId32 ")...", num_gener, index);

    env->build();
    env->apply_variation();
    indiv->reevaluate(env);

    // Check, and possibly update, the environment according to the backup files (update necessary if the env. was modified by aevol_modify at some point)
    if (ae_utils::mod( num_gener, backup_step ) == 0)
      {
        char env_file_name[255];
        sprintf( env_file_name, "./" ENV_FNAME_FORMAT, num_gener );
        gzFile env_file = gzopen( env_file_name, "r" );
        backup_env = new Environment();
        backup_env->load( env_file );

        if ( ! env->is_identical_to(*backup_env, tolerance) )
          {
            printf("Warning: At t=%" PRId32 ", the replayed environment is not the same\n", num_gener);
            printf("         as the one saved at generation %" PRId32 "... \n", num_gener );
            printf("         with tolerance of %lg\n", tolerance);
            printf("Replacing the replayed environment by the one stored in the backup.\n");
            delete env;
            env = new Environment(*backup_env);
          }
        delete backup_env;
      }


    // Warning: this portion of code won't work if the number of units changes
    // during the evolution, or if some translocations occurred between different genetic units

    genetic_unit_number = 0;
    std::list<DnaReplicReport*>::const_iterator dnareport = rep->get_dna_replic_reports().begin();
    std::list<ae_genetic_unit>::iterator unit = indiv->get_genetic_unit_list_std_nonconst().begin();

    if ( check_now && ae_utils::mod(num_gener, backup_step) == 0)
    {

      exp_manager_backup = new ae_exp_manager();
      exp_manager_backup->load( num_gener, false, true, false );
      stored_indiv = new ae_individual( * (ae_individual *)exp_manager_backup->get_indiv_by_id( index ), false );
      stored_unit = stored_indiv->get_genetic_unit_list_std().begin();
    }

    while (dnareport != rep->get_dna_replic_reports().end())
    {
      assert(unit != indiv->get_genetic_unit_list_std().end());

      unit->get_dna()->set_replic_report(*dnareport);

      // ***************************************
      //             Transfer events
      // ***************************************

      for (const auto& mutation: (*dnareport)->get_HT()) {
        metabolic_error_before = indiv->get_dist_to_target_by_feature( METABOLISM );
        unitlen_before = unit->get_dna()->get_length();
        unit->compute_nb_of_affected_genes(&mutation, nb_genes_at_breakpoints, nb_genes_in_segment, nb_genes_in_replaced_segment);


        unit->get_dna()->undergo_this_mutation(&mutation);
        indiv->reevaluate(env);


        metabolic_error_after = indiv->get_dist_to_target_by_feature( METABOLISM );
        impact_on_metabolic_error = metabolic_error_after - metabolic_error_before;


        mutation.get_generic_description_string( mut_descr_string );
        fprintf( output, "%" PRId32 " %" PRId32 " %s %" PRId32 " %.15f  %" PRId32 " %" PRId32 " %" PRId32 " \n",\
                 num_gener, genetic_unit_number, \
                 mut_descr_string, unitlen_before, \
                 impact_on_metabolic_error, nb_genes_at_breakpoints, nb_genes_in_segment, nb_genes_in_replaced_segment );
      }


      // ***************************************
      //           Rearrangement events
      // ***************************************

      for (const auto& mutation: (*dnareport)->get_rearrangements()) {
        metabolic_error_before = indiv->get_dist_to_target_by_feature( METABOLISM );
        unitlen_before = unit->get_dna()->get_length();
        unit->compute_nb_of_affected_genes(&mutation, nb_genes_at_breakpoints, nb_genes_in_segment,  nb_genes_in_replaced_segment);

        unit->get_dna()->undergo_this_mutation(&mutation);

        indiv->reevaluate(env);
        metabolic_error_after = indiv->get_dist_to_target_by_feature( METABOLISM );
        impact_on_metabolic_error = metabolic_error_after - metabolic_error_before;

        mutation.get_generic_description_string( mut_descr_string );
        fprintf( output, "%" PRId32 " %" PRId32 " %s %" PRId32 " %.15f %" PRId32 " %" PRId32 " %" PRId32 " \n",\
                 num_gener, genetic_unit_number, \
                 mut_descr_string, unitlen_before, \
                 impact_on_metabolic_error, nb_genes_at_breakpoints, nb_genes_in_segment,  nb_genes_in_replaced_segment );
      }


      // ***************************************
      // Local events (point mutations & small indels)
      // ***************************************

      for (const auto& mutation: (*dnareport)->get_mutations()) {
        metabolic_error_before = indiv->get_dist_to_target_by_feature( METABOLISM );
        unitlen_before = unit->get_dna()->get_length();
        unit->compute_nb_of_affected_genes(&mutation, nb_genes_at_breakpoints, nb_genes_in_segment, nb_genes_in_replaced_segment);

        unit->get_dna()->undergo_this_mutation(&mutation);

        indiv->reevaluate(env);
        metabolic_error_after = indiv->get_dist_to_target_by_feature( METABOLISM );
        impact_on_metabolic_error = metabolic_error_after - metabolic_error_before;

        mutation.get_generic_description_string( mut_descr_string );
        fprintf( output, "%" PRId32 " %" PRId32 " %s %" PRId32 " %.15f %" PRId32 " %" PRId32 " %" PRId32 " \n",\
                 num_gener, genetic_unit_number, \
                 mut_descr_string, unitlen_before, \
                 impact_on_metabolic_error, nb_genes_at_breakpoints, nb_genes_in_segment, nb_genes_in_replaced_segment );
      }

      if ( check_now && ae_utils::mod(num_gener, backup_step) == 0)
      {
        if ( verbose )
        {
          printf("Checking the sequence of the unit...");
          fflush(NULL);
        }

        assert(stored_unit != stored_indiv->get_genetic_unit_list_std().end());

        char * str1 = new char[unit->get_dna()->get_length() + 1];
        memcpy(str1, unit->get_dna()->get_data(), \
               unit->get_dna()->get_length()*sizeof(char));
        str1[unit->get_dna()->get_length()] = '\0';

        char * str2 = new char[(stored_unit->get_dna())->get_length() + 1];
        memcpy(str2, (stored_unit->get_dna())->get_data(), (stored_unit->get_dna())->get_length()*sizeof(char));
        str2[(stored_unit->get_dna())->get_length()] = '\0';

        if(strncmp(str1,str2, (stored_unit->get_dna())->get_length())==0)
        {
          if ( verbose ) printf(" OK\n");
        }
        else
        {
          if ( verbose ) printf( " ERROR !\n" );
          fprintf( stderr, "Error: the rebuilt unit is not the same as \n");
          fprintf( stderr, "the one saved at generation %" PRId32 "... ", begin_gener );
          fprintf( stderr, "Rebuilt unit : %zu bp\n %s\n", strlen(str1), str1 );
          fprintf( stderr, "Stored unit  : %zu bp\n %s\n", strlen(str2), str2 );
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

        ++stored_unit;
      }


      ++dnareport;
      ++unit;
      genetic_unit_number++;
    }

    assert(unit == indiv->get_genetic_unit_list_std().end());


    if ( verbose ) printf(" OK\n");

    delete rep;

    if ( check_now && ae_utils::mod(num_gener, backup_step) == 0 )
    {
      assert(stored_unit == stored_indiv->get_genetic_unit_list_std().end());
      delete stored_indiv;
      delete exp_manager_backup;
    }
  }

  gzclose(lineage_file);
  fclose(output);
  delete exp_manager;
  delete indiv;
  delete env;

  exit(EXIT_SUCCESS);

}




/*!
  \brief

*/
void print_help(char* prog_path)
{
  printf( "\n" );
  printf( "*********************** aevol - Artificial Evolution ******************* \n" );
  printf( "*                                                                      * \n" );
  printf( "*               Fixed mutations post-treatment program                 * \n" );
  printf( "*                                                                      * \n" );
  printf( "************************************************************************ \n" );
  printf( "\n\n" );
  printf( "This program is Free Software. No Warranty.\n" );
  printf( "Copyright (C) 2009  LIRIS.\n" );
  printf( "\n" );
#ifdef __REGUL
  printf( "Usage : rfixed_mutations -h\n");
  printf( "or :    rfixed_mutations [-vn] -f lineage_file \n" );
#else
  printf( "Usage : fixed_mutations -h\n");
  printf( "or :    fixed_mutations [-vn] -f lineage_file \n" );
#endif
  printf( "\n" );
  printf( "This program computes the fixed mutations of the individuals within the lineage\n" );
  printf( "of lineage_file\n" );
  printf( "\n" );
  printf( "WARNING: This program should not be used for simulations run with several genetic units\n" );
  printf( "(eg chromosome + plasmids). It works only for simulations with only a single chromosome.\n" );
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
  printf( "\t                       Compute the fixed mutations of the individuals within lineage_file.\n" );
  printf( "\n" );
  printf( "\t-t tolerance or --tolerance tolerance : \n");
  printf( "\t                       Tolerance used to compare the replayed environment to environment in backup\n");
  printf( "\n" );
}


/*!
  \brief Print aevol version number

*/
void print_version( void )
{
  printf( "aevol %s\n", VERSION );
}
