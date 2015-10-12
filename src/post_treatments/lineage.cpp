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
#include <errno.h>
#include <inttypes.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <sys/stat.h>

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"

using namespace aevol;

enum check_type
{
  FULL_CHECK  = 0,
  LIGHT_CHECK = 1,
  NO_CHECK    = 2
};

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);

int main(int argc, char** argv)
{
  // The output file (lineage.ae or lineage.rae) contains the following information:
  //
  // - common data                                                (ae_common::write_to_backup)
  // - begin gener                                                (int32_t)
  // - end gener                                                  (int32_t)
  // - final individual index                                     (int32_t)
  // - initial genome size                                        (int32_t)
  // - initial ancestor (nb genetic units + sequences)            (Individual::write_to_backup)
  // - replication report of ancestor at generation begin_gener+1 (ae_replic_report::write_to_backup)
  // - replication report of ancestor at generation begin_gener+2 (ae_replic_report::write_to_backup)
  // - replication report of ancestor at generation begin_gener+3 (ae_replic_report::write_to_backup)
  // - ...
  // - replication report of ancestor at generation end_gener     (ae_replic_report::write_to_backup)


  printf("\n  WARNING : Parameters' change in the middle of a simulation is not managed.\n");


  // =====================
  //  Parse command line
  // =====================

  // Default values
  check_type  check_genome      = LIGHT_CHECK;
  bool verbose = false;
  int64_t t0 = 0;
  int64_t t_end = -1;
  int32_t final_indiv_index = -1;
  int32_t final_indiv_rank  = -1;
  char tree_file_name[50];

  const char * short_options = "hVvncb:i:r:e:";
  static struct option long_options[] = {
    {"help",      no_argument,       NULL,  'h'},
    {"version",   no_argument,       NULL,  'V'},
    {"verbose",   no_argument,       NULL,  'v'},
    {"nocheck",   no_argument,       NULL,  'n'},
    {"fullcheck", no_argument,       NULL,  'c'},
    {"begin",     required_argument, NULL,  'b'},
    {"index",     required_argument, NULL,  'i'},
    {"rank",      required_argument, NULL,  'r'},
    {"end",       required_argument,  NULL, 'e' },
    {0, 0, 0, 0}
  };

  int option;
  while((option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
  {
    switch(option)
    {
      case 'h' :
      {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' :
      {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 'v' : verbose = true;                    break;
      case 'n' : check_genome = NO_CHECK;           break;
      case 'c' : check_genome = FULL_CHECK;         break;
      case 'b' : t0  = atol(optarg);                break;
      case 'i' : final_indiv_index  = atol(optarg); break;
      case 'r' : final_indiv_rank  = atol(optarg);  break;
      case 'e' :
      {
        if (strcmp(optarg, "") == 0)
        {
          printf("%s: error: Option -e or --end : missing argument.\n", argv[0]);
          exit(EXIT_FAILURE);
        }

        t_end = atol(optarg);

        break;
      }
    }
  }

  // Set undefined command line parameters to default values
  if (t_end == -1) {
    // Set t_end to the content of the LAST_GENER file if it exists.
    // If it doesn't, print help and exit
    FILE* lg_file = fopen(LAST_GENER_FNAME, "r");
    if (lg_file != NULL) {
      if (fscanf(lg_file, "%" PRId64, &t_end) == EOF) {
        printf("ERROR: failed to read last generation from file %s\n",
               LAST_GENER_FNAME);
        exit(EXIT_FAILURE);
      }
      fclose(lg_file);
    }
    else {
      printf("%s: error: You must provide a generation number.\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  // Load the simulation
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(t_end, true, false);

  // Check that the tree was recorded
  if (not exp_manager->get_record_tree()) {
    Utils::ExitWithUsrMsg("The phylogenetic tree wasn't recorded during "
                              "evolution, could not reconstruct the lineage");
  }

  int64_t tree_step = exp_manager->get_tree_step();

  //delete exp_manager;


  // The tree
  Tree* tree = NULL;

  // Indices, ranks and replication reports of the individuals in the lineage
  int32_t* indices = new int32_t[t_end - t0 + 1];
  //~ int32_t *                 ranks   = new int32_t[end_gener - begin_gener + 1];
  ReplicationReport** reports = new ReplicationReport*[t_end - t0];
  // NB: we do not need the report of the ancestor at generation begin_gener
  // (it might be the generation 0, for which we have no reports)
  // reports[0] = how ancestor at generation begin_gener + 1 was created
  // reports[i] = how ancestor at generation begin_gener + i + 1 was created
  // reports[end_gener - begin_gener - 1] = how the final individual was created
  //
  //            -----------------------------------------------------------------------------------------
  //  reports  | gener_0 => gener_1 | gener_1 => gener_2 | ... | gener_n-1 => gener_n | //////////////// |
  //            -----------------------------------------------------------------------------------------
  //  indices  |  index at gener_0  |  index at gener_1  | ... |  index at gener_n-1  | index at gener_n |
  //            -----------------------------------------------------------------------------------------



  // =========================
  //  Load the last tree file
  // =========================

  if (verbose)
  {
    printf("\n\n");
    printf("====================================\n");
    printf(" Loading the last tree file ... ");
    fflush(stdout);
  }


  // Example for ae_common::rec_params->get_tree_step() == 100 :
  //
  // tree_000100.ae ==>  generations   1 to 100.
  // tree_000200.ae ==>  generations 101 to 200.
  // tree_000300.ae ==>  generations 201 to 300.
  // etc.
  //
  // Thus, the information for generation end_gener are located
  // in the file called (end_gener/ae_common::rec_params->get_tree_step() + 1) * ae_common::rec_params->get_tree_step(),
  // except if end_gener%ae_common::rec_params->get_tree_step()==0.

  #ifdef __REGUL
    sprintf(tree_file_name,"tree/tree_%06" PRId64 ".rae", t_end);
  #else
    sprintf(tree_file_name,"tree/tree_%06" PRId64 ".ae", t_end);
  #endif

  tree = new Tree(exp_manager, tree_file_name);

  if (verbose)
  {
    printf("OK\n");
    printf("====================================\n");
  }


  // ============================================================================
  //  Find the index of the final individual and retrieve its replication report
  // ============================================================================
  if (final_indiv_index != -1)
  {
    // The index was directly provided, get the replication report and update the indices and ranks tables
    reports[t_end - t0 - 1] =
        new ReplicationReport(*(tree->get_report_by_index(t_end,
                                                          final_indiv_index)));
    final_indiv_rank = reports[t_end - t0 - 1]->rank();

    indices[t_end - t0]  = final_indiv_index;
  }
  else
  {
    if (final_indiv_rank == -1)
    {
      // No index nor rank was given in the command line.
      // By default, we construct the lineage of the best individual, the rank of which
      // is simply the number of individuals in the population.
      final_indiv_rank = exp_manager->get_nb_indivs();
    }

    // Retrieve the replication report of the individual of interest (at t_end)
    reports[t_end - t0 - 1] = new ReplicationReport(*(tree->get_report_by_rank(t_end, final_indiv_rank)));
    final_indiv_index = reports[t_end - t0 - 1]->id();

    indices[t_end - t0]  = final_indiv_index;
    //~ ranks[end_gener - begin_gener]    = final_indiv_rank;
  }

  if (verbose) printf("The final individual has the index %" PRId32 " (rank %" PRId32 ")\n", final_indiv_index, final_indiv_rank);


  // =======================
  //  Open the output file
  // =======================
  char output_file_name[101];

  #ifdef __REGUL
    snprintf(output_file_name, 100,
        "lineage-b%06" PRId64 "-e%06" PRId64 "-i%" PRId32 "-r%" PRId32 ".rae",
        t0, t_end, final_indiv_index, final_indiv_rank);
  #else
    snprintf(output_file_name, 100,
        "lineage-b%06" PRId64 "-e%06" PRId64 "-i%" PRId32 "-r%" PRId32 ".ae",
        t0, t_end, final_indiv_index, final_indiv_rank);
  #endif

  gzFile lineage_file = gzopen(output_file_name, "w");
  if (lineage_file == NULL)
  {
    fprintf(stderr, "File %s could not be created, exiting.\n", output_file_name);
    fprintf(stderr, "Please check your permissions in this directory.\n");
    exit(EXIT_FAILURE);
  }




  // ===================================================
  //  Retrieve the replication reports of the ancestors
  // ===================================================

  if (verbose)
  {
    printf("\n\n\n");
    printf("======================================================================\n");
    printf(" Parsing tree files to retrieve the ancestors' replication reports... \n");
    printf("======================================================================\n");
  }


  // Retrieve the index of the first ancestor from the last replication report
  indices[t_end - t0 -1] = reports[t_end - t0 - 1]->get_parent_id();

  // For each generation (going backwards), retrieve the index of the parent and
  // the corresponding replication report
  for (int64_t i = t_end - t0 - 2 ; i >= 0 ; i--)
  {
    int64_t t = t0 + i + 1;

    // We want to fill reports[i], that is to say, how the ancestor
    // at generation begin_gener + i + 1  was created
    if (verbose)
      printf("Getting the replication report for the ancestor at generation %" PRId64 "\n", t);

    // If we've exhausted the current tree file, load the next one
    if (Utils::mod(t, tree_step) == 0)
    {
      // Change the tree file
      delete tree;

      #ifdef __REGUL
        sprintf(tree_file_name,"tree/tree_%06" PRId64 ".rae", t);
      #else
        sprintf(tree_file_name,"tree/tree_%06" PRId64 ".ae", t);
      #endif

      tree = new Tree(exp_manager, tree_file_name);
    }

    // Copy the replication report of the ancestor
    reports[i] = new ReplicationReport(*(tree->get_report_by_index(t, indices[i + 1])));

    // Retreive the index and rank of the next ancestor from the report
    indices[i] = reports[i]->get_parent_id();
  }
  delete exp_manager;


  if (verbose)  printf("OK\n");


  // =============================================================================
  //  Get the initial genome from the backup file and write it in the output file
  // =============================================================================

  if (verbose)
  {
    printf("\n\n\n");
    printf("=============================================== \n");
    printf(" Getting the initial genome sequence... ");
    fflush(NULL);
  }

  // Load the simulation
  exp_manager = new ExpManager();
  exp_manager->load(t0, true, false);

  // Copy the initial ancestor
  // NB : The list of individuals is sorted according to the index
  const Individual* const initial_ancestor = exp_manager->get_indiv_by_id(indices[0]);

  gzwrite(lineage_file, &t0, sizeof(t0));
  gzwrite(lineage_file, &t_end, sizeof(t_end));
  gzwrite(lineage_file, &final_indiv_index, sizeof(final_indiv_index));
  gzwrite(lineage_file, &final_indiv_rank, sizeof(final_indiv_rank));

  initial_ancestor->get_grid_cell()->save(lineage_file);


  if (verbose)
  {
    printf("OK\n");
    printf("=============================================== \n");
  }


  // ===============================================================================
  //  Write the replication reports of the successive ancestors in the output file
  //  (and, optionally, check that the rebuilt genome is correct each time a backup
  //  is available)
  // ===============================================================================

  if (verbose)
  {
    printf("\n\n\n");
    printf("============================================================ \n");
    printf(" Write the replication reports in the output file... \n");
    printf("============================================================ \n");
  }

  std::list<GeneticUnit>::const_iterator unit;

  Individual* stored_indiv = nullptr;
  std::list<GeneticUnit>::const_iterator stored_gen_unit;

  ExpManager* exp_manager_backup = NULL;

  // NB: I must keep the genome encapsulated inside an Individual, because
  // replaying the mutations has side effects on the list of promoters,
  // which is stored in the individual
  bool check_genome_now = false;

  for (int64_t i = 0 ; i < t_end - t0 ; i++)
  {
    // Where are we in time...
    int64_t t = t0 + i + 1;

    // Do we need to check the genome now?
    check_genome_now = ((check_genome == FULL_CHECK &&
        Utils::mod(t, exp_manager->get_backup_step()) == 0) ||
        (check_genome == LIGHT_CHECK && t == t_end));

    // Write the replication report of the ancestor for current generation
    if (verbose)
    {
      printf("Writing the replication report for t= %" PRId64 " (built from indiv %" PRId32 " at t= %" PRId64 ")\n",
             t, indices[i], t-1);
    }
    reports[i]->write_to_tree_file(lineage_file);
    if (verbose) printf(" OK\n");


    if (check_genome_now)
    {
      // Load the simulation
      exp_manager_backup = new ExpManager();
      exp_manager_backup->load(t, true, false);

      // Copy the ancestor from the backup
      stored_indiv = exp_manager_backup->get_indiv_by_id(indices[i+1]);
      stored_gen_unit = stored_indiv->get_genetic_unit_list().cbegin();
    }


    // Warning: this portion of code won't work if the number of units changes
    // during the evolution

    // Replay the mutations stored in the current replication report on the
    // current genome
    unit = initial_ancestor->get_genetic_unit_list().cbegin();
    for (const auto& mut: reports[i]->dna_replic_report().get_HT()) {
      (unit->get_dna())->undergo_this_mutation(*mut);
    }
    for (const auto& mut: reports[i]->dna_replic_report().get_rearrangements()) {
      (unit->get_dna())->undergo_this_mutation(*mut);
    }
    for (const auto& mut: reports[i]->dna_replic_report().get_mutations()) {
      unit->get_dna()->undergo_this_mutation(*mut);
    }

    if (check_genome_now)
    {
      if (verbose)
      {
        printf("Checking the sequence of the unit...");
        fflush(stdout);
      }
      assert(stored_gen_unit != stored_indiv->get_genetic_unit_list().cend());

      char * str1 = new char[unit->get_dna()->length() + 1];
      memcpy(str1, unit->get_dna()->data(),
             unit->get_dna()->length() * sizeof(char));
      str1[unit->get_dna()->length()] = '\0';

      char * str2 = new char[stored_gen_unit->get_dna()->length() + 1];
      memcpy(str2, stored_gen_unit->get_dna()->data(),
             stored_gen_unit->get_dna()->length() * sizeof(char));
      str2[stored_gen_unit->get_dna()->length()] = '\0';

      if (strncmp(str1, str2, stored_gen_unit->get_dna()->length()) == 0)
      {
        if (verbose) printf(" OK\n");
      }
      else
      {
        if (verbose) printf(" ERROR !\n");
        fprintf(stderr, "Error: the rebuilt unit is not the same as \n");
        fprintf(stderr, "the one stored in backup file at %" PRId64 "\n", t);
        fprintf(stderr, "Rebuilt unit : %" PRId32 " bp\n %s\n",
                (int32_t)strlen(str1), str1);
        fprintf(stderr, "Stored unit  : %" PRId32 " bp\n %s\n",
                (int32_t)strlen(str2), str2);
        delete [] str1;
        delete [] str2;
        gzclose( lineage_file );
        delete exp_manager_backup;
        delete exp_manager;
        delete [] reports;
        fflush(stdout);
        exit(EXIT_FAILURE);
      }

      delete [] str1;
      delete [] str2;

      ++stored_gen_unit;
    }
    ++unit;

    assert(unit == initial_ancestor->get_genetic_unit_list().cend());
    if (check_genome_now)
    {
      assert(stored_gen_unit == stored_indiv->get_genetic_unit_list().cend());
      delete exp_manager_backup;
    }
  }


  gzclose(lineage_file);
  delete [] reports;
  delete exp_manager;

  exit(EXIT_SUCCESS);
}

/*!
  \brief

*/
void print_help(char* prog_path)
{
  // default values :
  // begin_gener = 0
  // indiv  = best individual at generation end_gener

  // there must be a genome backup file for begin_gener

  // not relevant if crossover

  printf("\n");
  printf("*********************** aevol - Artificial Evolution ******************* \n");
  printf("*                                                                      * \n");
  printf("*                      Lineage post-treatment program                  * \n");
  printf("*                                                                      * \n");
  printf("************************************************************************ \n");
  printf("\n\n");
  printf("This program is Free Software. No Warranty.\n");
  printf("Copyright (C) 2009  LIRIS.\n");
  printf("\n");
#ifdef __REGUL
  printf("Usage : rlineage -h\n");
  printf("or :    rlineage [-vn] [-i index | -r rank] [-b gener1] -e end_gener \n");
#else
  printf("Usage : lineage -h\n");
  printf("or :    lineage [-vn] [-i index | -r rank] [-b gener1] -e end_gener \n");
#endif
  printf("\n");
#ifdef __REGUL
  printf("This program retrieves the ancestral lineage of an individual and writes \n");
  printf("it in an output file called lineage.rae. Specifically, it retrieves the \n");
  printf("lineage of the individual of end_gener whose index is index, going \n");
  printf("back in time up to gener1. This program requires at least one population backup\n");
  printf("file (for the generation gener1), one environment backup file (for the generation gener1)\n");
  printf("and all tree files for generations gener1 to end_gener.\n");
#else
  printf("This program retrieves the ancestral lineage of an individual and writes \n");
  printf("it in an output file called lineage.ae. Specifically, it retrieves the \n");
  printf("lineage of the individual of end_gener whose index is index, going \n");
  printf("back in time up to gener1. This program requires at least one population backup\n");
  printf("file (for the generation gener1), one environment backup file (for the generation gener1)\n");
  printf("and all tree files for generations gener1 to end_gener.\n");
#endif
  printf("\n");
  printf("WARNING: This program should not be used for simulations run with lateral\n");
  printf("transfer. When an individual has more than one parent, the notion of lineage\n");
  printf("used here is not relevant.\n");
  printf("\n");
  printf("\t-h or --help    : Display this help.\n");
  printf("\n");
  printf("\t-v or --verbose : Be verbose, listing generations as they are \n");
  printf("\t                  treated.\n");
  printf("\n");
  printf("\t-n or --nocheck    : Disable genome sequence checking. Makes the \n");
  printf("\t                       program faster, but it is not recommended. \n");
  printf("\t                       It is better to let the program check that \n");
  printf("\t                       when we rebuild the genomes of the ancestors\n");
  printf("\t                       from the lineage file, we get the same sequences\n");
  printf("\t                       as those stored in the backup files.\n");
  printf("\n");
  printf("\t-c or --fullcheck  : Will perform the genome checks every <BACKUP_STEP>\n");
  printf("\t                       generations. Default behaviour is lighter as it\n");
  printf("\t                       only performs these checks at the ending generation.\n");
  printf("\n");
  printf("\t-i index or --index index : \n");
  printf("\t                  Retrieve the lineage of the individual whose\n");
  printf("\t                  index is index. The index must be comprised \n");
  printf("\t                  between 0 and N-1, with N the size of the \n");
  printf("\t                  population at the ending generation. If neither\n");
  printf("\t                  index nor rank are specified, the program computes \n");
  printf("\t                  the lineage of the best individual of the ending \n");
  printf("\t                  generation.\n");
  printf("\n");
  printf("\t-r rank or --rank rank : \n");
  printf("\t                  Retrieve the lineage of the individual whose\n");
  printf("\t                  rank is rank. The rank must be comprised \n");
  printf("\t                  between 1 and N, with N the size of the \n");
  printf("\t                  population at the endind generation. If neither\n");
  printf("\t                  index nor rank are specified, the program computes \n");
  printf("\t                  the lineage of the best individual of the ending \n");
  printf("\t                  generation.\n");
  printf("\n");
  printf("\t-b gener1 or --begin gener1 : \n");
  printf("\t                  Retrieve the lineage up to generation gener1.\n");
  printf("\t                  There must be a genome backup file for this\n");
  printf("\t                  generation. If not specified, the program \n");
  printf("\t                  retrieves the lineage up to generation 0.\n");
  printf("\n");
  printf("\t-e end_gener or --end end_gener : \n");
  printf("\t                  Retrieve the lineage of the individual of end_gener \n");
  printf("\t                  (default: that contained in file last_gener.txt, if any)\n");
  printf("\n");
}
