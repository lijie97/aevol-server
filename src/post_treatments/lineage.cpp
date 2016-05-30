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

// ============================================================================
//                                   Includes
// ============================================================================
#include <cerrno>
#include <cinttypes>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <zlib.h>
#include <sys/stat.h>
#include <getopt.h>

#include <list>

#include "aevol.h"

using namespace aevol;

// Helper functions
void print_help(char* prog_path);
void interpret_cmd_line_options(int argc, char* argv[]);

// Command-line option variables
static bool full_check = false;
static bool verbose = false;
static int64_t t0 = 0;
static int64_t t_end = -1;
static int32_t final_indiv_index = -1;
static int32_t final_indiv_rank  = -1;
static char tree_file_name[255]; // TODO(dpa) remove magic number

int main(int argc, char** argv) {
  // The output file (lineage.ae) contains the following information:
  // You may check that this information is up-to-date by searching
  // "lineage_file" in this source file
  //
  // - t0
  // - t_end
  // - final individual index
  // - final individual rank
  // - initial_ancestor information (including its genome)
  // - replication report of ancestor at generation t0+1
  // - replication report of ancestor at generation t0+2
  // - replication report of ancestor at generation t0+3
  // - ...
  // - replication report of ancestor at generation t_end

  interpret_cmd_line_options(argc, argv);

  printf("\n  WARNING : Parameter change in the middle of a simulation is not managed.\n");

  // Load the simulation
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(t_end, true, false);

  // Check that the tree was recorded
  if (not exp_manager->record_tree()) {
    Utils::ExitWithUsrMsg("The phylogenetic tree wasn't recorded during "
                          "evolution, could not reconstruct the lineage");
  }

  int64_t tree_step = exp_manager->tree_step();

  //delete exp_manager;


  // The tree
  Tree* tree = nullptr;

  // Indices, ranks and replication reports of the individuals in the lineage
  int32_t* indices = new int32_t[t_end - t0 + 1];
  ReplicationReport** reports = new ReplicationReport*[t_end - t0];
  // NB: we do not need the report of the ancestor at time t0 since we have
  // retrieved the individual itself from the initial backup
  // (plus it might be the generation 0, for which we have no reports)
  // reports[0] = how ancestor at t0 + 1 was created
  // reports[i] = how ancestor at t0 + i + 1 was created
  // reports[t_end - t0 - 1] = how the final individual was created
  //
  //           ---------------------------------------------------------------
  //  reports |  t0 => t1   |  t1 => t2   |...| t_n-1 => t_n   | XXXXXXXXXXXX |
  //           ---------------------------------------------------------------
  //  indices | index at t0 | index at t1 |...| index at t_n-1 | index at t_n |
  //           ---------------------------------------------------------------



  // =========================
  //  Load the last tree file
  // =========================

  if (verbose) {
    printf("\n\n");
    printf("====================================\n");
    printf(" Loading the last tree file ... ");
    fflush(stdout);
  }


  // Example for ae_common::rec_params->tree_step() == 100 :
  //
  // tree_000100.ae ==>  timesteps 1 to 100.
  // tree_000200.ae ==>  timesteps 101 to 200.
  // tree_000300.ae ==>  timesteps 201 to 300.
  // etc.
  //

  #ifdef __REGUL
    sprintf(tree_file_name,"tree/tree_" TIMESTEP_FORMAT ".rae", t_end);
  #else
    sprintf(tree_file_name,"tree/tree_" TIMESTEP_FORMAT ".ae", t_end);
  #endif

  tree = new Tree(exp_manager, tree_file_name);

  if (verbose) {
    printf("OK\n");
    printf("====================================\n");
  }


  // ============================================================================
  //  Find the index of the final individual and retrieve its replication report
  // ============================================================================
  if (final_indiv_index != -1) {
    // The index was directly provided, get the replication report and update the indices and ranks tables
    reports[t_end - t0 - 1] =
        new ReplicationReport(*(tree->report_by_index(t_end,
                                                      final_indiv_index)));
    final_indiv_rank = reports[t_end - t0 - 1]->rank();

    indices[t_end - t0]  = final_indiv_index;
  }
  else {
    if (final_indiv_rank == -1) {
      // No index nor rank was given in the command line.
      // By default, we construct the lineage of the best individual, the rank of which
      // is simply the number of individuals in the population.
      final_indiv_rank = exp_manager->nb_indivs();
    }

    // Retrieve the replication report of the individual of interest (at t_end)
    reports[t_end - t0 - 1] =
        new ReplicationReport(*(tree->report_by_rank(t_end, final_indiv_rank)));
    final_indiv_index = reports[t_end - t0 - 1]->id();

    indices[t_end - t0]  = final_indiv_index;
    //~ ranks[end_gener - begin_gener]    = final_indiv_rank;
  }

  if (verbose) {
    printf("The final individual has index %" PRId32
           " (rank %" PRId32 ")\n", final_indiv_index, final_indiv_rank);
  }


  // =======================
  //  Open the output file
  // =======================
  char output_file_name[101];

  #ifdef __REGUL
    snprintf(output_file_name, 100,
        "lineage-b" TIMESTEP_FORMAT "-e" TIMESTEP_FORMAT "-i%" PRId32 "-r%" PRId32 ".rae",
        t0, t_end, final_indiv_index, final_indiv_rank);
  #else
    snprintf(output_file_name, 100,
        "lineage-b" TIMESTEP_FORMAT "-e" TIMESTEP_FORMAT "-i%" PRId32 "-r%" PRId32 ".ae",
        t0, t_end, final_indiv_index, final_indiv_rank);
  #endif

  gzFile lineage_file = gzopen(output_file_name, "w");
  if (lineage_file == nullptr) {
    fprintf(stderr, "File %s could not be created.\n", output_file_name);
    fprintf(stderr, "Please check your permissions in this directory.\n");
    exit(EXIT_FAILURE);
  }




  // ===================================================
  //  Retrieve the replication reports of the ancestors
  // ===================================================

  if (verbose) {
    printf("\n\n\n");
    printf("======================================================================\n");
    printf(" Parsing tree files to retrieve the ancestors' replication reports... \n");
    printf("======================================================================\n");
  }


  // Retrieve the index of the first ancestor from the last replication report
  indices[t_end - t0 -1] = reports[t_end - t0 - 1]->parent_id();

  // For each generation (going backwards), retrieve the index of the parent and
  // the corresponding replication report
  for (int64_t i = t_end - t0 - 2 ; i >= 0 ; i--) {
    int64_t t = t0 + i + 1;

    // We want to fill reports[i], that is to say, how the ancestor
    // at generation begin_gener + i + 1  was created
    if (verbose)
      printf("Getting the replication report for the ancestor at generation %" PRId64 "\n", t);

    // If we've exhausted the current tree file, load the next one
    if (Utils::mod(t, tree_step) == 0) {
      // Change the tree file
      delete tree;

      #ifdef __REGUL
        sprintf(tree_file_name,"tree/tree_" TIMESTEP_FORMAT ".rae", t);
      #else
        sprintf(tree_file_name,"tree/tree_" TIMESTEP_FORMAT ".ae", t);
      #endif

      tree = new Tree(exp_manager, tree_file_name);
    }

    // Copy the replication report of the ancestor
    reports[i] =
        new ReplicationReport(*(tree->report_by_index(t, indices[i + 1])));

    // Retreive the index and rank of the next ancestor from the report
    indices[i] = reports[i]->parent_id();
  }
  delete exp_manager;


  if (verbose) printf("OK\n");


  // =============================================================================
  //  Get the initial genome from the backup file and write it in the output file
  // =============================================================================

  if (verbose) {
    printf("\n\n\n");
    printf("=============================================== \n");
    printf(" Getting the initial genome sequence... ");
    fflush(nullptr);
  }

  // Load the simulation
  exp_manager = new ExpManager();
  exp_manager->load(t0, true, false);

  // Copy the initial ancestor
  // NB : The list of individuals is sorted according to the index
  const Individual& initial_ancestor = *(exp_manager->indiv_by_id(indices[0]));

  // Write file "header"
  gzwrite(lineage_file, &t0, sizeof(t0));
  gzwrite(lineage_file, &t_end, sizeof(t_end));
  gzwrite(lineage_file, &final_indiv_index, sizeof(final_indiv_index));
  gzwrite(lineage_file, &final_indiv_rank, sizeof(final_indiv_rank));

  initial_ancestor.grid_cell()->save(lineage_file);


  if (verbose) {
    printf("OK\n");
    printf("=============================================== \n");
  }


  // ===============================================================================
  //  Write the replication reports of the successive ancestors in the output file
  //  (and, optionally, check that the rebuilt genome is correct each time a backup
  //  is available)
  // ===============================================================================

  if (verbose) {
    printf("\n\n\n");
    printf("============================================================ \n");
    printf(" Write the replication reports in the output file... \n");
    printf("============================================================ \n");
  }

  std::list<GeneticUnit>::const_iterator unit;

  Individual* stored_indiv = nullptr;
  std::list<GeneticUnit>::const_iterator stored_gen_unit;

  ExpManager* exp_manager_backup = nullptr;

  // NB: I must keep the genome encapsulated inside an Individual, because
  // replaying the mutations has side effects on the list of promoters,
  // which is stored in the individual
  bool check_genome_now = false;

  for (int64_t i = 0 ; i < t_end - t0 ; i++) {
    // Where are we in time...
    int64_t t = t0 + i + 1;

    // Do we need to check the genome now?
    check_genome_now = t == t_end ||
        (full_check && Utils::mod(t, exp_manager->backup_step()) == 0);

    // Write the replication report of the ancestor for current generation
    if (verbose) {
      printf("Writing the replication report for t= %" PRId64
             " (built from indiv %" PRId32 " at t= %" PRId64 ")\n",
             t, indices[i], t-1);
    }
    reports[i]->write_to_tree_file(lineage_file);
    if (verbose) printf(" OK\n");


    if (check_genome_now) {
      // Load the simulation
      exp_manager_backup = new ExpManager();
      exp_manager_backup->load(t, true, false);

      // Copy the ancestor from the backup
      stored_indiv = exp_manager_backup->indiv_by_id(indices[i + 1]);
      stored_gen_unit = stored_indiv->genetic_unit_list().cbegin();
    }


    // Warning: this portion of code won't work if the number of units changes
    // during the evolution

    // Replay the mutations stored in the current replication report on the
    // current genome
    unit = initial_ancestor.genetic_unit_list().cbegin();
    for (const auto& mut: reports[i]->dna_replic_report().HT()) {
      (unit->dna())->undergo_this_mutation(*mut);
    }
    for (const auto& mut: reports[i]->dna_replic_report().rearrangements()) {
      (unit->dna())->undergo_this_mutation(*mut);
    }
    for (const auto& mut: reports[i]->dna_replic_report().mutations()) {
      unit->dna()->undergo_this_mutation(*mut);
    }

    if (check_genome_now) {
      if (verbose) {
        printf("Checking the sequence of the unit...");
        fflush(stdout);
      }
      assert(stored_gen_unit != stored_indiv->genetic_unit_list().cend());

      char * str1 = new char[unit->dna()->length() + 1];
      memcpy(str1, unit->dna()->data(),
             unit->dna()->length() * sizeof(char));
      str1[unit->dna()->length()] = '\0';

      char * str2 = new char[stored_gen_unit->dna()->length() + 1];
      memcpy(str2, stored_gen_unit->dna()->data(),
             stored_gen_unit->dna()->length() * sizeof(char));
      str2[stored_gen_unit->dna()->length()] = '\0';

      if (strncmp(str1, str2, stored_gen_unit->dna()->length()) == 0) {
        if (verbose) printf(" OK\n");
      }
      else {
        if (verbose) printf(" ERROR !\n");
        fprintf(stderr, "Error: the rebuilt unit is not the same as \n");
        fprintf(stderr, "the one stored in backup file at %" PRId64 "\n", t);
        fprintf(stderr, "Rebuilt unit : %" PRId32 " bp\n %s\n",
                (int32_t)strlen(str1), str1);
        fprintf(stderr, "Stored unit  : %" PRId32 " bp\n %s\n",
                (int32_t)strlen(str2), str2);
        delete [] str1;
        delete [] str2;
        gzclose(lineage_file);
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

    assert(unit == initial_ancestor.genetic_unit_list().cend());
    if (check_genome_now) {
      assert(stored_gen_unit == stored_indiv->genetic_unit_list().cend());
      delete exp_manager_backup;
    }
  }


  gzclose(lineage_file);
  delete [] reports;
  delete exp_manager;

  exit(EXIT_SUCCESS);
}

/**
 * \brief print help and exist
 */
void print_help(char* prog_path) {
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) {
    prog_name++;
  }
  else {
    prog_name = prog_path;
  }

  printf("******************************************************************************\n");
  printf("*                                                                            *\n");
  printf("*                        aevol - Artificial Evolution                        *\n");
  printf("*                                                                            *\n");
  printf("* Aevol is a simulation platform that allows one to let populations of       *\n");
  printf("* digital organisms evolve in different conditions and study experimentally  *\n");
  printf("* the mechanisms responsible for the structuration of the genome and the     *\n");
  printf("* transcriptome.                                                             *\n");
  printf("*                                                                            *\n");
  printf("******************************************************************************\n");
  printf("\n");
  printf("%s:\n", prog_name);
  printf("\tReconstruct the lineage of a given individual from the tree files\n");
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s [-b TIMESTEP] [-e TIMESTEP] [-I INDEX | -R RANK] [-F] [-v]\n",
         prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n");
  printf("  -V, --version\n\tprint version number, then exit\n");
  printf("  -b, --begin TIMESTEP\n");
  printf("\tspecify time t0 up to which to reconstruct the lineage\n");
  printf("  -e, --end TIMESTEP\n");
  printf("\tspecify time t_end of the indiv whose lineage is to be reconstructed\n");
  printf("  -I, --index INDEX\n");
  printf("\tspecify the index of the indiv whose lineage is to be reconstructed\n");
  printf("  -R, --rank RANK\n");
  printf("\tspecify the rank of the indiv whose lineage is to be reconstructed\n");
  printf("  -F, --full-check\n");
  printf("\tperform genome checks whenever possible\n");
  printf("  -v, --verbose\n\tbe verbose\n");
}

void interpret_cmd_line_options(int argc, char* argv[]) {
  // Define allowed options
  const char * short_options = "hVb:e:FI:R:v";
  static struct option long_options[] = {
      {"help",      no_argument,       nullptr, 'h'},
      {"version",   no_argument,       nullptr, 'V'},
      {"begin",     required_argument, nullptr, 'b'},
      {"end",       required_argument, nullptr, 'e'},
      {"fullcheck", no_argument,       nullptr, 'F'},
      {"index",     required_argument, nullptr, 'I'},
      {"rank",      required_argument, nullptr, 'R'},
      {"verbose",   no_argument,       nullptr, 'v'},
      {0, 0, 0, 0}
  };

  // Get actual values of the command-line options
  int option;
  while((option = getopt_long(argc, argv, short_options,
                              long_options, nullptr)) != -1) {
    switch(option) {
      case 'h' : {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' : {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 'b' : {
        if (strcmp(optarg, "") == 0) {
          printf("%s: error: Option -b or --begin : missing argument.\n",
                 argv[0]);
          exit(EXIT_FAILURE);
        }
        t0  = atol(optarg);
        break;
      }
      case 'e' : {
        if (strcmp(optarg, "") == 0) {
          printf("%s: error: Option -e or --end : missing argument.\n",
                 argv[0]);
          exit(EXIT_FAILURE);
        }
        t_end = atol(optarg);
        break;
      }
      case 'F' : {
        full_check = true;
        break;
      }
      case 'I' : {
        final_indiv_index  = atoi(optarg);
        break;
      }
      case 'R' : {
        final_indiv_rank  = atoi(optarg);
        break;
      }
      case 'v' : {
        verbose = true;
        break;
      }
      default : {
        // An error message is printed in getopt_long, we just need to exit
        exit(EXIT_FAILURE);
      }
    }
  }

  // If t_end wasn't provided, use default
  if (t_end < 0) {
    t_end = OutputManager::last_gener();
  }
}
