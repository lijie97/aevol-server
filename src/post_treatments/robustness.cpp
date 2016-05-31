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




// =================================================================
//                              Libraries
// =================================================================
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <inttypes.h>
#include <getopt.h>
#include <math.h>
#include <sys/stat.h>

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"
#include "IndivAnalysis.h"

using std::list;
using namespace aevol;

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);





// =====================================================================
//                         Main Function
// =====================================================================


int main(int argc, char* argv[]) {
  // ----------------------------------------
  //     command-line option parsing
  // ----------------------------------------
  int32_t nb_children = 1000;
  int32_t wanted_rank = -1;
  int32_t wanted_index = -1;
  int64_t timestep = 0;

  const char* options_list = "hVg:n:r:i:";
  static struct option long_options_list[] = {
      {"help",        no_argument,       NULL, 'h'},
      {"version",     no_argument,       NULL, 'V'},
      {"gener",       required_argument, NULL, 'g'},
      {"nb-children", required_argument, NULL, 'n'},
      {"rank",        required_argument, NULL, 'r'},
      {"index",       required_argument, NULL, 'i'},
      {0,             0,                 0,    0}
  };

  int option = -1;
  bool rank_already_set = false;
  bool index_already_set = false;
  while ((option = getopt_long(argc, argv, options_list, long_options_list,
                               NULL)) != -1) {
    switch (option) {
      case 'h' : {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' : {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 'g' : {
        if (strcmp(optarg, "") == 0) {
          printf("%s: error: Option -g or --gener : missing argument.\n",
                 argv[0]);
          exit(EXIT_FAILURE);
        }

        timestep = atol(optarg);
        break;
      }
      case 'n' :
        nb_children = atol(optarg);
        break;
      case 'r' :
        if (index_already_set) {
          fprintf(stderr,
                  "%s: error: Options -r and -i are incompatible. Please choose one of them only.\n",
                  argv[0]);
          exit(EXIT_FAILURE);
        }
        wanted_rank = atol(optarg);
        rank_already_set = true;
        break;
      case 'i' :
        if (rank_already_set) {
          fprintf(stderr,
                  "%s: error: Options -r and -i are incompatible. Please choose one of them only.\n",
                  argv[0]);
          fprintf(stderr, "           Use %s --help for more information.\n",
                  argv[0]);
          exit(EXIT_FAILURE);
        }
        wanted_index = atol(optarg);
        index_already_set = true;
        break;
    }
  }

  
  // ----------------------
  //  Prepare the outputs
  // ----------------------
  char directory_name[255];
  snprintf(directory_name, 255, "analysis-generation_" TIMESTEP_FORMAT,
      timestep);
  
  // Check whether the directory already exists and is writable
  if (access(directory_name, F_OK) == 0) {
    if (access(directory_name, X_OK | W_OK) != 0) {
      fprintf(stderr, "Error: cannot enter or write in directory %s.\n",
              directory_name);
      exit(EXIT_FAILURE);
    }
  }
  else {
    // Create the directory with permissions : rwx r-x r-x
    if (mkdir(directory_name,
              S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0) {
      fprintf(stderr, "Error: cannot create directory %s.\n", directory_name);
      exit(EXIT_FAILURE);
    }
  }


  // -----------------------------------------------------------------------------------
  //  Load the backup and get the individual for which detailed information is desired
  // -----------------------------------------------------------------------------------
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(timestep, true, false);

  Individual* indiv_tmp = nullptr;

  // If neither a rank nor an index was provided, consider the best indiv
  if ((wanted_rank == -1) && (wanted_index == -1)) {
    indiv_tmp = exp_manager->best_indiv();
  }
  else if (wanted_index != -1) {
    indiv_tmp = exp_manager->indiv_by_id(wanted_index);
  }
  else {
    indiv_tmp = exp_manager->indiv_by_rank(wanted_rank);
  }
  // Update id and rank
  wanted_index = indiv_tmp->id();
  wanted_rank = indiv_tmp->rank();

  IndivAnalysis wanted_indiv(*indiv_tmp);
  indiv_tmp = nullptr;

  // Now that we have the index and rank of the indiv of interest, we can
  // generate the output file name and hence open that file
  char filename[255];
  snprintf(filename, 255, "%s/robustness-%" PRId64 "-i%" PRId32 "-r%" PRId32,
      directory_name, timestep, wanted_index, wanted_rank);
  FILE* output = fopen(filename, "w");
  if (output == nullptr) {
    Utils::ExitWithUsrMsg(std::string("Could not open file ") + filename);
  }

  fprintf(output, "###############################################################################\n");
  fprintf(output, "#  Summary of the mutants generated from individual with rank %" PRId32
  " and index %" PRId32 " at timestep %" PRId64 " \n",
      wanted_rank, wanted_index, timestep);
  fprintf(output, "###############################################################################\n");
  fprintf(output, "#  1.  Timestep\n");
  fprintf(output, "#  2.  Proportion of mutants that are better than their parent\n");
  fprintf(output, "#  3.  Proportion of mutants that are as good as their parent was\n");
  fprintf(output, "#  4.  Proportion of mutants that are worse than their parent\n");
  fprintf(output, "#  5.  Average difference in metabolic error btw good mutants and their parent\n");
  fprintf(output, "#  6.  Average difference in metabolic error btw bad mutants and their parent\n");
  fprintf(output, "#  7.  Maximum gain in metabolic error among good mutants\n");
  fprintf(output, "#  8.  Maximum loss of metabolic error among bad mutants\n");
  fprintf(output, "###############################################################################\n");

  /*fprintf(output, "###############################################################################\n");
  fprintf(output, "#  Mutants generated from individual with rank %" PRId32
      " and index %" PRId32 " at timestep %" PRId64 " \n",
      wanted_rank, wanted_index, timestep);
  fprintf(output, "###############################################################################\n");
  fprintf(output, "#  1.  Parent id\n");
  fprintf(output, "#  2.  Parent metabolic error\n");
  fprintf(output, "#  3.  Parent secretion\n");
  fprintf(output, "#  4.  Mutant id\n");
  fprintf(output, "#  5.  Mutant metabolic error\n");
  fprintf(output, "#  6.  Mutant secretion\n");
  fprintf(output, "###############################################################################\n");*/

  wanted_indiv.Evaluate();
  wanted_indiv.compute_statistical_data();
  wanted_indiv.compute_non_coding();

  double reproduction_statistics[3];
  double offsprings_statistics[6];
  double th_fv = wanted_indiv.compute_theoritical_f_nu();

  wanted_indiv.compute_experimental_f_nu(
      nb_children,
      std::make_shared<JumpingMT>(time(nullptr)),
      output);

  fclose(output);
  delete exp_manager;
  return EXIT_SUCCESS;
}


/*!
  \brief

*/
void print_help(char* prog_path) {
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) {
    prog_name++;
  }
  else { prog_name = prog_path; }

  printf("\n");
  printf(
      "*********************** aevol - Artificial Evolution ******************* \n");
  printf(
      "*                                                                      * \n");
  printf(
      "*                     Robustness post-treatment program                * \n");
  printf(
      "*                                                                      * \n");
  printf(
      "************************************************************************ \n");
  printf("\n\n");
  printf("This program is Free Software. No Warranty.\n");
  printf("\n");
  printf("Usage : %s -h\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s -g GENER [-n NBCHILDREN] [-r RANK | -i INDEX]\n",
         prog_name);
  printf("\n");
  printf(
      "This program computes the replication statistics of all the individuals of a given generation,\n");
  printf(
      "like the proportion of neutral, beneficial and deleterious offsprings. This is done by simulating\n");
  printf(
      "NBCHILDREN replications for each individual, with its mutation, rearrangement and transfer rates.\n");
  printf(
      "Depending on those rates and genome size, there can be several events per replication.\n");
  printf(
      "These statistics are written in analysis-generation_TIMESTEP/robustness-allindivs.out, with one \n");
  printf("line per individual in the specified generation. \n\n");
  printf(
      "The program also outputs detailed statistics for one of the individuals (the best one by default). \n");
  printf(
      "The detailed statistics for this individual are written in the file called \n");
  printf(
      "analysis-generation_TIMESTEP/robustness-singleindiv-details-iINDEX-rRANK.out, with one line per simulated\n");
  printf("child of this particular individual.\n");
  printf("\n");
  printf("\n");
  printf("\t-h or --help    : Display this help, then exit\n");
  printf("\n");
  printf("\t-V or --version : Print version number, then exit\n");
  printf("\n");
  printf("\t-g GENER or --gener GENER : \n");
  printf(
      "\t                  Generation at which the statistics are computed\n");
  printf("\n");
  printf("\t-i INDEX or --index INDEX : \n");
  printf(
      "\t                  Index of individual of interest. Should be comprised between 0 and N-1, where\n");
  printf("\t                  N is the size of the population.\n");
  printf("\n");
  printf("\t-r RANK or --rank RANK : \n");
  printf(
      "\t                  Rank of individual of interest. Should be comprised between 1 and N, where\n");
  printf(
      "\t                  N is the size of the population. Default = N (fittest individual).\n");
  printf("\n");
  printf("\t-n NBCHILDREN or --nb-children NBCHILDREN : \n");
  printf(
      "\t                  Perform NBCHILDREN replications per individual to compute its statistics. \n");
  printf("\t                  Default = 1000.\n");
  printf("\n");


  printf("\n");
}
