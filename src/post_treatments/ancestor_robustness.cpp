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

/*
 This post-treatment loads a population and create a large number of mutants of
 each individual in the same way they would be created during reproduction of
 the individual. For each mutant we report the phenotypic effect of the mutation
 on metabolism and on secretion. This allows to infer the mutational landscape
 of a population.
*/

#include <list>
#include <getopt.h>

#include "aevol.h"
#include "IndivAnalysis.h"

using namespace aevol;

void print_help(char* prog_name);

int main(int argc, char* argv[]) {
  // Load parameters from command line
  int32_t nb_mutants = 1000;  // Default number of mutants per individual
  int32_t begin = 0;        // Default starting generation
  int32_t end = -1;         // Default ending generation (-1 for last generation stored in lineage file)
  int32_t period = 1;       // Period of analyze
  char* output_file_name = "robustness_summary.txt";
  char* lineage_file_name = nullptr;
  bool verbose = false;

  const char* short_options = "hVvn:b:e:p:o:l";
  static struct option long_options[] =
      {
          {"help",    no_argument,       nullptr, 'h'},
          {"version", no_argument,       nullptr, 'V'},
          {"verbose", no_argument,       nullptr, 'v'},
          {"number",  required_argument, nullptr, 'n'},
          {"begin",   required_argument, nullptr, 'b'},
          {"end",     required_argument, nullptr, 'e'},
          {"period",  required_argument, nullptr, 'p'},
          {"output",  required_argument, nullptr, 'o'},
          {0, 0,                         0,       0}
      };

  int option;
  while ((option = getopt_long(argc, argv, short_options, long_options,
                               nullptr)) != -1) {
    switch (option) {
      case 'h' :
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      case 'V' :
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      case 'v' :
        verbose = true;
        break;
      case 'n' :
        nb_mutants = atol(optarg);
        break;
      case 'b' :
        begin = atol(optarg);
        break;
      case 'e' :
        end = atol(optarg);
        break;
      case 'p' :
        period = atol(optarg);
        break;
      case 'o' :
        output_file_name = new char[strlen(optarg) + 1];
        sprintf(output_file_name, "%s", optarg);
        break;
    }
  }

  // There should be only one remaining arg: the lineage file
  if (optind != argc - 1) {
    Utils::ExitWithUsrMsg("please specify a lineage file");
  }

  lineage_file_name = new char[strlen(argv[optind]) + 1];
  sprintf(lineage_file_name, "%s", argv[optind]);


  // =======================
  //  Open the lineage file
  // =======================
  gzFile lineage_file = gzopen(lineage_file_name, "r");
  if (lineage_file == Z_NULL) {
    Utils::ExitWithUsrMsg(std::string("Could not read lineage file ") +
                          lineage_file_name + "\n");
  }

  int64_t t0 = 0;
  int64_t t_end = 0;
  int32_t final_indiv_index = 0;
  int32_t final_indiv_rank = 0;

  gzread(lineage_file, &t0, sizeof(t0));
  gzread(lineage_file, &t_end, sizeof(t_end));
  gzread(lineage_file, &final_indiv_index, sizeof(final_indiv_index));
  gzread(lineage_file, &final_indiv_rank, sizeof(final_indiv_rank));

  if (verbose) {
    printf("\n\n");
    printf(
        "===============================================================================\n");
    printf(" Robustness of the ancestors of indiv. %" PRId32
               " (rank %" PRId32 ") from time %" PRId64 " to %" PRId64 "\n",
           final_indiv_index, final_indiv_rank, t0, t_end);
    printf(
        "================================================================================\n");
  }

  // =============================
  //  Open the experience manager
  // =============================
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(t0, true, false);

  // The current version doesn't allow for phenotypic variation nor for
  // different phenotypic targets among the grid
  if (not exp_manager->world()->phenotypic_target_shared()) {
    Utils::ExitWithUsrMsg("sorry, ancestor stats has not yet been implemented "
                              "for per grid-cell phenotypic target\n");
  }
  auto phenotypicTargetHandler =
      exp_manager->world()->phenotypic_target_handler();
  if (phenotypicTargetHandler->var_method() != NO_VAR) {
    Utils::ExitWithUsrMsg("sorry, ancestor stats has not yet been implemented "
                              "for variable phenotypic targets\n");
  }

  int64_t backup_step = exp_manager->backup_step();

  // =========================
  //  Open the output file(s)
  // =========================

  // // Create missing directories
  // int status;
  // status = mkdir("stats/ancestor_stats/", 0755);
  // if ((status == -1) && (errno != EEXIST))
  //   err(EXIT_FAILURE, "stats/ancestor_stats/");

  FILE* output_summary = fopen(output_file_name, "w");
  if (output_summary == nullptr) {
    Utils::ExitWithUsrMsg(std::string("Could not create ") + output_file_name);
  }

  std::shared_ptr<JumpingMT> prng = std::make_shared<JumpingMT>(9695);

  // ==============================
  //  Prepare the initial ancestor 
  // ==============================
  GridCell* grid_cell = new GridCell(lineage_file, exp_manager, nullptr);
  IndivAnalysis indiv(*(grid_cell->individual()));
  indiv.Evaluate();
  //  indiv->compute_statistical_data();
  //  indiv->compute_non_coding();

  // ==============================
  //  Compute robustness of the initial ancestor 
  // ==============================

  if (begin == 0) {
    indiv.compute_experimental_f_nu(nb_mutants, prng, output_summary, nullptr,
                                    verbose);
  }

  // ==========================================================================
  //  Replay the mutations to get the successive ancestors and analyze them
  // ==========================================================================
  ReplicationReport* rep = nullptr;

  int32_t index;

  aevol::AeTime::plusplus();
  while ((time() <= t_end) && (((time() < end) || (end == -1)))) {
    rep = new ReplicationReport(lineage_file, &indiv);
    index = rep->id(); // who we are building...
    indiv.Reevaluate();

    if (verbose) {
      printf("Ancestor at generation %" PRId64
                 " has index %" PRId32 "\n", time(), index);
    }


    // 2) Replay replication (create current individual's child)
    GeneticUnit& gen_unit = indiv.genetic_unit_nonconst(0);


    // For each genetic unit, replay the replication (undergo all mutations)
    // TODO <david.parsons@inria.fr> disabled for multiple GUs
    const auto& dnarep = rep->dna_replic_report();

    for (const auto& mut: dnarep.HT())
      gen_unit.dna()->undergo_this_mutation(*mut);
    for (const auto& mut: dnarep.rearrangements())
      gen_unit.dna()->undergo_this_mutation(*mut);
    for (const auto& mut: dnarep.mutations())
      gen_unit.dna()->undergo_this_mutation(*mut);

    // 3) All the mutations have been replayed, we can now evaluate the new individual
    indiv.Reevaluate();

    // if we are between "begin" and "end" and at the correct period, compute robustness

    if ((time() >= begin) && ((time() < end) || (end == -1)) &&
        (((time() - begin) % period) == 0)) {
      indiv.compute_experimental_f_nu(nb_mutants, prng, output_summary, nullptr,
                                      verbose);
    }
    delete rep;

    aevol::AeTime::plusplus();
  }

  gzclose(lineage_file);
  fclose(output_summary);
  delete exp_manager;
  return EXIT_SUCCESS;
}


// Print help
void print_help(char* prog_name) {
  printf("\n\
%s is a post-treatment that generates and analyses a large quantity of mutants for a lineage of ancestors.\
For each mutant we record the phenotypic effect on metabolism.\n\n\
Usage: %s [-h] -i input_file_name -o output_file_name [-b start_at_generation] [-e end_at_generation] [-p period] [-n num_mutants] [-r] [-h bin_size] [-v verbose] [-s mutation_seed]\n\
\t-h: display this screen\n\
\t-f input_file_name: lineage file to be analyzed\n\
\t-o output_file_name: name of the output file (to be written in ./stats/ancestor_stats). In case of histogram output (-h) one file will be produced for each histogram and output_file_name will be postfixed with the generation number\n\
\t-b start_at_generation: first generation of the lineage to be analyzed (default: 0)\n\
\t-e end_at_generation: last generation of the lineage to be analyzed (default: last generation stored in the input file)\n\
\t-p period: temporal resolution of the analyze (default: 1)\n\
\t-n nb_mutants : generate and analyse nb_mutants per individual (default: 1000)\n\
Example:\n\t%s -i lineage_file -o toto.out -b 4000 -e 5000 -p 10 -n 100000 -s 19769\n",
         prog_name, prog_name, prog_name);

// \t-r: raw output; store the difference of metabolic error for each mutant generated (warning: the output file may quickly grow)\n
// \t-h bin_size: store the histogram with a bin_size resolution. One output file is generated for each histogram (postfixed with the generation number)\n
// \t-s mutation_seed: specify the seed to be used for the mutation random generator\n\n

}
