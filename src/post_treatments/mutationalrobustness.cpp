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

using namespace aevol;

void print_help(char* prog_name);
void analyse_indiv(ExpManager* exp,
                   Individual* initial_indiv,
                   FILE* output,
                   int32_t ndiv);

int main(int argc, char* argv[]) {
  // Load parameters from command line
  int32_t ndiv = 100000; // Default number of mutants per individual
  int32_t gener = -1; // What generation to load
  char* output_file_name = NULL;
  bool best_only = false; // Treat only the best individual?

  // TODO <david.parsons@inria.fr> version
  const char * options_list = "hn:r:o:b";
  static struct option long_options_list[] = {
    { "help", 0, NULL, 'h' },
    { "number", 1, NULL, 'n' },
    { "generation", 1, NULL, 'r' },
    { "output", 1, NULL, 'o' },
    { "best", 0, NULL, 'b' },
    { 0, 0, 0, 0 }
  };

  int option;
  while ((option = getopt_long(argc, argv, options_list, long_options_list, NULL)) != -1)
  {
    switch (option) {
      case 'h' : print_help(basename(argv[0])); exit(EXIT_SUCCESS);
      case 'n' : ndiv = atol(optarg); break;
      case 'r' : gener = atol(optarg); break;
      case 'o' : {
        output_file_name = new char[strlen(optarg) + 1];
        sprintf(output_file_name, "%s", optarg);
        break;
      }
      case 'b' : best_only = true; break;
    }
  }

  // Load the population from the backup file
  if (gener == -1) {
    printf("You must specify a generation number. Please use the option -r or --generation.\n");
    exit(EXIT_FAILURE);
  }

  // Open output file and write the header
  FILE * output = fopen(output_file_name, "w");
  if (output == NULL) {
    fprintf(stderr, "ERROR : Could not create the output file %s\n", output_file_name);
    exit(EXIT_FAILURE);
  }

  // Load experiment
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(gener, true, false);

  // Positive impact means
  fprintf(output,
      "# #################################################################\n"
      "#              Mutations produced by mutationalrobustness\n"
      "# #################################################################\n"
      "# Number of replicate per individual : %" PRId32 "\n"
      "# Impact on metabolism SPACE impact on secretion\n"
      "#\n", ndiv);

  // Parse and treat the individuals
  if (!best_only) {
    for (Individual* indiv: exp_manager->world()->indivs()) {
      analyse_indiv(exp_manager, indiv, output, ndiv);
    }
  }
  else {
    Individual* indiv = exp_manager->world()->best_indiv();
    analyse_indiv(exp_manager, indiv, output, ndiv);
  }

  // Clean memory and exit
  delete exp_manager;
  delete [] output_file_name;
  return EXIT_SUCCESS;
}

// Treatment of one individual
void analyse_indiv(ExpManager* exp,
                   Individual* initial_indiv,
                   FILE* output,
                   int32_t ndiv) {
  double initial_metabolic_error = initial_indiv->
      dist_to_target_by_feature(METABOLISM);
  double initial_secretion_error = initial_indiv->
      dist_to_target_by_feature(SECRETION);
  double final_metabolic_error      = 0.0;
  double impact_on_metabolic_error  = 0.0;
  double final_secretion_error      = 0.0;
  double impact_on_secretion_error  = 0.0;

  const Habitat& habitat = initial_indiv->habitat();

  Individual* indiv = nullptr;
  int32_t i;

  // Perform ndiv reproductions with mutations
  for (i = 0; i < ndiv; i++) {
    if (i % 1000 == 0) {
      printf("*");
      fflush(stdout);
    }

    indiv = exp->sel()->do_replication(initial_indiv);

    indiv->EvaluateInContext(habitat);
    final_metabolic_error     = indiv->dist_to_target_by_feature(METABOLISM);
    impact_on_metabolic_error = final_metabolic_error - initial_metabolic_error;
    final_secretion_error     = indiv->dist_to_target_by_feature(SECRETION);
    impact_on_secretion_error = final_secretion_error - initial_secretion_error;

    fprintf(output, "%+.15f %+.15f \n",
            impact_on_metabolic_error,
            impact_on_secretion_error);

    delete indiv;
  }
  fprintf(output, "\n");
}

// Print help
void print_help(char* prog_name) {
  printf("\n\
%s is a post-treatment that generates and analyses a large quantity of mutants for all individuals in a backup.\
For each mutant we record the phenotypic effect on metabolism and on secretion.\n\n\
Usage: %s [-h] -r num_generation -o output_file_name -n num_mutants [-b] \n\
\t-h : display this screen\n\
\t-r num_generation  : read the generation num_generation from a full aevol backup\n\
\t-o output_file_name : write the results in file output_file_name\n\
\t-n nb_mutants : generate and analyse nb_mutants per individual\n\
\t-b : only treat the best individual\n\n\
Example:\n\t%s -r 20000 -n 1000 -o toto.out\n",prog_name,prog_name,prog_name);
}
