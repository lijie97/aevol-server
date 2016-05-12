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
#include <libgen.h>

#include "aevol.h"

using namespace aevol;


enum check_type {
  FULL_CHECK  = 0,
  LIGHT_CHECK = 1,
  ENV_CHECK   = 2,
  NO_CHECK    = 3
};

void print_help(char* prog_name);
void analyze_indiv(Individual* initial_indiv,
                   FILE* output,
                   int32_t ndiv,
		   std::shared_ptr<JumpingMT> prng,
		   bool verbose);

int main(int argc, char* argv[]) {
  // Load parameters from command line
  int32_t nb_indiv = 1000;  // Default number of mutants per individual
  int32_t begin = 0;        // Default starting generation
  int32_t end = -1;         // Default ending generation (-1 for last generation stored in lineage file)
  int32_t period = 1;       // Period of analyze
  char*   output_file_name    = NULL;
  char*   lineage_file_name   = NULL;
  bool    verbose             = false;

  const char* short_options = "hVvn:b:e:p:f:o:l";
  static struct option long_options[] =
  {
    {"help",        no_argument,       NULL, 'h'},
    {"version",     no_argument,       NULL, 'V'},
    {"verbose",     no_argument,       NULL, 'v'},
    {"number",      required_argument, NULL, 'n'},
    {"begin",       required_argument, NULL, 'b'},
    {"end",         required_argument, NULL, 'e'},
    {"period",      required_argument, NULL, 'p'},
    {"file",        required_argument, NULL, 'f'},
    {"output",      required_argument, NULL, 'o'},
    {0, 0, 0, 0}
  };

  int option;
  while ((option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
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
        nb_indiv = atol(optarg);
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
      case 'f' :
        lineage_file_name = new char[strlen(optarg) + 1];
        sprintf(lineage_file_name, "%s", optarg);
        break;
    }
  }

  // =======================
  //  Open the lineage file
  // =======================
  gzFile lineage_file = gzopen(lineage_file_name, "r");
  if (lineage_file == Z_NULL) {
    fprintf(stderr, "ERROR : Could not read lineage file %s\n", lineage_file_name);
    exit(EXIT_FAILURE);
  }

  int64_t t0 = 0;
  int64_t t_end = 0;
  int32_t final_indiv_index = 0;
  int32_t final_indiv_rank  = 0;

  gzread(lineage_file, &t0, sizeof(t0));
  gzread(lineage_file, &t_end, sizeof(t_end));
  gzread(lineage_file, &final_indiv_index, sizeof(final_indiv_index));
  gzread(lineage_file, &final_indiv_rank,  sizeof(final_indiv_rank));

  if (verbose) {
    printf("\n\n");
    printf("===============================================================================\n");
    printf(" Robustness of the ancestors of indiv. %" PRId32
           " (rank %" PRId32 ") from time %" PRId64 " to %" PRId64 "\n",
           final_indiv_index, final_indiv_rank, t0, t_end);
    printf("================================================================================\n");
  }

  // =============================
  //  Open the experience manager
  // =============================
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(t0, true, false);

  // The current version doesn't allow for phenotypic variation nor for
  // different phenotypic targets among the grid
  if (not exp_manager->world()->phenotypic_target_shared())
    Utils::ExitWithUsrMsg("sorry, ancestor stats has not yet been implemented "
                          "for per grid-cell phenotypic target\n");
  auto phenotypicTargetHandler =
      exp_manager->world()->phenotypic_target_handler();
  if (phenotypicTargetHandler->var_method() != NO_VAR)
    Utils::ExitWithUsrMsg("sorry, ancestor stats has not yet been implemented "
                          "for variable phenotypic targets\n");

  int64_t backup_step = exp_manager->backup_step();

  // =========================
  //  Open the output file(s)
  // =========================

  // // Create missing directories
  // int status;
  // status = mkdir("stats/ancstats/", 0755);
  // if ((status == -1) && (errno != EEXIST))
  //   err(EXIT_FAILURE, "stats/ancstats/");

  FILE * output_summary = fopen(output_file_name, "w");

  std::shared_ptr<JumpingMT> prng = std::make_shared<JumpingMT>(9695);
    
  // ==============================
  //  Prepare the initial ancestor 
  // ==============================
  GridCell* grid_cell = new GridCell(lineage_file, exp_manager, nullptr);
  auto* indiv = grid_cell->individual();
  indiv->Evaluate();
  //  indiv->compute_statistical_data();
  //  indiv->compute_non_coding();

  // ==============================
  //  Compute robustness of the initial ancestor 
  // ==============================

  if (begin == 0)
    analyze_indiv(indiv,output_summary,nb_indiv,prng,verbose);

  // ==========================================================================
  //  Replay the mutations to get the successive ancestors and analyze them
  // ==========================================================================
  ReplicationReport* rep = nullptr;

  int32_t index;

  ExpManager* exp_manager_backup = nullptr;
  Habitat *backup_habitat = nullptr;

  aevol::AeTime::plusplus();
  while ((time() <= t_end) && (((time() < end)||(end == -1))))
  {
    rep = new ReplicationReport(lineage_file, indiv);
    index = rep->id(); // who we are building...
    indiv->Reevaluate();

    if (verbose)
      printf("Ancestor at generation %" PRId64
             " has index %" PRId32 "\n", time(), index);
    
    
    // 2) Replay replication (create current individual's child)
    GeneticUnit& gen_unit = indiv->genetic_unit_nonconst(0);
    GeneticUnit* stored_gen_unit = nullptr;
    Individual* stored_indiv = nullptr;


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
    indiv->Reevaluate();

    // if we are between "begin" and "end" and at the correct period, compute robustness
    
    if ((time() >= begin)&&((time() < end)||(end == -1)))
      {
	if (((time() - begin) % period) == 0)
	  {
	    analyze_indiv(indiv,output_summary,nb_indiv,prng,verbose);
	  }
      }
    delete rep;

    aevol::AeTime::plusplus();
  }

  gzclose(lineage_file);
  fclose(output_summary);
  delete exp_manager;
  // delete mystats;
  delete indiv;
  return EXIT_SUCCESS;
}

void analyze_indiv(Individual* indiv,
                   FILE* output_summary,
                   int32_t nb_indiv,
		   std::shared_ptr<JumpingMT> prng,
		   bool verbose)
{
  double fracplus = 0;
  double meanplus = 0;
  double fracminus = 0;
  double meanminus = 0;
  double maxplus = 0;
  double maxminus = 0;
  double fracnull = 0;
  int32_t nb_events = 0;
  
  double indiv_metabolic_error = indiv->dist_to_target_by_feature(METABOLISM);    const Habitat& habitat = indiv->habitat();
  
  for (int32_t i = 0; i < nb_indiv; i++)
    {
      Individual* new_indiv = new Individual(indiv, 0, prng, prng);
      // Perform transfer, rearrangements and mutations
      if (not new_indiv->allow_plasmids()) {
	const GeneticUnit* chromosome = &new_indiv->genetic_unit_list().front();
	nb_events = chromosome->dna()->perform_mutations(indiv->id());
      }
      else {
	printf("WARNING: Mutational Robustness does not handle multiple Genetic Units\n");
      }
	if (nb_events == 0)
	  {
	    fracnull++;
	  }
	else
	  {
	    new_indiv->EvaluateInContext(habitat);
	    double new_metabolic_error    = new_indiv->dist_to_target_by_feature(METABOLISM);
	    
	    if (new_metabolic_error == indiv_metabolic_error)
	      fracnull++;
	    if (new_metabolic_error  > indiv_metabolic_error)
	      {
		fracminus++;
		if ((new_metabolic_error - indiv_metabolic_error) > maxminus) maxminus = new_metabolic_error - indiv_metabolic_error;
		meanminus += new_metabolic_error - indiv_metabolic_error;
	      }
	    if (new_metabolic_error  < indiv_metabolic_error)
	      {
		fracplus++;
			if ((new_metabolic_error - indiv_metabolic_error) < maxplus) maxplus = new_metabolic_error - indiv_metabolic_error;
		meanplus += new_metabolic_error - indiv_metabolic_error;
	      }
	  }
	
	delete new_indiv;
    }
  
  
  if (verbose)  printf("f+: %f   f0: %f    f-:%f\n",fracplus,fracnull,fracminus);
  
   if ( output_summary == NULL ){
    fprintf( stderr, "ERROR : Could not create robustness_summary.txt\n");
    exit( EXIT_FAILURE );
  }
  fprintf(output_summary, "%lld %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n",time(),fracplus/nb_indiv,fracnull/nb_indiv,fracminus/nb_indiv, meanplus/fracplus,meanminus/fracminus,maxplus,maxminus);
}


// Print help
void print_help(char* prog_name) {
  printf("\n\
%s is a post-treatment that generates and analyses a large quantity of mutants for a lineage of ancestors.\
For each mutant we record the phenotypic effect on metabolism.\n\n\
Usage: %s [-h] -i input_file_name -o output_file_name [-b start_at_generation] [-e end_at_generation] [-p period] [-n num_mutants] [-r] [-h bin_size] [-v verbose] [-s mutation_seed]\n\
\t-h: display this screen\n\
\t-f input_file_name: lineage file to be analyzed\n\
\t-o output_file_name: name of the output file (to be written in ./stats/ancstats). In case of histogram output (-h) one file will be produced for each histogram and output_file_name will be postfixed with the generation number\n\
\t-b start_at_generation: first generation of the lineage to be analyzed (default: 0)\n\
\t-e end_at_generation: last generation of the lineage to be analyzed (default: last generation stored in the input file)\n\
\t-p period: temporal resolution of the analyze (default: 1)\n\
\t-n nb_mutants : generate and analyse nb_mutants per individual (default: 1000)\n\
Example:\n\t%s -i lineage_file -o toto.out -b 4000 -e 5000 -p 10 -n 100000 -s 19769\n", prog_name, prog_name, prog_name);

// \t-r: raw output; store the difference of metabolic error for each mutant generated (warning: the output file may quickly grow)\n
// \t-h bin_size: store the histogram with a bin_size resolution. One output file is generated for each histogram (postfixed with the generation number)\n
// \t-s mutation_seed: specify the seed to be used for the mutation random generator\n\n

}
