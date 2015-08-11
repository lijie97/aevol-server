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

 
const char* DEFAULT_PARAM_FILE_NAME = "param.in";
 
 
// =================================================================
//                              Includes
// =================================================================
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <getopt.h>

#include <list>

#if __cplusplus == 201103L
#include "make_unique.h"
#endif

#include "ParameterLine.h"
#ifdef __X11
#include "ExpManager_X11.h"
#else
#include "ExpManager.h"
#include "libaevol/ParameterLine.h"
#include "libaevol/JumpingMT.h"

#endif

using namespace aevol;

// =================================================================
//                         Function declarations
// =================================================================
enum population_change_type
{
  SUBPOPULATIONS_BASED_ON_NON_CODING_BASES = 3,
  REMOVE_NON_CODING_BASES_BEST_IND = 4,
  REMOVE_NON_CODING_BASES_POPULATION = 5,
  DOUBLE_NON_CODING_BASES_BEST_IND = 6,
  DOUBLE_NON_CODING_BASES_POPULATION = 7
};

void print_help(char* prog_path);

ParameterLine* get_line(FILE* param_file);
void format_line(ParameterLine* formated_line, char* line, bool* line_is_interpretable);
// void change_by_cloning_best(ae_population* pop, ae_exp_manager* exp_m);
// void change_based_on_non_coding_bases_of_best_individual(ae_population* pop, ae_exp_manager* exp_m, population_change_type type);
// void change_based_on_non_coding_bases_in_population(ae_population* pop, ae_exp_manager* exp_m, population_change_type type);



int main(int argc, char* argv[])
{
  // 1) Initialize command-line option variables with default values
  char* param_file_name = NULL;
  bool verbose          = false;
  int64_t num_gener = -1;
  
  // 2) Define allowed options
  const char * options_list = "hf:g:V";
  static struct option long_options_list[] = {
    { "help",     no_argument,        NULL, 'h' },
    { "file",     required_argument,  NULL, 'f' }, // Provide file with parameters to change
    { "gener",    required_argument,  NULL, 'g' },
    { "version",  no_argument,        NULL, 'V' },
    { 0, 0, 0, 0 }
  };
      
  // 3) Get actual values of the command-line options
  int option;
  while ((option = getopt_long(argc, argv, options_list, long_options_list, NULL)) != -1) 
    {
      switch (option) 
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
        case 'f' :
          {
            if (strcmp(optarg, "") == 0)
              {
                printf("%s: error: Option -f or --file : missing argument.\n", argv[0]);
                exit(EXIT_FAILURE);
              }
        
            param_file_name = optarg;
            break;
          }
        case 'g' :
          {
            num_gener = atoi(optarg);
            break;
          }
        default :
          {
            // An error message is printed in getopt_long, we just need to exit
            exit(EXIT_FAILURE);
          }
        }
    }  
  
  // 4) Set undefined command line parameters to default values
  if (param_file_name == NULL) {
    param_file_name = new char[strlen(DEFAULT_PARAM_FILE_NAME)+1];
    sprintf(param_file_name, "%s", DEFAULT_PARAM_FILE_NAME);
  }
  if (num_gener == -1) {
    // Set num_gener to the content of the LAST_GENER file if it exists.
    // If it doesn't, print help and exit
    FILE* lg_file = fopen(LAST_GENER_FNAME, "r");
    if (lg_file != NULL) {
      if (fscanf(lg_file, "%" PRId64 "\n", &num_gener) == EOF) {
        Utils::ExitWithUsrMsg(
            std::string("failed to read last generation from file ") +
            LAST_GENER_FNAME);
        exit(EXIT_FAILURE);
      }
      fclose(lg_file);
    }
    else
      Utils::ExitWithUsrMsg("you must provide a generation number");
  }
  
  // 5) Check the consistancy of the command-line options
  if (num_gener == -1)
    {
      printf("%s: error: You must provide a generation number.\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  
  // 6) Initialize the experiment manager
#ifndef __NO_X
  ExpManager* exp_manager = new ExpManager_X11();
#else
  ExpManager* exp_manager = new ExpManager();
#endif
  exp_manager->load(num_gener, false, verbose);

  // 7) Define syntaxic sugars for the population, the environment, the selection...  
//  Environment* env = exp_manager->get_env();
  Selection* sel = exp_manager->get_sel();
  World* world = exp_manager->world();


  // If relevant, load the tree information 
  char tree_file_name[50];
  Tree* tree = nullptr;
  bool take_care_of_the_tree = exp_manager->get_record_tree() &&
                               get_time() > 0;

  if (take_care_of_the_tree)
  {
    // If a tree is available, assign the replication reports to the individuals
    #ifdef __REGUL
      sprintf(tree_file_name,"tree/tree_%06" PRId64 ".rae", num_gener);
    #else
      sprintf(tree_file_name,"tree/tree_%06" PRId64 ".ae", num_gener);
    #endif
    
    tree = new Tree(exp_manager, tree_file_name);
  }

 




  // 8) Interpret and apply changes
  printf("Interpret and apply changes\n");
  FILE* param_file  = fopen(param_file_name,  "r");
  if (param_file == NULL)
  {
    printf("%s:%d: error: could not open parameter file %s\n", __FILE__, __LINE__, param_file_name);
    exit(EXIT_FAILURE);
  }

  std::list<Gaussian> new_gaussians;
  bool phen_target_change = false;
  bool start_to_record_tree = false;
  bool set_tree_step = false;
  int32_t tree_step = 100;
  
  ParameterLine* line;
  int32_t cur_line = 0;
  while ((line = get_line(param_file)) != NULL) 
  {
    cur_line++;
    if (strcmp(line->words[0], "ENV_AXIS_FEATURES") == 0)
    {
      // TODO <david.parsons@inria.fr> adapt to new organization
      printf("%s:%d: error: ENV_AXIS_FEATURES has to be adapted to the new organization.\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
//      int16_t env_axis_nb_segments = line->nb_words / 2;
//      double* env_axis_segment_boundaries = new double [env_axis_nb_segments + 1];
//      env_axis_segment_boundaries[0] = X_MIN;
//      for (int16_t i = 1 ; i < env_axis_nb_segments ; i++)
//      {
//        env_axis_segment_boundaries[i] = atof(line->words[2*i]);
//      }
//      env_axis_segment_boundaries[env_axis_nb_segments] = X_MAX;
//
//      // Set segment features
//      PhenotypicFeature* env_axis_features = new PhenotypicFeature[env_axis_nb_segments];
//      for (int16_t i = 0 ; i < env_axis_nb_segments ; i++)
//      {
//        if (strcmp(line->words[2*i+1], "NEUTRAL") == 0)
//        {
//          env_axis_features[i] = NEUTRAL;
//        }
//        else if (strcmp(line->words[2*i+1], "METABOLISM") == 0)
//        {
//          env_axis_features[i] = METABOLISM;
//        }
//        else if (strcmp(line->words[2*i+1], "SECRETION") == 0)
//        {
//          exp_manager->get_exp_s()->set_with_secretion(true);
//          env_axis_features[i] = SECRETION;
//        }
//        else if (strcmp(line->words[2*i+1], "DONOR") == 0)
//        {
//          env_axis_features[i] = DONOR;
//        }
//        else if (strcmp(line->words[2*i+1], "RECIPIENT") == 0)
//        {
//          env_axis_features[i] = RECIPIENT;
//        }
//        else
//        {
//          printf("ERROR in param file \"%s\" on line %" PRId32 " : unknown axis feature \"%s\".\n",
//                  param_file_name, cur_line, line->words[2*i+1]);
//          exit(EXIT_FAILURE);
//        }
//      }
//      env->set_segmentation(env_axis_nb_segments,
//                             env_axis_segment_boundaries,
//                             env_axis_features);
//      env_hasbeenmodified = true;
//      delete env_axis_segment_boundaries;
//      delete env_axis_features;
    }
    else if (strcmp(line->words[0], "RECORD_TREE") == 0)
    {
      if (strncmp(line->words[1], "true", 4) == 0)
      {
        start_to_record_tree = true;
      }
      else if (strncmp(line->words[1], "false", 5) == 0)
      {
        printf("ERROR stop recording tree is not implemented yet.\n");
        exit(EXIT_FAILURE);
      }
      else
      {
        printf("ERROR in param file \"%s\" on line %" PRId32" : unknown tree recording option (use true/false).\n",
               param_file_name, cur_line);
        exit(EXIT_FAILURE);
      }
      if (exp_manager->get_output_m()->get_record_tree())
      {
        printf("ERROR modification of already existing tree not impemented yet\n");
        exit(EXIT_FAILURE);
      }
    }
    else if (strcmp(line->words[0], "TREE_STEP") == 0)
    {
      tree_step = atol(line->words[1]);
      set_tree_step = true;
    }
    else if (strcmp(line->words[0], "TREE_MODE") == 0)
    {
      printf("ERROR : Tree mode management has been removed.\n");
      exit(EXIT_FAILURE);
    }
    else if (strcmp(line->words[0], "DUMP_STEP") == 0)
    {
      int step = atoi(line->words[1]);
      if (step>0)
      {
        exp_manager->get_output_m()->set_dump_step(step);
      }
    }
    else if (strcmp(line->words[0], "BACKUP_STEP") == 0)
    {
      exp_manager->get_output_m()->set_backup_step(atol(line->words[1]));
    }
    else if (strcmp(line->words[0], "BIG_BACKUP_STEP") == 0)
    {
      exp_manager->get_output_m()->set_big_backup_step(atol(line->words[1]));
    }
    else if (strcmp(line->words[0], "POPULATION_SIZE") == 0)
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
          ": the change of population size is not implemented yet\n for spatially structured populations",
              param_file_name, cur_line);
      exit(EXIT_FAILURE);
    }
    else if (strcmp(line->words[0], "SELECTION_SCHEME") == 0)
    {
      if (strncmp(line->words[1], "lin", 3) == 0)
      {
        if (line->nb_words != 3)
        {
          printf("ERROR in param file \"%s\" on line %" PRId32 " : selection pressure parameter is missing.\n",
                  param_file_name, cur_line);
          exit(EXIT_FAILURE);
        }
        sel->set_selection_scheme(RANK_LINEAR);
        sel->set_selection_pressure(atof(line->words[2]));
      }
      else if (strncmp(line->words[1], "exp", 3) == 0)
      {
        if (line->nb_words != 3)
        {
          printf("ERROR in param file \"%s\" on line %" PRId32 " : selection pressure parameter is missing.\n",
                  param_file_name, cur_line);
          exit(EXIT_FAILURE);
        }
        sel->set_selection_scheme(RANK_EXPONENTIAL);
        sel->set_selection_pressure(atof(line->words[2]));
      }
      else if (strncmp(line->words[1], "fitness", 7) == 0)
      {
        if (line->nb_words != 3)
        {
          printf("ERROR in param file \"%s\" on line %" PRId32 " : selection pressure parameter is missing.\n",
                  param_file_name, cur_line);
          exit(EXIT_FAILURE);
        }
        sel->set_selection_scheme(FITNESS_PROPORTIONATE);
        sel->set_selection_pressure(atof(line->words[2]));
      }
      else if (strcmp(line->words[1], "fittest") == 0)
      {
        sel->set_selection_scheme(FITTEST);
      }
      else
      {
        printf("ERROR in param file \"%s\" on line %" PRId32" : unknown selection scheme \"%s\".\n",
            param_file_name, cur_line, line->words[1]);
        exit(EXIT_FAILURE);
      }
    }
    else if (strcmp(line->words[0], "SELECTION_PRESSURE") == 0)
    {
      printf("WARNING: SELECTION_PRESSURE keyword is outdated, you should specify a value for selection pressure using SELECTION_SCHEME\n");
      sel->set_selection_pressure(atof(line->words[1]));
      printf("\tChange of selection pressure to %f\n",atof(line->words[1]));
    }
    else if (strcmp(line->words[0], "POINT_MUTATION_RATE") == 0)
    {
      double point_mutation_rate = atof(line->words[1]);
      for (auto& indiv: exp_manager->get_indivs())
        indiv->set_point_mutation_rate(point_mutation_rate);
      printf("\tChange of overall point mutation rate to %f\n", point_mutation_rate);
    }
    else if (strcmp(line->words[0], "SMALL_INSERTION_RATE") == 0)
    {
      double small_insertion_rate = atof(line->words[1]);
      for (auto& indiv: exp_manager->get_indivs())
        indiv->set_small_insertion_rate(small_insertion_rate);
      printf("\tChange of overall small insertion rate to %f\n", small_insertion_rate);
    }
    else if (strcmp(line->words[0], "SMALL_DELETION_RATE") == 0)
    {
      double small_deletion_rate = atof(line->words[1]);
      for (auto& indiv: exp_manager->get_indivs())
        indiv->set_small_deletion_rate(small_deletion_rate);
      printf("\tChange of overall small deletion rate to %f\n", small_deletion_rate);
    }
    else if (strcmp(line->words[0], "MAX_INDEL_SIZE") == 0)
    {
      int16_t max_indel_size = atol(line->words[1]);
      for (auto& indiv: exp_manager->get_indivs())
        indiv->set_max_indel_size(max_indel_size);
      printf("\tChange of overall maximum indel size to %d\n", max_indel_size);
    }
    else if (strcmp(line->words[0], "DUPLICATION_RATE") == 0)
    {
      double duplication_rate = atof(line->words[1]);
      for (auto& indiv: exp_manager->get_indivs())
        indiv->set_duplication_rate(duplication_rate);
      printf("\tChange of overall duplication rate to %f\n", duplication_rate);
    }
    else if (strcmp(line->words[0], "DELETION_RATE") == 0)
    {
      double deletion_rate = atof(line->words[1]);
      for (auto& indiv: exp_manager->get_indivs())
        indiv->set_deletion_rate(deletion_rate);
      printf("\tChange of overall deletion rate to %f\n", deletion_rate);
    }
    else if (strcmp(line->words[0], "TRANSLOCATION_RATE") == 0)
    {
      double translocation_rate = atof(line->words[1]);
      for (auto& indiv: exp_manager->get_indivs())
        indiv->set_translocation_rate(translocation_rate);
      printf("\tChange of overall translocation rate to %f\n", translocation_rate);
    }
    else if (strcmp(line->words[0], "INVERSION_RATE") == 0)
    {
      double inversion_rate = atof(line->words[1]);
      for (auto& indiv: exp_manager->get_indivs())
        indiv->set_inversion_rate(inversion_rate);
      printf("\tChange of overall inversion to %f\n", inversion_rate);
    }
    else if (strcmp(line->words[0], "TRANSFER_INS_RATE") == 0)
    {
      double transfer_ins_rate = atof(line->words[1]);
      for (auto& indiv: exp_manager->get_indivs())
        indiv->set_HT_ins_rate(transfer_ins_rate);
      exp_manager->set_HT_ins_rate(transfer_ins_rate);
      printf("\tChange of overall transfer insertion rate to %f\n", transfer_ins_rate);
    }
    else if (strcmp(line->words[0], "TRANSFER_REPL_RATE") == 0)
    {
      double transfer_repl_rate = atof(line->words[1]);
      for (auto& indiv: exp_manager->get_indivs())
        indiv->set_HT_repl_rate(transfer_repl_rate);
      exp_manager->set_HT_repl_rate(transfer_repl_rate);
      printf("\tChange of overall transfer replacement rate to %f\n", transfer_repl_rate);
    }
    else if ((strcmp(line->words[0], "ENV_ADD_GAUSSIAN") == 0) || (strcmp(line->words[0], "ENV_GAUSSIAN") == 0))
    {
      // TODO <david.parsons@inria.fr> adapt to new organization
//      printf("%s:%d: error: ENV_ADD_GAUSSIAN has to be adapted to the new organization.\n", __FILE__, __LINE__);
//      exit(EXIT_FAILURE);

      new_gaussians.emplace_back(atof(line->words[1]),
                                 atof(line->words[2]),
                                 atof(line->words[3]));
      printf("\tAdding a gaussian with %f, %f, %f \n",
             atof(line->words[1]), atof(line->words[2]), atof(line->words[3]));
      phen_target_change = true;
    }
    else if (strcmp(line->words[0], "ENV_ADD_POINT") == 0) 
    {
      // custom_points
      printf("%s:%d: error: Custom points management has been removed.\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    else if (strcmp(line->words[0], "ENV_VARIATION") == 0)
    {
      // TODO <david.parsons@inria.fr> adapt to new organization
//      printf("%s:%d: error: ENV_VARIATION has to be adapted to the new organization.\n", __FILE__, __LINE__);
//      exit(EXIT_FAILURE);

      static bool env_var_already_set = false;
      if (env_var_already_set)
      {
        printf("%s:%d: ERROR in param file : duplicate entry for %s.\n", __FILE__, __LINE__, line->words[0]);
        exit(EXIT_FAILURE);
      }
      env_var_already_set = true;

      if (strcmp(line->words[1], "none") == 0) {
        assert(line->nb_words == 2);
          exp_manager->world()->phenotypic_target_handler()->set_var_method(NO_VAR);
        printf("\tNo more environmental variation\n");
      }
      else if (strcmp(line->words[1], "autoregressive_mean_variation") == 0) {
        assert(line->nb_words == 5);
        auto pt_handler = exp_manager->world()->phenotypic_target_handler();
        pt_handler->set_var_method(AUTOREGRESSIVE_MEAN_VAR);
        pt_handler->set_var_sigma_tau(atof(line->words[2]),
                                      atol(line->words[3]));
        pt_handler->set_var_prng(
            std::make_shared<JumpingMT>(atoi(line->words[4])));
        printf("\tChange of environmental variation to a autoregressive mean variation with sigma=%f, tau=%ld and seed=%d\n", atof(line->words[2]),atol(line->words[3]),atoi(line->words[4]));
      }
      else if (strcmp(line->words[1], "autoregressive_height_variation") == 0) {
        assert(line->nb_words == 5);
        auto pt_handler = exp_manager->world()->phenotypic_target_handler();
        pt_handler->set_var_method(AUTOREGRESSIVE_HEIGHT_VAR);
        pt_handler->set_var_sigma_tau(atof(line->words[2]),
                                      atol(line->words[3]));
        pt_handler->set_var_prng(
            std::make_shared<JumpingMT>(atoi(line->words[4])));
        printf("\tChange of environmental variation to a autoregressive height variation with sigma=%f, tau=%ld and seed=%d\n", atof(line->words[2]),atol(line->words[3]),atoi(line->words[4]));
      }
      else if (strcmp(line->words[1], "add_local_gaussians") == 0) {
        assert(line->nb_words == 3);
        auto pt_handler = exp_manager->world()->phenotypic_target_handler();
        pt_handler->set_var_method(LOCAL_GAUSSIANS_VAR);
        pt_handler->set_var_prng(std::make_shared<JumpingMT>(atoi(line->words[2])));
        printf("\tChange of environmental variation to a local gaussians variation with seed=%d\n", atoi(line->words[2]));
      }
      else {
        Utils::ExitWithUsrMsg("unknown environment variation method");
      }
    }
    else if (strcmp(line->words[0], "SECRETION_CONTRIB_TO_FITNESS") == 0)
    {
      exp_manager->get_exp_s()->set_secretion_contrib_to_fitness(atof(line->words[1]));
    }
    else if (strcmp(line->words[0], "SECRETION_COST") == 0)
    {
      exp_manager->get_exp_s()->set_secretion_cost(atof(line->words[1]));
    }
    else if (strcmp(line->words[0], "PLASMID_MINIMAL_LENGTH") == 0)
    {
      if (not exp_manager->get_with_plasmids())
      {
        printf("ERROR: option PLASMID_MINIMAL_LENGTH has no sense because there are no plasmids in this population.\n");
        exit(EXIT_FAILURE);
      }
      int32_t plasmid_minimal_length = atoi(line->words[1]);
      for (const auto& indiv: exp_manager->get_indivs())
      {
        if (indiv->get_genetic_unit(1).get_seq_length()<plasmid_minimal_length)
        {
          printf("ERROR: there is one genetic unit with a smaller length than the new minimum.\n");
          exit(EXIT_FAILURE);
        }
        indiv->get_genetic_unit_nonconst(1).set_min_gu_length(plasmid_minimal_length);
      }
    }
    else if (strcmp(line->words[0], "PLASMID_MAXIMAL_LENGTH") == 0)
    {
      if (!exp_manager->get_with_plasmids())
      {
        printf("ERROR: option PLASMID_MAXIMAL_LENGTH has no sense because there are no plasmids in this population.\n");
        exit(EXIT_FAILURE);
      }
      int32_t plasmid_maximal_length = atoi(line->words[1]);
      for (const auto& indiv: exp_manager->get_indivs())
      {
        if (indiv->get_genetic_unit_nonconst(1).get_seq_length()>plasmid_maximal_length)
        {
          printf("ERROR: there is one genetic unit with a higher length than the new maximum.\n");
          exit(EXIT_FAILURE);
        }
        indiv->get_genetic_unit_nonconst(1).set_max_gu_length(plasmid_maximal_length);
      }
    }
    else if (strcmp(line->words[0], "CHROMOSOME_MINIMAL_LENGTH") == 0)
    {
      int32_t chromosome_minimal_length = atoi(line->words[1]);
      for (const auto& indiv: exp_manager->get_indivs())
      {
        if (indiv->get_genetic_unit_nonconst(0).get_seq_length()<chromosome_minimal_length)
        {
          printf("ERROR: there is one genetic unit with a smaller length than the new minimum.\n");
          exit(EXIT_FAILURE);
        }
        indiv->get_genetic_unit_nonconst(0).set_min_gu_length(chromosome_minimal_length);
      }
    }
    else if (strcmp(line->words[0], "CHROMOSOME_MAXIMAL_LENGTH") == 0)
    {
      int32_t chromosome_maximal_length = atoi(line->words[1]);
      for (const auto& indiv: exp_manager->get_indivs()) {
        if (indiv->get_genetic_unit_nonconst(0).get_seq_length()>chromosome_maximal_length)
        {
          printf("ERROR: there is one genetic unit with a higher length than the new maximum.\n");
          exit(EXIT_FAILURE);
        }
        indiv->get_genetic_unit_nonconst(0).set_max_gu_length(chromosome_maximal_length);
      }
    }
    else if (strcmp(line->words[0], "SEED") == 0)
    {
      int32_t seed = atoi(line->words[1]);
  
      // Change prngs
#if __cplusplus == 201103L
      sel->set_prng(make_unique<JumpingMT>(seed));
      world->set_prng(make_unique<JumpingMT>(seed));
#else
      sel->set_prng(std::make_unique<JumpingMT>(seed));
      world->set_prng(std::make_unique<JumpingMT>(seed));
#endif
  
      printf("\tChange of the seed to %d in selection and world \n",atoi(line->words[1]));
    }
    else if (strcmp(line->words[0], "MUT_SEED") == 0)
    {
      int32_t mut_seed = atoi(line->words[1]);
  
      // Change mutation prng
      world->set_mut_prng(std::make_shared<JumpingMT>(mut_seed));
      printf("\tChange of the seed to %d in mutations \n",atoi(line->words[1]));
    }
    else if (strcmp(line->words[0], "STOCH_SEED") == 0)
    {
      int32_t stoch_seed = atoi(line->words[1]);
  
      // Change stochasticity prng
      world->set_stoch_prng(std::make_shared<JumpingMT>(stoch_seed));
      printf("\tChange of the seed to %d in individuals' stochasticity \n",atoi(line->words[1]));
    }
    else if (strcmp(line->words[0], "CLONE_BEST") == 0)
    {
      exp_manager->FillGridWithClones(*(exp_manager->get_best_indiv()));
      printf("\tChange of the population for a population with %" PRId32
          " individuals, all clones of the best one\n",
          exp_manager->get_nb_indivs());
    }
    // TODO: re-enable these options
    // else if (strcmp(line->words[0], "CREATE_3_SUBPOPULATIONS_BASED_ON_NON_CODING_BASES") == 0)
    // {
    //   change_based_on_non_coding_bases_of_best_individual(pop, exp_manager, SUBPOPULATIONS_BASED_ON_NON_CODING_BASES);
    //   printf("\tChange of the population for a population with %" PRId32 " individuals in 3 equal subpopulations (A: clones of the previous best individual, B: clones of the previous best individual without any non coding bases, C: clones of the previous best individual with twice non bases\n",pop->get_nb_indivs());
    //   printf("WARNING: lineage will not work properly if called with \n");
    //   printf("         a begin generation anterior to this modification \n");
    // }
    // else if (strcmp(line->words[0], "REMOVE_NON_CODING_BASES_BEST") == 0)
    // {
    //   change_based_on_non_coding_bases_of_best_individual(pop, exp_manager, REMOVE_NON_CODING_BASES_BEST_IND);
    //   printf("\tChange of the population for a population with %" PRId32 " clones of the best individual ancestor without any non coding bases\n",pop->get_nb_indivs());
    //   printf("WARNING: lineage will not work properly if called with \n");
    //   printf("         a begin generation anterior to this modification \n");
    // }
    // else if (strcmp(line->words[0], "REMOVE_NON_CODING_BASES_POP") == 0)
    // {
    //   change_based_on_non_coding_bases_in_population(pop, exp_manager,  REMOVE_NON_CODING_BASES_POPULATION);
    //   printf("\tChange of the population for a population with %" PRId32 " individuals without any non coding bases\n",pop->get_nb_indivs());
    //   printf("WARNING: lineage will not work properly if called with \n");
    //   printf("         a begin generation anterior to this modification \n");
    // }
    // else if (strcmp(line->words[0], "DOUBLE_NON_CODING_BASES_BEST") == 0)
    // {
    //   change_based_on_non_coding_bases_of_best_individual(pop, exp_manager,  DOUBLE_NON_CODING_BASES_BEST_IND);
    //   printf("\tChange of the population for a population with %" PRId32 " clones of the best individual ancestor with twice the non coding bases number \n",pop->get_nb_indivs());
    //   printf("WARNING: lineage will not work properly if called with \n");
    //   printf("         a begin generation anterior to this modification \n");
    // }
    // else if (strcmp(line->words[0], "DOUBLE_NON_CODING_BASES_POP") == 0)
    // {
    //   change_based_on_non_coding_bases_in_population(pop, exp_manager, DOUBLE_NON_CODING_BASES_POPULATION);
    //   printf("\tChange of the population for a population with %" PRId32 " individuals with twice the non coding bases number\n",pop->get_nb_indivs());
    //   printf("WARNING: lineage will not work properly if called with \n");
    //   printf("         a begin generation anterior to this modification \n");
    // }
    else
    {
      printf("%s:%d: error: the change %s is not implemented yet \n", __FILE__, __LINE__, line->words[0]);
      exit(EXIT_FAILURE);
    }

    delete line;
  }
  fclose(param_file);

  printf("OK\n");

  if (phen_target_change) {
      // The current version doesn't allow for phenotypic variation nor for
      // different phenotypic targets among the grid
      if (not exp_manager->world()->phenotypic_target_shared())
          Utils::ExitWithUsrMsg("sorry, aevol_modify has not yet been implemented "
                                        "for per grid-cell phenotypic target");
      auto phenotypicTargetHandler =
              exp_manager->world()->phenotypic_target_handler();
      phenotypicTargetHandler->set_gaussians(new_gaussians);
      phenotypicTargetHandler->BuildPhenotypicTarget();
  }

  // 9) Save the modified experiment
  if (start_to_record_tree)
  {
    if (!set_tree_step)
    {
      printf("WARNING: you modifed parameter RECORD_TREE without specifying "
          "TREE_STEP in the same parameter modification file. TREE_STEP will "
          "be set to its default value even if you previously gave another "
          "value.\n");
    }
    exp_manager->get_output_m()->init_tree(exp_manager, tree_step);
  }

  if (take_care_of_the_tree)
  {
    printf("Save the modified replication reports into tree...\t");

    #ifdef __REGUL
      sprintf(tree_file_name,"tree/tree_%06" PRId64 ".rae", num_gener);
    #else
      sprintf(tree_file_name,"tree/tree_%06" PRId64 ".ae", num_gener);
    #endif
    gzFile tree_file = gzopen(tree_file_name, "w");
    tree->write_to_tree_file(tree_file);
    gzclose(tree_file);
    printf("OK\n");
  }
  printf("Save the modified experiment into backup...\t");
  exp_manager->WriteSetupFiles();
  exp_manager->WriteDynamicFiles();
  printf("OK\n");

  delete exp_manager;
}


/*!
  \brief Get a line in a file and format it
  
  \param param_file file with param in which a line is reading
  \return line (pointer)
  
  \see format_line(ParameterLine* formated_line, char* line, bool* line_is_interpretable)
*/
ParameterLine* get_line(FILE* param_file)
{
  char line[255];
  ParameterLine* formated_line = new ParameterLine();

  bool found_interpretable_line = false; 

  while (!feof(param_file) && !found_interpretable_line)
    {
      if (!fgets(line, 255, param_file))
        {
          delete formated_line;
          return NULL;
        }
      format_line(formated_line, line, &found_interpretable_line);
    }

  if (found_interpretable_line)
    {
      return formated_line;
    }
  else
    {
      delete formated_line;
      return NULL;
    }
}

/*!
  \brief Format a line by parsing it and the words inside
  
  \param formated_line the resulted formated line
  \param line original line in char*
  \param line_is_interpretable boolean with about the possible intrepretation of the line
*/
void format_line(ParameterLine* formated_line, char* line, bool* line_is_interpretable)
{
  int16_t i = 0;
  int16_t j;

  // Parse line
  while (line[i] != '\n' && line[i] != '\0' && line[i] != '\r')
    {
      j = 0;
    
      // Flush white spaces and tabs
      while (line[i] == ' ' || line[i] == 0x09) i++; // 0x09 is the ASCII code for TAB
    
      // Check comments
      if (line[i] == '#') break;

      // If we got this far, there is content in the line
      *line_is_interpretable = true;

      // Parse word
      while (line[i] != ' '  && line[i] != '\n' && line[i] != '\0' && line[i] != '\r')
        {
          formated_line->words[formated_line->nb_words][j++] = line[i++];
        }

      // Add '\0' at end of word if it's not empty (line ending with space or tab)
      if (j != 0)
        {
          formated_line->words[formated_line->nb_words++][j] = '\0';
        }
    }
}






// /*!
//   \brief Change in the population based on non coding bases on the best individual. 3 types of changes
  
//   SUBPOPULATIONS_BASED_ON_NON_CODING_BASES:
//   Create the 3 subpopulations in the population. The definition of 3 subpopulations is based on non coding bases.
    
//   The subpopulation are clonal and based on the ancestor of best individual of pop at begin.
//   The individuals in first subpopulation are clones of the best individual. 
//   The individuals in second subpopulation are clones of the best individual without any bases that are not in coding RNA.  
//   The individuals in third subpopulation are clones of the best individual with addition of bases that are not in coding RNA to double them.
    
//   pop is changed into the new population with the 3 subpopulations
    
//   REMOVE_NON_CODING_BASES_BEST_IND: 
//   The individuals of the new population are clones of the best individual but without any bases that are not in coding RNA.  
    
//   DOUBLE_NON_CODING_BASES_BEST_IND:
//   The individuals of the new population are clones of the best individual but with addition of bases that are not in coding RNA to double them.
    
//   \param pop population to change
//   \param exp_m global exp_manager
//   \param type type of change in the population
// */
// void change_based_on_non_coding_bases_of_best_individual(ae_population* pop, ExpManager* exp_m, population_change_type type)
// {
//   if(type == SUBPOPULATIONS_BASED_ON_NON_CODING_BASES || type == REMOVE_NON_CODING_BASES_BEST_IND || type == DOUBLE_NON_CODING_BASES_BEST_IND)
//     {
//       // 1) Compute the population size
//       int32_t subpopulation_size = (int)floor(pop->get_nb_indivs()/3);
  
//       // 2) Get the best individual
//       ae_individual* best_indiv = exp_m->get_best_indiv();

    
//       // 3) Create the new population 

   
//       std::list<ae_individual*> new_generation;

//       ae_individual* indiv = create_clone(best_indiv, -1);
            
//       ae_individual* only_coding_indiv = create_clone(best_indiv, -1); //one individual being the clone of the chosen individual but without any non coding bases
//       only_coding_indiv->remove_non_coding_bases();
    
//       ae_individual* twice_non_coding_indiv = create_clone(best_indiv, -1); //one individual being the clone of the chosen individual but without any non coding bases
//       twice_non_coding_indiv->double_non_coding_bases();
    
    
//       int32_t* probe_A = new int32_t[5];
//       int32_t* probe_B = new int32_t[5];
//       int32_t* probe_C = new int32_t[5];
//       for(int32_t i = 0 ; i<5; i++)
//         {
//           probe_A[i] = 1;
//           probe_B[i] = 10;
//           probe_C[i] = 100;
//         }
//       indiv->set_int_probes(probe_A);
//       only_coding_indiv->set_int_probes(probe_B);
//       twice_non_coding_indiv->set_int_probes(probe_C);
    
//       double* probe_double_A = new double[5];
//       double* probe_double_B = new double[5];
//       double* probe_double_C = new double[5];
//       for(int32_t i = 0 ; i<5; i++)
//         {
//           probe_double_A[i] = 1;
//           probe_double_B[i] = 10;
//           probe_double_C[i] = 100;
//         }
//       indiv->set_double_probes(probe_double_A);
//       only_coding_indiv->set_double_probes(probe_double_B);
//       twice_non_coding_indiv->set_double_probes(probe_double_C);
      

//       switch(type)
//         {
//         case SUBPOPULATIONS_BASED_ON_NON_CODING_BASES:
//           {
//             int32_t  index_new_indiv = 0;
//             for (int32_t i = 0 ; i < subpopulation_size ; i++) // clones of the 3 individuals
//               {
//                 new_generation.push_back(create_clone(indiv, index_new_indiv++));
//                 new_generation.push_back(create_clone(only_coding_indiv, index_new_indiv++));
//                 new_generation.push_back(create_clone(twice_non_coding_indiv, index_new_indiv++));
//               }
//             break;
//           }
//         case REMOVE_NON_CODING_BASES_BEST_IND:
//           {
//             for (int32_t i = 0 ; i < pop->get_nb_indivs() ; i++)
//               {
//                 new_generation.push_back(create_clone(only_coding_indiv, i));
//               }
//             break;
//           }
//         case DOUBLE_NON_CODING_BASES_BEST_IND:
//           {
//             for (int32_t i = 0 ; i < pop->get_nb_indivs() ; i++)
//               {
//                 new_generation.push_back(create_clone(twice_non_coding_indiv, i));
//               }
//             break;
//           }
//         default:
//           {
//             fprintf(stderr, "%s:%d: error: wrong population_change_type %d\n", __FILE__, __LINE__, type);
//             exit(EXIT_FAILURE);
//             break;
//           }
//         }
      
//       //  4) Replace the current population by the new one
//       //     -> Useless since it is done by replace_population.
//       pop->replace_population(std::move(new_generation));
    
    

//       // TODO
//       // If the population is spatially structured, set each individual's position
//       // There will be a problem however for the "3 subpopulations" type of change,
//       // if the population size has changed (which is likely given that we do not 
//       // generally used population size that are multiple of 3)

//       pop->evaluate_individuals(exp_m->get_env());
//       pop->sort_individuals();
//     }
//   else
//     {
//       printf("%s:%d: error: wrong population_change_type %d\n", __FILE__, __LINE__, type);
//       exit(EXIT_FAILURE);
//     }
// }

// /*!
//   \brief Change in the population based on non coding bases. 2 types of changes
  
//   REMOVE_NON_CODING_BASES_POPULATION:
//   The individual of the new population are the individuals without any bases that are not in coding RNA.
    
//   DOUBLE_NON_CODING_BASES_POPULATION:
//   The individual of the new population are the individuals with addition of bases that are not in coding RNA to double them.
    
//   \param pop population to change
//   \param exp_m global exp_manager
//   \param type type of change in the population
// */
// void change_based_on_non_coding_bases_in_population(ae_population* pop, ExpManager* exp_m, population_change_type type)
// {
//   if(type == REMOVE_NON_CODING_BASES_POPULATION || type == DOUBLE_NON_CODING_BASES_POPULATION)
//     {
//       for (auto& indiv: pop->get_indivs())
//         if (type == REMOVE_NON_CODING_BASES_POPULATION)
//           indiv->remove_non_coding_bases();
//         else
//           indiv->double_non_coding_bases();
//     }
//   else
//     {
//       printf("%s:%d: error: wrong population_change_type %d\n", __FILE__, __LINE__, type);
//       exit(EXIT_FAILURE);
//     }
// }







/*!
  \brief 
  
*/
void print_help(char* prog_path) 
{
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) prog_name++;
  else prog_name = prog_path;
  
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
  printf("%s: modify an experiment as specified in param_file.\n", prog_name);
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s -g GENER [-f param_file]\n", prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n\n");
  printf("  -V, --version\n\tprint version number, then exit\n\n");
  printf("  -g, --gener GENER\n\tspecify generation number\n\n");
  printf("  -f, --file param_file\n\tspecify parameter file (default: param.in)\n");
}
