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


using std::list;
using namespace aevol;

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);





// =====================================================================
//                         Main Function
// =====================================================================


int main(int argc, char* argv[])
{
  // ----------------------------------------
  //     command-line option parsing
  // ----------------------------------------
  int32_t nb_children     = 1000;
  int32_t wanted_rank     = -1;
  int32_t wanted_index    = -1;
  int64_t num_gener       = 0;

  const char * options_list = "hVg:n:r:i:";
  static struct option long_options_list[] = {
    {"help",          no_argument,        NULL, 'h'},
    {"version",       no_argument,        NULL, 'V'},
    {"gener",         required_argument,  NULL, 'g'},
    {"nb-children",   required_argument,  NULL, 'n'},
    {"rank",          required_argument,  NULL, 'r'},
    {"index",         required_argument,  NULL, 'i'},
    {0, 0, 0, 0}
  };

  int option = -1;
  bool rank_already_set = false;
  bool index_already_set = false;
  while((option=getopt_long(argc,argv,options_list,long_options_list,NULL))!=-1)
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
        case 'g' :
          {
            if (strcmp(optarg, "") == 0)
              {
                printf("%s: error: Option -g or --gener : missing argument.\n", argv[0]);
                exit(EXIT_FAILURE);
              }

            num_gener = atol(optarg);
            break;
          }
        case 'n' :
          nb_children = atol(optarg);
          break;
        case 'r' :
          if (index_already_set) 
            {
              fprintf(stderr, "%s: error: Options -r and -i are incompatible. Please choose one of them only.\n", argv[0]);
              exit(EXIT_FAILURE);
            }
          wanted_rank = atol(optarg);
          rank_already_set = true;
          break;
        case 'i' :
          if (rank_already_set) 
            {
              fprintf(stderr, "%s: error: Options -r and -i are incompatible. Please choose one of them only.\n", argv[0]);
              fprintf(stderr, "           Use %s --help for more information.\n", argv[0]);
              exit(EXIT_FAILURE);
            }
          wanted_index = atol(optarg);
          index_already_set = true;
          break;
        }
    }



  // -----------------------------------------------------------------------------------
  //  Load the backup and get the individual for which detailed information is desired
  // -----------------------------------------------------------------------------------

  ae_exp_manager* exp_manager = new ae_exp_manager();
  exp_manager->load(num_gener, true, false);

  if ((wanted_rank == -1) && (wanted_index == -1)) 
    {
      wanted_rank = exp_manager->get_nb_indivs();  // the best one has rank N
    }

  // TODO: factor with duplicated code in mutagenesis.cpp
  ae_individual* wanted_indiv = nullptr;
  { // (local scope for `indivs` used as a shorthand)
    bool found = false;
    int32_t current_rank = -1;
    int32_t current_index = -1;
    std::list<ae_individual*> indivs = exp_manager->get_indivs();
    for (auto indiv = indivs.rbegin(); not found and indiv != indivs.rend(); ++indiv) {
      current_index = (*indiv)->get_id();
      current_rank = (*indiv)->get_rank();

      if (wanted_index != -1 and current_index == wanted_index) {
        found = true;
        wanted_indiv = *indiv;
        wanted_rank = current_rank;
      }
      else if (current_rank == wanted_rank) {
        // no index was specified, we use the desired rank
        found = true;
        wanted_indiv = *indiv;
        wanted_index = current_index;
      }
    }
    assert(found);
  }

  wanted_indiv->evaluate();
  wanted_indiv->compute_statistical_data();
  wanted_indiv->compute_non_coding();



  // ----------------------
  //  Prepare the outputs
  // ----------------------

  char directory_name[64];
  sprintf(directory_name, "analysis-generation%06" PRId64, num_gener);
  
  // Check whether the directory already exists and is writable
  if (access(directory_name, F_OK) == 0)
    {
      if (access(directory_name, X_OK | W_OK) != 0)
        {
          fprintf(stderr, "Error: cannot enter or write in directory %s.\n", directory_name);
          exit(EXIT_FAILURE);
        }
    }
  else
    {
      // Create the directory with permissions : rwx r-x r-x
      if (mkdir(directory_name, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0)
        {
          fprintf(stderr, "Error: cannot create directory %s.\n", directory_name);
          exit(EXIT_FAILURE);
        }
    }
  

  char filename[256];
  snprintf( filename, 255, "%s/robustness-allindivs-g%06" PRId64 ".out", directory_name, num_gener);
  FILE * outputfile_wholepop = fopen(filename, "w");
  
  snprintf(filename, 255,
           "%s/robustness-singleindiv-details-g%06" PRId64
               "-i%" PRId32 "-r%" PRId32 ".out",
           directory_name, num_gener, wanted_index, wanted_rank);
  FILE * outputfile_details = fopen(filename, "w");


  //  Write headers

  fprintf(outputfile_wholepop, "# ######################################################################\n");
  fprintf(outputfile_wholepop, "# Robustness data of individuals at generation %" PRId64 "\n",num_gener);
  fprintf(outputfile_wholepop, "# ######################################################################\n");
  fprintf(outputfile_wholepop, "#  1.  Rank\n");
  fprintf(outputfile_wholepop, "#  2.  Index\n");
  fprintf(outputfile_wholepop, "#  3.  Fitness\n");
  fprintf(outputfile_wholepop, "#  4.  Metabolic error\n");
  fprintf(outputfile_wholepop, "#  5.  Genome size\n");
  fprintf(outputfile_wholepop, "#  6.  Functional gene number\n");
  fprintf(outputfile_wholepop, "#  7.  Reproduction probability\n");
  fprintf(outputfile_wholepop, "#  8.  Proportion of neutral offsprings\n");
  fprintf(outputfile_wholepop, "#  9.  Proportion of beneficial offsprings\n");
  fprintf(outputfile_wholepop, "#  10. Proportion of deleterious offsprings\n");
  fprintf(outputfile_wholepop, "#  11. Theoretical proportion of neutral offsprings\n");
  fprintf(outputfile_wholepop, "#  12. Fitness mean of offsprings\n");
  fprintf(outputfile_wholepop, "#  13. Fitness variance of offsprings\n");
  fprintf(outputfile_wholepop, "#  14. Genome size mean of offsprings\n");
  fprintf(outputfile_wholepop, "#  15. Genome size variance of offsprings\n");
  fprintf(outputfile_wholepop, "#  16. Functional gene number mean of offsprings\n");
  fprintf(outputfile_wholepop, "#  17. Functional gene number variance of offsprings\n");
  fprintf(outputfile_wholepop, "# ######################################################################\n");
  
  fprintf(outputfile_details, "# #######################################################################################################\n");
  fprintf(outputfile_details, "#  Offspring details of individual with rank %" PRId32 " and index %" PRId32 " at generation %" PRId64 " \n", \
          wanted_rank, wanted_index, num_gener);
  fprintf(outputfile_details, "# #######################################################################################################\n");
  fprintf(outputfile_details, "#  1.  Fitness\n");
  fprintf(outputfile_details, "#  2.  Metabolic error\n");
  fprintf(outputfile_details, "#  3.  Genome size\n");
  fprintf(outputfile_details, "#  4.  Functional gene number\n");
  fprintf(outputfile_details, "#  5.  Number of coding bases\n");
  fprintf(outputfile_details, "#  6.  Number of transcribed but not translated bases\n");
  fprintf(outputfile_details, "#  7.  Number of non transcribed bases\n");
  fprintf(outputfile_details, "# #######################################################################################################\n");



  // ---------------------------------------------------------------------------
  //  Force each individual of the population to produce nb_children offspring 
  // ---------------------------------------------------------------------------

  
  double reproduction_statistics[3];
  double offsprings_statistics[6];
  double th_fv;
  
  // exp_manager->get_exp_s()->get_sel()->compute_prob_reprod();
  // double* tmp_reprod = exp_manager->get_exp_s()->get_sel()->prob_reprod();

  { // (local scope for `indivs` used as a shorthand)
    list<ae_individual*> indivs = exp_manager->get_indivs();
    for (auto indiv : indivs)
    {
      int32_t current_index = indiv->get_id();
      int32_t current_rank = indiv->get_rank();

      // Compute Fv ----------------------------------------------------------------
      th_fv = indiv->compute_theoritical_f_nu();

      if (indiv == wanted_indiv)
        indiv->compute_experimental_f_nu(nb_children,
                                         reproduction_statistics,
                                         offsprings_statistics,
                                         outputfile_details);
      else
        indiv->compute_experimental_f_nu(nb_children,
                                         reproduction_statistics,
                                         offsprings_statistics);


      // Write to file -------------------------------------------------------------
      fprintf(outputfile_wholepop,
              "%" PRId32 " %" PRId32 " %le %le %" PRId32 " %" PRId32
              " %le %le %le %le %le %le %le %le %le %le\n",
              current_rank,
              current_index,
              indiv->get_fitness(),
              indiv->dist_to_target_by_feature(METABOLISM),
              indiv->get_total_genome_size(),
              indiv->get_nb_functional_genes(),
              reproduction_statistics[0],
              reproduction_statistics[1],
              reproduction_statistics[2],
              th_fv,
              offsprings_statistics[0],
              offsprings_statistics[1],
              offsprings_statistics[2],
              offsprings_statistics[3],
              offsprings_statistics[4],
              offsprings_statistics[5]
             );
    }
  }



  fclose(outputfile_wholepop);
  fclose(outputfile_details); 
  delete exp_manager;

  
  return EXIT_SUCCESS;
}




/*!
  \brief

*/
void print_help(char* prog_path)
{
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) prog_name++;
  else prog_name = prog_path;

  printf("\n");
  printf("*********************** aevol - Artificial Evolution ******************* \n");
  printf("*                                                                      * \n");
  printf("*                     Robustness post-treatment program                * \n");
  printf("*                                                                      * \n");
  printf("************************************************************************ \n");
  printf("\n\n");
  printf("This program is Free Software. No Warranty.\n");
  printf("\n");
  printf("Usage : %s -h\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s -g GENER [-n NBCHILDREN] [-r RANK | -i INDEX]\n", prog_name);
  printf("\n");
  printf("This program computes the replication statistics of all the individuals of a given generation,\n");
  printf("like the proportion of neutral, beneficial, deleterious offsprings. This is done by simulating\n");
  printf("NBCHILDREN replications for each individual, with its mutation, rearrangement and transfer rates.\n");
  printf("Depending on those rates and genome size, there can be several events per replication.\n");
  printf("Those statistics are written in analysis-generationGENER/robustness-allindivs-gGENER.out, with one \n");
  printf("line per individual in the specified generation. \n\n");
  printf("The program also outputs detailed statistics for one of the individuals (the best one by default). \n");
  printf("The detailed statistics for this individual are written in the file called \n");
  printf("analysis-generationGENER/robustness-singleindiv-details-gGENER-iINDEX-rRANK.out, with one line per simulated\n");
  printf("child of this particular individual.\n");
  printf("\n");
  printf("\n");
  printf("\t-h or --help    : Display this help, then exit\n");
  printf("\n");
  printf("\t-V or --version : Print version number, then exit\n");
  printf("\n");
  printf("\t-g GENER or --gener GENER : \n");
  printf("\t                  Generation at which the statistics are computed\n");
  printf("\n");
  printf("\t-i INDEX or --index INDEX : \n");
  printf("\t                  Index of individual of interest. Should be comprised between 0 and N-1, where\n");
  printf("\t                  N is the size of the population.\n");
  printf("\n");
  printf("\t-r RANK or --rank RANK : \n");
  printf("\t                  Rank of individual of interest. Should be comprised between 1 and N, where\n");
  printf("\t                  N is the size of the population. Default = N (fittest individual).\n");
  printf("\n");
  printf("\t-n NBCHILDREN or --nb-children NBCHILDREN : \n");
  printf("\t                  Perform NBCHILDREN replications per individual to compute its statistics. \n");
  printf("\t                  Default = 1000.\n");
  printf("\n");


  printf("\n");
}
