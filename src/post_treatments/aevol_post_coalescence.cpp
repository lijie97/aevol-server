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
#include <inttypes.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <sys/stat.h>

#include <list>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <fstream>
// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"

using namespace aevol;

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);

int main(int argc, char** argv)
{
  printf("\n  WARNING : Parameters' change in the middle of a simulation is not managed.\n");


  // =====================
  //  Parse command line
  // =====================

  // Default values
  bool verbose = false;

  char tree_file_name[50];

  const char * short_options = "hVv::";
  static struct option long_options[] = {
    {"help",      no_argument,       NULL,  'h'},
    {"version",   no_argument,       NULL,  'V'},
    {"verbose",   no_argument,       NULL,  'v'},
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
      //case 'n' : check_genome = NO_CHECK;           break;
      //case 'c' : check_genome = FULL_CHECK;         break;
      //case 'b' : t0  = atol(optarg);                break;
      //case 'i' : final_indiv_index  = atol(optarg); break;
      //case 'r' : final_indiv_rank  = atol(optarg);  break;
    }
  }

  // There should be only one remaining arg: the lineage file
  if (optind != argc - 1) {
    Utils::ExitWithUsrMsg("please specify a lineage file");
  }




    char* lineage_file_name = new char[strlen(argv[optind]) + 1];
    sprintf(lineage_file_name, "%s", argv[optind]);

  // =======================
  //  Open the lineage file
  // =======================
  gzFile lineage_file = gzopen(lineage_file_name, "r");
  if (lineage_file == Z_NULL) {
    fprintf(stderr, "ERROR : Could not read the lineage file %s\n", lineage_file_name);
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




    // Load the simulation
    ExpManager* exp_manager = new ExpManager();
    exp_manager->load(t0, true, false);

    // Check that the tree was recorded
    if (not exp_manager->record_tree()) {
        Utils::ExitWithUsrMsg("The phylogenetic tree wasn't recorded during "
                              "evolution, could not reconstruct the lineage");
    }

    int64_t tree_step = exp_manager->tree_step();

    //delete exp_manager;


    // The tree
    Tree* tree = NULL;


    // =========================
    //  Load the last tree file
    // =========================

    GridCell* grid_cell = new GridCell(lineage_file, exp_manager, nullptr);
    auto* indiv = grid_cell->individual();
    int32_t index = indiv->id();


    int nb_muts = 0;

  ReplicationReport* rep_f = nullptr;

  World* world = exp_manager->world();
  int16_t grid_width  = world->width();
  int16_t grid_height = world->height();
  unsigned int pop_size = grid_height * grid_width;


  std::vector<int> coalescence_time;
  coalescence_time.resize(t_end);

  aevol::AeTime::set_time(t0);
    std::ofstream coalescence_file;
    coalescence_file.open("coalescence.csv",std::ofstream::trunc);
    coalescence_file<<"Generation,"<<"Coalescence"<<std::endl;

    map<int,Tree*> map_tree;

  while (time() < t_end)
  {
    if (verbose)
      printf("Computing Coalescence at generation %" PRId64
    " for the lineage (index %" PRId32 ")...", time(), index);

    if (time() != t0) {
#ifdef __REGUL
      printf("Coalescence is not supported yet for RAevol\n");
      exit(-1);
#else
      rep_f = new ReplicationReport(lineage_file, indiv);
#endif
      index = rep_f->id(); // who we are building...
              //printf("Update index %d\n",index);


      // For each genetic unit, replay the replication (undergo all mutations)
      const auto &dnarep = rep_f->dna_replic_report();

      nb_muts = dnarep.rearrangements().size() + dnarep.mutations().size();
    }

    if (nb_muts >= 1) {
      // Search the coalescence time for this individual

      std::vector<int> previous;
      previous.push_back(index);
      std::vector<int> current;

      bool coal_found = false;
      int coal_time = 1;
      int64_t local_time = time()+1;

      for (auto t : map_tree) {
          if (t.first <= time()) {
              delete t.second;

              map_tree.erase(t.first);
          }
      }

      delete map_tree[time()];

      if (map_tree.find(((int) ((local_time - 1) / tree_step) + 1) * tree_step) == map_tree.end()) {
              sprintf(tree_file_name, "tree/tree_"
              TIMESTEP_FORMAT
              ".ae", ((int) ((local_time - 1) / tree_step) + 1) * tree_step);
              map_tree[((int) ((local_time - 1) / tree_step) + 1) * tree_step] = new Tree(exp_manager, tree_file_name);
              tree = map_tree[((int) ((local_time - 1) / tree_step) + 1) * tree_step];

              printf("Loading tree %ld\n",((int) ((local_time - 1) / tree_step) + 1) * tree_step);
      } else
          tree = map_tree[((int) ((local_time - 1) / tree_step) + 1) * tree_step];

      while (!coal_found) {
          if (local_time >= t_end)
            break;


              if (Utils::mod(local_time-1, tree_step) == 0) {

                  if (map_tree.find(((int) ((local_time - 1) / tree_step) + 1) * tree_step) == map_tree.end()) {
                          sprintf(tree_file_name, "tree/tree_"
                          TIMESTEP_FORMAT
                          ".ae", ((int) ((local_time - 1) / tree_step) + 1) * tree_step);
                          map_tree[((int) ((local_time - 1) / tree_step) + 1) * tree_step] = new Tree(exp_manager, tree_file_name);
                          tree = map_tree[((int) ((local_time - 1) / tree_step) + 1) * tree_step];
                          printf("Loading tree %ld\n",((int) ((local_time - 1) / tree_step) + 1) * tree_step);
                  } else {
                      tree = map_tree[((int) ((local_time - 1) / tree_step) + 1) * tree_step];
                  }
              }


        ReplicationReport** reports = tree->reports(local_time);

          #pragma omp parallel for
          for (int i = 0; i < pop_size; i++) {
            ReplicationReport* rep =  new ReplicationReport(*(reports[i]));

            auto foundPrevious = std::find(previous.begin(),previous.end(),rep->parent_id());

            if ( foundPrevious != previous.end() ) {
              #pragma omp critical
              {
                current.push_back(rep->id());
              }
            }

            delete rep;
          }


        if (current.size() == pop_size) {
          coalescence_time[time()] = coal_time;
          coal_found = true;
        } else {
          local_time++;
          coal_time++;
          previous.swap(current);
          current.clear();
        }


        }



    delete rep_f;

      coalescence_file<<AeTime::time()<<","<<coalescence_time[AeTime::time()]<<std::endl;
    }

    aevol::AeTime::plusplus();
    if (verbose) printf(" OK\n");
  }





//  for (int gen = 0; gen < t_end; gen++) {
//
//  }

  coalescence_file.flush();
  coalescence_file.close();

  //delete exp_manager;

    free(lineage_file_name);

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
