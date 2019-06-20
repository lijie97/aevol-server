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

class Node {
 public:
    Node(unsigned long long lid) { id = lid; };

    unsigned long long id;
    std::unordered_map<unsigned long long, Node*> next_nodes;
    Node* root = nullptr;
    int dist_to_parent = 0;
    bool to_delete = false;
    bool is_last = false;
    std::string nhx = "";
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
  //check_type  check_genome      = LIGHT_CHECK;
  bool verbose = false;
  int64_t t0 = 0;
  int64_t t_end = -1;

  std::unordered_map<unsigned long long, Node*> allNodes;
  std::vector<unsigned long long> last;

  char tree_file_name[50];

  const char * short_options = "hVv:e:";
  static struct option long_options[] = {
    {"help",      no_argument,       NULL,  'h'},
    {"version",   no_argument,       NULL,  'V'},
    {"verbose",   no_argument,       NULL,  'v'},
    {"end",       required_argument,  NULL, 'e'},
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
      case 'b' : t0  = atol(optarg);                break;
      //case 'i' : final_indiv_index  = atol(optarg); break;
      //case 'r' : final_indiv_rank  = atol(optarg);  break;
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

  if (verbose)
  {
    printf("\n\n");
    printf("====================================\n");
    printf(" Loading the last tree file ... ");
    fflush(stdout);
  }


  // Example for ae_common::rec_params->tree_step() == 100 :
  //
  // tree_000100.ae ==>  generations   1 to 100.
  // tree_000200.ae ==>  generations 101 to 200.
  // tree_000300.ae ==>  generations 201 to 300.
  // etc.
  //
  // Thus, the information for generation end_gener are located
  // in the file called (end_gener/ae_common::rec_params->tree_step() + 1) * ae_common::rec_params->tree_step(),
  // except if end_gener%ae_common::rec_params->tree_step()==0.

  sprintf(tree_file_name,"tree/tree_%06" PRId64 ".ae", t_end);

  tree = new Tree(exp_manager, tree_file_name);

  if (verbose)
  {
    printf("OK\n");
    printf("====================================\n");
  }


  World* world = exp_manager->world();
  int16_t grid_width  = world->width();
  int16_t grid_height = world->height();
  int32_t pop_size = grid_height * grid_width;

  // ============================================================================
  //  Find the index of the final individual and retrieve its replication report
  // ============================================================================
  std::vector<unsigned long long> previous;
  std::vector<unsigned long long> current;

  for (int16_t x = 0 ; x < grid_width ; x++)
    for (int16_t y = 0 ; y < grid_height ; y++) {
      ReplicationReport* rep =  new ReplicationReport(*(tree->report_by_index(t_end,
                                                                              x * grid_height + y)));

      auto foundCurrent = std::find(current.begin(),current.end(),rep->id());

      if ( foundCurrent == current.end() ) {

        auto found = allNodes.find(rep->id());

        Node* child;
        if (found == allNodes.end()) {
          child = new Node(rep->id());
          allNodes[child->id] = child;
          last.push_back(child->id);
          child->is_last = true;

          current.push_back(child->id);
        } else {
          child = found->second;
        }

        auto foundFather_prev = std::find(previous.begin(),previous.end(),rep->parent_id());

        if (foundFather_prev == previous.end()) {
          previous.push_back(rep->parent_id());
        }


        if (rep->id() != rep->parent_id()) {
          auto foundFather = allNodes.find(rep->parent_id());

          if (foundFather == allNodes.end()) {
            Node* father = new Node(rep->parent_id());
            child->dist_to_parent = 1;
            child->root = father;
            father->next_nodes[child->id] = child;
            allNodes[father->id] = father;
          } else {
            Node* father = foundFather->second;
            auto foundChild = father->next_nodes.find(rep->id());

            if (foundChild == allNodes.end()) {
              father->next_nodes[child->id] = child;
              child->root = father;
            }
          }
        } else {
          child->dist_to_parent++;
        }

      }

      delete rep;

    }

  // =======================
  //  Open the output file
  // =======================



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

  // For each generation (going backwards), retrieve the index of the parent and
  // the corresponding replication report
  for (int64_t t = t_end - 1 ; t > 0 ; t--)
  {

    if (verbose)
      printf("Getting the replication report for the ancestor at generation %" PRId64 "\n", t);

    // If we've exhausted the current tree file, load the next one
    if (Utils::mod(t, tree_step) == 0)
    {
      // Change the tree file
      delete tree;


      sprintf(tree_file_name,"tree/tree_%06" PRId64 ".ae", t);

      tree = new Tree(exp_manager, tree_file_name);
    }

    std::vector<unsigned long long> filter = previous;
    previous.clear();
    current.clear();

    for (int16_t x = 0 ; x < grid_width ; x++)
      for (int16_t y = 0 ; y < grid_height ; y++) {
        ReplicationReport* rep =  new ReplicationReport(*(tree->report_by_index(t,
                                                                                x * grid_height + y)));


        auto foundFilter = std::find(filter.begin(),filter.end(),rep->id());

        if ( foundFilter != filter.end() ) {

          auto foundCurrent = find(current.begin(),current.end(),rep->id());

          if (foundCurrent == current.end()) {

            auto found = allNodes.find(rep->id());

            Node* child;

            if (found == allNodes.end()) {
              child = new Node(rep->id());
              allNodes[child->id] = child;
              current.push_back(child->id);
            } else {
              child = found->second;
            }

            auto foundFather_prev = std::find(previous.begin(),previous.end(),rep->parent_id());

            if (foundFather_prev == previous.end()) {
              previous.push_back(rep->parent_id());
            }


            if (rep->id() != rep->parent_id()) {
              auto foundFather = allNodes.find(rep->parent_id());

              if (foundFather == allNodes.end()) {
                Node* father = new Node(rep->parent_id());
                child->dist_to_parent = 1;
                child->root = father;
                father->next_nodes[child->id] = child;
                allNodes[father->id] = father;
              } else {
                Node* father = foundFather->second;
                auto foundChild = father->next_nodes.find(
                    rep->id());

                if (foundChild == allNodes.end()) {
                  father->next_nodes[child->id] = child;
                  child->root = father;
                }
              }
            } else {
              child->dist_to_parent++;
            }

          }
        }

        delete rep;

      }

  }


  printf("Loading data complete, Building Phylogenetic Tree...\n");

  std::unordered_map<unsigned long long,Node*> cleaned;

  printf("All nodes %ld\n",allNodes.size());

  bool stillOne=true;
  while (stillOne) {
    bool haveOne = false;

    for (auto nodX : allNodes) {
      if (nodX.second->next_nodes.size() == 1 && nodX.second->root != nullptr
            && !nodX.second->to_delete)
        haveOne = true;

    }

    if (!haveOne)
      break;

    for (auto nodX : allNodes) {

      Node* nod = nodX.second;

      //printf("Node %d : childs %ld\n", nod->id, nod->next_nodes.size());

      if (nod != nullptr && nod->root != nullptr) {
        if (nod->next_nodes.size() == 1) {
          //printf("Processing\n");

          Node* father = nod->root;//
          Node* child = nod->next_nodes.begin().operator->()->second;

          child->dist_to_parent += nod->dist_to_parent;

          //printf("Removing %d\n",nod->id);

          child->root = father;

          //printf("Father child %ld BEFORE\n",father->next_flushnodes.size());
          father->next_nodes[child->id] = child;
          //printf("Father child %ld DURING\n",father->next_nodes.size());
          father->next_nodes.erase(nodX.first);
          //printf("Father child %ld AFTER\n",father->next_nodes.size());

          nod->to_delete = true;
//        cleaned[child->id] = child;
        }
      }/* else
        cleaned[nod->id] = nod;*/
    }
  }




  for (auto node : allNodes) {
    if (! node.second->to_delete) {
      printf("Node %d : childs %ld distance %d\n", node.second->id,
             node.second->next_nodes.size(), node.second->dist_to_parent);
      cleaned[node.second->id] = node.second;
    }


    if (node.second->is_last && node.second->to_delete) {
      auto foundCurrent = find(last.begin(), last.end(), node.second->id);

      if (foundCurrent != last.end()) {
        last.erase(foundCurrent);
      }
    }
  }

  allNodes.clear();
  allNodes = cleaned;

  printf("All nodes %ld\n",allNodes.size());


  printf("Phylogenetic Tree : Distance computed\n");

  delete exp_manager;

  std::unordered_map<unsigned long long,Node*> beforeLast;

  for (unsigned long long lnode_id : last) {
    Node* nod = allNodes[lnode_id];
    beforeLast[nod->root->id] = nod->root;
  }


  std::unordered_map<unsigned long long,Node*> before = beforeLast;
  std::unordered_map<unsigned long long,Node*> next;
  std::string global = "";
  Node* root;
  while (before.size() > 0) {

    for (auto cell : before) {
      std::string local = "(";
      bool first = true;
      for (auto child_cell : cell.second->next_nodes) {
        //std::cout<<"Processing "<<child_cell.second->id<<std::endl;

        if (!first)
          local=local+",";
        else
          first = false;
        //std::cout<<"CHILD NHX "<<child_cell.second->nhx<<" LOCAL "<<local<<" Distance "<<std::to_string(child_cell.second->dist_to_parent)<<std::endl;

        if (child_cell.second->nhx == "")
          local=local+std::to_string(child_cell.second->id)+":"+std::to_string(child_cell.second->dist_to_parent);
        else {
          local = local + child_cell.second->nhx + ":" +
                  std::to_string(child_cell.second->dist_to_parent);
        }

//        std::cout<<local<<std::endl;
      }
      local=local+")";
      cell.second->nhx = local;

      printf("-----------------------------------------------------------------------------------------------\n");
        std::cout<<local<<std::endl;
      printf("-----------------------------------------------------------------------------------------------\n");

      if (cell.second->root != nullptr)
        next[cell.second->root->id] = cell.second->root;
    }
    if (next.size() == 0) {
      for (auto cell : before) {
        root = cell.second;

        if (root->next_nodes.size() == 1) {
          root = root->next_nodes.begin().operator->()->second;
        }
        std::cout<<root->nhx<<":"<<root->dist_to_parent<<std::endl;

        std::ofstream myfile;
        myfile.open ("output.nh");
        myfile << root->nhx<<":"<<root->dist_to_parent<<";"<<std::endl;
        myfile.close();

      }
      break;
    }
    before = next;
    next.clear();
  }

  printf("Root is %ld\n",root->id);
  bool notMoreChild=false;
  std::vector<Node*> nextChilds;
  for (auto child : root->next_nodes) {
    nextChilds.push_back(child.second);
  }

  while(!notMoreChild) {
    printf("-----------------------------------------------------------------------------------------------\n");
    printf("Next tree level:\n");
    std::vector<Node*> nextChilds2;

    for (auto child : nextChilds) {
      printf("%ld ",child->id);
      for (auto gchild : child->next_nodes) {
        nextChilds2.push_back(gchild.second);
      }
    }
    printf("\n");

    nextChilds = nextChilds2;

    if (nextChilds.size()==0)
      break;
  }

  if (verbose)  printf("OK\n");

  // Dump the tre into NHX format

  //delete exp_manager;

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
