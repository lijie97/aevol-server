//
// Created by arrouan on 14/11/16.
//
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

int main(int argc, char* argv[])
{
  bool verbose = false;

  int generation_to_dump;

  const char * short_options = "Vv:g:";
  static struct option long_options[] = {
      {"version",   no_argument,       NULL,  'V'},
      {"verbose",   no_argument,       NULL,  'v'},
      {"generation",       required_argument,  NULL, 'g'},
      {0, 0, 0, 0}
  };

  int option;
  while((option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
  {
    switch(option)
    {
      case 'V' :
      {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 'v' : verbose = true;                    break;
      case 'g' :
      {
        if (strcmp(optarg, "") == 0)
        {
          printf("%s: error: Option -g or --generation : missing argument.\n", argv[0]);
          exit(EXIT_FAILURE);
        }

        generation_to_dump = atoi(optarg);

        break;
      }
    }
  }

  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(generation_to_dump, true, false);


  World* world = exp_manager->world();
  int16_t grid_width  = world->width();
  int16_t grid_height = world->height();

  std::ofstream protein_concentration_file("protein_concentration.csv");
  std::ofstream rna_basal_concentration_file("rna_basal_concentration.csv");
  std::ofstream rna_produce_protein_file("rna_produce_protein.csv");
  std::ofstream rna_influence_enhancing_coef("rna_influence_enhancing_coef.csv");
  std::ofstream rna_influence_operating_coef("rna_influence_operating_coef.csv");
  std::ofstream nb_protein_signal("nb_protein_signal.csv");
  int nb_signal = 0;


  for (int16_t x = 0 ; x < grid_width ; x++)
    for (int16_t y = 0 ; y < grid_height ; y++) {
      Individual_R* indiv = dynamic_cast<Individual_R*>(world->indiv_at(x,y));
      int prot_id = 0;
      int indiv_id = x*grid_height+y;

      indiv->init_indiv();
      for (auto prot : indiv->protein_list_) {
        if (x==0 && y==0 ) {
          if (dynamic_cast<Protein_R*>(prot)->is_signal()) {
            nb_signal++;
          }
        }
        protein_concentration_file<<indiv_id<<","<<prot_id<<","<<
                  prot->concentration_<<std::endl;

        for (auto rna : dynamic_cast<Protein_R*>(prot)->_rna_R_list) {
          rna_produce_protein_file<<indiv_id<<","<<prot_id<<","<<
                            rna->get_local_id()<<std::endl;
        }
        prot_id++;
      }

      int rna_id = 0;
      for(auto& rna : indiv->_rna_list_coding) {
        rna_basal_concentration_file<<indiv_id<<","<<rna_id<<","<<
            rna->basal_level()<<std::endl;

        for (int i = 0; i < rna->_nb_influences; i++) {
          rna_influence_enhancing_coef<<indiv_id<<","<<rna_id<<","<<i<<","<<
              rna->_enhancing_coef_list[i]<<std::endl;

          rna_influence_operating_coef<<indiv_id<<","<<rna_id<<","<<i<<","<<
              rna->_operating_coef_list[i]<<std::endl;
        }

        rna_id++;
      }
    }

  nb_protein_signal<<nb_signal<<std::endl;

  protein_concentration_file.close();
  rna_basal_concentration_file.close();
  rna_influence_enhancing_coef.close();
  rna_influence_operating_coef.close();
  nb_protein_signal.close();
}
